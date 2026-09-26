"""
Observation sources for the cycle-time analysis -- one adapter per source.

    from sources import collect
    results = collect(datetime(2026, 9, 21, 12))
    print(availability_table(results))

WHAT "RELIABLE" MEANS HERE (prompt 96, 2026-09-22)

The model takes in every reliable source that is available at the cycle time,
and a source that is missing is skipped, not waited for and not faked. Each
adapter therefore returns a `SourceResult` with a STATUS, and the four
statuses are kept apart on purpose:

    ok           observations returned
    unavailable  the source has nothing for this time (no 06Z soundings; a
                 5-day buoy file that does not reach back this far)
    failed       the request or the parse raised -- a defect or an outage
    disabled     the adapter exists but is not trusted yet, and says why

"unavailable" and "failed" need different fixes, and "all rejected by QC" is
a third thing again (see observations.run_qc). Merging any two of them is how
an empty analysis gets mistaken for a quiet day.

    source   what                                   error_std (T / wind / p)
    asos     ASOS/AWOS, US states + Canada (IEM)    1.5 K / 2.0 m/s / 150 Pa
    raob     radiosondes, domain + ring (IEM)       1.0 K / 2.5 m/s / --
    ndbc     moored buoys and C-MAN (NDBC 5-day)    1.0 K / 2.0 m/s / 100 Pa
    mrms     radar mosaic -- ARCHIVED, not assimilated: the model is dry,
             and reflectivity has nothing to go into yet
    vad      NEXRAD VAD wind profiles -- DISABLED until a decoder has been
             checked against real files with the server's existing stack

Left out, and why: aircraft (AMDAR/ACARS) are access-restricted; satellite
retrievals are derived products whose height assignment uses a model
background; state RWIS road-weather stations are sited for roads, not
meteorology. Any of these can be added as an adapter later.

INTERFACE FACTS MEASURED 2026-09-22 (see P-54)

  * IEM raob.py wants `sts`/`ets`, ONE station per request, `K` prefix.
  * IEM asos.py accepts `sts`/`ets` and `mslp`/`alti`.
  * IEM's network list has look-alike codes for other countries
    (`DE__ASOS` is not Delaware); networks are named exactly, never by prefix.
"""

import csv
import io
import json
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src" / "verification"))

import numpy as np

from observations import Observation
from fetchers import (f_to_k, c_to_k, knots_to_ms, wind_to_uv, _num,
                      rh_from_dewpoint)

# Sources: Iowa Environmental Mesonet (2026a, b) for ASOS and RAOB; NOAA
# National Data Buoy Center (2026); MRMS (Zhang et al. 2016; NOAA 2026b).
# Full entries in docs/REFERENCES.md.
IEM = "https://mesonet.agron.iastate.edu"
IEM_ASOS = f"{IEM}/cgi-bin/request/asos.py"
IEM_RAOB = f"{IEM}/cgi-bin/request/raob.py"
IEM_RAOB_TABLE = f"{IEM}/geojson/network/RAOB.geojson"
NDBC_ACTIVE = "https://www.ndbc.noaa.gov/activestations.xml"
NDBC_5DAY = "https://www.ndbc.noaa.gov/data/5day2/{sid}_5day.txt"
MRMS_BUCKET = "https://noaa-mrms-pds.s3.amazonaws.com"

# States and provinces touching the analysis box. Exact network codes.
US_STATES = ["ME", "NH", "VT", "MA", "RI", "CT", "NY", "NJ", "PA", "DE",
             "MD", "VA", "WV", "OH", "MI", "IN", "KY", "NC", "TN"]
CA_PROVINCES = ["ON", "QC", "NB", "NS", "PE"]
ASOS_NETWORKS = ([f"{s}_ASOS" for s in US_STATES]
                 + [f"CA_{p}_ASOS" for p in CA_PROVINCES])

# Observation + representativeness error, per source and variable.
ERROR_STD = {
    "asos": {"TMP": 1.5, "RH": 12.0, "UGRD": 2.0, "VGRD": 2.0, "PMSL": 150.0},
    "ndbc": {"TMP": 1.0, "RH": 10.0, "UGRD": 2.0, "VGRD": 2.0, "PMSL": 100.0},
    "raob": {"TMP": 1.0, "RH": 10.0, "UGRD": 2.5, "VGRD": 2.5, "HGT": 10.0},
}

# The window ENDS at the cycle time (prompt 99: "the 00z run uses 00z
# conditions"). Nothing observed after it may enter a run, so a report 5 min
# late is not used even though it would be closer. Hourly METARs are taken at
# :51-:56, so the latest report at or before the hour is the one used.
WINDOW_MIN = 60


def in_window(t, cycle, window_min=WINDOW_MIN):
    """Minutes before the cycle time, or None if outside (cycle-window, cycle]."""
    lag = (cycle - t).total_seconds() / 60.0
    return lag if 0.0 <= lag <= window_min else None


@contextmanager
def single_try(fetcher):
    """
    One attempt per request, for sources where a missing file is NORMAL.

    PoliteFetcher retries with exponential backoff, which is right for a
    flaky server and wrong for a buoy that is simply off station: four
    attempts cost ~14 s each, and a dozen dead buoys would eat minutes of the
    1.5 h budget waiting for files that do not exist.
    """
    old = getattr(fetcher, "retries", None)
    if old is not None:
        fetcher.retries = 1
    try:
        yield fetcher
    finally:
        if old is not None:
            fetcher.retries = old


def utcnow():
    return datetime.now(timezone.utc).replace(tzinfo=None)


# ---------------------------------------------------------------------------
# Result record
# ---------------------------------------------------------------------------

@dataclass
class SourceResult:
    name: str
    status: str                        # ok | unavailable | failed | disabled
    obs: list = field(default_factory=list)
    raw: dict = field(default_factory=dict)   # name -> payload, for the archive
    message: str = ""
    seconds: float = 0.0
    platforms: int = 0

    def summary(self):
        n = len(self.obs)
        return (f"{self.name:5s} {self.status:11s} {n:6d} obs "
                f"{self.platforms:4d} platforms {self.seconds:6.1f}s  "
                f"{self.message}")


def analysis_box(domain, pad_lat=3.0, pad_lon=4.0):
    """The domain plus a ring, so observed upstream air sits at the edges."""
    return {"lat_min": domain["lat_min"] - pad_lat,
            "lat_max": domain["lat_max"] + pad_lat,
            "lon_min": domain["lon_min"] - pad_lon,
            "lon_max": domain["lon_max"] + pad_lon}


def in_box(lat, lon, box):
    return (box["lat_min"] <= lat <= box["lat_max"]
            and box["lon_min"] <= lon <= box["lon_max"])


def _obs(t, lat, lon, var, value, source, station, elev, pressure=None):
    return Observation(time=t, lat=float(lat), lon=float(lon), variable=var,
                       value=float(value), source=source, station=station,
                       elevation=None if elev is None else float(elev),
                       pressure=pressure,
                       error_std=ERROR_STD[source].get(var, 2.0))


def station_pressure_from_altimeter(alti_inhg, elev_m):
    """Station pressure (Pa) from an altimeter setting (NWS formula)."""
    a = float(alti_inhg) * 33.8639                       # hPa
    p = (a ** 0.190284 - 8.4228807e-5 * float(elev_m)) ** (1.0 / 0.190284)
    return p * 100.0


def slp_from_station_pressure(p_stn_pa, elev_m, t_k):
    """Reduce station pressure to sea level with a standard lapse rate."""
    g, rd, gamma = 9.80665, 287.05, 0.0065
    t_mean = float(t_k) + 0.5 * gamma * float(elev_m)
    return float(p_stn_pa) * np.exp(g * float(elev_m) / (rd * t_mean))


# ---------------------------------------------------------------------------
# ASOS / AWOS
# ---------------------------------------------------------------------------

def asos_url(cycle, networks=None, window_min=WINDOW_MIN):
    sts = cycle - timedelta(minutes=window_min)
    ets = cycle
    params = [("sts", sts.strftime("%Y-%m-%dT%H:%MZ")),
              ("ets", ets.strftime("%Y-%m-%dT%H:%MZ"))]
    params += [("data", d) for d in
               ("tmpf", "dwpf", "relh", "drct", "sknt", "mslp", "alti")]
    params += [("tz", "UTC"), ("format", "onlycomma"), ("latlon", "yes"),
               ("elev", "yes"), ("missing", "M"), ("trace", "T"),
               ("direct", "no"), ("report_type", "3")]
    params += [("network", n) for n in (networks or ASOS_NETWORKS)]
    return f"{IEM_ASOS}?{urlencode(params)}"


def parse_asos(text, box, cycle=None, window_min=WINDOW_MIN):
    """
    IEM ASOS CSV -> Observations: TMP, RH, UGRD, VGRD, PMSL.

    Sea-level pressure is the reported MSLP where present. Otherwise it is
    reduced from the altimeter setting: altimeter -> station pressure (the
    measured quantity) -> sea level with the observed temperature. Only the
    report nearest the cycle time is kept per station, so a station reporting
    specials every 10 minutes does not outweigh its neighbours, and nothing
    after the cycle time is kept.
    """
    best = {}
    for row in csv.DictReader(io.StringIO(text)):
        stn = (row.get("station") or "").strip()
        try:
            t = datetime.strptime(row["valid"].strip(), "%Y-%m-%d %H:%M")
        except (KeyError, ValueError, AttributeError):
            continue
        lat, lon = _num(row.get("lat")), _num(row.get("lon"))
        if not stn or lat is None or lon is None or not in_box(lat, lon, box):
            continue
        if cycle is not None:
            dt = in_window(t, cycle, window_min)
            if dt is None:
                continue
        else:
            dt = 0.0
        if stn not in best or dt < best[stn][0]:
            best[stn] = (dt, t, lat, lon, row)

    out = []
    for stn, (_, t, lat, lon, row) in best.items():
        elev = _num(row.get("elevation"))
        tmpf = _num(row.get("tmpf"))
        t_k = f_to_k(tmpf) if tmpf is not None else None
        if t_k is not None:
            out.append(_obs(t, lat, lon, "TMP", t_k, "asos", stn, elev))
        rh = _num(row.get("relh"))
        if rh is not None:
            out.append(_obs(t, lat, lon, "RH", rh, "asos", stn, elev))
        drct, sknt = _num(row.get("drct")), _num(row.get("sknt"))
        if drct is not None and sknt is not None:
            u, v = wind_to_uv(drct, knots_to_ms(sknt))
            out.append(_obs(t, lat, lon, "UGRD", u, "asos", stn, elev))
            out.append(_obs(t, lat, lon, "VGRD", v, "asos", stn, elev))
        mslp, alti = _num(row.get("mslp")), _num(row.get("alti"))
        slp = None
        if mslp is not None:
            slp = mslp * 100.0
        elif alti is not None and elev is not None and t_k is not None:
            slp = slp_from_station_pressure(
                station_pressure_from_altimeter(alti, elev), elev, t_k)
        if slp is not None:
            out.append(_obs(t, lat, lon, "PMSL", slp, "asos", stn, elev))
    return out


def fetch_asos(cycle, box, fetcher):
    url = asos_url(cycle)
    text = fetcher.get_text(url, timeout=300, use_cache=False)
    obs = parse_asos(text, box, cycle)
    if not obs:
        return SourceResult("asos", "unavailable", raw={"asos": text},
                            message="no reports in the window")
    return SourceResult("asos", "ok", obs, {"asos": text},
                        platforms=len({o.station for o in obs}))


# ---------------------------------------------------------------------------
# Radiosondes
# ---------------------------------------------------------------------------

def parse_raob_table(geojson_text, box):
    """Active radiosonde sites in the box, from IEM's RAOB network table."""
    g = json.loads(geojson_text)
    sites = {}
    for f in g.get("features", []):
        p = f.get("properties", {})
        lon, lat = f["geometry"]["coordinates"][:2]
        sid = p.get("sid", "")
        # Composite "_XXX" area ids and closed sites are skipped: the one
        # returns nothing, the other returns nothing new.
        if sid.startswith("_") or p.get("archive_end"):
            continue
        if in_box(lat, lon, box):
            sites[sid] = {"lat": float(lat), "lon": float(lon),
                          "elev": _num(p.get("elevation"))}
    return sites


def raob_url(station, cycle):
    """One station per request -- the service rejects lists (P-54)."""
    params = [("sts", (cycle - timedelta(hours=1)).strftime("%Y-%m-%dT%H:%MZ")),
              ("ets", (cycle + timedelta(minutes=1)).strftime("%Y-%m-%dT%H:%MZ")),
              ("station", station)]
    return f"{IEM_RAOB}?{urlencode(params)}"


def parse_raob(text, sites, cycle=None):
    """
    IEM RAOB CSV -> Observations on pressure levels.

    Header measured 2026-09-22: station, validUTC, levelcode, pressure_mb,
    height_m, tmpc, dwpc, drct, speed_kts, bearing, range_sm.
    Levels with no pressure are dropped; everything else is kept, so the
    analysis sees the significant levels as well as the mandatory ones.
    """
    out = []
    for row in csv.DictReader(io.StringIO(text)):
        stn = (row.get("station") or "").strip()
        meta = sites.get(stn)
        if meta is None:
            continue
        raw = (row.get("validUTC") or "").strip()
        t = None
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"):
            try:
                t = datetime.strptime(raw, fmt)
                break
            except ValueError:
                pass
        # Soundings are stamped with the nominal time (12:00) and launched
        # about an hour before it, so the whole profile is at or before the
        # cycle time.
        if t is None or (cycle is not None and in_window(t, cycle) is None):
            continue
        p = _num(row.get("pressure_mb"))
        if p is None or p <= 0:
            continue
        pa = p * 100.0
        lat, lon, elev = meta["lat"], meta["lon"], meta.get("elev")

        def add(var, value):
            if value is not None:
                out.append(_obs(t, lat, lon, var, value, "raob", stn, elev, pa))

        tc, td = _num(row.get("tmpc")), _num(row.get("dwpc"))
        add("TMP", c_to_k(tc) if tc is not None else None)
        add("HGT", _num(row.get("height_m")))
        if tc is not None and td is not None:
            add("RH", rh_from_dewpoint(tc, td))
        d, s = _num(row.get("drct")), _num(row.get("speed_kts"))
        if d is not None and s is not None:
            u, v = wind_to_uv(d, knots_to_ms(s))
            add("UGRD", u)
            add("VGRD", v)
    return out


def fetch_raob(cycle, box, fetcher):
    if cycle.hour not in (0, 12):
        return SourceResult("raob", "unavailable",
                            message=f"no launches at {cycle:%H}Z")
    table = fetcher.get_text(IEM_RAOB_TABLE, timeout=60)
    sites = parse_raob_table(table, box)
    obs, raw, missing = [], {"raob_table": table}, []
    for sid in sorted(sites):
        text = fetcher.get_text(raob_url(sid, cycle), timeout=120,
                                use_cache=False)
        got = parse_raob(text, sites, cycle)
        if got:
            obs.extend(got)
            raw[f"raob_{sid}"] = text
        else:
            missing.append(sid)
    n = len({o.station for o in obs})
    msg = f"{n}/{len(sites)} sites" + (f"; none from {','.join(missing)}"
                                        if missing else "")
    if not obs:
        return SourceResult("raob", "unavailable", raw=raw, message=msg)
    return SourceResult("raob", "ok", obs, raw, msg, platforms=n)


# ---------------------------------------------------------------------------
# NDBC buoys and C-MAN
# ---------------------------------------------------------------------------

def parse_ndbc_active(xml_text, box):
    sites = {}
    for s in ET.fromstring(xml_text):
        a = s.attrib
        if a.get("met") != "y":
            continue
        lat, lon = _num(a.get("lat")), _num(a.get("lon"))
        if lat is None or lon is None or not in_box(lat, lon, box):
            continue
        sites[a["id"].upper()] = {"lat": lat, "lon": lon,
                                  "elev": _num(a.get("elev")) or 0.0,
                                  "type": a.get("type", "")}
    return sites


def parse_ndbc_5day(text, sid, meta, cycle, window_min=WINDOW_MIN):
    """
    NDBC standard-meteorology text -> the report nearest the cycle time.

    Columns: YY MM DD hh mm WDIR WSPD GST WVHT DPD APD MWD PRES ATMP WTMP DEWP
    ... with 'MM' for missing. PRES is sea-level pressure in hPa.
    """
    best = None
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 16:
            continue
        try:
            t = datetime(int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]))
        except ValueError:
            continue
        dt = in_window(t, cycle, window_min)
        if dt is not None and (best is None or dt < best[0]):
            best = (dt, t, f)
    if best is None:
        return []
    _, t, f = best
    val = lambda i: None if f[i] == "MM" else _num(f[i])
    lat, lon, elev = meta["lat"], meta["lon"], meta["elev"]
    out = []
    wdir, wspd = val(5), val(6)
    if wdir is not None and wspd is not None:
        u, v = wind_to_uv(wdir, wspd)
        out += [_obs(t, lat, lon, "UGRD", u, "ndbc", sid, elev),
                _obs(t, lat, lon, "VGRD", v, "ndbc", sid, elev)]
    pres, atmp, dewp = val(12), val(13), val(15)
    if pres is not None:
        out.append(_obs(t, lat, lon, "PMSL", pres * 100.0, "ndbc", sid, elev))
    if atmp is not None:
        out.append(_obs(t, lat, lon, "TMP", c_to_k(atmp), "ndbc", sid, elev))
        if dewp is not None:
            out.append(_obs(t, lat, lon, "RH", rh_from_dewpoint(atmp, dewp),
                            "ndbc", sid, elev))
    return out


def fetch_ndbc(cycle, box, fetcher, now=None):
    now = now or utcnow()
    if now - cycle > timedelta(days=4, hours=12):
        return SourceResult("ndbc", "unavailable",
                            message="5-day files do not reach back this far")
    active = fetcher.get_text(NDBC_ACTIVE, timeout=60)
    sites = parse_ndbc_active(active, box)
    obs, raw, failed = [], {"ndbc_active": active}, 0
    with single_try(fetcher):
        for sid in sorted(sites):
            try:
                text = fetcher.get_text(NDBC_5DAY.format(sid=sid), timeout=60,
                                        use_cache=False)
            except RuntimeError:
                failed += 1                   # many listed stations are down
                continue
            got = parse_ndbc_5day(text, sid, sites[sid], cycle)
            if got:
                obs.extend(got)
                raw[f"ndbc_{sid}"] = text
    n = len({o.station for o in obs})
    msg = f"{n}/{len(sites)} stations reported" + (
        f", {failed} unreachable" if failed else "")
    if not obs:
        return SourceResult("ndbc", "unavailable", raw=raw, message=msg)
    return SourceResult("ndbc", "ok", obs, raw, msg, platforms=n)


# ---------------------------------------------------------------------------
# Radar: archived, not assimilated
# ---------------------------------------------------------------------------

MRMS_PRODUCT = "CONUS/MergedReflectivityQCComposite_00.50"


def mrms_list_url(cycle):
    """S3 listing of that hour's files. File times carry seconds (every ~2 min),
    so the key cannot be built from the cycle time -- it has to be looked up."""
    prev = cycle - timedelta(hours=1)          # the file at or just before HH:00
    prefix = (f"{MRMS_PRODUCT}/{prev:%Y%m%d}/MRMS_MergedReflectivityQCComposite"
              f"_00.50_{prev:%Y%m%d-%H}")
    return f"{MRMS_BUCKET}/?list-type=2&prefix={prefix}"


def nearest_mrms_key(listing_xml, cycle):
    keys = [k.text for k in ET.fromstring(listing_xml).iter()
            if k.tag.endswith("Key") and k.text]
    best = None
    for k in keys:
        try:
            stamp = k.rsplit("_", 1)[-1].split(".")[0]      # YYYYmmdd-HHMMSS
            t = datetime.strptime(stamp, "%Y%m%d-%H%M%S")
        except ValueError:
            continue
        dt = (cycle - t).total_seconds()
        if 0 <= dt and (best is None or dt < best[0]):
            best = (dt, k)
    return None if best is None or best[0] > 600 else best[1]


def fetch_mrms(cycle, box, fetcher):
    """
    The composite reflectivity at the cycle time, kept verbatim for later.

    Not assimilated: the model carries no water, so there is no state
    variable reflectivity could correct. Archiving it now costs ~1 MB per run
    and means the history exists when moisture does.
    """
    with single_try(fetcher):
        key = nearest_mrms_key(
            fetcher.get_text(mrms_list_url(cycle), timeout=60, use_cache=False),
            cycle)
        if key is None:
            return SourceResult("mrms", "unavailable",
                                message=f"no file within 10 min of {cycle:%H%M}Z")
        data = fetcher.get(f"{MRMS_BUCKET}/{key}", timeout=120, use_cache=False)
    return SourceResult("mrms", "ok", raw={"mrms_refl_composite": data},
                        message=f"{len(data)/1e6:.1f} MB archived, "
                                f"not assimilated (dry model)")


def fetch_vad(cycle, box, fetcher):
    return SourceResult(
        "vad", "disabled",
        message="decoder not yet checked against real NEXRAD Level III "
                "files with the server's stack (P-53)")


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

SOURCES = [("asos", fetch_asos), ("raob", fetch_raob), ("ndbc", fetch_ndbc),
           ("vad", fetch_vad), ("mrms", fetch_mrms)]


def collect(cycle, domain=None, fetcher=None, sources=None, skip=(),
            verbose=True):
    """
    Run every adapter, never letting one source stop the others.

    `skip` removes sources by name -- used to measure what a source
    contributes (prediction P6), and to run offline.
    """
    import config
    from netpolicy import PoliteFetcher
    fetcher = fetcher or PoliteFetcher()
    box = analysis_box(domain or config.DOMAIN)
    results = []
    for name, fn in (sources or SOURCES):
        t0 = time.time()
        if name in skip:
            r = SourceResult(name, "disabled", message="skipped on request")
        else:
            try:
                r = fn(cycle, box, fetcher)
            except Exception as e:                       # noqa: BLE001
                r = SourceResult(name, "failed",
                                 message=f"{type(e).__name__}: {e}"[:300])
        r.seconds = time.time() - t0
        results.append(r)
        if verbose:
            print("  " + r.summary(), flush=True)
    return results


def reparse(raw, cycle, box):
    """
    Observations from archived raw payloads -- the same parsers, no network.

    The archive keeps what each service sent, verbatim, so an analysis can be
    rebuilt after a parser fix; this is the function that does the rebuilding.
    """
    obs = []
    if "asos" in raw:
        obs += parse_asos(raw["asos"], box, cycle)
    if "raob_table" in raw:
        sites = parse_raob_table(raw["raob_table"], box)
        for k, text in raw.items():
            if k.startswith("raob_") and k != "raob_table":
                obs += parse_raob(text, sites, cycle)
    if "ndbc_active" in raw:
        sites = parse_ndbc_active(raw["ndbc_active"], box)
        for k, text in raw.items():
            if k.startswith("ndbc_") and k != "ndbc_active":
                sid = k[5:]
                if sid in sites:
                    obs += parse_ndbc_5day(text, sid, sites[sid], cycle)
    return obs


def availability_table(results):
    return "\n".join(r.summary() for r in results)


def all_observations(results):
    return [o for r in results if r.status == "ok" for o in r.obs]
