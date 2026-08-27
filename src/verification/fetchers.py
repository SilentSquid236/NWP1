"""
Observation fetchers -- real measurements only.

SOURCES, and why these

  ASOS/AWOS + state mesonets  Iowa Environmental Mesonet (IEM) archive.
      Surface: temperature, dewpoint, wind, pressure. IEM aggregates the
      federal ASOS network AND most state mesonets behind one interface, so
      NY/PA/NE mesonets come through the same call. Plain CSV, no auth, and
      decades of history.

  Radiosondes (RAOB)          IEM's RAOB archive.
      The only true upper-air OBSERVATIONS. ~8 sites in the Northeast domain
      (OKX, ALB, BUF, GYX, CAR, PIT, IAD, WAL), launching 00Z and 12Z.
      Note CHH (Chatham MA) was decommissioned in 2021 and appears in older
      station lists.

  MRMS radar                  AWS noaa-mrms-pds.
      Already mosaicked to a national grid, so it drops onto a model grid
      without writing a radar compositor. NEXRAD Level II is per-site
      volumetric and enormous -- hundreds of GB per CONUS day -- and needs
      compositing before it is usable for verification.

WHAT IS DELIBERATELY NOT HERE

  SPC mesoanalysis is a gridded MODEL ANALYSIS, not observations. Verifying
  against it would be circular in the same way as verifying against HRRR:
  it would measure agreement with someone else's model, not with the
  atmosphere. It is useful for context, never as truth.

  Any model output as a source of "observations". Every value produced by
  this module is a measurement made by an instrument.

Parsing is separated from fetching so the parsers can be unit-tested offline
against captured samples -- which is how the tests here work, since network
access is not guaranteed in every environment.
"""

import io
import csv
from datetime import datetime, timedelta
from urllib.parse import urlencode

import numpy as np

from observations import Observation, default_error_std

IEM_ASOS = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"
IEM_RAOB = "https://mesonet.agron.iastate.edu/cgi-bin/request/raob.py"

MISSING = {"M", "", "None", "null", "NA", "-99", "-999", "T"}


# ---------------------------------------------------------------------------
# Unit conversion. Observation networks report in mixed units; everything
# downstream assumes SI, and a missed conversion is a silent bias.
# ---------------------------------------------------------------------------

def f_to_k(f):
    return (float(f) - 32.0) * 5.0 / 9.0 + 273.15


def c_to_k(c):
    return float(c) + 273.15


def knots_to_ms(kt):
    return float(kt) * 0.514444


def wind_to_uv(direction_deg, speed_ms):
    """
    Meteorological wind direction is the direction the wind comes FROM,
    measured clockwise from north. The conversion to u/v components carries
    a minus sign for exactly that reason -- a wind FROM the west (270 deg)
    blows TOWARD the east, so u is positive.
    """
    d = np.radians(float(direction_deg))
    u = -float(speed_ms) * np.sin(d)
    v = -float(speed_ms) * np.cos(d)
    return float(u), float(v)


def _num(x):
    """Parse a numeric field, returning None for the many spellings of missing."""
    if x is None:
        return None
    s = str(x).strip()
    if s in MISSING:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v if np.isfinite(v) else None


# ---------------------------------------------------------------------------
# ASOS / mesonet surface observations
# ---------------------------------------------------------------------------

def parse_asos_csv(text, station_meta=None, source="asos"):
    """
    Parse IEM ASOS CSV into Observation records.

    Expected columns: station, valid, tmpf, dwpf, relh, drct, sknt, elevation
    (elevation is supplied via station_meta when the CSV lacks it).

    Rows with missing values for a variable simply do not produce an
    Observation for that variable -- a partial row still yields the fields it
    does have, rather than being discarded wholesale.
    """
    station_meta = station_meta or {}
    out = []

    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        stn = (row.get("station") or "").strip()
        if not stn:
            continue
        try:
            t = datetime.strptime(row["valid"].strip(), "%Y-%m-%d %H:%M")
        except (KeyError, ValueError):
            continue

        meta = station_meta.get(stn, {})
        lat = _num(row.get("lat")) if row.get("lat") else meta.get("lat")
        lon = _num(row.get("lon")) if row.get("lon") else meta.get("lon")
        if lat is None or lon is None:
            continue
        elev = _num(row.get("elevation")) if row.get("elevation") else meta.get("elev")

        def add(var, value):
            if value is None:
                return
            out.append(Observation(
                time=t, lat=float(lat), lon=float(lon), variable=var,
                value=float(value), source=source, station=stn,
                elevation=elev, pressure=None,
                error_std=default_error_std(source, var)))

        tmpf = _num(row.get("tmpf"))
        add("TMP", f_to_k(tmpf) if tmpf is not None else None)
        add("RH", _num(row.get("relh")))

        drct, sknt = _num(row.get("drct")), _num(row.get("sknt"))
        if drct is not None and sknt is not None:
            u, v = wind_to_uv(drct, knots_to_ms(sknt))
            add("UGRD", u)
            add("VGRD", v)

    return out


def asos_url(networks, start, end, stations=None):
    """
    Build the IEM request URL.

    networks : e.g. ["NY_ASOS", "PA_ASOS", "MA_ASOS"] -- state mesonets use
               the same pattern, e.g. "NY_RWIS".
    """
    params = [("data", "tmpf"), ("data", "dwpf"), ("data", "relh"),
              ("data", "drct"), ("data", "sknt"),
              ("year1", start.year), ("month1", start.month), ("day1", start.day),
              ("year2", end.year), ("month2", end.month), ("day2", end.day),
              ("tz", "UTC"), ("format", "onlycomma"), ("latlon", "yes"),
              ("elev", "yes"), ("missing", "M"), ("trace", "T"),
              ("direct", "no"), ("report_type", "3")]
    for n in networks:
        params.append(("network", n))
    for s in (stations or []):
        params.append(("station", s))
    return f"{IEM_ASOS}?{urlencode(params)}"


def fetch_asos(networks, start, end, stations=None, timeout=120):
    """
    Download and parse surface observations.

    Network access required. Kept as a thin wrapper so the parser above can be
    tested without it.
    """
    import urllib.request
    url = asos_url(networks, start, end, stations)
    with urllib.request.urlopen(url, timeout=timeout) as r:
        text = r.read().decode("utf-8", errors="replace")
    return parse_asos_csv(text)


# ---------------------------------------------------------------------------
# Radiosondes
# ---------------------------------------------------------------------------

# Upper-air sites inside the Northeast domain (37-47.5N, 82-66W).
NORTHEAST_RAOB = {
    "OKX": {"name": "Upton NY", "lat": 40.87, "lon": -72.86, "elev": 20.0},
    "ALB": {"name": "Albany NY", "lat": 42.69, "lon": -73.83, "elev": 96.0},
    "BUF": {"name": "Buffalo NY", "lat": 42.94, "lon": -78.73, "elev": 218.0},
    "GYX": {"name": "Gray ME", "lat": 43.89, "lon": -70.25, "elev": 125.0},
    "CAR": {"name": "Caribou ME", "lat": 46.87, "lon": -68.01, "elev": 191.0},
    "PIT": {"name": "Pittsburgh PA", "lat": 40.53, "lon": -80.23, "elev": 360.0},
    "IAD": {"name": "Sterling VA", "lat": 38.98, "lon": -77.47, "elev": 85.0},
    "WAL": {"name": "Wallops Island VA", "lat": 37.85, "lon": -75.48, "elev": 13.0},
    # CHH (Chatham MA) decommissioned April 2021 -- deliberately absent.
}


def parse_raob_csv(text, station_meta=None, source="raob"):
    """
    Parse IEM RAOB CSV into Observation records.

    Expected columns: station, validUTC, levelcode, pressure_mb, tmpc, dwpc,
    drct, speed_kts, height_m

    Only mandatory and significant levels with a real pressure are used;
    surface rows are kept but flagged with pressure so the vertical operator
    handles them consistently.
    """
    meta_table = dict(NORTHEAST_RAOB)
    if station_meta:
        meta_table.update(station_meta)
    out = []

    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        stn = (row.get("station") or "").strip().upper().lstrip("K")
        meta = meta_table.get(stn)
        if meta is None:
            continue

        raw_time = (row.get("validUTC") or row.get("valid") or "").strip()
        t = None
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M"):
            try:
                t = datetime.strptime(raw_time, fmt)
                break
            except ValueError:
                continue
        if t is None:
            continue

        p_mb = _num(row.get("pressure_mb"))
        if p_mb is None or p_mb <= 0:
            continue
        p_pa = p_mb * 100.0

        def add(var, value):
            if value is None:
                return
            out.append(Observation(
                time=t, lat=meta["lat"], lon=meta["lon"], variable=var,
                value=float(value), source=source, station=stn,
                elevation=meta.get("elev"), pressure=p_pa,
                error_std=default_error_std(source, var)))

        tmpc = _num(row.get("tmpc"))
        add("TMP", c_to_k(tmpc) if tmpc is not None else None)
        add("HGT", _num(row.get("height_m")))

        # RH from temperature and dewpoint (Magnus formula).
        dwpc = _num(row.get("dwpc"))
        if tmpc is not None and dwpc is not None:
            add("RH", rh_from_dewpoint(tmpc, dwpc))

        drct, kts = _num(row.get("drct")), _num(row.get("speed_kts"))
        if drct is not None and kts is not None:
            u, v = wind_to_uv(drct, knots_to_ms(kts))
            add("UGRD", u)
            add("VGRD", v)

    return out


def rh_from_dewpoint(t_c, td_c):
    """Relative humidity (%) from temperature and dewpoint in Celsius (Magnus)."""
    a, b = 17.625, 243.04
    num = np.exp(a * td_c / (b + td_c))
    den = np.exp(a * t_c / (b + t_c))
    return float(np.clip(100.0 * num / den, 0.0, 100.0))


def raob_url(stations, start, end):
    params = [("ts1", start.strftime("%Y%m%d%H%M")),
              ("ts2", end.strftime("%Y%m%d%H%M")),
              ("format", "comma"), ("dl", "1")]
    for s in stations:
        params.append(("station", s))
    return f"{IEM_RAOB}?{urlencode(params)}"


def fetch_raob(stations, start, end, timeout=180):
    import urllib.request
    url = raob_url(stations, start, end)
    with urllib.request.urlopen(url, timeout=timeout) as r:
        text = r.read().decode("utf-8", errors="replace")
    return parse_raob_csv(text)


# ---------------------------------------------------------------------------
# MRMS radar
# ---------------------------------------------------------------------------

MRMS_BUCKET = "https://noaa-mrms-pds.s3.amazonaws.com"

MRMS_PRODUCTS = {
    "reflectivity_composite": "CONUS/MergedReflectivityQCComposite_00.50",
    "precip_rate": "CONUS/PrecipRate_00.00",
    "precip_1h": "CONUS/RadarOnly_QPE_01H_00.00",
    "echo_top_18dbz": "CONUS/EchoTop_18_00.50",
}


def mrms_url(product, when):
    """
    Build an MRMS object URL. Files appear every ~2 minutes.

    MRMS is a gridded product, so verification against it is grid-to-grid
    rather than grid-to-point: regrid MRMS to the model grid and compare
    fields. That is a different code path from the point-observation operator
    and is not yet implemented -- noted rather than half-built.
    """
    if product not in MRMS_PRODUCTS:
        raise ValueError(f"unknown MRMS product {product!r}; "
                         f"have {sorted(MRMS_PRODUCTS)}")
    path = MRMS_PRODUCTS[product]
    return (f"{MRMS_BUCKET}/{path}/{when:%Y%m%d}/"
            f"MRMS_{path.split('/')[-1]}_{when:%Y%m%d-%H%M%S}.grib2.gz")
