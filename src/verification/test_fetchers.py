"""
Tests for observation parsing and unit conversion.

Parsers are tested against captured sample text, not live network calls, so
these run anywhere. Unit conversion is where silent biases come from: a
Fahrenheit-to-Kelvin slip or a wind-direction sign error produces plausible
numbers that are consistently wrong, which is the hardest kind of bug to see.

Run:  python test_fetchers.py
"""

from datetime import datetime

import numpy as np

from fetchers import (parse_asos_csv, parse_raob_csv, f_to_k, c_to_k,
                      knots_to_ms, wind_to_uv, rh_from_dewpoint,
                      asos_url, raob_url, mrms_url, NORTHEAST_RAOB, _num)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


ASOS_SAMPLE = """station,valid,lat,lon,elevation,tmpf,dwpf,relh,drct,sknt
ALB,2026-08-01 12:00,42.7483,-73.8017,84.0,68.0,55.0,63.12,270.00,10.00
ALB,2026-08-01 13:00,42.7483,-73.8017,84.0,71.6,M,M,M,M
BUF,2026-08-01 12:00,42.9408,-78.7358,215.0,32.0,30.0,92.00,0.00,20.00
BAD,2026-08-01 12:00,M,M,M,70.0,50.0,50.0,180.0,5.0
"""

RAOB_SAMPLE = """station,validUTC,levelcode,pressure_mb,height_m,tmpc,dwpc,drct,speed_kts
KALB,2026-08-01 12:00:00,4,1000.00,111.0,20.0,10.0,180.00,15.00
KALB,2026-08-01 12:00:00,4,850.00,1500.0,10.0,5.0,225.00,25.00
KALB,2026-08-01 12:00:00,4,500.00,5600.0,-15.0,-25.0,270.00,50.00
KOKX,2026-08-01 12:00:00,4,1000.00,90.0,22.0,18.0,90.00,10.00
KXXX,2026-08-01 12:00:00,4,1000.00,50.0,20.0,10.0,180.00,15.00
"""


# ---------------------------------------------------------------------------
def test_unit_conversions():
    """Check against values that are exactly known."""
    checks = [
        ("32F = 273.15K", abs(f_to_k(32.0) - 273.15) < 1e-9),
        ("212F = 373.15K", abs(f_to_k(212.0) - 373.15) < 1e-9),
        ("-40F = -40C", abs(f_to_k(-40.0) - c_to_k(-40.0)) < 1e-9),
        ("0C = 273.15K", abs(c_to_k(0.0) - 273.15) < 1e-9),
        ("1kt = 0.514444 m/s", abs(knots_to_ms(1.0) - 0.514444) < 1e-6),
    ]
    ok = all(c[1] for c in checks)
    report("unit conversions exact at known points", ok,
           "; ".join(f"{n} {'ok' if v else 'FAIL'}" for n, v in checks))


# ---------------------------------------------------------------------------
def test_wind_direction_convention():
    """
    Meteorological direction is where wind comes FROM. A westerly wind
    (270 deg) blows toward the east, so u > 0 and v ~ 0. Getting this backwards
    reverses every wind observation -- and the model would look badly wrong in
    a way that no amount of physics debugging would fix.
    """
    u_w, v_w = wind_to_uv(270.0, 10.0)      # from the west
    u_n, v_n = wind_to_uv(0.0, 10.0)        # from the north
    u_s, v_s = wind_to_uv(180.0, 10.0)      # from the south

    ok = (u_w > 9.9 and abs(v_w) < 1e-9 and
          abs(u_n) < 1e-9 and v_n < -9.9 and
          abs(u_s) < 1e-9 and v_s > 9.9)
    report("wind direction convention correct (FROM, not TOWARD)", ok,
           f"270deg -> u={u_w:+.2f} v={v_w:+.2f}; "
           f"0deg -> u={u_n:+.2f} v={v_n:+.2f}; "
           f"180deg -> u={u_s:+.2f} v={v_s:+.2f}")


# ---------------------------------------------------------------------------
def test_asos_parsing():
    """Parse a sample, including missing values and an unusable row."""
    obs = parse_asos_csv(ASOS_SAMPLE)

    alb_t = [o for o in obs if o.station == "ALB" and o.variable == "TMP"]
    buf_t = [o for o in obs if o.station == "BUF" and o.variable == "TMP"]
    bad = [o for o in obs if o.station == "BAD"]
    winds = [o for o in obs if o.variable in ("UGRD", "VGRD")]

    ok = (len(alb_t) == 2 and abs(alb_t[0].value - f_to_k(68.0)) < 1e-9
          and abs(buf_t[0].value - 273.15) < 1e-9
          and len(bad) == 0 and len(winds) == 4)
    report("ASOS CSV parses, converts, and skips unusable rows", ok,
           f"{len(obs)} obs; ALB TMP {alb_t[0].value:.2f}K (68F); "
           f"BUF {buf_t[0].value:.2f}K (32F); rows without lat/lon dropped: "
           f"{len(bad) == 0}; wind components {len(winds)}")


# ---------------------------------------------------------------------------
def test_asos_partial_row_keeps_what_it_has():
    """
    A row with temperature but no wind should still yield a temperature
    observation. Discarding whole rows for one missing field throws away a
    large fraction of real surface data.
    """
    obs = parse_asos_csv(ASOS_SAMPLE)
    hour13 = [o for o in obs
              if o.station == "ALB" and o.time == datetime(2026, 8, 1, 13)]
    kinds = {o.variable for o in hour13}

    ok = kinds == {"TMP"}
    report("partial rows keep the fields they do have", ok,
           f"13Z row had tmpf only -> variables {sorted(kinds)}")


# ---------------------------------------------------------------------------
def test_raob_parsing_and_station_filter():
    """
    Radiosonde rows parse with pressure attached, and unknown stations
    (outside the domain) are dropped rather than given invented coordinates.
    """
    obs = parse_raob_csv(RAOB_SAMPLE)

    alb = [o for o in obs if o.station == "ALB"]
    unknown = [o for o in obs if o.station == "XXX"]
    t850 = [o for o in alb if o.variable == "TMP" and abs(o.pressure - 85000) < 1]
    heights = [o for o in obs if o.variable == "HGT"]

    ok = (len(unknown) == 0 and len(t850) == 1
          and abs(t850[0].value - 283.15) < 1e-9
          and t850[0].lat == NORTHEAST_RAOB["ALB"]["lat"]
          and len(heights) == 4)
    report("RAOB parses with pressure, drops unknown stations", ok,
           f"{len(obs)} obs; ALB 850hPa T = {t850[0].value:.2f}K (10C); "
           f"unknown station KXXX dropped: {len(unknown) == 0}")


# ---------------------------------------------------------------------------
def test_rh_from_dewpoint():
    """RH must be 100% when dewpoint equals temperature, and fall as it drops."""
    same = rh_from_dewpoint(20.0, 20.0)
    dry = rh_from_dewpoint(20.0, 5.0)
    very_dry = rh_from_dewpoint(20.0, -10.0)

    ok = (abs(same - 100.0) < 1e-6 and 30.0 < dry < 45.0
          and very_dry < dry and very_dry > 0)
    report("RH from dewpoint behaves correctly", ok,
           f"Td=T -> {same:.1f}%; T=20 Td=5 -> {dry:.1f}%; "
           f"T=20 Td=-10 -> {very_dry:.1f}%")


# ---------------------------------------------------------------------------
def test_missing_value_spellings():
    """Observation networks spell 'missing' many ways; all must become None."""
    spellings = ["M", "", "None", "NA", "-999", "null", "  M  ", "abc"]
    got = [_num(s) for s in spellings]
    real = [_num("12.5"), _num("-3"), _num("0")]

    ok = all(g is None for g in got) and real == [12.5, -3.0, 0.0]
    report("all missing-value spellings parse to None", ok,
           f"{spellings} -> all None: {all(g is None for g in got)}; "
           f"real values preserved: {real}")


# ---------------------------------------------------------------------------
def test_no_model_sources():
    """
    Guard against regression: every source this module produces must be a
    measurement. If a model-derived source is ever added here, verification
    stops being independent and this test should fail loudly.
    """
    obs = parse_asos_csv(ASOS_SAMPLE) + parse_raob_csv(RAOB_SAMPLE)
    sources = {o.source for o in obs}
    model_words = {"hrrr", "gfs", "rap", "nam", "mesoanalysis", "analysis",
                   "model", "forecast"}
    contaminated = {s for s in sources
                    if any(w in s.lower() for w in model_words)}

    ok = not contaminated and sources == {"asos", "raob"}
    report("no model-derived sources in the observation stream", ok,
           f"sources present: {sorted(sources)}; model-derived: "
           f"{sorted(contaminated) if contaminated else 'none'}")


# ---------------------------------------------------------------------------
def test_urls_well_formed():
    """URL builders should produce something plausible without network access."""
    s, e = datetime(2026, 8, 1), datetime(2026, 8, 2)
    a = asos_url(["NY_ASOS", "PA_ASOS"], s, e)
    r = raob_url(["OKX", "ALB"], s, e)
    m = mrms_url("reflectivity_composite", datetime(2026, 8, 1, 12, 30))

    ok = ("network=NY_ASOS" in a and "network=PA_ASOS" in a and "tz=UTC" in a
          and "station=OKX" in r and "noaa-mrms-pds" in m
          and "20260801-123000" in m)
    report("request URLs are well formed", ok,
           f"ASOS has both networks and UTC; RAOB has stations; "
           f"MRMS path: ...{m[-52:]}")


if __name__ == "__main__":
    print("\nObservation fetchers and parsing\n" + "=" * 62)
    for fn in (test_unit_conversions,
               test_wind_direction_convention,
               test_asos_parsing,
               test_asos_partial_row_keeps_what_it_has,
               test_raob_parsing_and_station_filter,
               test_rh_from_dewpoint,
               test_missing_value_spellings,
               test_no_model_sources,
               test_urls_well_formed):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)
