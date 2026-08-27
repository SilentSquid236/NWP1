"""
Verification: compare model output with real observations.

Feeds two consumers -- data assimilation (correct the state using the
difference) and post-processing (learn the systematic part of the difference).
Both need the same observation operator, so it lives here once.
"""
from observations import (Observation, run_qc, range_check, gross_error_check,
                          buddy_check, default_error_std, RANGE_LIMITS)
from obs_operator import GridInterpolator, elevation_correct_temperature
from fetchers import (parse_asos_csv, parse_raob_csv, fetch_asos, fetch_raob,
                      asos_url, raob_url, mrms_url, wind_to_uv, f_to_k, c_to_k,
                      knots_to_ms, rh_from_dewpoint, NORTHEAST_RAOB,
                      MRMS_PRODUCTS)
from scoring import (scores, scores_by, skill_vs_reference,
                     match_forecast_to_obs, ForecastArchive, report)

__all__ = ["Observation", "run_qc", "range_check", "gross_error_check",
           "buddy_check", "default_error_std", "RANGE_LIMITS",
           "GridInterpolator", "elevation_correct_temperature",
           "scores", "scores_by", "skill_vs_reference",
           "match_forecast_to_obs", "ForecastArchive", "report",
           "parse_asos_csv", "parse_raob_csv", "fetch_asos", "fetch_raob",
           "asos_url", "raob_url", "mrms_url", "wind_to_uv", "f_to_k",
           "c_to_k", "knots_to_ms", "rh_from_dewpoint", "NORTHEAST_RAOB",
           "MRMS_PRODUCTS"]
