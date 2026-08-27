"""
Adaptive bias correction: adjust forecasts using recent forecast-observation
error, per station and lead time.

WHY THIS IS NOT CIRCULAR
------------------------
The corrections are learned against OBSERVATIONS -- radiosondes, METARs --
never against another model. Training a correction toward HRRR would make
this an HRRR emulator with extra steps, which is the failure mode that
motivated building a physics core in the first place. Every error here is
(forecast - observation), so the target is the real atmosphere.

METHOD
------
A Kalman filter tracks a slowly varying bias for each (station, lead time,
variable). After each verification:

    innovation = (forecast - observation) - bias_estimate
    bias      <- bias + K * innovation
    K          = P / (P + R)

P is the uncertainty in our bias estimate, R the observation-plus-noise
variance. K adapts automatically: uncertain at first, so it learns fast, then
settles as evidence accumulates. This is why it works with days of data
instead of the years a trained network needs.

WHAT IT CANNOT DO
-----------------
It removes SYSTEMATIC error -- a station that always runs 1.5 K warm because
its grid cell is at the wrong elevation. It cannot fix random error, and it
cannot fix a missed weather event. Applied too aggressively it drags the
forecast toward recent conditions, which destroys genuine signal in exactly
the situations that matter most: regime changes.

Guardrails, all deliberate:
  * corrections are capped (max_correction)
  * a minimum sample count is required before any correction is applied
  * corrections apply to OUTPUT ONLY and are never fed back into the model
    state, which would violate the conservation properties of the core
"""

import numpy as np


class KalmanBiasCorrector:
    """
    Tracks one bias estimate per (station, lead_time, variable) key.

    Parameters
    ----------
    process_var  : how fast the bias is allowed to drift (units^2 per cycle).
                   Larger = adapts faster but noisier.
    obs_var      : observation + representativeness error variance (units^2).
    max_correction : hard cap on the applied correction (units).
    min_samples  : refuse to correct until this many verifications exist.
    """

    def __init__(self, process_var=0.01, obs_var=1.0, max_correction=5.0,
                 min_samples=5):
        self.process_var = float(process_var)
        self.obs_var = float(obs_var)
        self.max_correction = float(max_correction)
        self.min_samples = int(min_samples)

        self.bias = {}          # key -> current bias estimate
        self.P = {}             # key -> estimate variance
        self.n = {}             # key -> number of updates
        self.history = {}       # key -> list of raw errors (diagnostics)

    # --- core ---------------------------------------------------------------

    def update(self, key, forecast, observation):
        """
        Fold one verified forecast into the bias estimate.

        Returns the updated bias.
        """
        if not np.isfinite(forecast) or not np.isfinite(observation):
            return self.bias.get(key, 0.0)

        error = float(forecast) - float(observation)

        if key not in self.bias:
            self.bias[key] = 0.0
            self.P[key] = self.obs_var          # start uncertain
            self.n[key] = 0
            self.history[key] = []

        # Predict: the bias may have drifted since last time.
        P = self.P[key] + self.process_var

        # Update.
        K = P / (P + self.obs_var)
        self.bias[key] += K * (error - self.bias[key])
        self.P[key] = (1.0 - K) * P
        self.n[key] += 1
        self.history[key].append(error)

        return self.bias[key]

    def correction(self, key):
        """
        Correction to SUBTRACT from a forecast. Zero until enough evidence.
        """
        if self.n.get(key, 0) < self.min_samples:
            return 0.0
        return float(np.clip(self.bias[key],
                             -self.max_correction, self.max_correction))

    def apply(self, key, forecast):
        return forecast - self.correction(key)

    # --- diagnostics --------------------------------------------------------

    def gain(self, key):
        """Current Kalman gain -- how much a new observation moves the estimate."""
        if key not in self.P:
            return 1.0
        P = self.P[key] + self.process_var
        return P / (P + self.obs_var)

    def stats(self, key):
        h = self.history.get(key, [])
        if not h:
            return None
        h = np.asarray(h)
        return {
            "n": len(h),
            "bias_estimate": self.bias[key],
            "raw_mean_error": float(h.mean()),
            "raw_rmse": float(np.sqrt((h**2).mean())),
            "gain": self.gain(key),
            "applied_correction": self.correction(key),
        }

    def summary(self):
        rows = []
        for key in sorted(self.bias, key=str):
            s = self.stats(key)
            if s:
                rows.append((key, s))
        return rows

    def __repr__(self):
        return (f"KalmanBiasCorrector({len(self.bias)} keys, "
                f"obs_var={self.obs_var}, cap={self.max_correction})")


def verify(forecasts, observations):
    """
    Standard scores for a set of matched forecast/observation pairs.

    bias : mean error, the part bias correction can remove
    rmse : total error
    mae  : less sensitive to outliers than rmse
    """
    f = np.asarray(forecasts, dtype=float)
    o = np.asarray(observations, dtype=float)
    m = np.isfinite(f) & np.isfinite(o)
    f, o = f[m], o[m]
    if len(f) == 0:
        return {"n": 0}

    e = f - o
    return {
        "n": len(f),
        "bias": float(e.mean()),
        "rmse": float(np.sqrt((e**2).mean())),
        "mae": float(np.abs(e).mean()),
        "corr": float(np.corrcoef(f, o)[0, 1]) if len(f) > 1 else np.nan,
    }


def skill_score(corrected, raw, observations):
    """
    Fractional RMSE improvement of corrected over raw.

    Positive = the correction helped. NEGATIVE MEANS IT HURT, and that is a
    real possibility worth checking every time -- an over-aggressive corrector
    drags forecasts toward recent conditions and loses more in regime changes
    than it gains on quiet days.
    """
    vc = verify(corrected, observations)
    vr = verify(raw, observations)
    if not vc.get("n") or not vr.get("n") or vr["rmse"] == 0:
        return np.nan
    return 1.0 - vc["rmse"] / vr["rmse"]
