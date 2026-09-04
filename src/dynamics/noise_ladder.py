"""
Noise-amplitude threshold on flat ground with a properly balanced 41 m/s jet.

The question the earlier matrix could not answer: is the model's failure on
"analysis-like noise" a property of the model, or of the noise amplitude that
was chosen? A ladder locates the threshold instead of asserting one.
"""
import numpy as np
np.seterr(all="ignore")
import probe_failure as P

print("flat ground, 41 m/s balanced jet, mixing+drag on", flush=True)
for noise in (0.0, 0.15, 0.3, 0.6, 1.2):
    m = P.build(0.0, noise, dT=1.5, clip=None)
    dt = m.max_dt()
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    fin = np.isfinite(m.u).all()
    print("  noise=%.2f m/s -> %2d/12 h, max|u|=%.1f"
          % (noise, done, np.abs(m.u).max() if fin else float("nan")),
          flush=True)
