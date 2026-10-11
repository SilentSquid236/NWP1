"""
Physical constants for every part of NWP1: the one place they are defined.

Until 2026-10-07 three values of the dry-air gas constant were in use (287.058
in the hydrostatic core, 287.04 in the non-hydrostatic core, the Kessler scheme
and the nest, 287.05 in the analysis and the verification) and two of cp
(1004.6 and 1004.5); docs/EQUATIONS.md section 9 found them. No reason was
recorded for any of the differences. The values below are the hydrostatic
production core's, so that the production model's numbers do not move:
- RD = R* / M_d = 8.314462618 / 0.0289647 = 287.058 J/(kg K);
- CP = 1004.6 J/(kg K) (dry air near 0 C; 3.5 RD would be 1004.7);
- G0 = 9.80665 m/s^2 (standard gravity), P0 = 1000 hPa.
"""
RD = 287.058          # J/(kg K), dry air
CP = 1004.6           # J/(kg K), dry air at constant pressure
CV = CP - RD          # J/(kg K)
KAPPA = RD / CP
G0 = 9.80665          # m/s^2
P0 = 100_000.0        # Pa
