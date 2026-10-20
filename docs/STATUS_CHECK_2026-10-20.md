# NWP1 status check — what needs to be improved (2026-10-20)

Measured 2026-10-20, 12:50–17:10Z, from the server's live data. Asked for in
prompt 224. The progress report of 2026-10-18 describes the work done; this
document ranks what is still wrong.

## 1. Where the model stands

**Live verification**, scored with `tools/scoreboard.py`. This is run in
scratch, because the scoreboard is not deployed yet (§4). It uses 29 cycles
ending 10-19 00Z. Skill = 1 − MSE(model)/MSE(best of P0, P24); above 0 the
model beats both persistence forecasts.

| 12 km | 1–6 h | 7–12 h | 13–18 h | 19–24 h |
|---|---|---|---|---|
| 2 m T skill (bias, K) | +0.35 (+0.08) | +0.48 (+0.39) | +0.22 (+0.68) | −0.04 (+0.72) |
| 10 m wind skill (speed bias, m/s) | −0.82 (−0.72) | −0.05 (−0.49) | +0.16 (−0.46) | +0.20 (−0.45) |
| MSLP skill (bias, hPa) | −1.19 (−0.90) | −0.22 (−1.80) | −0.09 (−2.60) | −0.05 (−3.20) |

| 3 km | 1–6 h | 7–12 h | 13–18 h | 19–24 h |
|---|---|---|---|---|
| 2 m T skill (bias, K) | +0.36 (−0.03) | +0.48 (+0.53) | +0.23 (+0.92) | −0.00 (+0.23) |
| 2 m dewpoint skill (bias, K) | −0.95 (−1.01) | −0.02 (−0.74) | +0.18 (−0.18) | +0.03 (−0.98) |
| 10 m wind skill (speed bias, m/s) | −0.86 (−0.41) | −0.10 (−0.17) | +0.10 (−0.16) | |

What the model does well:
- 2 m temperature from 1 to 18 h: it beats both persistence forecasts at 12 km and 3 km.
- 10 m wind beyond 12 h.

What it does badly:
- Temperature skill is gone by 19–24 h.
- Sea-level pressure never beats persistence.
- In the first 6 h, 10 m wind and dewpoint are worse than holding the start-time observation.
- The three newest verified 12 km cycles (10-18 12Z to 10-19 00Z) lost to P24, which scored 3.6–4.3 K in that quiet weather.

**3 km runs.**
- 8 of the last 10 cycles stopped at the deadline, at 20.2–23.6 h. The last fully completed run was 10-19 18Z.
- The new start-up (`exact-free`, in cron since 10-19 00Z) did not change this: 6 of its 7 cycles hit the deadline.

## 2. The new finding: the 3 km time step is set by a wind maximum that no observation supports

The 3 km step is the advective limit, dt = 0.7·dx/max|u| (`nh3d.max_dt`). With the initial
max|u| printed in each cycle's log, it matches the logged step to within 0.1 s:
- 2100/115.5 = 18.2 s (10-20 12Z; logged 18.2 s);
- 2100/109.7 = 19.1 s (10-19 18Z; logged 19.2 s);
- 2100/123.6 = 17.0 s (10-18 18Z; logged 16.9 s).

The printed max|u| is that of the 12 km driver state. The 3 km step uses the
3 km state, which accounts for the 0.1 s differences.

The log line says "external wave ~290 m/s sets this", but that is a fixed
string in `forecast.py`. It is wrong for the 3 km core.

Where the maximum is, in the 12 km analysis at or above 400 hPa:

| cycle | analysed max (m/s) | where | strongest in-domain sounding (m/s) | points > 60 m/s | median (m/s) |
|---|---|---|---|---|---|
| 10-13 00Z | 108.9 | 47.45 N 80.33 W (north edge, Ontario gap) | | 8.9 % | 32.8 |
| 10-17 12Z | 124.5 | 37.49 N 66.95 W (SE corner, ocean) | 73.6 | 33.5 % | 54.1 |
| 10-19 12Z | 119.0 | 37.70 N 70.44 W (south edge, ocean) | 23.1 | 11.2 % | 20.0 |
| 10-20 12Z | 125.3 | 37.81 N 70.00 W (south edge, ocean) | 30.4 | 11.5 % | 24.8 |

The 3 km analysis has the same maximum at the same place, because it is
interpolated from the 12 km analysis.

**Interpretation.**
- The core sits on the southern boundary over the Atlantic, where no sounding reaches.
- The first guess there is the previous forecast. The boundaries are frozen (P-53), so the edge holds whatever the analysis put there for 24 h, and the next first guess carries it on.
- Nothing in the cycle corrects it. The gross-error check is not active in production.
- On 10-19 the sounding maximum was 23 m/s and the analysis maximum 119 m/s. A 119 m/s core is not plausible on a day whose strongest observed wind is 23 m/s, but there is no observation at that point to prove it wrong.

**Consequence.** If the spurious maximum were removed, the 3 km step would rise
to the 26 s acoustic cap, 1.43× faster. The run that completed on 10-19 18Z
took 123 min at 19.2 s. That gain is likely enough to clear the deadline every
cycle. The 12 km step would also gain about 18 % (12.0 to 14.2 s from its own limit at 125 against 60 m/s).

## 3. What needs to be improved, in order

1. **P-82: the data-void jet.**
   - It costs 3 km forecast hours every cycle, and it fills the frozen boundary with a wind that is probably spurious.
   - Proposed test J5: where the nearest upper-air observation is far away, relax the first guess toward the observation-only analysis. One arm uses a fixed radius; the other an error that grows with cycling.
   - Measures: analysed max against soundings, 3 km dt and hours reached, aloft error against persistence, surface scores unchanged.
   - Before J5 runs, correct the 3 km time-step log label to name the actual limiter.
2. **P-53: frozen lateral boundaries.**
   - This is the largest error after 18 h: temperature skill falls to ≤ 0 at 19–24 h, and the 16–17 Oct cold outbreak was missed by up to +12 K at 24 h.
   - Observation-only options to test: upstream advection of the boundary from interior tendencies, or open (radiation) boundaries for the inflow edges.
3. **Sea-level pressure.**
   - It never beats persistence.
   - Over the last 29 cycles the 12 km bias drifts from −0.9 hPa at 1–6 h to −3.2 hPa at 24 h. In the earlier 29-cycle window it was −0.2 to +0.3 hPa, so this drift is new.
   - Check the domain mass budget, i.e. surface pressure lost through the boundaries, before blaming P-74.
4. **The first six hours** (spin-up).
   - 10 m wind (skill −0.82 / −0.86) and 3 km dewpoint (−0.95, dry by 1.0 K) are worse than persistence at 1–6 h.
   - 2 m T is already about 2 K wrong at +1 h, against 1.1 K for P0.
   - Diagnose what the first model hour does to the analysed surface layer.
5. **Boundary layer and vertical grid (P-84).**
   - L1 showed that lowering the lowest level alone makes the surface forecast worse.
   - The next candidate, L2 (height-scaled mixing length with a low first level), waits until items 1–4 move, because those dominate the error.
6. **P-81: 12 km rain.** No new cases since the tropical-remnant event; it stays open.

## 4. Operations and record

- **The last push has not reached the server.** origin/main is still 00d8512. Until it does:
  - the scoreboard and the status-page score section are not live;
  - verification runs once a day at 03:30Z and the newest verified run is 41 h old;
  - the deployment waiter (deadline 10-22) is still waiting.
- **P-83 confirmation is pending.** Only one `exact-free` cycle (10-19 00Z) has been verified so far, so the confirmation needs about 8 more verified cycles.
- **The problem register needs a sweep.**
  - Six entries marked FIXED are filed under OPEN: P-57, P-58, P-76, P-77, P-78, P-80.
  - Nine OPEN entries appear to describe states that no longer exist, for example P-06 "the fetchers have never touched the network" and P-07 "the archive has no data". The other seven are P-55, P-59, P-64, P-65, P-68, P-69 and P-75; each needs its evidence checked before it is closed.
  - `problem.py check` reports 9 format issues.
- **Prompt log:** prompt 40 has no row, 24 is untagged and 143 carries a tag outside the scheme.
- **Token ledger:** no rows since 09-22, so the cost ratio in TOKEN_COST is a lower bound.
