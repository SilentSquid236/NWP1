# NWP1 status — 2026-11-14 21:40Z

Rebuilt every 10 min on the school server from NWP1's own logs. Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.

## Attention

- cron 2026-11-14T00 cam_cycle.sh status 2
- cron 2026-11-14T06 cam_cycle.sh status 2
- cron 2026-11-14T12 cam_cycle.sh status 2

## Forecast skill against live observations

Scored from ASOS observations by `tools/scoreboard.py`, 2026-11-14 17:16Z. Each forecast is compared with two persistence forecasts built from the same observations: **P0** holds each station's value at the start time; **P24** uses the value 24 h before the valid time. Skill = 1 − MSE(model) / MSE(best of P0, P24): above 0 the model beats both.

### 12 km — 29 cycles to 11-13 12Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 57097 | 2.13 | -0.08 | 3.28 | 3.99 | 0.58 |
| 2 m temperature (K) | 7–12 h | 55701 | 2.47 | -0.18 | 5.22 | 4.02 | 0.62 |
| 2 m temperature (K) | 13–18 h | 54296 | 2.46 | -0.31 | 5.29 | 3.98 | 0.62 |
| 2 m temperature (K) | 19–24 h | 54470 | 2.27 | -0.41 | 4.22 | 3.89 | 0.66 |
| 2 m dewpoint (K) | 1–6 h | 57009 | 2.19 | -0.90 | 1.69 | 4.39 | -0.68 |
| 2 m dewpoint (K) | 7–12 h | 55627 | 2.18 | -0.54 | 2.80 | 4.43 | 0.39 |
| 2 m dewpoint (K) | 13–18 h | 54232 | 2.09 | -0.24 | 3.56 | 4.40 | 0.65 |
| 2 m dewpoint (K) | 19–24 h | 54412 | 2.08 | -0.16 | 4.10 | 4.34 | 0.74 |
| 10 m wind (vector) (m/s) | 1–6 h | 49863 | 3.06 | -0.12 | 2.52 | 4.26 | -0.47 |
| 10 m wind (vector) (m/s) | 7–12 h | 48455 | 2.70 | -0.05 | 3.45 | 4.25 | 0.39 |
| 10 m wind (vector) (m/s) | 13–18 h | 47200 | 2.64 | -0.04 | 3.95 | 4.23 | 0.56 |
| 10 m wind (vector) (m/s) | 19–24 h | 47716 | 2.67 | 0.00 | 4.15 | 4.21 | 0.58 |
| sea-level pressure (hPa) | 1–6 h | 39231 | 3.06 | -1.10 | 1.44 | 5.43 | -3.50 |
| sea-level pressure (hPa) | 7–12 h | 38587 | 1.70 | -1.02 | 3.20 | 5.52 | 0.72 |
| sea-level pressure (hPa) | 13–18 h | 37989 | 1.87 | -1.13 | 4.49 | 5.59 | 0.83 |
| sea-level pressure (hPa) | 19–24 h | 38238 | 2.24 | -1.31 | 5.42 | 5.64 | 0.83 |

2 m temperature by cycle, 1–24 h (12 km): model RMSE / P24 RMSE (K)

11-10 18Z 3.6/2.9 **worse**; 11-11 00Z 3.1/2.6 **worse**; 11-11 06Z 2.6/2.8 better; 11-11 12Z 3.5/3.1 **worse**; 11-11 18Z 2.8/3.4 better; 11-12 00Z 2.2/4.0 better; 11-12 06Z 1.9/4.1 better; 11-12 12Z 2.0/4.1 better; 11-12 18Z 2.0/4.4 better; 11-13 00Z 1.6/4.4 better; 11-13 06Z 1.9/4.5 better; 11-13 12Z 2.5/4.5 better

### 3 km — 29 cycles to 11-13 12Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 58521 | 2.30 | -0.25 | 3.29 | 4.00 | 0.51 |
| 2 m temperature (K) | 7–12 h | 57097 | 2.85 | -0.38 | 5.23 | 4.03 | 0.50 |
| 2 m temperature (K) | 13–18 h | 55810 | 2.80 | -0.48 | 5.30 | 3.98 | 0.50 |
| 2 m temperature (K) | 19–24 h | 41406 | 2.45 | -0.51 | 4.32 | 3.89 | 0.60 |
| 2 m dewpoint (K) | 1–6 h | 58437 | 2.34 | -0.93 | 1.69 | 4.38 | -0.91 |
| 2 m dewpoint (K) | 7–12 h | 57023 | 2.43 | -0.44 | 2.80 | 4.42 | 0.25 |
| 2 m dewpoint (K) | 13–18 h | 55747 | 2.27 | -0.00 | 3.55 | 4.39 | 0.59 |
| 2 m dewpoint (K) | 19–24 h | 41379 | 2.24 | 0.08 | 4.11 | 4.42 | 0.70 |
| 10 m wind (vector) (m/s) | 1–6 h | 51189 | 3.04 | -0.11 | 2.51 | 4.24 | -0.46 |
| 10 m wind (vector) (m/s) | 7–12 h | 49729 | 2.67 | -0.07 | 3.44 | 4.23 | 0.40 |
| 10 m wind (vector) (m/s) | 13–18 h | 48577 | 2.58 | -0.08 | 3.94 | 4.21 | 0.57 |
| 10 m wind (vector) (m/s) | 19–24 h | 36290 | 2.51 | -0.04 | 4.09 | 4.14 | 0.62 |
| sea-level pressure (hPa) | 1–6 h | 40228 | 3.21 | -1.31 | 1.44 | 5.43 | -3.95 |
| sea-level pressure (hPa) | 7–12 h | 39555 | 1.94 | -1.31 | 3.20 | 5.53 | 0.63 |
| sea-level pressure (hPa) | 13–18 h | 38987 | 2.08 | -1.40 | 4.50 | 5.61 | 0.79 |
| sea-level pressure (hPa) | 19–24 h | 29018 | 2.00 | -1.30 | 4.92 | 5.18 | 0.84 |

2 m temperature by cycle, 1–24 h (3 km): model RMSE / P24 RMSE (K)

11-10 18Z 4.1/3.0 **worse**; 11-11 00Z 3.2/2.6 **worse**; 11-11 06Z 3.1/2.8 **worse**; 11-11 12Z 4.4/2.7 **worse**; 11-11 18Z 3.2/3.3 better; 11-12 00Z 2.2/4.0 better; 11-12 06Z 2.1/4.2 better; 11-12 12Z 2.1/4.1 better; 11-12 18Z 2.0/4.2 better; 11-13 00Z 1.7/4.4 better; 11-13 06Z 2.0/4.5 better; 11-13 12Z 2.4/4.6 better

### 12 km dry run (overlap) — 29 cycles to 11-13 12Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 57097 | 2.57 | -0.05 | 3.28 | 3.99 | 0.39 |
| 2 m temperature (K) | 7–12 h | 55701 | 3.57 | 0.19 | 5.22 | 4.02 | 0.21 |
| 2 m temperature (K) | 13–18 h | 54296 | 3.96 | 0.48 | 5.29 | 3.98 | 0.01 |
| 2 m temperature (K) | 19–24 h | 54470 | 4.06 | 0.69 | 4.22 | 3.89 | -0.09 |
| 10 m wind (vector) (m/s) | 1–6 h | 49863 | 3.44 | -0.47 | 2.52 | 4.26 | -0.87 |
| 10 m wind (vector) (m/s) | 7–12 h | 48455 | 3.94 | -0.29 | 3.45 | 4.25 | -0.30 |
| 10 m wind (vector) (m/s) | 13–18 h | 47200 | 4.22 | -0.30 | 3.95 | 4.23 | -0.14 |
| 10 m wind (vector) (m/s) | 19–24 h | 47716 | 4.38 | -0.30 | 4.15 | 4.21 | -0.12 |
| sea-level pressure (hPa) | 1–6 h | 39231 | 3.05 | 0.64 | 1.44 | 5.43 | -3.47 |
| sea-level pressure (hPa) | 7–12 h | 38587 | 4.12 | 0.74 | 3.20 | 5.52 | -0.65 |
| sea-level pressure (hPa) | 13–18 h | 37989 | 5.05 | 0.89 | 4.49 | 5.59 | -0.26 |
| sea-level pressure (hPa) | 19–24 h | 38238 | 5.70 | 1.06 | 5.42 | 5.64 | -0.11 |

2 m temperature by cycle, 1–24 h (12 km dry run (overlap)): model RMSE / P24 RMSE (K)

11-10 18Z 4.4/2.9 **worse**; 11-11 00Z 4.8/2.6 **worse**; 11-11 06Z 4.3/2.8 **worse**; 11-11 12Z 4.2/3.1 **worse**; 11-11 18Z 3.7/3.4 **worse**; 11-12 00Z 3.7/4.0 better; 11-12 06Z 3.3/4.1 better; 11-12 12Z 3.1/4.1 better; 11-12 18Z 3.9/4.4 better; 11-13 00Z 3.0/4.4 better; 11-13 06Z 3.1/4.5 better; 11-13 12Z 3.6/4.5 better

## Notes from the AI collaborator (newest first)

- **2026-10-25T18:11:02Z — First-guess void blend in production.** The 12 km analysis now blends the previous forecast halfway toward the mean sounding far from any sounding, at 00Z and 12Z (test J5b: wind RMSE -3.2 %, temperature -1.6 %, false jet maximum 110 to 88 m/s; fixes the false low north of NY). Commit 3354e64. Crontab backup: data/logs/crontab_backup_20261025.txt
- 2026-10-22T01:17:04Z: deployed 98222db. 12 km product is now the moist NH run (forecast_nhm.npz); click-to-sounding viewer; moist map products; 3 km maps after each 3 km run (NWP_12KM_MOIST=1, NWP_CAM_MAPS=1). Crontab backup: data/logs/crontab_backup_20261021.txt
- 2026-10-20 17:25Z: commit cbe50b9 deployed by scratch/sb_sw/ff_sb.sh. tools/scoreboard.py is live and the status page has a scores section. Verify cron changed from 30 3 to 10 5,11,17,23 (backup data/logs/crontab_backup_20261019.txt). Tests 11/11. First scoreboard: 1180 rows, 0 alerts.
- **2026-10-20 05:15Z — test L1: a lower lowest level does not help on its own.** The 12 km model's lowest level is 235 m above ground (P-84). Rerunning six forecasts with it at 51, 38 and 42 m (20, 30 and 40 levels) made nights about 0.3 K less warm but days 0.25 K colder: 2 m temperature error rose 0.16-0.22 K overall and 10 m wind error 0.12-0.18 m/s. Winds aloft were unaffected. Production keeps its levels. The surface physics was tuned around the 235 m level, so a low level would need mixing that scales with height near the ground.
- **2026-10-18 23:20Z — the new start-up is in cron.** After your push (00d8512) the tests passed and both switches were set: the 12 km and the 3 km forecasts start with the new start-up from the 2026-10-19 00Z cycle (run at 01:45Z). The next cycles will show whether the gains in tests J4 and K1 hold in daily runs.
- **2026-10-18 00:55Z — test K1: the new start-up also works at 3 km.** One 3 km cycle (16 October 12Z) rerun with it: stable, starts closer to the analysis aloft (8.8 against 11.2 m/s), the fastest wind at the start 8 m/s lower so the time step is 7 % longer, and the surface scores and edge ring unchanged. Both start-up switches go into cron automatically after your next commit and push.
- **2026-10-17 17:40Z — test J4: a start-up fix that works.** A new start-up that lets the wind change on all four edges brings the 12 km forecast above 400 hPa to the level of persistence at +12 h (22.0 against 26.0 m/s for the current start-up), keeps the jet, and keeps the clean edges and pressure of the current one; 2 m temperature is 0.06 K better, 10 m wind mixed (u 0.19 better, v 0.12 worse). It is opt-in and not yet in cron.
- **2026-10-17 06:10Z — test J3: the start-up is the main source of the error aloft.** With the old (pre-E4) start-up the 12 km forecast comes within 3-6 % of persistence above 400 hPa (the current one is 18-24 % worse) and the 10 m wind improves by 0.24 m/s, but the pressure ring at the edges comes back (3.5 against 2.3) with a +2.4 hPa pressure bias. A variant that only removes the net inflow evenly changed little. Neither start-up is right yet; production is unchanged.
- **2026-10-15 18:30Z — correction: the sponge is NOT the cause.** Test J2 changed the 12 km upper boundary four ways (control, a thinner sponge, a sponge that follows the wind, a radiative lid); none helped aloft (all 1-8 % worse than the control), and the radiative lid ran clean on all six forecasts. The forecast's error aloft is made in its first hour: the start-up's divergence removal pushes each level's net inflow out through the north and east edges, adding about +15 m/s of southerly wind at 250-300 hPa (new problem P-83). The A3 analysis setting is in cron from the 19:45Z cycle.
- **2026-10-15 ~14:00Z — cause found, test J2 queued.** The 12 km model's wind sponge holds 206-299 hPa at the starting winds, so its forecasts cannot move the jet; above 400 hPa they are worse than persistence (24.2 against 21.1 m/s at +12 h). J2 compares four upper-boundary settings on six forecasts; it starts when the 13:45Z cycle finishes and should be done before 19:45Z.
- **2026-10-15 04:05Z — test J1.** Leaving the soundings outside the domain out of the analysis lowers the jet by about 20 m/s but does not stop it growing. You chose to adopt it anyway (with a 300 km length scale) as a mitigation: it goes into cron after your next push.
- **2026-10-15 — the analysed jet is too strong (P-82).** The 12 km analysis had 138 m/s at 200 hPa on 15 October 00Z; the strongest sounding was 69 m/s. Because the 3 km time step follows the fastest wind, it fell from 25.7 to 16.8 s, and every 3 km cycle since 14 October 12Z has hit its time limit, at 18.9-23.8 h. The 15 October 00Z cycle stopped at 18.9 h because my test J1 ran alongside it; tests now wait for gaps between cycles.
- **2026-10-15 04:27Z — simple ice is in cron.** On your commit 0e17989 the C tests passed 32/32 and `NWP_CAM_ICE=1` went into the crontab; the 06Z cycle was the first with ice. On the 13 October case (test R3) it raised hourly rain skill (ETS at 0.25 mm 0.077 -> 0.116) and cut the Bancroft maximum from 411 to 307 mm, still far above the 28 mm the nearest gauge measured (P-81 stays open).
- **2026-10-14 07:15Z — R2 (simple ice): September rain skill +52 % (ETS 0.023 -> 0.035), no cost to T/wind/dewpoint; the Bancroft extreme only -12 %.** Gauge check: CWRK near Bancroft measured 28 mm where the cron 3 km forecast 333 mm. Next suspect: no convection scheme at 12 km. Asking you whether to switch ice on in cron.
- **2026-10-14 01:30Z — R1b mixed, not adopted.** Fading the edge humidity toward the domain mean cut the Bancroft case (217 -> 93 mm at 12 km) but doubled another (10-12 18Z: 222 -> 605 mm) by moistening dry edges. Side finding: September dewpoint error improved 2.72 -> 2.52 K. R2 (simple ice) is running now.

## Experiments (active in the last 24 h)

_None active._

## Operational cron cycles (last 3 days)

| cycle | started (UTC) | script | wall (min) | status |
|---|---|---|---|---|
| 2026-11-13T18Z | 11-13 19:45 | cam_cycle.sh | 157.8 | 2 ⚠ |
| 2026-11-13T18Z | 11-13 22:22 | maps_3km | 9.2 | 0 |
| 2026-11-14T00Z | 11-14 01:45 | daily.sh | 20.8 | 0 |
| 2026-11-14T00Z | 11-14 01:45 | cam_cycle.sh | 159.8 | 2 ⚠ |
| 2026-11-14T00Z | 11-14 04:24 | maps_3km | 8.4 | 0 |
| 2026-11-14T06Z | 11-14 07:45 | daily.sh | 20.1 | 0 |
| 2026-11-14T06Z | 11-14 07:45 | cam_cycle.sh | 157.8 | 2 ⚠ |
| 2026-11-14T06Z | 11-14 10:22 | maps_3km | 8.8 | 0 |
| 2026-11-14T12Z | 11-14 13:45 | daily.sh | 21.2 | 0 |
| 2026-11-14T12Z | 11-14 13:45 | cam_cycle.sh | 158.9 | 2 ⚠ |
| 2026-11-14T12Z | 11-14 16:23 | maps_3km | 8.0 | 0 |
| 2026-11-14T18Z | 11-14 19:45 | daily.sh | 20.5 | 0 |

## Server

- load average 23.14 / 20.47 / 19.34 (104 logical CPUs)
- memory available 367 of 394 GB
- next cron cycle 01:45Z (in 245 min)

| running forecast | grid | latest progress |
|---|---|---|
| 20261114_18 | 3 km | step 2944/5424 t+13.03 h 0.6 steps/s elapsed 82.3 min ETA 69.3 min |

