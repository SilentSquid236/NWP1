# NWP1 status — 2026-11-28 17:20Z

Rebuilt every 10 min on the school server from NWP1's own logs. Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.

## Forecast skill against live observations

Scored from ASOS observations by `tools/scoreboard.py`, 2026-11-28 11:19Z. Each forecast is compared with two persistence forecasts built from the same observations: **P0** holds each station's value at the start time; **P24** uses the value 24 h before the valid time. Skill = 1 − MSE(model) / MSE(best of P0, P24): above 0 the model beats both.

### 12 km — 27 cycles to 11-27 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 52667 | 2.68 | -0.05 | 3.84 | 5.25 | 0.51 |
| 2 m temperature (K) | 7–12 h | 51227 | 3.33 | -0.31 | 6.40 | 5.25 | 0.60 |
| 2 m temperature (K) | 13–18 h | 50054 | 3.20 | -0.26 | 6.73 | 5.20 | 0.62 |
| 2 m temperature (K) | 19–24 h | 50176 | 2.75 | -0.17 | 5.66 | 5.26 | 0.73 |
| 2 m dewpoint (K) | 1–6 h | 52187 | 2.42 | -0.43 | 2.26 | 6.93 | -0.15 |
| 2 m dewpoint (K) | 7–12 h | 50759 | 2.74 | -0.30 | 4.16 | 7.05 | 0.56 |
| 2 m dewpoint (K) | 13–18 h | 49643 | 2.79 | 0.01 | 5.54 | 7.12 | 0.75 |
| 2 m dewpoint (K) | 19–24 h | 49770 | 2.84 | 0.07 | 6.72 | 7.16 | 0.82 |
| 10 m wind (vector) (m/s) | 1–6 h | 46542 | 2.79 | -0.04 | 2.44 | 4.30 | -0.31 |
| 10 m wind (vector) (m/s) | 7–12 h | 44968 | 2.67 | 0.06 | 3.71 | 4.79 | 0.48 |
| 10 m wind (vector) (m/s) | 13–18 h | 43939 | 2.68 | 0.13 | 4.63 | 5.20 | 0.67 |
| 10 m wind (vector) (m/s) | 19–24 h | 44511 | 2.73 | 0.20 | 5.26 | 5.43 | 0.73 |
| sea-level pressure (hPa) | 1–6 h | 36962 | 2.32 | -0.04 | 3.02 | 13.87 | 0.41 |
| sea-level pressure (hPa) | 7–12 h | 36366 | 1.21 | -0.08 | 7.20 | 14.93 | 0.97 |
| sea-level pressure (hPa) | 13–18 h | 35936 | 1.20 | -0.14 | 11.12 | 15.63 | 0.99 |
| sea-level pressure (hPa) | 19–24 h | 36198 | 1.37 | -0.37 | 14.66 | 15.91 | 0.99 |

2 m temperature by cycle, 1–24 h (12 km): model RMSE / P24 RMSE (K)

11-24 12Z 2.2/8.0 better; 11-24 18Z 3.2/6.8 better; 11-25 00Z 2.1/5.3 better; 11-25 06Z 2.6/4.2 better; 11-25 12Z 3.2/4.0 better; 11-25 18Z 2.6/4.1 better; 11-26 00Z 2.2/4.2 better; 11-26 06Z 2.3/4.5 better; 11-26 12Z 2.7/4.7 better; 11-26 18Z 2.3/5.0 better; 11-27 00Z 2.1/5.2 better; 11-27 06Z 1.7/4.8 better

### 3 km — 27 cycles to 11-27 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 53431 | 2.80 | -0.37 | 3.77 | 5.27 | 0.45 |
| 2 m temperature (K) | 7–12 h | 52111 | 3.74 | -0.74 | 6.34 | 5.28 | 0.50 |
| 2 m temperature (K) | 13–18 h | 50981 | 3.60 | -0.64 | 6.70 | 5.23 | 0.53 |
| 2 m temperature (K) | 19–24 h | 44195 | 2.91 | -0.49 | 5.77 | 5.48 | 0.72 |
| 2 m dewpoint (K) | 1–6 h | 53004 | 2.44 | -0.58 | 2.25 | 6.96 | -0.17 |
| 2 m dewpoint (K) | 7–12 h | 51691 | 2.87 | -0.30 | 4.16 | 7.07 | 0.52 |
| 2 m dewpoint (K) | 13–18 h | 50611 | 2.89 | 0.17 | 5.56 | 7.15 | 0.73 |
| 2 m dewpoint (K) | 19–24 h | 43853 | 2.84 | -0.06 | 7.14 | 7.63 | 0.84 |
| 10 m wind (vector) (m/s) | 1–6 h | 47367 | 2.79 | -0.06 | 2.45 | 4.33 | -0.30 |
| 10 m wind (vector) (m/s) | 7–12 h | 45878 | 2.65 | 0.04 | 3.72 | 4.81 | 0.49 |
| 10 m wind (vector) (m/s) | 13–18 h | 44883 | 2.60 | 0.07 | 4.63 | 5.20 | 0.69 |
| 10 m wind (vector) (m/s) | 19–24 h | 39135 | 2.76 | 0.13 | 5.56 | 5.73 | 0.75 |
| sea-level pressure (hPa) | 1–6 h | 37706 | 2.40 | -0.10 | 3.02 | 13.88 | 0.37 |
| sea-level pressure (hPa) | 7–12 h | 37154 | 1.33 | -0.16 | 7.20 | 14.93 | 0.97 |
| sea-level pressure (hPa) | 13–18 h | 36624 | 1.34 | -0.22 | 11.14 | 15.65 | 0.99 |
| sea-level pressure (hPa) | 19–24 h | 31969 | 1.61 | -0.49 | 15.71 | 17.06 | 0.99 |

2 m temperature by cycle, 1–24 h (3 km): model RMSE / P24 RMSE (K)

11-24 12Z 2.2/8.0 better; 11-24 18Z 3.3/6.8 better; 11-25 00Z 2.1/5.3 better; 11-25 06Z 2.9/4.2 better; 11-25 12Z 3.8/4.0 better; 11-25 18Z 2.8/4.1 better; 11-26 00Z 2.4/4.2 better; 11-26 06Z 3.0/4.5 better; 11-26 12Z 3.4/4.7 better; 11-26 18Z 2.5/5.0 better; 11-27 00Z 2.4/5.2 better; 11-27 06Z 1.8/4.8 better

### 12 km dry run (overlap) — 27 cycles to 11-27 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 52667 | 2.90 | -0.11 | 3.84 | 5.25 | 0.43 |
| 2 m temperature (K) | 7–12 h | 51227 | 4.16 | -0.38 | 6.40 | 5.25 | 0.37 |
| 2 m temperature (K) | 13–18 h | 50054 | 4.82 | -0.21 | 6.73 | 5.20 | 0.14 |
| 2 m temperature (K) | 19–24 h | 50176 | 5.06 | 0.00 | 5.66 | 5.26 | 0.07 |
| 10 m wind (vector) (m/s) | 1–6 h | 46542 | 3.21 | -0.77 | 2.44 | 4.30 | -0.73 |
| 10 m wind (vector) (m/s) | 7–12 h | 44968 | 3.92 | -0.70 | 3.71 | 4.79 | -0.12 |
| 10 m wind (vector) (m/s) | 13–18 h | 43939 | 4.54 | -0.83 | 4.63 | 5.20 | 0.04 |
| 10 m wind (vector) (m/s) | 19–24 h | 44511 | 5.01 | -0.91 | 5.26 | 5.43 | 0.09 |
| sea-level pressure (hPa) | 1–6 h | 36962 | 4.06 | 0.84 | 3.02 | 13.87 | -0.81 |
| sea-level pressure (hPa) | 7–12 h | 36366 | 7.74 | 1.63 | 7.20 | 14.93 | -0.16 |
| sea-level pressure (hPa) | 13–18 h | 35936 | 11.44 | 2.49 | 11.12 | 15.63 | -0.06 |
| sea-level pressure (hPa) | 19–24 h | 36198 | 14.82 | 3.36 | 14.66 | 15.91 | -0.02 |

2 m temperature by cycle, 1–24 h (12 km dry run (overlap)): model RMSE / P24 RMSE (K)

11-24 12Z 3.3/8.0 better; 11-24 18Z 3.6/6.8 better; 11-25 00Z 2.7/5.3 better; 11-25 06Z 4.1/4.2 better; 11-25 12Z 5.5/4.0 **worse**; 11-25 18Z 2.9/4.1 better; 11-26 00Z 2.7/4.2 better; 11-26 06Z 4.3/4.5 better; 11-26 12Z 5.2/4.7 **worse**; 11-26 18Z 3.8/5.0 better; 11-27 00Z 3.2/5.2 better; 11-27 06Z 2.7/4.8 better

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
| 2026-11-27T18Z | 11-27 19:45 | daily.sh | 17.2 | 0 |
| 2026-11-27T18Z | 11-27 19:45 | cam_cycle.sh | 129.1 | 0 |
| 2026-11-27T18Z | 11-27 21:54 | maps_3km | 10.1 | 0 |
| 2026-11-28T00Z | 11-28 01:45 | daily.sh | 17.6 | 0 |
| 2026-11-28T00Z | 11-28 01:45 | cam_cycle.sh | 130.4 | 0 |
| 2026-11-28T00Z | 11-28 03:55 | maps_3km | 10.0 | 0 |
| 2026-11-28T06Z | 11-28 07:45 | daily.sh | 17.4 | 0 |
| 2026-11-28T06Z | 11-28 07:45 | cam_cycle.sh | 131.6 | 0 |
| 2026-11-28T06Z | 11-28 09:56 | maps_3km | 9.6 | 0 |
| 2026-11-28T12Z | 11-28 13:45 | daily.sh | 17.2 | 0 |
| 2026-11-28T12Z | 11-28 13:45 | cam_cycle.sh | 133.1 | 0 |
| 2026-11-28T12Z | 11-28 15:58 | maps_3km | 9.8 | 0 |

## Server

- load average 1.02 / 0.92 / 0.60 (104 logical CPUs)
- memory available 386 of 394 GB
- next cron cycle 19:45Z (in 145 min)
- no forecast running

