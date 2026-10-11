# NWP1 status — 2026-11-07 20:10Z

Rebuilt every 10 min on the school server from NWP1's own logs. Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.

## Forecast skill against live observations

Scored from ASOS observations by `tools/scoreboard.py`, 2026-11-07 17:20Z. Each forecast is compared with two persistence forecasts built from the same observations: **P0** holds each station's value at the start time; **P24** uses the value 24 h before the valid time. Skill = 1 − MSE(model) / MSE(best of P0, P24): above 0 the model beats both.

### 12 km — 29 cycles to 11-06 12Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 57836 | 3.13 | -0.63 | 4.43 | 4.68 | 0.50 |
| 2 m temperature (K) | 7–12 h | 56544 | 4.33 | -0.65 | 6.95 | 4.69 | 0.15 |
| 2 m temperature (K) | 13–18 h | 55256 | 4.64 | -0.54 | 6.87 | 4.68 | 0.02 |
| 2 m temperature (K) | 19–24 h | 55252 | 4.60 | -0.70 | 5.24 | 4.65 | 0.02 |
| 2 m dewpoint (K) | 1–6 h | 57543 | 2.86 | -0.93 | 2.02 | 4.97 | -1.01 |
| 2 m dewpoint (K) | 7–12 h | 56181 | 3.46 | -0.62 | 3.23 | 4.97 | -0.15 |
| 2 m dewpoint (K) | 13–18 h | 54940 | 3.92 | -0.47 | 3.87 | 4.87 | -0.02 |
| 2 m dewpoint (K) | 19–24 h | 54957 | 4.30 | -0.65 | 4.52 | 4.79 | 0.09 |
| 10 m wind (vector) (m/s) | 1–6 h | 50252 | 2.78 | -0.19 | 2.28 | 3.48 | -0.49 |
| 10 m wind (vector) (m/s) | 7–12 h | 49141 | 3.07 | -0.12 | 2.99 | 3.41 | -0.06 |
| 10 m wind (vector) (m/s) | 13–18 h | 48045 | 3.23 | -0.02 | 3.19 | 3.28 | -0.02 |
| 10 m wind (vector) (m/s) | 19–24 h | 48590 | 3.39 | 0.05 | 3.23 | 3.27 | -0.10 |
| sea-level pressure (hPa) | 1–6 h | 39661 | 2.21 | -0.29 | 1.34 | 6.48 | -1.71 |
| sea-level pressure (hPa) | 7–12 h | 39029 | 3.31 | -0.21 | 3.01 | 6.40 | -0.21 |
| sea-level pressure (hPa) | 13–18 h | 38563 | 4.43 | 0.01 | 4.39 | 6.32 | -0.01 |
| sea-level pressure (hPa) | 19–24 h | 38725 | 5.43 | 0.41 | 5.58 | 6.08 | 0.05 |

2 m temperature by cycle, 1–24 h (12 km): model RMSE / P24 RMSE (K)

11-03 18Z 3.3/3.2 **worse**; 11-04 00Z 3.3/3.3 better; 11-04 06Z 4.9/3.5 **worse**; 11-04 12Z 6.0/4.3 **worse**; 11-04 18Z 3.1/4.4 better; 11-05 00Z 2.9/4.2 better; 11-05 06Z 3.3/4.4 better; 11-05 12Z 3.1/4.0 better; 11-05 18Z 3.4/4.1 better; 11-06 00Z 3.3/4.6 better; 11-06 06Z 3.1/4.9 better; 11-06 12Z 2.1/5.2 better

### 3 km — 29 cycles to 11-06 12Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 59384 | 3.12 | -0.61 | 4.40 | 4.68 | 0.50 |
| 2 m temperature (K) | 7–12 h | 58063 | 4.25 | -0.42 | 6.91 | 4.69 | 0.18 |
| 2 m temperature (K) | 13–18 h | 56747 | 4.57 | -0.17 | 6.83 | 4.68 | 0.05 |
| 2 m temperature (K) | 19–24 h | 56738 | 4.63 | -0.24 | 5.23 | 4.64 | 0.01 |
| 2 m dewpoint (K) | 1–6 h | 59099 | 3.07 | -1.23 | 2.01 | 4.96 | -1.32 |
| 2 m dewpoint (K) | 7–12 h | 57714 | 3.83 | -0.88 | 3.23 | 4.97 | -0.40 |
| 2 m dewpoint (K) | 13–18 h | 56446 | 4.28 | -0.70 | 3.88 | 4.87 | -0.22 |
| 2 m dewpoint (K) | 19–24 h | 56455 | 4.65 | -0.90 | 4.53 | 4.80 | -0.05 |
| 10 m wind (vector) (m/s) | 1–6 h | 51681 | 2.90 | -0.08 | 2.27 | 3.49 | -0.62 |
| 10 m wind (vector) (m/s) | 7–12 h | 50539 | 3.24 | -0.05 | 2.99 | 3.42 | -0.18 |
| 10 m wind (vector) (m/s) | 13–18 h | 49420 | 3.38 | 0.02 | 3.19 | 3.29 | -0.13 |
| 10 m wind (vector) (m/s) | 19–24 h | 49978 | 3.55 | 0.08 | 3.23 | 3.27 | -0.20 |
| sea-level pressure (hPa) | 1–6 h | 40670 | 2.47 | 0.29 | 1.34 | 6.49 | -2.38 |
| sea-level pressure (hPa) | 7–12 h | 40014 | 3.60 | 0.25 | 3.02 | 6.40 | -0.42 |
| sea-level pressure (hPa) | 13–18 h | 39523 | 4.73 | 0.39 | 4.41 | 6.33 | -0.15 |
| sea-level pressure (hPa) | 19–24 h | 39677 | 5.73 | 0.76 | 5.59 | 6.09 | -0.05 |

2 m temperature by cycle, 1–24 h (3 km): model RMSE / P24 RMSE (K)

11-03 18Z 3.5/3.2 **worse**; 11-04 00Z 3.0/3.3 better; 11-04 06Z 4.7/3.5 **worse**; 11-04 12Z 5.9/4.3 **worse**; 11-04 18Z 3.3/4.3 better; 11-05 00Z 2.8/4.2 better; 11-05 06Z 3.1/4.3 better; 11-05 12Z 2.9/4.0 better; 11-05 18Z 3.6/4.1 better; 11-06 00Z 4.0/4.6 better; 11-06 06Z 3.8/4.9 better; 11-06 12Z 2.3/5.2 better

### 12 km dry run (overlap) — 29 cycles to 11-06 12Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 57836 | 3.04 | -0.32 | 4.43 | 4.68 | 0.53 |
| 2 m temperature (K) | 7–12 h | 56544 | 4.16 | -0.39 | 6.95 | 4.69 | 0.21 |
| 2 m temperature (K) | 13–18 h | 55256 | 4.47 | -0.30 | 6.87 | 4.68 | 0.09 |
| 2 m temperature (K) | 19–24 h | 55252 | 4.50 | -0.38 | 5.24 | 4.65 | 0.06 |
| 10 m wind (vector) (m/s) | 1–6 h | 50252 | 2.75 | -0.34 | 2.28 | 3.48 | -0.46 |
| 10 m wind (vector) (m/s) | 7–12 h | 49141 | 3.04 | -0.26 | 2.99 | 3.41 | -0.04 |
| 10 m wind (vector) (m/s) | 13–18 h | 48045 | 3.17 | -0.22 | 3.19 | 3.28 | 0.01 |
| 10 m wind (vector) (m/s) | 19–24 h | 48590 | 3.33 | -0.20 | 3.23 | 3.27 | -0.07 |
| sea-level pressure (hPa) | 1–6 h | 39661 | 2.23 | -0.35 | 1.34 | 6.48 | -1.76 |
| sea-level pressure (hPa) | 7–12 h | 39029 | 3.33 | -0.25 | 3.01 | 6.40 | -0.22 |
| sea-level pressure (hPa) | 13–18 h | 38563 | 4.45 | -0.07 | 4.39 | 6.32 | -0.02 |
| sea-level pressure (hPa) | 19–24 h | 38725 | 5.46 | 0.34 | 5.58 | 6.08 | 0.04 |

2 m temperature by cycle, 1–24 h (12 km dry run (overlap)): model RMSE / P24 RMSE (K)

11-03 18Z 3.1/3.2 better; 11-04 00Z 3.2/3.3 better; 11-04 06Z 4.6/3.5 **worse**; 11-04 12Z 5.7/4.3 **worse**; 11-04 18Z 3.2/4.4 better; 11-05 00Z 3.0/4.2 better; 11-05 06Z 3.2/4.4 better; 11-05 12Z 3.0/4.0 better; 11-05 18Z 3.6/4.1 better; 11-06 00Z 3.3/4.6 better; 11-06 06Z 3.2/4.9 better; 11-06 12Z 3.1/5.2 better

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
| 2026-11-06T18Z | 11-06 19:45 | cam_cycle.sh | 131.4 | 0 |
| 2026-11-06T18Z | 11-06 21:56 | maps_3km | 10.2 | 0 |
| 2026-11-07T00Z | 11-07 01:45 | daily.sh | 18.2 | 0 |
| 2026-11-07T00Z | 11-07 01:45 | cam_cycle.sh | 128.0 | 0 |
| 2026-11-07T00Z | 11-07 03:52 | maps_3km | 10.2 | 0 |
| 2026-11-07T06Z | 11-07 07:45 | daily.sh | 17.9 | 0 |
| 2026-11-07T06Z | 11-07 07:45 | cam_cycle.sh | 125.7 | 0 |
| 2026-11-07T06Z | 11-07 09:50 | maps_3km | 10.2 | 0 |
| 2026-11-07T12Z | 11-07 13:45 | daily.sh | 18.2 | 0 |
| 2026-11-07T12Z | 11-07 13:45 | cam_cycle.sh | 125.2 | 0 |
| 2026-11-07T12Z | 11-07 15:50 | maps_3km | 10.2 | 0 |
| 2026-11-07T18Z | 11-07 19:45 | daily.sh | 18.2 | 0 |

## Server

- load average 1.04 / 3.09 / 4.17 (104 logical CPUs)
- memory available 379 of 394 GB
- next cron cycle 01:45Z (in 335 min)

| running forecast | grid | latest progress |
|---|---|---|
| 20261107_18 | 3 km | (no step line yet) |

