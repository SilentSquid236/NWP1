# NWP1 status — 2026-10-24 13:30Z

Rebuilt every 10 min on the school server from NWP1's own logs. Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.

## Attention

- cron 2026-10-23T12 cam_cycle.sh status 2
- cron 2026-10-23T18 cam_cycle.sh status 2
- cron 2026-10-24T00 cam_cycle.sh status 2
- cron 2026-10-24T06 cam_cycle.sh status 2

## Forecast skill against live observations

Scored from ASOS observations by `tools/scoreboard.py`, 2026-10-24 11:22Z. Each forecast is compared with two persistence forecasts built from the same observations: **P0** holds each station's value at the start time; **P24** uses the value 24 h before the valid time. Skill = 1 − MSE(model) / MSE(best of P0, P24): above 0 the model beats both.

### 12 km — 29 cycles to 10-23 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 56690 | 3.54 | 0.08 | 5.03 | 4.94 | 0.49 |
| 2 m temperature (K) | 7–12 h | 55408 | 5.19 | 0.23 | 8.46 | 4.89 | -0.13 |
| 2 m temperature (K) | 13–18 h | 54380 | 5.47 | 0.57 | 8.39 | 4.86 | -0.27 |
| 2 m temperature (K) | 19–24 h | 54162 | 5.15 | 0.80 | 5.87 | 4.92 | -0.09 |
| 2 m dewpoint (K) | 1–6 h | 12300 | 3.94 | -1.23 | 2.40 | 3.05 | -1.70 |
| 2 m dewpoint (K) | 7–12 h | 13966 | 4.52 | -1.04 | 3.57 | 3.29 | -0.89 |
| 2 m dewpoint (K) | 13–18 h | 15691 | 4.14 | -0.06 | 3.89 | 3.64 | -0.29 |
| 2 m dewpoint (K) | 19–24 h | 17288 | 4.20 | -0.16 | 3.72 | 3.65 | -0.33 |
| 10 m wind (vector) (m/s) | 1–6 h | 48647 | 3.18 | -0.54 | 2.45 | 4.59 | -0.68 |
| 10 m wind (vector) (m/s) | 7–12 h | 47305 | 3.45 | -0.44 | 3.44 | 4.49 | -0.01 |
| 10 m wind (vector) (m/s) | 13–18 h | 46503 | 3.69 | -0.40 | 3.93 | 4.32 | 0.12 |
| 10 m wind (vector) (m/s) | 19–24 h | 47093 | 3.93 | -0.25 | 4.06 | 4.20 | 0.07 |
| sea-level pressure (hPa) | 1–6 h | 39557 | 3.24 | -0.03 | 2.23 | 10.11 | -1.11 |
| sea-level pressure (hPa) | 7–12 h | 38970 | 5.37 | -0.84 | 5.04 | 10.04 | -0.14 |
| sea-level pressure (hPa) | 13–18 h | 38490 | 7.35 | -1.75 | 7.16 | 9.75 | -0.05 |
| sea-level pressure (hPa) | 19–24 h | 38541 | 9.01 | -2.64 | 8.80 | 9.47 | -0.05 |

2 m temperature by cycle, 1–24 h (12 km): model RMSE / P24 RMSE (K)

10-20 12Z 4.7/3.2 **worse**; 10-20 18Z 5.3/3.2 **worse**; 10-21 00Z 3.6/3.1 **worse**; 10-21 06Z 4.0/2.9 **worse**; 10-21 12Z 5.0/2.6 **worse**; 10-21 18Z 5.3/2.4 **worse**; 10-22 00Z 3.5/2.4 **worse**; 10-22 06Z 4.6/2.8 **worse**; 10-22 12Z 5.5/2.9 **worse**; 10-22 18Z 5.2/3.1 **worse**; 10-23 00Z 3.1/3.4 better; 10-23 06Z 3.8/4.1 better

### 3 km — 29 cycles to 10-23 12Z (newest verified run started 23 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 57830 | 3.60 | -0.22 | 5.22 | 4.87 | 0.45 |
| 2 m temperature (K) | 7–12 h | 56504 | 5.25 | 0.14 | 8.57 | 4.85 | -0.17 |
| 2 m temperature (K) | 13–18 h | 55305 | 5.33 | 0.74 | 8.32 | 4.90 | -0.18 |
| 2 m temperature (K) | 19–24 h | 32985 | 4.37 | 0.21 | 5.95 | 4.19 | -0.08 |
| 2 m dewpoint (K) | 1–6 h | 57288 | 3.91 | -1.35 | 2.43 | 6.42 | -1.59 |
| 2 m dewpoint (K) | 7–12 h | 55931 | 4.96 | -0.62 | 4.02 | 6.28 | -0.52 |
| 2 m dewpoint (K) | 13–18 h | 54791 | 5.33 | 0.29 | 5.00 | 6.16 | -0.14 |
| 2 m dewpoint (K) | 19–24 h | 32741 | 5.29 | -0.14 | 4.64 | 4.92 | -0.30 |
| 10 m wind (vector) (m/s) | 1–6 h | 49604 | 3.24 | -0.28 | 2.45 | 4.50 | -0.75 |
| 10 m wind (vector) (m/s) | 7–12 h | 48394 | 3.43 | -0.17 | 3.35 | 4.31 | -0.05 |
| 10 m wind (vector) (m/s) | 13–18 h | 47580 | 3.62 | -0.14 | 3.79 | 4.20 | 0.09 |
| 10 m wind (vector) (m/s) | 19–24 h | 28868 | 3.56 | -0.09 | 3.49 | 3.71 | -0.04 |
| sea-level pressure (hPa) | 1–6 h | 40265 | 3.25 | 0.65 | 2.12 | 10.04 | -1.35 |
| sea-level pressure (hPa) | 7–12 h | 39648 | 5.03 | -0.48 | 4.74 | 9.75 | -0.12 |
| sea-level pressure (hPa) | 13–18 h | 39125 | 6.85 | -1.50 | 6.82 | 9.47 | -0.01 |
| sea-level pressure (hPa) | 19–24 h | 23584 | 6.92 | -0.75 | 6.90 | 7.81 | -0.00 |

2 m temperature by cycle, 1–24 h (3 km): model RMSE / P24 RMSE (K)

10-20 18Z 6.1/3.2 **worse**; 10-21 00Z 3.4/3.1 **worse**; 10-21 06Z 3.8/2.9 **worse**; 10-21 12Z 4.7/2.6 **worse**; 10-21 18Z 5.8/2.5 **worse**; 10-22 00Z 3.4/2.3 **worse**; 10-22 06Z 4.3/2.7 **worse**; 10-22 12Z 5.0/2.9 **worse**; 10-22 18Z 5.7/3.1 **worse**; 10-23 00Z 3.0/3.2 better; 10-23 06Z 3.4/3.9 better; 10-23 12Z 4.2/4.4 better

### 12 km dry run (overlap) — 16 cycles to 10-23 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 26059 | 3.57 | 0.21 | 5.52 | 3.03 | -0.39 |
| 2 m temperature (K) | 7–12 h | 27408 | 5.05 | 0.26 | 8.95 | 3.01 | -1.81 |
| 2 m temperature (K) | 13–18 h | 28906 | 4.97 | 0.60 | 8.81 | 3.09 | -1.59 |
| 2 m temperature (K) | 19–24 h | 30397 | 3.69 | 0.62 | 5.06 | 3.35 | -0.21 |
| 10 m wind (vector) (m/s) | 1–6 h | 22834 | 2.69 | -0.53 | 2.23 | 3.07 | -0.46 |
| 10 m wind (vector) (m/s) | 7–12 h | 23809 | 2.93 | -0.42 | 2.94 | 3.13 | 0.01 |
| 10 m wind (vector) (m/s) | 13–18 h | 24891 | 3.22 | -0.44 | 3.32 | 3.21 | -0.01 |
| 10 m wind (vector) (m/s) | 19–24 h | 26643 | 3.37 | -0.35 | 3.19 | 3.19 | -0.12 |
| sea-level pressure (hPa) | 1–6 h | 17817 | 2.32 | 0.58 | 1.43 | 5.21 | -1.62 |
| sea-level pressure (hPa) | 7–12 h | 18927 | 3.61 | 0.21 | 3.31 | 5.77 | -0.19 |
| sea-level pressure (hPa) | 13–18 h | 20062 | 4.90 | -0.03 | 4.71 | 6.22 | -0.08 |
| sea-level pressure (hPa) | 19–24 h | 21356 | 6.11 | -0.09 | 6.10 | 6.45 | -0.00 |

2 m temperature by cycle, 1–24 h (12 km dry run (overlap)): model RMSE / P24 RMSE (K)

10-20 12Z 4.7/3.2 **worse**; 10-20 18Z 5.3/3.2 **worse**; 10-21 00Z 3.8/3.1 **worse**; 10-21 06Z 3.8/2.9 **worse**; 10-21 12Z 4.5/2.6 **worse**; 10-21 18Z 5.0/2.4 **worse**; 10-22 00Z 3.6/2.4 **worse**; 10-22 06Z 4.3/2.8 **worse**; 10-22 12Z 4.8/2.9 **worse**; 10-22 18Z 4.7/3.1 **worse**; 10-23 00Z 3.2/3.4 better; 10-23 06Z 3.5/4.1 better

## Notes from the AI collaborator (newest first)

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
- **2026-10-13 19:25Z — simple ice added (opt-in --ice), all tests pass** (C 32/32, moisture 12/12, forecast 17/17). Below 0 C cloud water becomes ice and rain becomes snow (slower fall, melting at the freezing level). Queued: R1b (RH-space edge fade, ~23Z) then R2 (ice, and ice + fade, ~00:30Z), both after the 19:45Z cron cycle.

## Experiments (active in the last 24 h)

_None active._

## Operational cron cycles (last 3 days)

| cycle | started (UTC) | script | wall (min) | status |
|---|---|---|---|---|
| 2026-10-23T12Z | 10-23 13:45 | daily.sh | 18.9 | 0 |
| 2026-10-23T12Z | 10-23 13:45 | cam_cycle.sh | 161.5 | 2 ⚠ |
| 2026-10-23T12Z | 10-23 16:26 | maps_3km | 9.2 | 0 |
| 2026-10-23T18Z | 10-23 19:45 | daily.sh | 18.5 | 0 |
| 2026-10-23T18Z | 10-23 19:45 | cam_cycle.sh | 159.4 | 2 ⚠ |
| 2026-10-23T18Z | 10-23 22:24 | maps_3km | 9.0 | 0 |
| 2026-10-24T00Z | 10-24 01:45 | daily.sh | 18.8 | 0 |
| 2026-10-24T00Z | 10-24 01:45 | cam_cycle.sh | 163.7 | 2 ⚠ |
| 2026-10-24T00Z | 10-24 04:28 | maps_3km | 8.8 | 0 |
| 2026-10-24T06Z | 10-24 07:45 | daily.sh | 17.9 | 0 |
| 2026-10-24T06Z | 10-24 07:45 | cam_cycle.sh | 157.8 | 2 ⚠ |
| 2026-10-24T06Z | 10-24 10:22 | maps_3km | 9.8 | 0 |

## Server

- load average 0.00 / 0.00 / 0.00 (104 logical CPUs)
- memory available 386 of 394 GB
- next cron cycle 13:45Z (in 15 min)
- no forecast running

