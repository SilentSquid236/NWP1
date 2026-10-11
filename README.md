# NWP1 status — 2026-10-31 16:50Z

Rebuilt every 10 min on the school server from NWP1's own logs. Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.

## Forecast skill against live observations

Scored from ASOS observations by `tools/scoreboard.py`, 2026-10-31 11:20Z. Each forecast is compared with two persistence forecasts built from the same observations: **P0** holds each station's value at the start time; **P24** uses the value 24 h before the valid time. Skill = 1 − MSE(model) / MSE(best of P0, P24): above 0 the model beats both.

### 12 km — 29 cycles to 10-30 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 60308 | 3.76 | -0.45 | 5.36 | 4.50 | 0.30 |
| 2 m temperature (K) | 7–12 h | 59018 | 5.52 | -0.66 | 8.84 | 4.54 | -0.48 |
| 2 m temperature (K) | 13–18 h | 57708 | 5.68 | -0.51 | 8.69 | 4.65 | -0.49 |
| 2 m temperature (K) | 19–24 h | 56438 | 5.01 | -0.48 | 5.99 | 4.78 | -0.10 |
| 2 m dewpoint (K) | 1–6 h | 59865 | 3.39 | -1.04 | 2.36 | 4.48 | -1.06 |
| 2 m dewpoint (K) | 7–12 h | 58618 | 4.29 | -0.94 | 3.67 | 4.53 | -0.37 |
| 2 m dewpoint (K) | 13–18 h | 57324 | 4.95 | -0.60 | 4.29 | 4.69 | -0.33 |
| 2 m dewpoint (K) | 19–24 h | 56176 | 5.48 | -0.56 | 4.79 | 4.88 | -0.31 |
| 10 m wind (vector) (m/s) | 1–6 h | 51339 | 2.73 | 0.25 | 1.98 | 3.04 | -0.90 |
| 10 m wind (vector) (m/s) | 7–12 h | 49767 | 3.24 | 0.31 | 2.78 | 3.23 | -0.35 |
| 10 m wind (vector) (m/s) | 13–18 h | 48814 | 3.54 | 0.24 | 3.19 | 3.42 | -0.23 |
| 10 m wind (vector) (m/s) | 19–24 h | 48771 | 3.79 | 0.29 | 3.32 | 3.41 | -0.31 |
| sea-level pressure (hPa) | 1–6 h | 40526 | 2.37 | 0.56 | 1.51 | 7.12 | -1.47 |
| sea-level pressure (hPa) | 7–12 h | 40048 | 3.76 | 0.61 | 3.46 | 7.10 | -0.18 |
| sea-level pressure (hPa) | 13–18 h | 39547 | 5.15 | 0.72 | 4.98 | 6.95 | -0.07 |
| sea-level pressure (hPa) | 19–24 h | 39226 | 6.53 | 0.83 | 6.45 | 6.90 | -0.02 |

2 m temperature by cycle, 1–24 h (12 km): model RMSE / P24 RMSE (K)

10-27 12Z 4.8/5.4 better; 10-27 18Z 4.2/5.0 better; 10-28 00Z 3.3/4.6 better; 10-28 06Z 4.8/3.9 **worse**; 10-28 12Z 4.9/3.7 **worse**; 10-28 18Z 4.2/4.4 better; 10-29 00Z 4.1/5.0 better; 10-29 06Z 5.8/5.6 **worse**; 10-29 12Z 6.6/6.2 **worse**; 10-29 18Z 4.2/6.1 better; 10-30 00Z 4.8/6.0 better; 10-30 06Z 5.4/6.1 better

### 3 km — 29 cycles to 10-30 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 60512 | 3.78 | -0.53 | 5.39 | 4.49 | 0.29 |
| 2 m temperature (K) | 7–12 h | 59365 | 5.51 | -0.36 | 8.84 | 4.54 | -0.47 |
| 2 m temperature (K) | 13–18 h | 58012 | 5.62 | 0.03 | 8.67 | 4.66 | -0.45 |
| 2 m temperature (K) | 19–24 h | 46292 | 5.18 | 0.13 | 6.28 | 5.01 | -0.07 |
| 2 m dewpoint (K) | 1–6 h | 60109 | 3.63 | -1.45 | 2.37 | 4.50 | -1.35 |
| 2 m dewpoint (K) | 7–12 h | 59015 | 4.63 | -1.02 | 3.67 | 4.54 | -0.59 |
| 2 m dewpoint (K) | 13–18 h | 57672 | 5.32 | -0.54 | 4.30 | 4.69 | -0.53 |
| 2 m dewpoint (K) | 19–24 h | 46111 | 5.95 | -0.24 | 5.08 | 5.19 | -0.37 |
| 10 m wind (vector) (m/s) | 1–6 h | 51622 | 2.87 | 0.34 | 1.99 | 3.06 | -1.08 |
| 10 m wind (vector) (m/s) | 7–12 h | 50231 | 3.32 | 0.37 | 2.79 | 3.24 | -0.42 |
| 10 m wind (vector) (m/s) | 13–18 h | 49254 | 3.63 | 0.30 | 3.20 | 3.42 | -0.28 |
| 10 m wind (vector) (m/s) | 19–24 h | 39951 | 3.99 | 0.25 | 3.48 | 3.58 | -0.31 |
| sea-level pressure (hPa) | 1–6 h | 41033 | 2.67 | 0.92 | 1.51 | 7.12 | -2.13 |
| sea-level pressure (hPa) | 7–12 h | 40609 | 3.91 | 0.90 | 3.46 | 7.09 | -0.28 |
| sea-level pressure (hPa) | 13–18 h | 40074 | 5.23 | 0.96 | 4.99 | 6.95 | -0.10 |
| sea-level pressure (hPa) | 19–24 h | 32223 | 6.86 | 0.92 | 6.62 | 7.14 | -0.08 |

2 m temperature by cycle, 1–24 h (3 km): model RMSE / P24 RMSE (K)

10-27 12Z 4.4/5.4 better; 10-27 18Z 4.7/5.0 better; 10-28 00Z 3.2/4.6 better; 10-28 06Z 4.6/3.9 **worse**; 10-28 12Z 4.8/3.7 **worse**; 10-28 18Z 4.6/4.4 **worse**; 10-29 00Z 4.1/5.0 better; 10-29 06Z 5.5/5.6 better; 10-29 12Z 6.1/6.1 better; 10-29 18Z 4.6/6.1 better; 10-30 00Z 5.1/6.0 better; 10-30 06Z 5.5/6.0 better

### 12 km dry run (overlap) — 29 cycles to 10-30 06Z (newest verified run started 29 h ago)

| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |
|---|---|---|---|---|---|---|---|
| 2 m temperature (K) | 1–6 h | 60308 | 3.66 | -0.05 | 5.36 | 4.50 | 0.34 |
| 2 m temperature (K) | 7–12 h | 59018 | 5.27 | -0.29 | 8.84 | 4.54 | -0.35 |
| 2 m temperature (K) | 13–18 h | 57708 | 5.37 | -0.20 | 8.69 | 4.65 | -0.33 |
| 2 m temperature (K) | 19–24 h | 56438 | 4.84 | -0.13 | 5.99 | 4.78 | -0.02 |
| 10 m wind (vector) (m/s) | 1–6 h | 51339 | 2.62 | 0.04 | 1.98 | 3.04 | -0.75 |
| 10 m wind (vector) (m/s) | 7–12 h | 49767 | 3.15 | 0.18 | 2.78 | 3.23 | -0.28 |
| 10 m wind (vector) (m/s) | 13–18 h | 48814 | 3.47 | 0.11 | 3.19 | 3.42 | -0.18 |
| 10 m wind (vector) (m/s) | 19–24 h | 48771 | 3.65 | 0.10 | 3.32 | 3.41 | -0.21 |
| sea-level pressure (hPa) | 1–6 h | 40526 | 2.33 | 0.42 | 1.51 | 7.12 | -1.38 |
| sea-level pressure (hPa) | 7–12 h | 40048 | 3.77 | 0.51 | 3.46 | 7.10 | -0.19 |
| sea-level pressure (hPa) | 13–18 h | 39547 | 5.15 | 0.65 | 4.98 | 6.95 | -0.07 |
| sea-level pressure (hPa) | 19–24 h | 39226 | 6.50 | 0.76 | 6.45 | 6.90 | -0.02 |

2 m temperature by cycle, 1–24 h (12 km dry run (overlap)): model RMSE / P24 RMSE (K)

10-27 12Z 4.4/5.4 better; 10-27 18Z 3.9/5.0 better; 10-28 00Z 3.3/4.6 better; 10-28 06Z 4.6/3.9 **worse**; 10-28 12Z 4.8/3.7 **worse**; 10-28 18Z 4.1/4.4 better; 10-29 00Z 4.0/5.0 better; 10-29 06Z 5.3/5.6 better; 10-29 12Z 5.7/6.2 better; 10-29 18Z 4.0/6.1 better; 10-30 00Z 4.8/6.0 better; 10-30 06Z 5.3/6.1 better

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
| 2026-10-30T18Z | 10-30 19:45 | daily.sh | 16.8 | 0 |
| 2026-10-30T18Z | 10-30 19:45 | cam_cycle.sh | 132.0 | 0 |
| 2026-10-30T18Z | 10-30 21:57 | maps_3km | 10.1 | 0 |
| 2026-10-31T00Z | 10-31 01:45 | daily.sh | 16.7 | 0 |
| 2026-10-31T00Z | 10-31 01:45 | cam_cycle.sh | 134.0 | 0 |
| 2026-10-31T00Z | 10-31 03:58 | maps_3km | 10.1 | 0 |
| 2026-10-31T06Z | 10-31 07:45 | daily.sh | 16.1 | 0 |
| 2026-10-31T06Z | 10-31 07:45 | cam_cycle.sh | 127.4 | 0 |
| 2026-10-31T06Z | 10-31 09:52 | maps_3km | 10.0 | 0 |
| 2026-10-31T12Z | 10-31 13:45 | daily.sh | 16.6 | 0 |
| 2026-10-31T12Z | 10-31 13:45 | cam_cycle.sh | 141.6 | 0 |
| 2026-10-31T12Z | 10-31 16:06 | maps_3km | 10.1 | 0 |

## Server

- load average 0.00 / 0.00 / 1.08 (104 logical CPUs)
- memory available 387 of 394 GB
- next cron cycle 19:45Z (in 175 min)
- no forecast running

