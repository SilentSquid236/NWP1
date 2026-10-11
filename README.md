# NWP1 status — 2026-10-18 01:00Z

Rebuilt every 10 min on the school server from NWP1's own logs. Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.

## Attention

- cron 2026-10-17T00 cam_cycle.sh status 2
- cron 2026-10-17T06 cam_cycle.sh status 2
- cron 2026-10-17T12 cam_cycle.sh status 2
- cron 2026-10-17T18 cam_cycle.sh status 2

## Notes from the AI collaborator (newest first)

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
- **2026-10-13 19:20Z — R1 refuted (my design flaw).** Fading the edge vapour as a mixing ratio supersaturated the cold edges; rain went from 217 to 1306 mm. Rewritten to fade relative humidity instead (cannot exceed 100 %); R1b tests it after the 19:45Z cron cycle, with a sanity gate first. Result ~23Z.
- **2026-10-13 17:15Z — P-81 cause found, test R1 running.** A new water counter shows the heavy rain is fed through the domain's outer faces: the edges hold the 00Z tropical humidity and winds all day (no outside forecast), importing 2.5x the domain's rain; no conservation bug. R1 tests fading the edge humidity (6 h / 12 h) on 4 October and 20 September cycles at 12 km (~75 min). Also: a 4-day-old stuck E4 submission started by itself at 13:06Z; I stopped it after 3 min and rebuilt the 6 archives it deleted.
- **2026-10-13 04:30Z — first moist cron cycle ran cleanly (124 min), but its rain is far too heavy.** 411 mm in 24 h near Bancroft, Ontario, where Environment Canada called for 40-60 mm (tropical remnants). Right place, 7-10x too much; the moist 12 km parent does it too (217 mm). Not an edge or saturation bug; likely the warm-rain-only physics plus explicit convection at 12 km (P-81). Treat moist 3 km rain amounts as unreliable for now. Waiting for your call on keeping moist on.
- **2026-10-12 21:30Z — baseline for scoring the cron 3 km.** Re-verified with one operator so old and new runs compare fairly. The old-code dry 3 km (3 cycles, 10-10 12Z to 10-11 00Z) had 2 m T RMSE 5.9 K against 4.4 K for the 12 km (it had no land surface; daytime bias -1.8 K); winds about equal. The new dry and moist cycles get the same treatment as their windows close.

## Experiments (active in the last 24 h)

_None active._

## Operational cron cycles (last 3 days)

| cycle | started (UTC) | script | wall (min) | status |
|---|---|---|---|---|
| 2026-10-16T12Z | 10-16 13:45 | daily.sh | 16.9 | 0 |
| 2026-10-16T12Z | 10-16 13:45 | cam_cycle.sh | 161.1 | 2 ⚠ |
| 2026-10-16T18Z | 10-16 19:45 | daily.sh | 16.5 | 0 |
| 2026-10-16T18Z | 10-16 19:45 | cam_cycle.sh | 158.6 | 2 ⚠ |
| 2026-10-17T00Z | 10-17 01:45 | daily.sh | 17.1 | 0 |
| 2026-10-17T00Z | 10-17 01:45 | cam_cycle.sh | 156.8 | 2 ⚠ |
| 2026-10-17T06Z | 10-17 07:45 | daily.sh | 16.6 | 0 |
| 2026-10-17T06Z | 10-17 07:45 | cam_cycle.sh | 156.1 | 2 ⚠ |
| 2026-10-17T12Z | 10-17 13:45 | daily.sh | 17.1 | 0 |
| 2026-10-17T12Z | 10-17 13:45 | cam_cycle.sh | 157.8 | 2 ⚠ |
| 2026-10-17T18Z | 10-17 19:45 | daily.sh | 16.8 | 0 |
| 2026-10-17T18Z | 10-17 19:45 | cam_cycle.sh | 159.0 | 2 ⚠ |

## Server

- load average 0.00 / 0.86 / 6.94 (104 logical CPUs)
- memory available 386 of 394 GB
- next cron cycle 01:45Z (in 45 min)
- no forecast running

