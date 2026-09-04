# Prompt log

Every instruction given by the human collaborator, in order, verbatim, with
what it caused. This is the **input record** for the AI-to-build study: the
research log records what the AI did, and this records what it was asked.

Verbatim means verbatim — typos and all. Prompts are not cleaned up, because
one of the questions the study can actually answer is how much precision a
direction needs in order to land.

## Why keep this

Three things are only visible from the input side:

1. **Where direction came from.** The AI proposed nearly all of the technical
   detail. Almost none of the *turns* were its idea. Separating the two is the
   study's central measurement, and it cannot be reconstructed from the code.
2. **Cost of a wrong direction.** A prompt that sent work down a dead end is
   as informative as one that worked, and the dead ends disappear from a
   repository unless they are written down.
3. **Prompt length against effect.** Several of the highest-leverage prompts
   here are under twelve words.

## Classification

| tag | meaning |
|---|---|
| **DIR** | direction-setting: chooses what to build or in what order |
| **CON** | constraint: a fact about the environment that bounds every design |
| **COR** | correction: the AI was wrong or off-track and was redirected |
| **MET** | methodological: changes *how* the work is done, not what |
| **OBS** | observation: the human reporting what they saw, usually pasted output |
| **ADM** | administrative: tooling, transfer, logistics |

---

## Phase 1 — Environment and access (prompts 1-16)

| # | prompt | tag | effect |
|---|---|---|---|
| 1 | "Can you familiarize yourself with the file structure in the files and get the background of the project from the pdf" | DIR | established the project as a resumed effort, not a fresh start |
| 2 | "Lets get WSL done out1 should be fixed the server will handle the processing not my PC" | DIR | moved execution off Windows |
| 3 | "I also cant install any new packages" | **CON** | the single most binding constraint in the project |
| 4 | "since I share this server with 30 other people I cannot make environments I have no sudo powers" | **CON** | ruled out venvs, conda, containers — everything downstream assumes this |
| 5 | "can I integrate you with termius to streamline this task" | ADM | explored; terminals proved click-only |
| 6 | "lets let you take control and I will approve" | ADM | computer-use attempt, later abandoned |
| 7 | "you should have access the 242 ip is where you need to go" | ADM | — |
| 8 | "try again" | ADM | — |
| 9 | "so the project was originally built to run on my local pc when I had no access to the server but now can be made for linux" | DIR | reframed the port as the task |
| 10 | "how about we use github" | ADM | transfer mechanism |
| 11 | "can we make it so you scp the files over and I run them" | ADM | settled the working loop: AI writes, human runs |
| 12 | "https://github.com/SilentSquid236/NWP1 / Continue from where you left off." | ADM | — |
| 13 | "my main pc is windows" | CON | — |
| 14 | "password authentication not supported" | OBS | — |
| 15 | "sadly git is not on the server" | **CON** | killed the git-based transfer plan outright |
| 16 | "its on the server" | OBS | — |

**Observation.** Sixteen prompts before a line of physics. Four of them
(3, 4, 13, 15) are environment facts that no amount of AI reasoning could have
produced, and each invalidated a plan the AI had already proposed.

---

## Phase 2 — Resource discipline (prompts 17-21, 37)

| # | prompt | tag | effect |
|---|---|---|---|
| 17 | "make sure we only ever use 50% of the resources on the server unless specified to use more" | **CON** | produced `resources.py` |
| 18 | "Does that count the fact we only have access to the resources of data5" | COR | caught that the ceiling was being computed against the wrong pool |
| 19 | "I need a summary for github" | ADM | — |
| 20 | "Everything is on the xenon" | OBS | — |
| 21 | "torch is there this server is meant for meteorology research" | OBS | — |
| 37 | "ah we must make sure we dont use the entire bandwith" | **CON** | produced `netpolicy.py`: token bucket, cache, sequential fetch |

**Observation.** Prompt 18 is a correction disguised as a question, and it is
the pattern that recurs: the human notices a *scope* error — right computation,
wrong domain — that internal review did not catch. Prompt 37, nine words,
generated an entire module.

---

## Phase 3 — First real data (prompts 22-40)

| # | prompt | tag | effect |
|---|---|---|---|
| 22 | "I want to do the NE region of the CONUS like how the hrrr has it NE region" | DIR | fixed the domain |
| 23 | "Honestly we should also increase the levels to 20" | DIR | 20 pressure levels |
| 24 | "Sounds good" | — | — |
| 25 | "Lets get ready to do a test run" | DIR | first contact with real data |
| 26 | [pasted ingest failure output] | OBS | GRIB regex bug |
| 27 | "there are 2 src folders" | **COR** | human caught `src/src` — the AI had not |
| 28-30 | [pasted `ls`, diagnostic, and CF-name output] | OBS | three interface bugs found by real data |
| 31 | "It worked with nothing written" | OBS | silent write failure |
| 32, 34 | [pasted forecast divergence output] | OBS | the instability, first seen |
| 33 | "we are under 5" | OBS | — |
| 35 | "With what we have done so far what can I expect from this" | **MET** | forced `docs/CAPABILITIES.md` — an honest statement of limits, written before there was pressure to overclaim |
| 36 | "Yeah. We should also do a run with this on the xenon how long would it take to run 12 hours out" | DIR | — |
| 38 | "passed" | OBS | — |
| 39 | "lets keep moving" | DIR | — |

**Observation.** Prompt 35 is the highest-value prompt of this phase and asks
for nothing to be built. Asking "what can I actually expect" before the model
worked produced a document that later kept the project honest about failures.

---

## Phase 4 — The pivot (prompts 41-50)

| # | prompt | tag | effect |
|---|---|---|---|
| 41 | "this isnt about the NN anymore we are going full weather prediction model using the things learned from the NN" | **DIR** | *the* decision of the project — abandoned the neural emulator for real physics |
| 42 | "We must start from the basics of 2D its fundamental" | **MET** | forced shallow water first. Every C-grid operator bug was found in 2D, where it was cheap |
| 43 | "lets keep moving do visualization at the end" | DIR | — |
| 44 | "Lets do dry conditions before we touch wet at all" | **MET** | deferred moisture; the dry core is still not finished, which vindicates it |
| 45 | "We need to factor in no lab conditions there needs to be variance in the equations that allow evolving conditions" | DIR | produced stochastic physics (SPPT) |
| 46 | "Im thinking we should also integrate the NN to take a predicted output compare it to recent patterns and adjust the output" | DIR | the NN returns as post-processing, not as the model |
| 47 | "lets keep moving" | DIR | — |
| 48 | "make sure our data isnt coming from the hrrr but observations, asos, mesonets, and radar products" | **CON** | verification may never use model output. Shapes `src/verification/` entirely |
| 49 | "Lets start documenting everything we do I am thinking of doing a research project of AI and Weather Modeling" | **MET** | created the research log |
| 50 | "AI to build" | **MET** | two words. Redefined the study's subject from *AI as forecaster* to *AI as builder* — a completely different paper |

**Observation.** Prompts 42 and 44 are ordering constraints, and both were
right. The AI's own instinct was to attempt the full 3D moist system early;
being forced down to 2D dry is why the operator bugs were found at all.
Prompt 48 is a scientific-integrity constraint the AI would not have imposed
on itself.

---

## Phase 5 — Debugging the instability (prompts 51-60)

| # | prompt | tag | effect |
|---|---|---|---|
| 51 | "lets get back into it" | DIR | — |
| 52 | "yes" | DIR | approved the sigma-coordinate port |
| 53 | **"I think taking a step back and probing the error is a better idea then guess checking"** | **MET** | the pivotal prompt. Came after five failed patch cycles and ended the patch-and-pray loop |
| 54 | "Try using really tall terrain vs 0 terrain to get a better idea of extreme values that could be affecting output" | **MET** | the experimental design that isolated terrain as a variable |
| 55 | "Lets put that idea into effect reformulate the integration paying attention to how coriolis and pgf interact with the boundary layer that could throw problems in" | **DIR** | correctly predicted the PGF residual would be largest at the *bottom*. Measured: 2.1e-03 at the surface against 3.0e-05 at the top |
| 56 | "can I get a visualization of whats happening" | MET | produced the cross-sections that showed lid reflection |
| 57 | "look into interactions with elevated terrain such as orographic effects" | DIR | — |
| 58 | "Yeah lets continue debugging the flat terrain" | DIR | — |
| 59 | "Yes lets work on surface interactions" | DIR | produced `surface.py`, and the sweep that exposed the initialization defects |
| 60 | "We should log prompts given and create a tree diagram of the project structure" | MET | this file and `docs/STRUCTURE.md` |
| 61 | "lets continue fixing the failures at higher levels" | DIR | traced the 2500 m failure to mountain-wave overturning; produced `convection.py` |
| 62 | "Lets keep going" | DIR | probed 4000 m; eliminated timestep and diffusivity ceiling; identified Nh/U ~ 1 as the regime boundary |
| 63 | "Lets take problems and document them and when one is fixed add an explanation of what was done to fix it" | **MET** | produced `docs/PROBLEMS.md` and `tools/problem.py`; the audit caught 13 fixes asserted without a measurement |
| 64 | "I agree with the mountain point as it stands if we can hit 2km of terrain thats a good start and we can move on to other problems" | **DIR** | closed P-01 at an agreed target rather than an open-ended one; freed the effort that ported the driver to the sigma core |
| 65 | "Lets continue" | DIR | took P-46; three hypotheses, two wrong, and the defect turned out to be the test |
| 66 | "lets tackle p-07" | **DIR** | built the archive machinery; found a 74 K surface-observation trap in the old operator |
| 67 | "What is the command to pull from the github again I think files got installed in the wrong folders" | ADM | produced `tools/checklayout.py` |
| 68 | "I would like to use the curl method" | ADM | `tools/pull.sh` — git is not installable on the server |
| 69 | "even after pull.sh I am missing files" | **OBS** | the files had never been pushed; produced `tools/manifest.py` |
| 70 | "All files should be in the NWP_deployment_package folder" | COR | the archive had a wrapper directory and nested the package inside itself |
| 71 | "whe nrunning the archive everything freezes" | **OBS** | three buffering/progress defects (P-47) and a raw-observation fetch bug (P-48) |
| 72 | "Log that the model effort was switched to high I would also like all files here C:\\Users\\Epier\\Desktop\\NWP\\NWP_Deployment_Package" | **MET**, ADM | recorded the instrument change; synced the working copy |

**Observation.** Prompt 53 is 18 words and is the most consequential
instruction in the project. Before it, nine candidate causes had been patched
and measured one at a time, all negative. After it, the failure was traced in
a single session to three defects in the *test setup* — a super-geostrophic
clipped jet, a one-term geostrophic wind over sigma terrain, and unfiltered
white noise. The AI had been debugging a model that was not broken.

Prompt 55 is the clearest case of domain intuition beating the AI's search
order. The AI was examining the model top; the human said look at the boundary
layer; the residual was at the bottom, by a factor of 70.

---

## Aggregate

| tag | count | share |
|---|---|---|
| DIR direction-setting | 22 | 37% |
| OBS observation / pasted output | 13 | 22% |
| ADM administrative | 10 | 17% |
| MET methodological | 9 | 15% |
| CON constraint | 8 | 13% |
| COR correction | 3 | 5% |

(Counts exceed the prompt total slightly: a few prompts carry two tags.)

Prompt 61 is worth noting for length against effect. Eight words, no
technique named, and the ambiguity in "higher levels" — terrain height or
model level — turned out not to matter, because the answer was at a specific
model level over specific terrain and the measurement found it either way.

**Median prompt length: 11 words.** The four prompts that changed the project
most — 41, 50, 53, 55 — average 19 words. None specifies a technique. Every
one specifies a *direction* or a *method*, and leaves the technique to the AI.

**The AI never once proposed:** starting from 2D, deferring moisture, banning
HRRR from verification, capping server usage, or probing instead of patching.
It proposed nearly every equation, discretisation and test in the repository.

That split is the study's main finding so far, and it is the reason this file
exists.

---

## Maintaining this log

Append each new prompt with its tag and effect as it happens, not in a batch
afterwards. Reconstructing intent later is exactly the kind of self-report the
methodology section warns about.
