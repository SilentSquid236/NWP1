# Bringing NWP1 into Claude Science

Claude Science is Anthropic's research workbench: a desktop app (macOS,
Windows, Linux) in which each project keeps its own memory, data sources,
installed skills and artifacts. Skills are `SKILL.md` files and can be imported
from GitHub. It can also reach remote machines over SSH.

**What was confirmed and what was not** (checked 2026-09-22 from Anthropic's
announcement and an independent review): projects, per-project memory, and
`SKILL.md` skills importable from GitHub are confirmed. How an existing
repository or a `CLAUDE.md` context file is picked up is **not** documented in
either source. The steps below are written to work whichever way it turns out
to behave, and step 2 is the one to check first.

## Status after the first Claude Science session (2026-09-22)

| step | outcome |
|---|---|
| 1. project | exists |
| 2. context | **not picked up automatically** — project memory was empty at start. `CLAUDE.md` was found on request at `Desktop\NWP\NWP1\CLAUDE.md`, the five constraints were listed back, and they, the method and the known pitfalls were written into project memory |
| 3. skills | `nwp-debug`, `nwp-record-session`, `nwp-sync` imported as personal skills; `nwp-sync` first revised for git on the server |
| 4. research record | not copied in: the session reads `docs/` in place from the granted folder `Desktop\NWP\NWP1` |
| 5. record the move | row in "Instrument changes"; research-log entry of the same date |

The working copy on the desktop is `NWP1\`. `NWP_Deployment_Package\` is the
git clone; it had not been fetched and still showed `main` at `b97bc06`. GitHub
itself had already merged `package/claude-science` (#2) and
`p40/ceiling-ladder` (#3): the server's `main` is at `44068f2`.

## Do not install it on the shared server

The Xeon has an admin policy against installing new packages, and that covers
this app. Install it on the **Windows desktop**. If working directly on the
server is wanted later, that is a question for the server's admin, not a
workaround to find.

## Steps

1. **Create a project** in Claude Science named `NWP1 — AI-built weather model`.

2. **Give it the repository.** Either connect GitHub
   (`SilentSquid236/NWP1`, branch `main` once the package branch is merged), or
   point it at the local clone `Desktop\NWP\NWP_Deployment_Package`, or upload
   the zip. Then open a session and ask it to read `CLAUDE.md` and summarise the
   constraints back. **If it cannot list the five constraints in
   `CLAUDE.md`, it has not picked up the context, and nothing else will go
   right** — paste `CLAUDE.md` into the project's instructions or memory by hand.

3. **Import the three skills** from `skills/`:
   - `nwp-debug` — probe-first failure diagnosis
   - `nwp-record-session` — end-of-session bookkeeping
   - `nwp-sync` — moving work between GitHub, the desktop and the server

   These carry the project's working habits. The learning log found that most
   lessons in this project survived only as habits, which a fresh session does
   not have; a skill is how a habit survives the move.

4. **Add the research record as project knowledge** if the app does not index
   the repository on its own: `docs/RESEARCH_LOG.md`, `docs/PROBLEMS.md`,
   `docs/LEARNING_LOG.md`, `docs/PROMPT_LOG.md`, `docs/AI_COLLABORATION.md`,
   `docs/TOKEN_COST.md`.

5. **Record the move.** It is a change of instrument for the study: add a row
   to the "Instrument changes" table in `docs/AI_COLLABORATION.md` (new
   environment, and the model in use), and a research-log entry. Counts in the
   study should not be pooled across it without saying so.

## What does not come across

- **This conversation.** Everything a new session needs was written into
  `CLAUDE.md`, the skills and `docs/`. If something important is missing from
  those, it is lost — that is a reason to check step 2 carefully.
- **The claude.ai project docs** that have no repository counterpart:
  `claude/sync-state.md` (operational notes on the device bridge; its procedure
  is carried by the `nwp-sync` skill) and `claude/initialization-findings.md`
  (an excerpt of a research-log entry that is already in the repository).
- **`data/`**, which lives only on the server and should stay there.
- **Past token usage.** Sessions before 2026-09-16 are estimates in
  `docs/token_ledger.csv`. Record Claude Science sessions as measured rows from
  the start.

## First things to do in the new environment

In priority order, matching the problem register:

1. **P-07 / P-06** — run `bash tools/daily.sh` on the server once (it is a bash script; `python tools/daily.sh` fails with a SyntaxError), by hand, and
   bring back the log. The verification archive is the only item that gets
   permanently more expensive each day it stays empty, and the live fetch paths
   have never been exercised.
2. **P-52** — detect a blowing-up run inside a forecast hour, not at its end.
3. **P-50** — the fourth cause is shear at the model top that the sponge was
   masking; P-49 is solved in principle by the radiative lid if P-50 can be.
