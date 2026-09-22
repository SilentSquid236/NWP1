---
name: nwp-record-session
description: Use at the end of any NWP1 working session, or when a problem is fixed, to update the research record — research log, problem register, prompt log, learning log and token ledger — and run the checks that must pass before work is called done.
---

# Recording an NWP1 session

The research record is half this project's deliverable: the subject of study is
the AI collaboration itself. Some of it can only be captured on the day.

## Always

1. **Research log** — `docs/RESEARCH_LOG.md`, a dated entry above the
   "Recording for the AI-collaboration study" heading. Shape:
   Context · Hypothesis (stated before the result) · Method · Result (numbers,
   tables) · Interpretation (what it means and what it does *not* mean) ·
   Status. Failed hypotheses are recorded, not dropped.

2. **Problem register** — `docs/PROBLEMS.md`, edited in place. Statuses: OPEN,
   FIXED, ELIMINATED, REVERTED, ACCEPTED. A FIXED entry must carry a
   **Confirmed by** line with a measurement. `python tools/problem.py new "<title>"`
   appends a skeleton.

3. **Prompt log** — `docs/PROMPT_LOG.md`, every human prompt verbatim with a tag
   (DIR, CON, COR, MET, OBS, ADM) and its effect.

4. **Token ledger** — do this before the session ends; the counts do not exist
   afterwards:

       python tools/tokens.py --add <date> --in N --cache-read N \
           --cache-write N --out N --basis measured --note "<what it did>"

## When it applies

5. **Learning log** — `docs/LEARNING_LOG.md`. Add a callback only when an
   earlier lesson visibly changed what happened next; mark it H (human
   prompted), A (AI unprompted) or T (a tool enforced it). A new lesson gets its
   own entry with its origin.

6. **Configuration changes** to the collaborator (model, effort level) go in
   the "Instrument changes" table in `docs/AI_COLLABORATION.md`. They are a seam
   in the study and counts should not be pooled across them.

## Before calling anything done

    python tools/problem.py check
    python tools/checklayout.py
    python tools/tree.py
    python tools/manifest.py
    python tools/manifest.py --check

and every test suite the change could have touched. A suite that was not run is
not evidence.
