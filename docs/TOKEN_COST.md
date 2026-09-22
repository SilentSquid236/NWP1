# What this project costs

Two numbers, and conflating them would be the easiest way to get the write-up
wrong.

| | what it is | current figure |
|---|---|---|
| **BILLED** | money that left the account — the subscription, prorated over the months the project has run | **$400** (2 months × $200) |
| **API-EQUIV** | what the same tokens would cost at list API prices | **$79**, uncertain by ~2× ($40–$158) |

`python tools/tokens.py --report` produces both from `docs/token_ledger.csv`.

## Why both

The project is built through a Claude subscription, so the marginal cost of a
token is zero: the cost is the plan, whatever is used. That is the number that
answers "what did you pay."

It is not the number a paper can use. *"What does it cost to build a regional
NWP model with an AI collaborator"* has to be asked in units someone else can
reproduce, and the reproducible unit is API pricing. So the ledger carries a
shadow price alongside the real one.

**API-EQUIV must never be described as money spent.** It is a conversion, for
comparison only.

## The result so far, which runs the wrong way

At list prices the tokens are worth about **a fifth** of what the subscription
costs — 0.20× per dollar billed. The obvious expectation was the opposite, that
an agentic project like this would be extracting many times the plan price in
tokens.

It does not, and the reason is visible in the research log. The expensive parts
of this project have not been tokens:

- A 12-hour forecast sweep costs **wall-clock**, not tokens. Several sessions
  spent most of their duration waiting on runs of 10–60 minutes each.
- The costly mistakes cost **rework**, not context. Nine failed patch cycles
  against a model that was not broken, a K_MAX ladder that varied nothing, four
  test-design defects — each burned sessions, and sessions are mostly waiting.
- Two-thirds of the turns are short: a measurement comes back, it is read, the
  next one is launched.

So the honest framing for the write-up is that the AI-assisted cost of this
project is dominated by **elapsed time and human attention**, and the token bill
is close to a rounding error against the subscription. That is a more
interesting claim than a large token number would have been, and it is
falsifiable: measured rows will either hold it up or not.

## What is measured and what is not

Every row is tagged. The four rows dated before 2026-09-16 are **estimates**,
not recovered measurements — there is no way to read a past session's token
counts from inside a later one.

The estimate is turn count (a hard number from `docs/PROMPT_LOG.md`) times an
assumed number of tool-calling steps per turn, times an assumed average
conversation size re-sent per step. That last parameter is the weak one: cache
reads dominate the token count, and the average context is uncertain by roughly
a factor of two, which propagates straight through. Hence the ±2× band. Every
parameter is written into the `note` field of each row, so the arithmetic can
be redone with different assumptions.

`python tools/tokens.py --estimate` prints the method.

## Recording a session — do it on the day

This is the same shape of problem as the verification archive (P-07): the
number exists at the end of the session and not afterwards.

    python tools/tokens.py --add 2026-09-16 --in 120000 --cache-read 8400000 \
        --cache-write 950000 --out 145000 --basis measured \
        --note "radiative boundary, P-50 probes"

Read the counts from the session's usage breakdown before the session ends.
One measured session is worth more than all four estimated rows, because it
calibrates the step-count and context assumptions that the estimate rests on.

## Prices

`python tools/tokens.py --prices`. They are dated in the source
(`PRICES_ASOF`), because an undated number cannot expire — L7 in the learning
log, learned the hard way when a constraint about the sponge outlived the thing
it described by four days.

Checked 2026-09-16 against a pricing aggregator, **not** against Anthropic
directly. Re-verify at `docs.claude.com/en/docs/about-claude/pricing` before
any figure from this tool is published.

`PLAN_USD_PER_MONTH` and `PLAN_NAME` at the top of `tools/tokens.py` are
assumptions about the subscription and should be corrected to whatever is
actually being paid.
