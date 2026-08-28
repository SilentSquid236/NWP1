"""
Tests for network politeness.

The rate-limiter arithmetic is tested WITHOUT sleeping -- consume() returns
the wait it would take, so the token-bucket maths is verified in milliseconds
instead of minutes. Only one test actually sleeps, to confirm the throttle
does take real time.

Run:  python test_netpolicy.py
"""

import tempfile
import time
from pathlib import Path

from netpolicy import (RateLimiter, DownloadCache, PoliteFetcher,
                       estimate_ingest_mb, max_mbps)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# ---------------------------------------------------------------------------
def test_burst_allowed_then_throttled():
    """
    One second's worth of transfer should pass without waiting; beyond that
    the bucket empties and the caller must wait. This is what protects
    interactive users -- a short burst is fine, a sustained stream is not.
    """
    rl = RateLimiter(mbps=10.0, burst_seconds=1.0)

    first = rl.consume(10_000_000, sleep=False)      # exactly the burst
    second = rl.consume(10_000_000, sleep=False)     # another 10 MB

    ok = first == 0.0 and 0.9 < second < 1.1
    report("burst allowed, sustained transfer throttled", ok,
           f"first 10 MB waited {first:.3f}s, next 10 MB waited {second:.3f}s "
           f"(expect ~1.0s at 10 MB/s)")


# ---------------------------------------------------------------------------
def test_wait_scales_with_size():
    """Wait time must be proportional to bytes over the cap."""
    waits = []
    for mb in (5, 10, 20):
        rl = RateLimiter(mbps=5.0, burst_seconds=0.0)
        waits.append(rl.consume(mb * 1_000_000, sleep=False))

    r1 = waits[1] / waits[0]
    r2 = waits[2] / waits[1]
    ok = abs(r1 - 2.0) < 0.05 and abs(r2 - 2.0) < 0.05
    report("throttle time scales linearly with transfer size", ok,
           f"5/10/20 MB at 5 MB/s -> {waits[0]:.1f}/{waits[1]:.1f}/"
           f"{waits[2]:.1f}s (ratios {r1:.2f}, {r2:.2f})")


# ---------------------------------------------------------------------------
def test_tokens_refill_over_time():
    """
    Waiting should restore capacity -- otherwise the limiter would penalise a
    job that politely paused, which is backwards.
    """
    rl = RateLimiter(mbps=10.0, burst_seconds=1.0)
    rl.consume(10_000_000, sleep=False)             # empty the bucket

    rl.last -= 0.5                                   # pretend half a second passed
    wait_after_pause = rl.consume(5_000_000, sleep=False)

    ok = wait_after_pause < 0.01
    report("bucket refills while idle", ok,
           f"after a 0.5 s pause, 5 MB waited {wait_after_pause:.4f}s "
           f"(refill covered it)")


# ---------------------------------------------------------------------------
def test_throttle_actually_sleeps():
    """One test that spends real time, confirming the sleep path works."""
    rl = RateLimiter(mbps=2.0, burst_seconds=0.0)
    t0 = time.monotonic()
    rl.consume(500_000, sleep=True)                  # 0.5 MB at 2 MB/s = 0.25 s
    elapsed = time.monotonic() - t0

    ok = 0.2 < elapsed < 0.45
    report("throttling spends real wall-clock time", ok,
           f"0.5 MB at 2 MB/s took {elapsed:.3f}s (expect ~0.25s)")


# ---------------------------------------------------------------------------
def test_cache_avoids_second_download():
    """The cheapest bandwidth saving is not downloading twice."""
    with tempfile.TemporaryDirectory() as d:
        cache = DownloadCache(d)
        url = "https://example.invalid/data.grib2"

        miss = cache.get(url)
        cache.put(url, b"payload-bytes")
        hit = cache.get(url)

        ok = (miss is None and hit == b"payload-bytes"
              and cache.hits == 1 and cache.misses == 1)
        report("cache returns stored content on a repeat request", ok,
               f"first {miss}, second {hit!r}, "
               f"{cache.hits} hit / {cache.misses} miss")


# ---------------------------------------------------------------------------
def test_cache_write_is_atomic():
    """
    A crash mid-write must not leave a truncated file that later reads as
    valid data -- silently corrupt input is worse than a re-download.
    """
    with tempfile.TemporaryDirectory() as d:
        cache = DownloadCache(d)
        url = "https://example.invalid/x"
        p = cache.put(url, b"complete")

        partials = list(Path(d).rglob("*.partial"))
        ok = p.exists() and p.read_bytes() == b"complete" and not partials
        report("cache writes atomically, leaving no partial files", ok,
               f"final file intact, {len(partials)} .partial left behind")


# ---------------------------------------------------------------------------
def test_fetcher_is_sequential_by_construction():
    """
    No connection pool, no thread pool, and a minimum gap between requests.
    Concurrency is the fastest way to saturate a shared link.
    """
    f = PoliteFetcher(mbps=5.0, gap=0.2)
    has_pool = any(hasattr(f, a) for a in
                   ("pool", "executor", "session", "threads", "workers"))

    t0 = time.monotonic()
    f._wait_gap()
    f._wait_gap()
    elapsed = time.monotonic() - t0

    ok = not has_pool and elapsed >= 0.2
    report("fetcher is sequential with an enforced request gap", ok,
           f"no pool/executor attribute: {not has_pool}; "
           f"two requests spaced {elapsed:.2f}s apart (gap 0.2s)")


# ---------------------------------------------------------------------------
def test_default_cap_is_conservative():
    """
    The default must be well below a typical link so it is safe unattended.
    A research server usually has 1 Gb/s (~125 MB/s); 8 MB/s is ~6% of that.
    """
    cap = max_mbps()
    ok = 1.0 <= cap <= 20.0
    report("default bandwidth cap is conservative", ok,
           f"{cap:.1f} MB/s = {cap*8:.0f} Mb/s, roughly "
           f"{cap/125*100:.0f}% of a 1 Gb/s link")


# ---------------------------------------------------------------------------
def test_ingest_estimate_is_sane():
    """Knowing the cost in advance beats discovering it from a complaint."""
    e = estimate_ingest_mb(hours=13)
    ok = (0 < e["per_hour_download_MB"] < 200
          and e["total_download_MB"] > e["per_hour_download_MB"]
          and e["minutes_at_cap"] > 0)
    report("ingest download size can be estimated up front", ok,
           f"13 h ingest ~ {e['total_download_MB']:.0f} MB, "
           f"{e['minutes_at_cap']:.1f} min at the cap "
           f"({e['per_hour_download_MB']:.0f} MB/hour)")


if __name__ == "__main__":
    print("\nNetwork politeness\n" + "=" * 62)
    for fn in (test_burst_allowed_then_throttled,
               test_wait_scales_with_size,
               test_tokens_refill_over_time,
               test_throttle_actually_sleeps,
               test_cache_avoids_second_download,
               test_cache_write_is_atomic,
               test_fetcher_is_sequential_by_construction,
               test_default_cap_is_conservative,
               test_ingest_estimate_is_sane):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)
