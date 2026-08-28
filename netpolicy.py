"""
Network politeness for a shared machine.

CPU is not the only shared resource. This box has ~30 users behind one
connection, and a download loop that saturates the link is far more disruptive
than one that saturates a core -- everyone's ssh session, file transfer, and
data fetch degrades at once, and it is not obvious to them why.

Policy, mirroring the CPU governor in resources.py:

  * a DEFAULT CAP on download rate, well below the link capacity
  * SEQUENTIAL requests -- no parallel connection pools
  * a pause between requests, so bursts do not queue up
  * exponential backoff on failure, so a struggling server is not hammered
  * a local CACHE, so nothing is ever downloaded twice

Override with NWP_MAX_MBPS when you know the link is quiet, exactly as
NWP_RESOURCE_FRACTION overrides the CPU cap.

The rate limiter is a token bucket: tokens accumulate at the target rate, a
transfer spends tokens equal to its byte count, and the caller sleeps when the
bucket runs dry. It smooths bursts rather than merely capping an average,
which is what actually protects interactive users.
"""

import hashlib
import os
import time
from pathlib import Path

DEFAULT_MAX_MBPS = 8.0          # megabytes/sec -- polite on a shared link
DEFAULT_GAP_S = 0.5             # minimum pause between requests
DEFAULT_RETRIES = 4


def max_mbps():
    return float(os.environ.get("NWP_MAX_MBPS", DEFAULT_MAX_MBPS))


class RateLimiter:
    """
    Token bucket. Capacity is one second's worth of transfer, so a short burst
    is allowed but a sustained stream is held at the target rate.
    """

    def __init__(self, mbps=None, burst_seconds=1.0):
        self.rate = float(mbps if mbps is not None else max_mbps()) * 1e6
        self.capacity = self.rate * burst_seconds
        self.tokens = self.capacity
        self.last = time.monotonic()
        self.total_bytes = 0
        self.total_wait = 0.0

    def _refill(self, now=None):
        now = now if now is not None else time.monotonic()
        elapsed = max(0.0, now - self.last)
        self.tokens = min(self.capacity, self.tokens + elapsed * self.rate)
        self.last = now

    def consume(self, n_bytes, sleep=True):
        """
        Account for n_bytes of transfer, returning the seconds to wait.

        Separating the calculation from the sleeping makes the arithmetic
        testable without spending real time.
        """
        self._refill()
        self.total_bytes += n_bytes
        self.tokens -= n_bytes

        wait = 0.0
        if self.tokens < 0:
            wait = -self.tokens / self.rate
            self.total_wait += wait
            if sleep:
                time.sleep(wait)
                self._refill()
        return wait

    @property
    def effective_mbps(self):
        elapsed = max(1e-9, time.monotonic() - (self.last - self.total_wait))
        return self.total_bytes / 1e6 / elapsed

    def __repr__(self):
        return (f"RateLimiter({self.rate/1e6:.1f} MB/s, "
                f"{self.total_bytes/1e6:.1f} MB moved, "
                f"{self.total_wait:.1f} s throttled)")


class DownloadCache:
    """
    Content-addressed cache keyed by URL.

    The cheapest bandwidth saving is not re-downloading. Re-running an ingest
    after a crash, or re-fetching observations while debugging, should cost
    nothing.
    """

    def __init__(self, root=None):
        self.root = Path(root or os.environ.get(
            "NWP_CACHE_DIR", Path.home() / ".cache" / "nwp"))
        self.hits = 0
        self.misses = 0

    def path_for(self, url):
        h = hashlib.sha256(url.encode()).hexdigest()[:24]
        return self.root / h[:2] / h

    def get(self, url):
        p = self.path_for(url)
        if p.exists():
            self.hits += 1
            return p.read_bytes()
        self.misses += 1
        return None

    def put(self, url, data):
        p = self.path_for(url)
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(".partial")
        tmp.write_bytes(data)
        tmp.replace(p)          # atomic: a crash never leaves a half file
        return p

    def __repr__(self):
        return f"DownloadCache({self.root}, {self.hits} hits, {self.misses} misses)"


class PoliteFetcher:
    """
    Sequential, rate-limited, cached HTTP fetching with backoff.

    Deliberately has no parallelism. Concurrent downloads are the single
    easiest way to saturate a shared link, and the speedup is not worth it
    for a job that runs in the background anyway.
    """

    def __init__(self, mbps=None, gap=DEFAULT_GAP_S, retries=DEFAULT_RETRIES,
                 cache=None, user_agent="NWP-research/1.0"):
        self.limiter = RateLimiter(mbps)
        self.gap = float(gap)
        self.retries = int(retries)
        self.cache = cache if cache is not None else DownloadCache()
        self.user_agent = user_agent
        self._last_request = 0.0

    def _wait_gap(self):
        elapsed = time.monotonic() - self._last_request
        if elapsed < self.gap:
            time.sleep(self.gap - elapsed)
        self._last_request = time.monotonic()

    def get(self, url, timeout=120, use_cache=True):
        if use_cache:
            hit = self.cache.get(url)
            if hit is not None:
                return hit

        import urllib.request
        import urllib.error

        last_err = None
        for attempt in range(self.retries):
            if attempt:
                # Exponential backoff: a server returning errors should be
                # given room, not retried harder.
                time.sleep(min(60.0, 2.0 ** attempt))
            try:
                self._wait_gap()
                req = urllib.request.Request(
                    url, headers={"User-Agent": self.user_agent})
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    data = r.read()
                self.limiter.consume(len(data))
                if use_cache:
                    self.cache.put(url, data)
                return data
            except (urllib.error.URLError, OSError) as e:
                last_err = e

        raise RuntimeError(f"failed after {self.retries} attempts: {url}\n"
                           f"  last error: {last_err}")

    def get_text(self, url, encoding="utf-8", **kw):
        return self.get(url, **kw).decode(encoding, errors="replace")

    def describe(self):
        return (f"  bandwidth cap  : {self.limiter.rate/1e6:.1f} MB/s "
                f"(NWP_MAX_MBPS to override)\n"
                f"  request gap    : {self.gap:.1f} s, sequential only\n"
                f"  cache          : {self.cache.root}\n"
                f"  moved          : {self.limiter.total_bytes/1e6:.1f} MB, "
                f"throttled {self.limiter.total_wait:.1f} s, "
                f"cache {self.cache.hits} hit / {self.cache.misses} miss")


def estimate_ingest_mb(hours, levels=20, channels=5, ny=388, nx=438):
    """
    Rough download size for an ingest run, so the cost is known in advance
    rather than discovered from a complaint.

    Herbie byte-range downloads only the matching GRIB messages, but GRIB is
    compressed and subset server-side, so the transfer is far smaller than the
    decoded array.
    """
    decoded_mb = channels * levels * ny * nx * 4 / 1e6
    grib_mb = decoded_mb * 0.25          # typical GRIB2 compression
    return {
        "per_hour_decoded_MB": decoded_mb,
        "per_hour_download_MB": grib_mb,
        "total_download_MB": grib_mb * hours,
        "minutes_at_cap": grib_mb * hours / max_mbps() / 60,
    }
