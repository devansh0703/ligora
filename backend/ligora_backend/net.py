"""Shared HTTP timeout policy for all outbound requests.

requests' plain ``timeout`` only bounds a single socket operation; a server
that trickles bytes can hold a connection open indefinitely. Every HTTP call
in Ligora therefore uses urllib3's ``Timeout(total=...)`` — a real total
deadline per request. When the deadline passes the request fails and the
caller reports the source as unavailable (honest unavailability), never a
fallback value.
"""

from __future__ import annotations

import urllib3


def http_timeout(total_seconds: int) -> urllib3.Timeout:
    """A total-deadline timeout for one HTTP request.

    connect is capped at 10 s so unreachable hosts fail fast; the total
    deadline is the authoritative cap.
    """
    total = max(1.0, float(total_seconds))
    return urllib3.Timeout(total=total, connect=min(10.0, total))
