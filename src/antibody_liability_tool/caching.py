"""TTL-based disk caching utilities built on *diskcache*."""

from __future__ import annotations

import functools
import hashlib
import logging
from pathlib import Path
from typing import Any, Callable, TypeVar

import diskcache

logger = logging.getLogger(__name__)

F = TypeVar("F", bound=Callable[..., Any])

_DEFAULT_DIR = Path(".cache/antibody_tool")
_DEFAULT_TTL = 604800  # 7 days in seconds

# Module-level singleton (lazily initialised).
_cache_instance: diskcache.Cache | None = None


def get_cache(
    directory: str | Path | None = None,
    ttl: int | None = None,
) -> diskcache.Cache:
    """Return the module-level :class:`diskcache.Cache`, creating it if needed.

    Parameters
    ----------
    directory:
        Filesystem path for the cache database. Defaults to
        ``.cache/antibody_tool`` relative to cwd.
    ttl:
        Default time-to-live in seconds. Only used on first creation.
    """
    global _cache_instance  # noqa: PLW0603
    if _cache_instance is None:
        cache_dir = Path(directory) if directory else _DEFAULT_DIR
        cache_dir.mkdir(parents=True, exist_ok=True)
        _cache_instance = diskcache.Cache(str(cache_dir))
        logger.debug("Disk cache initialised at %s", cache_dir)
    return _cache_instance


def _make_key(*args: Any, **kwargs: Any) -> str:
    """Build a deterministic cache key from positional and keyword arguments."""
    raw = repr((args, sorted(kwargs.items())))
    return hashlib.sha256(raw.encode()).hexdigest()


def cached(ttl: int | None = None) -> Callable[[F], F]:
    """Decorator that caches a function's return value on disk.

    Parameters
    ----------
    ttl:
        Time-to-live in seconds. Falls back to ``_DEFAULT_TTL``.

    Usage::

        @cached(ttl=3600)
        def expensive(x: int) -> int:
            ...
    """
    expire = ttl if ttl is not None else _DEFAULT_TTL

    def decorator(func: F) -> F:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            cache = get_cache()
            key = f"{func.__module__}.{func.__qualname__}:{_make_key(*args, **kwargs)}"
            result = cache.get(key, default=diskcache.UNKNOWN)
            if result is not diskcache.UNKNOWN:
                logger.debug("Cache hit for %s", key)
                return result
            result = func(*args, **kwargs)
            cache.set(key, result, expire=expire)
            return result

        return wrapper  # type: ignore[return-value]

    return decorator


def invalidate(func: Callable[..., Any] | None = None) -> int:
    """Remove cached entries.

    Parameters
    ----------
    func:
        If given, only entries whose key starts with the function's
        qualified name are deleted. Otherwise the entire cache is cleared.

    Returns
    -------
    int
        Number of entries removed.
    """
    cache = get_cache()
    if func is None:
        count = len(cache)
        cache.clear()
        logger.info("Cleared entire cache (%d entries)", count)
        return count

    prefix = f"{func.__module__}.{func.__qualname__}:"
    removed = 0
    for key in list(cache):
        if isinstance(key, str) and key.startswith(prefix):
            cache.delete(key)
            removed += 1
    logger.info("Invalidated %d entries for %s", removed, prefix)
    return removed
