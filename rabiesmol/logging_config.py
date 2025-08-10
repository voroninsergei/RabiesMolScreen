"""Centralised logging configuration for RabiesMol.

This module provides a helper to obtain a ``loguru`` logger bound to a
given module name.  It resets any default sinks to avoid duplicate
logging when RabiesMol is used as a library, and attaches a simple sink
that prints messages to stdout at the INFO level.
"""

from __future__ import annotations
# Prefer loguru for rich logging; fall back to the built-in logging module if
# loguru is unavailable.  This allows the package to operate in minimal
# environments where extra dependencies are not installed.
try:
    from loguru import logger as _logger  # type: ignore
except ImportError:  # pragma: no cover - fallback used only if loguru missing
    import logging
    _logger = logging.getLogger(__name__)  # type: ignore

def get_logger(name: str):
    """Return a preconfigured logger bound to ``name``.

    When loguru is available, this function returns a loguru logger with a
    single sink printing at INFO level and bound with the ``mod`` name.  If
    loguru is not installed, it returns a standard library logger for the
    given module name without additional configuration.

    Parameters
    ----------
    name:
        Name of the module requesting a logger.

    Returns
    -------
    Any
        A configured logger instance (either ``loguru.Logger`` or
        ``logging.Logger``).
    """
    # If using loguru, configure sinks and binding
    if hasattr(_logger, "bind"):
        # Remove existing sinks to prevent duplicate output when reused
        try:
            _logger.remove()
        except Exception:
            pass
        _logger.add(lambda m: print(m, end=""), level="INFO")
        return _logger.bind(mod=name)
    else:
        # Fallback: return a standard library logger without extra sinks
        import logging
        return logging.getLogger(name)