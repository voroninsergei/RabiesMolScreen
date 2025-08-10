from loguru import logger as _logger

def get_logger(name: str):
    """Return a preconfigured Loguru logger bound to module name."""
    # reset default sinks to avoid duplicates when used as a library
    try:
        _# logger.remove()  # do not remove external sinks
    except Exception:
        pass
    _logger.add(lambda m: print(m, end=""), level="INFO")
    return _logger.bind(mod=name)
