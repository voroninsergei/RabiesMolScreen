from __future__ import annotations
from typing import Callable, Dict, Optional

_TOOLS: Dict[str, Callable[[], str]] = {}

def register_tool(name: str):
    def _wrap(fn: Callable[[], str]):
        _TOOLS[name] = fn
        return fn
    return _wrap

def resolve(name: str) -> Optional[str]:
    fn = _TOOLS.get(name)
    if not fn:
        return None
    try:
        return fn()
    except Exception:
        return None

def available_tools() -> Dict[str, bool]:
    return {k: (resolve(k) is not None) for k in _TOOLS}
