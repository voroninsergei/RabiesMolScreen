from __future__ import annotations
from .vina import VinaEngine

class SminaEngine(VinaEngine):
    name = "smina"
    # In real life we'd call 'smina' binary; here we reuse VinaEngine logic for demo.
