"""PhysiCymatics Python reference implementations (Phase 1).

Tidy, test-pinned refactors of the oracles in ../oracles/. Import the physics
modules directly, e.g. `from reference import plate_modal, faraday_linear`.
"""
from . import plate_modal, faraday_linear, caustics, zhang_vinals  # noqa: F401
