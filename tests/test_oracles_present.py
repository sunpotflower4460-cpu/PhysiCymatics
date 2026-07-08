"""Every oracle referenced by CLAUDE.md must exist and carry its result JSON."""
import os
ROOT = os.path.dirname(os.path.dirname(__file__))
def _has(*parts): assert os.path.exists(os.path.join(ROOT, *parts)), "/".join(parts)
def test_c001(): _has("oracles","c001_plate","c001.py"); _has("oracles","c001_plate","c001_result.json")
def test_c002(): _has("oracles","c002_sand","c002b.py")
def test_c003(): _has("oracles","c003_faraday","c003.py"); _has("oracles","c003_faraday","c003_sweep.json")
def test_c005(): _has("oracles","c005_shape","c005.py")
def test_c006(): _has("oracles","c006_caustics","c006_caustic_oracle.json")
def test_pack(): _has("physics_pack","PhysiCymatics_physics_pack_v1.json")
