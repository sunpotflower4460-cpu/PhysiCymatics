"""Phase-2 fast guard: the Swift golden JSON must stay equal to what the pinned
Python reference + physics_pack produce. The Swift PhysicsCore tests compare
against these same files, so this keeps both sides honest even before a Swift
toolchain runs in CI. Regenerate with `python reference/gen_golden.py`."""
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(ROOT, "reference"))
import gen_golden as gg  # noqa: E402

RES = os.path.join(ROOT, "core", "PhysicsCore", "Tests",
                   "PhysicsCoreTests", "Resources")


def _load(name):
    with open(os.path.join(RES, name)) as f:
        return json.load(f)


def test_golden_files_exist():
    for name in ("golden_plate.json", "golden_faraday.json", "physics_pack.json"):
        assert os.path.exists(os.path.join(RES, name)), name


def test_plate_golden_matches_reference():
    committed = _load("golden_plate.json")
    fresh = gg.gen_plate()
    assert committed["eigenfrequencies_hz"] == fresh["eigenfrequencies_hz"]
    assert committed["center_drive_reachable_hz"] == fresh["center_drive_reachable_hz"]
    assert committed["n_modes"] == fresh["n_modes"] == 30
    for key in ("at_152hz", "at_300hz"):
        for a, b in zip(committed["forced_response"][key],
                        fresh["forced_response"][key]):
            assert abs(a - b) < 1e-12


def test_faraday_golden_matches_reference():
    committed = _load("golden_faraday.json")
    fresh = gg.gen_faraday()
    for key in ("kc", "Gc", "a_c_ms2", "lam_mm", "lam_disp_mm"):
        assert abs(committed[key] - fresh[key]) / abs(fresh[key]) < 1e-9
    assert committed["subharmonic"] == fresh["subharmonic"] is True
    assert committed["liouville_maxdev"] < 1e-12


def test_pack_resource_is_the_committed_pack():
    # The Swift test bundle must load the SAME pack the repo ships.
    bundle = _load("physics_pack.json")
    canonical = json.load(open(os.path.join(
        ROOT, "physics_pack", "PhysiCymatics_physics_pack_v1.json")))
    assert bundle["schema"] == canonical["schema"]
    assert bundle["plate_modal_shapes"]["modes"] == canonical["plate_modal_shapes"]["modes"]
