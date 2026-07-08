"""Phase-4 pinned tests: freeze the literature-anchor results (pack v1.4).

These guard the four原典照合 (external anchors) so a future edit cannot silently
drop a GREEN back to YELLOW. Fast (reads committed JSON only)."""
import json
import os

ROOT = os.path.dirname(os.path.dirname(__file__))
PACK = json.load(open(os.path.join(
    ROOT, "physics_pack", "PhysiCymatics_physics_pack_v1.json")))


def test_pack_schema_is_v1_4():
    assert PACK["schema"] == "physicymatics.physics_pack.v1_4"


def test_leissa_plate_green():
    a = PACK["external_anchors"]["leissa_plate"]
    assert a["status"] == "GREEN"
    assert a["max_dev_pct"] <= 0.5           # measured 0.34%


def test_bottom_boundary_layer_green():
    assert PACK["external_anchors"]["bottom_boundary_layer"]["status"] == "GREEN"


def test_chen_vinals_phase_diagram_green():
    assert PACK["external_anchors"]["chen_vinals_phase_diagram"]["status"] == "GREEN"


def test_water_wavelength_green():
    assert PACK["external_anchors"]["water_wavelength"]["status"] == "GREEN"


def test_all_yellows_cleared():
    assert PACK["anchor_summary"]["all_yellows_cleared"] is True


def test_capillary_regime_matches_chen_vinals():
    lm = PACK["pattern_shape"]["capillary_regime"]["literature_match"]
    assert "Chen-Vinals" in lm


def test_walker_oil_threshold_is_documented_floor():
    # The one HONEST floor left: lambda GREEN, threshold recorded (not passed),
    # needs Kumar-Tuckerman full-viscous (still YELLOW by design).
    a = PACK["external_anchors"]["walker_oil"]
    assert "GREEN" in a["status"]            # lambda side is GREEN
    w = json.load(open(os.path.join(
        ROOT, "oracles", "c003_faraday", "c003_walker_check.json")))
    assert w["wavelength"]["dev_pct"] <= 10.0
    assert w["threshold"]["verdict"] == "documented floor"


def test_anchor_oracle_files_present():
    for path in (
        ("oracles", "c001_plate", "c001_leissa_check.json"),
        ("oracles", "c003_faraday", "c003_bl_coefficient_check.json"),
        ("oracles", "c003_faraday", "c003_walker_check.json"),
        ("oracles", "c005_shape", "c005_chenvinals_check.json"),
    ):
        assert os.path.exists(os.path.join(ROOT, *path)), "/".join(path)
