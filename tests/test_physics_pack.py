"""Phase-0 pinned tests: the physics pack must keep its verified numbers.
These guard against silent drift. They read only committed JSON (fast)."""
import json, os, math, pytest
ROOT = os.path.dirname(os.path.dirname(__file__))
PACK = json.load(open(os.path.join(ROOT, "physics_pack", "PhysiCymatics_physics_pack_v1.json")))

def test_pack_schema_version():
    assert PACK["schema"].startswith("physicymatics.physics_pack.v1")

def test_water_constants_established():
    c = PACK["fluid_water"]["constants"]
    assert abs(c["rho_kgm3"] - 998.0) < 1.0
    assert abs(c["sigma_nm"] - 0.0728) < 1e-3
    assert abs(c["nu_m2s"] - 1.004e-6) < 1e-8

def test_faraday_response_is_subharmonic():
    assert "f/2" in PACK["fluid_water"]["response"]

def test_plate_modes_present_and_sorted():
    modes = PACK["plate_modal_shapes"]["modes"]
    assert len(modes) == 30
    fs = [m["f_hz"] for m in modes]
    assert fs == sorted(fs)
    # center-drive reachable rings must include the 152 Hz prediction
    reach = PACK["plate_detents"]["circular_free_steel_24cm"]["center_drive_reachable_hz"]
    assert any(abs(f - 152.0) < 1.0 for f in reach)

def test_modal_shape_grid_length():
    m0 = PACK["plate_modal_shapes"]["modes"][0]
    assert len(m0["R"]) == 129  # 129-point radial table

def test_sand_gamma_window_monotone():
    w = PACK["sand_on_plate"]["gamma_window"]
    assert w["frozen_at"] < w["barely"] < w["clear"] < w["strong"]

def test_lag_law_recorded():
    lag = PACK["sand_on_plate"]["lag_law"]
    assert "f^-3.5" in lag["transport_speed_scaling"] or "f^-3" in lag["transport_speed_scaling"]

def test_shape_capillary_square_earned():
    cap = PACK["pattern_shape"]["capillary_regime"]
    assert "SQUARE" in cap["result"].upper()

def test_zeta_is_flagged_placeholder():
    z = PACK["plate_modal_shapes"]["zeta"]
    assert "PLACEHOLDER" in z["status"].upper()

def test_temperature_robustness_tracks_viscosity():
    t = PACK["fluid_water"]["robustness_temperature"]
    # threshold ratio should track viscosity ratio within a few percent
    assert abs(t["a_c_ratio"]["10C"] - t["nu_ratio"]["10C"]) < 0.05
