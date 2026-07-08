#!/usr/bin/env python3
"""Generate webapp/src/data/{physics_pack.trimmed.json, c005_snapshot.json,
truth/*} from the committed physics_pack + oracles. Run from repo root:
    python3 webapp/scripts/gen_pack_subset.py
These are pack-derived artifacts; do not hand-edit."""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "webapp", "src", "data")
pack = json.load(open(os.path.join(ROOT, "physics_pack", "PhysiCymatics_physics_pack_v1.json")))
fw = pack["fluid_water"]
sub = {
    "schema": pack["schema"],
    "_provenance": "GENERATED subset of physics_pack/PhysiCymatics_physics_pack_v1.json (pack-derived, do not hand-edit).",
    "plate_modal_shapes": pack["plate_modal_shapes"],
    "plate_detents": {"circular_free_steel_24cm": pack["plate_detents"]["circular_free_steel_24cm"]},
    "fluid_water": {k: fw[k] for k in ("constants", "response", "threshold_curve_h20mm", "threshold_curve_h5mm", "app_default_440hz")},
    "pattern_shape": pack["pattern_shape"],
}
json.dump(sub, open(os.path.join(DATA, "physics_pack.trimmed.json"), "w"), ensure_ascii=False, separators=(",", ":"))
c5 = json.load(open(os.path.join(ROOT, "oracles", "c005_shape", "c005_result.json")))
json.dump({k: c5[k] for k in ("top2_pair_angle", "top2_share", "temporal_peak_over_drive", "angle_peaks", "rms_h")},
          open(os.path.join(DATA, "c005_snapshot.json"), "w"), ensure_ascii=False, separators=(",", ":"))
os.makedirs(os.path.join(DATA, "truth"), exist_ok=True)
c1 = json.load(open(os.path.join(ROOT, "oracles", "c001_plate", "c001_result.json")))
eig = sorted(d["f_hz"] for d in c1["eigenfrequencies_hz"])
json.dump({"_provenance": "oracles/c001_plate/c001_result.json eigenfrequencies_hz (sorted)", "eigenfrequencies_hz": eig},
          open(os.path.join(DATA, "truth", "c001_eigenfrequencies.json"), "w"), ensure_ascii=False)
sweep = json.load(open(os.path.join(ROOT, "oracles", "c003_faraday", "c003_sweep.json")))
pts = sorted(({"f_hz": r["f_d"], "a_c_ms2": r["Gc"] * 9.81, "lam_mm": r["lam_mm"]}
              for r in sweep if r["h"] == 0.020 and r["variant"] == "bulk"), key=lambda x: x["f_hz"])
json.dump({"_provenance": "oracles/c003_faraday/c003_sweep.json h=0.020 bulk (a_c=Gc*9.81)", "points": pts},
          open(os.path.join(DATA, "truth", "c003_threshold_h20_bulk.json"), "w"), ensure_ascii=False)
print("regenerated webapp data subset")
