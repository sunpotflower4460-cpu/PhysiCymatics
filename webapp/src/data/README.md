# webapp/src/data — pack由来の生成物（手編集禁止）

`scripts/gen_pack_subset.py` が `physics_pack/` と `oracles/` から生成する。
再生成: リポジトリ直下で `python3 webapp/scripts/gen_pack_subset.py`。

- `physics_pack.trimmed.json` — packの必要部分（plate_modal_shapes / plate_detents /
  fluid_water.threshold_curve* / pattern_shape）。schema=physicymatics.physics_pack.v1_4。
- `c005_snapshot.json` — c005_result.json の形スナップショット（正方格子）。
- `truth/*.json` — 照合テストの真値（c001固有値・c003閾値）。
