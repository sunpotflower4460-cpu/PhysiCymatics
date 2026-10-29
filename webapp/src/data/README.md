# webapp/src/data — pack由来の生成物（手編集禁止）

`scripts/gen_pack_subset.py` が `physics_pack/` と `oracles/` から生成する。
再生成: リポジトリ直下で `python3 webapp/scripts/gen_pack_subset.py`。

- `physics_pack.trimmed.json` — packの必要部分（plate_modal_shapes / plate_detents /
  fluid_water.threshold_curve* / pattern_shape）。schema=physicymatics.physics_pack.v1_4。
- `c005_snapshot.json` — c005_result.json の形スナップショット（正方格子）。
- `truth/*.json` — 照合テストの真値（c001固有値・c003閾値）。
- `sand_mode01.json` — c002単一モード4096点符号付き表（`scripts/gen_sand_web_data.py`）。
- `truth/sand_m1_short.json` — M1 fixtureのWeb照合用コピー。
- `plate_modes_transient.json` — M3過渡用モード基底（`scripts/gen_transient_web_data.py`）。
- `truth/transient_m6_parity.json` — M6-1 Python/Web照合用（Python plate_transient 出力）。
