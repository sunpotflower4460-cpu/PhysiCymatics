# reference/ — Python参照実装（Phase 1）

オラクル（`../oracles/`）の探索的コードを、テスト可能な整理版にした場所。
すべて `../tests/` のピン留めテストで数値を固定している。各モジュール先頭に
PUT-IN / EMERGED / claim-tier / floors（正直な床）のヘッダーを付す。

| ファイル | 出どころ | 合格線（tests/） | 走行 |
|---|---|---|---|
| `plate_modal.py` | c001 円板固有値＋強制応答 | 2ルート一致 <1e-6（実測1.58e-10）、中央駆動152Hz到達・90.5Hz不可、Chladni p∈[2.0,2.3] | fast |
| `faraday_linear.py` | c003 閾値/波長/f2 | kc・a_c(≈1.1 m/s²)・subharmonic・Liouville<1e-12・2ルート<0.1%・λ=disp(f/2) | fast |
| `caustics.py` | c006 光線追跡 | 近軸交差深度=nR/(n-1) <0.01%、負の球面収差（周辺光線が深く交差） | fast |
| `zhang_vinals.py` | c005 ZV方程式（線形オンセット骨組み） | 成長率が駆動で単調増加＋オラクルlinearモードにピン一致 | **slow（夜間）** |

## 使い方（例）

```python
from reference import plate_modal, faraday_linear, caustics
modes = plate_modal.mode_table(nu=0.30)
plate_modal.center_drive_reachable_hz(modes)[0]      # → 152.0 Hz
faraday_linear.sweep_f(60.0, 0.020, "bulk")["a_c_ms2"]  # → ≈1.105 m/s²
caustics.paraxial_deviation_pct(20e-6)               # → ≈0.000 %
```

## 分担（CLAUDE.md Phase 1↔4）

- 線形・閾値・光線追跡＝ここ（fast、Phase 1で固定）。
- c005 の**非線形フラッグシップ**（どの角度が生き残るか＝正方格子の獲得）と
  Kumar-Tuckerman 完全粘性化は重いので `../oracles/c005_shape/c005.py` に温存。
  `zhang_vinals.py` は**線形オンセットの route-2 自己検証**だけを担う骨組みで、
  夜間 slow ジョブに載せてある（Phase 4 で本走行へ育てる）。

## 罠（`../docs/working_ledger.md` の図鑑と対応）

- T1 k格子エイリアシング → `faraday_linear.sweep_f` の3段階ズームで回避。
- T2/T3 焦点の偏り・符号の鏡像 → `caustics` は空気→水（下向き）・山側集中で、
  焦点は**頂点(z=A)基準の近軸交差深度をu→0外挿**。argmax輝度は使わない。
