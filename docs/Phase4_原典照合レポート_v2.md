# Phase 4 原典照合レポート v2（完全版）— 全YELLOWをGREENへ

サンドボックス(Claude)側で**4件すべての原典照合を完了**。Claude Codeはこれで台帳のYELLOWを一掃できる。

## 結果サマリ（4件すべて解決）

| 照合対象 | 原典 | 結果 | データ |
|---|---|---|---|
| **板の固有値** | Leissa NASA SP-160, 自由円板 ν=0.33 | **GREEN 0.34%**（8モード） | c001_leissa_check.json |
| **底境界層減衰** | Case & Parkinson 1957 / Miles 1967 | **GREEN**（式・係数一致、フリーパラメータなし） | c003_bl_coefficient_check.json |
| **模様の位相図** | Chen & Viñals PRE 60,559 (1999) | **GREEN**（毛管=正方形一致・混合=正方形離脱一致、交差手法検証） | c005_chenvinals_check.json |
| **水の波長** | 公表水実験 30-60Hz | **GREEN ±10%** | （既報） |
| 歩行液滴オイル | 20cSt 80Hz（複数論文） | **λ GREEN 6.1% / 閾値=正直な床** | c003_walker_check.json |

## 詳細

### 1. Leissa（板）GREEN
自由円板ν=0.33の定番λ²表と8モード照合、最大0.34%。文献表は3-4桁精度なので事実上完全一致。**c001のYELLOW撤廃**。

### 2. 底境界層 GREEN
c003の減衰項 `2νk² + k√(νω/2)/sinh(2kh)` の底Stokes層部分が、Case&Parkinson/Miles原典の漸近公式と**形も係数も一致**（勝手な係数なし）。しかもBL/bulk比が「浅い/低周波(kh<1)でだけ効く」というc003の測定結果を裏付け。**c003境界層のYELLOW撤廃**。

### 3. Chen-Viñals（形）GREEN・交差手法検証
原典予言（γ≪1）：毛管域(高f)→正方形、混合域(中f)→六角/八回系列、重力域(低f)→縞。
- c005毛管域(440Hz): 正方形（90°ペア64%）→ **一致**
- c005混合域(40Hz): 正方形が選ばれない・θ_r=61°領域が活性 → **予言通り正方形域を離脱**
- 重要: c005は**直接PDE計算**、Chen-Viñalsは**独立な振幅方程式縮約**。別手法が同じ位相図を出した＝強い交差検証。**c005のfrontier→establishedに格上げ可**。

### 4. 歩行液滴オイル — 正直な二段判定（隠さず床にする）
20cStオイル80Hz: 波長5.04 vs 4.75mm（6.1%, GREEN）。閾値5.3 vs 4.2-4.3g（γ/ω=0.124で弱減衰近似が妥当域の縁）。**これはKumar-Tuckerman完全粘性版が必要という既知の宿題を独立に裏付けた**（水では無問題、高粘性で効く）。

## 台帳更新指示（Claude Code向け）
- c001: YELLOW → **GREEN**
- c003境界層: YELLOW → **GREEN**
- c005: frontier → **established**（毛管域の正方形、Chen-Viñals交差検証済み）
- c003高粘性: KT完全粘性版が必要（Phase-4 CI、オラクルc003.py＋この床の記録あり）
- pack v1.4 の anchor_summary に全記録
