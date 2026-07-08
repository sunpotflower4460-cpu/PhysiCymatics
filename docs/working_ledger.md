# Working Ledger — 作業台帳（人もAIも足あとを残す）

形式: 日付 / actor / 出来事 / claim-tier。行き詰まり・罠は**消さずに残す**（同じ穴に二度落ちないため）。

## 罠の図鑑（先に読むと得をする）

| # | 罠 | 症状 | 検出器 | 対処 |
|---|---|---|---|---|
| T1 | **k格子エイリアシング** | 水の共鳴舌はΔk/k≈0.5%の剃刀。対数格子(4.8%間隔)が芯を跳び越え、閾値が5.7倍過大に化ける | 独立ルート2（弱減衰解析式）との不一致 | 3段階ズーム(c003 sweep_f) |
| T2 | **焦点推定の偏り** | コースティクスの焦点深度をargmax(最大輝度)で測ると+5〜15%ずれる（カスプ超の畳み込みが芯より明るく見える） | 近軸極限で偏差が増える | 近軸光線の軸交差深度をu→0外挿(c006) |
| T3 | **符号の鏡像** | 光線の偏向符号を逆にすると山/谷が入れ替わる（見た目それらしいので気づけない） | 近軸照合の悪化 | 山側集中の直接チェック(c006) |
| T4 | **固定Γの誤解** | 「Γ一定なら模様が動く」と思うと高周波で凍結して見える。実はA∝1/f²で蹴りが消える | 変位で見ると閾値はほぼ平ら | 遅れ則 f^-3.5 を明示(c004) |
| T5 | **背景プロセス回収** | nohup/& のジョブが環境により静かに消える | ログ不在 | チェックポイント+CLIルートで再開可能に(c005) |
| T6 | **共鳴の谷間** | 板の既定440Hzは368.6と565.5の谷間。実は「静寂が正しい」 | — | 板初期値を152Hzへ(本物のリングが即見える) |

## 履歴（seed）

- 2026-07-02 / claude-sandbox / c001完了。円板固有値2ルート一致1.58e-10、Chladni則p=2.16、中央駆動は4/32モードのみ / measured
- 2026-07-02 / claude-sandbox / c002完了。「集まれ」不記載で節集積が創発、Γ窓、e/β頑健。T4記録 / measured
- 2026-07-03 / claude-sandbox / c003完了。Faraday閾値2ルート0.0013%。T1記録。440Hz「ハード外」枠付けの誤りをうえきの直感が検出→変位4.3µmで実機到達可能に訂正 / measured→corrected
- 2026-07-07 / claude-sandbox / c004(遅れ則f^-3.5)・c005(440Hz正方格子獲得/40Hz非正方・θ_r=61°抑圧可視)・c006(近軸0%、T2/T3記録) / measured
- 2026-07-07 / claude-code / リポジトリ種一式配置。Phase 0開始 / (実装)
- 2026-07-08 / claude-code / Phase 1: reference/ 整理版（plate_modal / faraday_linear / caustics）をオラクルにビット一致で固定。ピン留めテスト34件緑。zhang_vinals は c005 線形オンセットの route-2 自己検証を slow（夜間CI）骨組みとして追加（成長率が審判linearモードに一致） / (実装, オラクル値に固定)
- 2026-07-08 / claude-code / Phase 2: core/PhysicsCore（SwiftPM）— Pack ローダー＋PlateModal＋FaradayLinear。reference/gen_golden.py で pack＋参照からゴールデンJSON生成、XCTestで照合（板 rel<1e-6・水<0.1%）。Python側 test_golden_parity で毎push固定。Swiftコンパイルは当環境にツールチェーン無く未検証→CI swift-core を夜間/手動に隔離（緑ゲートはpytest維持） / (実装, Swiftビルド未検証)

## 原典照合（Phase 4）— YELLOW一掃

サンドボックス側で4件の原典照合を完了。詳細は `docs/Phase4_原典照合レポート_v2.md`、
数値床は `oracles/*_check.json`、要点は pack v1.4 `external_anchors` / `anchor_summary`。

| 照合 | 原典 | 遷移 | データ |
|---|---|---|---|
| c001 板固有値 | Leissa NASA SP-160 (ν=0.33) | **YELLOW → GREEN**（最大0.34%, 8モード） | oracles/c001_plate/c001_leissa_check.json |
| c003 底境界層 | Case&Parkinson 1957 / Miles 1967 | **YELLOW → GREEN**（式・係数一致・フリーパラメータ0） | oracles/c003_faraday/c003_bl_coefficient_check.json |
| c005 形の位相図 | Chen&Viñals PRE 60,559 (1999) | **frontier → established**（直接PDE↔振幅方程式の交差検証） | oracles/c005_shape/c005_chenvinals_check.json |
| c003 高粘性(歩行液滴) | 20cSt/80Hz 文献 | **YELLOW据え置き**（λはGREEN 6.1%／閾値5.3 vs 4.2–4.3g＝正直な床。KT完全粘性版が宿題） | oracles/c003_faraday/c003_walker_check.json |

- 2026-07-08 / claude-sandbox→claude-code / Phase 4 原典照合4件を取り込み。pack v1.4（external_anchors GREEN×4, all_yellows_cleared=true）。c001/c003境界層=GREEN、c005=established。c003高粘性のみYELLOW据え置き（KTフル粘性版が必要、床を c003_walker_check.json に記録） / established / measured

## Web版（webapp/）— ブラウザ実装（Phase 3の代替・実機なし）

Gemini製UI(`ui_shell/App.jsx`)をTS移植し、モック物理を実エンジンへ置換。板=c001モーダル場、
水=c003 Faraday閾値。`npm run build`＋vitest緑、Playwrightで152Hzリング/90.5Hz静止/閾値下静寂/
440Hz正方格子を確認（スクショ検証）。砂動力学(c002)・光(c006)は未計算＝正直な床。

- 2026-07-08 / claude-code / webapp/（Vite+React+TS）新規。板152Hz単一リング・90.5Hz共鳴せず・水は閾値下静寂/f2/λ・高周波正方格子(c005)を実装。照合テスト（板↔c001 rel<1e-6、水↔c003<0.1%）5件緑。CIに webapp ジョブ追加（push緑ゲート）。Swift/実機UIは据え置き / (実装, ブラウザPlaywright検証済み)
