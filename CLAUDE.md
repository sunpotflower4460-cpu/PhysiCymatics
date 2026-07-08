# CLAUDE.md — PhysiCymatics 作業憲法（Claude Code用）

あなたはこのリポジトリの実装担当。**物理優先サイマティクス**——全表示は実在の物理方程式の数値解であり、サンドボックス（Claude）で検証済みの物理を、iOS（Swift 6 / SwiftUI / Metal）とPython参照実装に忠実に移す。

## NO-COMPROMISE規律（違反はマージ不可）

1. **起きるなら表示、起きないなら起きていない**。文献知識・見た目の貼り付け禁止。模様の形は計算が選ぶまで出さない
2. 動力学コードの**禁止語**: 節/node目標・模様名・ターゲットパターン（測定・分析コードのみ可）
3. **語の予約**: LIVE / TIER1 / EMERGED は実測・実接続のみ。モックは常に紫MOCK表示（physics_coreフィールド駆動）
4. **432Hz等の周波数優遇禁止**
5. claim tiers必須: measured / observed / established / interpretive / frontier。「同じ数学≠同じもの」
6. モジュールヘッダー必須: PUT-IN / EMERGED / claim-tier / floors（正直な床）
7. **レンダリング正直条項**: 描画は物理状態を読むだけ。書き換え・誇張・イージング禁止。美は実光学から
8. うえきさんへの報告は**やさしい日本語**、専門用語には一行の説明、質問は一度に1つ
9. 行き詰まり・バグ・罠は削除せず `docs/working_ledger.md` に記録（罠の図鑑）

## オラクル地図（実装はこれらと一致しなければ緑にならない）

| oracle | 検証済みの事実 | 実装の合格線 |
|---|---|---|
| `oracles/c001_plate/` | 円板固有値：独立2法一致 **1.58e-10**、節トポロジー、Chladni則p=2.16 | Swift PlateModal の固有周波数が c001_result.json と相対 <1e-6 |
| `oracles/c002_sand/` | 節集積の創発、Γ窓(0.8凍結/15/40/80)、形成10-20s、e/β頑健 | 砂層はΓゲート＋前状態保持（履歴）。輸送遅れ則 ∝f^-3.5 |
| `oracles/c003_faraday/` | 閾値：2ルート一致 **0.0013%**、Liouville<1e-12、λ=f/2分散に4桁一致、水温頑健 | FaradayLinear が sweep/ext JSON と <0.1% |
| `oracles/c005_shape/` | ZV方程式PDE：線形オンセット=c003一致、**440Hz→正方格子獲得**(90.0°,64%)、40Hz→非正方(73.8°,θ_r=61°抑圧可視,未収束) | Tier2.5表示は毛管域=正方のみ。混合域=「未確定」表示 |
| `oracles/c006_caustics/` | 光線追跡：近軸交差深度 **0.000%** vs nR/(n-1)、負の球面収差測定済 | Metalコースティクスは profiles と照合 |

**再現性の定義**: 参照実装（倍精度）=ビット厳密／GPU=機種内決定的・機種間は許容誤差内／審判=Pythonオラクル。

## フェーズ計画（受入条件つき）

- **Phase 0 配置**: この種一式＋既存UI（Geminiコード→`ui_shell/`）を整理しコミット。CI（`pytest -m "not slow"`）緑。
- **Phase 1 参照実装**: `reference/`にPythonオラクルの整理版（plate_modal / faraday_linear / caustics）。ピン留めテスト全緑。c005/KT完全粘性は`-m slow`でCI夜間走行の骨組み。
- **Phase 2 Swiftコア**: `PhysicsCore`パッケージ（pack読込・PlateModal・FaradayLinear）。**golden tests**: Swift出力↔Pythonオラクル JSON照合。
- **Phase 3 アプリ**: `docs/PhysiCymatics_Gemini配線指示_v03final.md` 完全準拠でSwiftUI移植＋Metal（c006照合）。板初期152Hz、水初期60Hz。
- **Phase 4 CI重課題**: c005 seed束・f走査・混合域収束／Kumar-Tuckerman完全粘性化／**原典照合4件**（Leissa・底境界層係数・Chen-Viñals PRE 60,559 位相図・歩行液滴閾値）→ 台帳のYELLOWをGREEN化。

## 記録ルーム（ledger/）

JSONL、schema `physicymatics.ledger.v0`。actor ∈ {human-app, human-hand, claude-sandbox, claude-code-ci, gemini-shell}。
手動実験は Issue テンプレ「実験報告」→ あなたが ledger へ転記（tier=observed, actor=human-hand）。
CIの重走行結果も必ず ledger に1行残す。

## いま既知の罠（先に読め）

`docs/working_ledger.md` の「罠の図鑑」参照。特に：k格子は共鳴舌(Δk/k~0.5%)をエイリアスする／argmax系の焦点推定は偏る／符号の鏡像はパラキシャル照合が検出器／背景プロセスは環境により回収される。
