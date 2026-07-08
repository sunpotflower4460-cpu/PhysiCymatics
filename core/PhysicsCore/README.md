# PhysicsCore — Swiftコア（Phase 2）

`physics_pack` を読み込み、Python 参照実装（`../../reference/`）とオラクルを
**許容誤差内で再現**する Swift パッケージ。アプリ本体（Phase 3）はこの上に載る。

## 中身

| ファイル | 役割 | claim-tier |
|---|---|---|
| `Sources/PhysicsCore/Pack.swift` | physics_pack の型付きローダー（読むだけ） | established |
| `Sources/PhysicsCore/PlateModal.swift` | 板：固有周波数・中央駆動の到達モード・強制応答場 | measured（c001モデル内） |
| `Sources/PhysicsCore/FaradayLinear.swift` | 水：Faraday 閾値（monodromy＋2ルート照合、3段ズーム） | measured（線形/減衰モデル内） |

## ゴールデンテスト（Swift 出力 ↔ Python JSON）

`Tests/PhysicsCoreTests/` が、`reference/gen_golden.py` の生成した JSON と
Swift の出力を突き合わせる（CLAUDE.md Phase 2 の合格線）：

- 板：固有周波数30本・中央到達4本・強制応答場（152Hz/300Hz、半径5点）を rel<1e-6。
- 水：60Hz/20mm の閾値 kc・a_c(≈1.1 m/s²)・λ・subharmonic を **<0.1%**、
  Liouville<1e-12、2ルート一致<0.1%。

```bash
python reference/gen_golden.py        # ゴールデン再生成（pack＋参照が真実）
cd core/PhysicsCore && swift test      # Swift 側の照合
```

## 正直条項（重要・未検証事項）

- このパッケージは **Swift ツールチェーンの無い環境で執筆**したため、
  **コンパイル/`swift test` は著者の手元で未検証**です。数値ロジックは
  `reference/` の Python を1対1で移植しており、ズレの検出器として
  Python 側 `tests/test_golden_parity.py`（毎push・緑ゲート）がゴールデンJSONを
  参照実装に固定しています。
- CI の `swift-core` ジョブは **夜間＋手動（workflow_dispatch）** で走り、
  最初の実ツールチェーン走行がコンパイルと照合を検証します。push/PR の
  速い緑ゲート（pytest）はこれに影響されません。
- ζ（減衰比）は pack のプレースホルダ **0.002**（未測定）。手実験#1のリングダウンで測る予定。
