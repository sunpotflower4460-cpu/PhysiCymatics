# webapp — PhysiCymatics Web版（Vite + React + TypeScript）

ブラウザで動く実装。**板=c001**（円板モーダル場）と **水=c003**（Faraday閾値）の
本物の物理を、`physics_pack` と `oracles` の検証済み数値から計算して描く。
土台UIは `ui_shell/App.jsx`（Gemini製）をTSに移植し、モック物理を実エンジンへ置換。

## 起動
```bash
cd webapp && npm install && npm run dev
# → http://localhost:5173/  （板は初期152Hz、水は初期60Hz）
# 板: 中央のオレンジ点=駆動点をドラッグ。152Hz→単一リング／90.5Hz→共鳴せず（ほぼ静止）
# 水: Fluidタブ→振幅を上げると閾値で波が出現、f/2応答、閾値下は静寂
```
`npm run build`（tsc + vite）／`npm test`（vitest, 照合テスト）。

## 構成
- `src/physics/pack.ts` — 物理packの型付きローダー（読むだけ）
- `src/physics/plate.ts` — TIER1 板：強制モーダル応答（節線=砂／c001の式）
- `src/physics/faraday.ts` — TIER2 水：閾値の対数補間＋f/2定在波場（c003）
- `src/physics/shape.ts` — Tier2.5 形：毛管域=正方格子(c005)／混合域=未確定
- `src/physics/audio.ts` — 内部駆動トーン（Web Audio、実合成／物理ではない）
- `src/data/*` — pack由来の生成物（`scripts/gen_pack_subset.py`で再生成、手編集禁止）
- `src/physics/__tests__/*` — 板固有値↔c001(rel<1e-6)、水閾値↔c003(<0.1%)

## 正直な床（未実装・MOCK）
- 板の減衰ζは pack のPLACEHOLDER **0.002**（未測定）。
- 水の**振幅は線形理論では決まらない** → 表示は正規化。閾値の5倍超は「非線形域・計算範囲外」。
- 砂の形成動力学(c002)・光/コースティクス(c006)は**本Web版では未計算**（スコープ外）。
- 音声のマイク/ファイルFFTは未接続（ミュージック画面に紫のMOCK表示）。
- 形の混合域は c005 のCI収束待ち（frontier、未確定表示）。
