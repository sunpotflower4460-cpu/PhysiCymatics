# reference/ — Python参照実装（Phase 1で作る）

オラクル（`../oracles/`）の探索的コードを、テスト可能な整理版にする場所。
- `plate_modal.py` — c001の円板固有値＋強制応答（pack出力）
- `faraday_linear.py` — c003の閾値/波長/f2（2ルート自己検証）
- `caustics.py` — c006の光線追跡（近軸0%照合）
- `zhang_vinals.py` — c005のPDE（重い、`-m slow`）
すべて `../tests/` のピン留めテストで数値固定する。
