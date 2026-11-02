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
- 2026-10-29 / claude-code-ci / M0-M1 (砂指示書): reference/sand_bounce.py が c002b とビット一致（50粒0.3s）。Γ=80・500粒・24s で中央値|W| 0.415→0.144・リング2.96倍、Γ=0.8 は凍結。映像 renders/m1_gamma80（無音）、報告 docs/sand_implementation_report.md。M2以降未着手 / measured (within c002 model)
- 2026-10-29 / claude-code-ci / M2 (砂指示書): TS砂コアがM1 fixtureと1e-9m未満で一致。Worker+Sand表示をAppへ接続（再生/停止/再開/撒き直し/粒子数/Γ/作品表示）。ブラウザ実測 500粒×1.00・2000粒×0.42・5000粒×0.15（実時間超過なし）。既存Plate/Fluid回帰OK。報告 docs/sand_implementation_report.md / measured (within c002 model)
- 2026-10-29 / claude-code-ci / M3 (砂指示書): 多モード過渡応答（1300Hz迄10モード18自由度・モーダル質量・RK4+Hermite）を Python参照+TS に実装。定常振幅は解析式と2%以内、中央駆動対称性・切替連続性を確認。score から力とPCMを生成し演奏プリセット（152→649→368.6Hz）と自由周波数（吸着なし、173.3Hzで振幅低下を観察）をブラウザ確認。500粒で約×0.5。報告 docs/sand_implementation_report.md / measured (within model)
- 2026-10-29 / claude-code-ci / M4-M5 (砂指示書): 確定時刻CLIで音声入り縦長動画 renders/m4_performance（200粒・3区間・12s・A/V一致をffprobe確認）。規模ベンチ 100k粒で約1064万step/s、100k粒動画 renders/m5_scale100k（リング2.74倍・593s）まで到達。60万粒は単一モード約3.4h／多モード約100h規模と実測外挿し未達として記録（チェックポイント+高速backendが次段）。報告 docs/sand_implementation_report.md / measured (within model)
- 2026-10-29 / claude-code-ci / M5前進 (砂指示書): 表面評価をバッチ化（従来計算と1e-12未満一致）、render CLIにチェックポイント再開を追加（中断→再開が最終状態・フレームとも一致）。過渡backend大規模実測 600k粒=約34万step/s → 36秒演奏は約301時間と判明。60万粒ジョブ renders/m5_600k を開始し再開を証明（3フレーム+checkpoint）。完成動画は未達として記録 / measured (within model)
- 2026-10-29 / claude-code-ci / M6-1: advance_to厳密化、ステップ内イベント正確化（first-contact探索+残り時間処理、離脱はHermite線形性で厳密解、inelastic collapseは着座正則化として記録）。Python/TS照合 1e-6m未満・回数一致。音はAudioWorklet逐次生成へ（全長PCM廃止、位相連続の周波数変更）。実測 docs/m6_1_profile.md、設計 docs/m6_2_event_driven_design.md。Python過渡は大幅に遅いまま（審判用途） / measured (within model)
- 2026-10-29 / claude-code-ci / M6-2: 物理条件による休止を Python/TS 両方へ実装（着座は加速度下限、飛行は接触下限+縁到達、入力変更で全無効化）。on/off 一致（位置1e-6m未満・回数一致）を両言語で固定。TS実測: 強駆動で約1.3-1.5倍、弱駆動で約6倍（docs/m6_2_results.md）。空間共有は未実装 / measured (within model)
- 2026-10-30 / claude-code-ci / M6-2R: Codex反例（70Hz入力+649Hz残留で休止が離脱を見逃す）を確認し、休止判定をODEエネルギー上限から「同一RK4+Hermiteの板先読みによる厳密区間範囲」へ置換。Hermite高さ係数 0.105→4/27。c_roll>0 は着座休止を無効化。イベント上限時の強制着座を廃止（両言語とも上限後は区間分割、全測定で到達0回）。非弾性崩壊は分解能正則化（反発頂点1e-8m未満を粘着、根拠は結果文書）として分離。反例を両言語回帰テスト化。代償として強駆動の速度利益はほぼ消失（約±4%）、弱駆動で約25%短縮。on/off は演算順序差のためビット一致せず、許容値と根拠を docs/m6_2_results.md に記録 / measured (within model)
- 2026-10-30 / claude-code-ci / M6-3: TSにイベントスケジューラ（活動集合+起床ヒープ、旧走査と回数完全一致・位置差1e-9m未満）、領域単位の着座判定共有（実装Φの厳密なセル範囲+先読み、暧昧なら粒子単位へ復帰）、描画差分転送（再構成=snapshot一致をテスト固定）を実装。実測では支配項は活動粒子の表面評価のままで、改善は条件により数%〜20%程度。100,000粒弱駆動で物理1秒あたり約480秒相当、60万粒リアルタイムは現構造では未達と記録（docs/m6_3_results.md）。Python参照は無変更 / measured (within model)- 2026-10-31 / claude-code-ci / M6-4: surfAt最適化（時刻係数共有・割当除去・平坦化、参照版とビット一致、全体約35%短縮）、複数Worker並列（板複製+連続分割、統合結果は単一とビット一致、ブラウザで8Worker約8倍）、WebGL点描画（60万点で約10ms以下、Canvas2D約281ms）を実装。600,000粒弱駆動の初期0.05秒を実測（単一67.4秒、8分割最遅7.9秒）。リアルタイムには約2桁不足と記録（docs/m6_4_results.md） / measured (within model)
- 2026-10-31 / claude-code-ci / M6-4R: 並列同期をラウンド制（世代・バッチ番号・制御版・フレーム種別で照合）へ再設計し、制御は専用ラウンドで全Worker適用を確認してから次へ進む。偽Workerハーネスで順序競合7件を固定。初期化は乱数スキップによる区間生成へ変更（全生成とビット一致）。ブラウザで再生中周波数変更・停止再開・撒き直しを確認 / measured (within model)
- 2026-10-31 / claude-code-ci / M6-5: 保守的非接触判定（通過セル群×現在セグメントの厳密な板高さ上限と弾道下限の区間包含）をTSへ実装しデフォルトON化。ON/OFFで位置・速度・状態・衝突回数がビット一致。実測 5,000粒強駆動1.80倍・弱駆動1.12倍・形成後1.43倍・20,000粒強駆動1.56倍。負荷平準化はシャード偏り約3%で不採用、着座軽量キャッシュとGPUは不採用/調査のみ。開発用 `?count=N`（上限60万）経路を追加し、ブラウザで20,000粒実計算を確認。600,000粒初期0.05秒=単一57.6秒、8分割最遅6.8秒（理想並列見積もり）。リアルタイムまで約2桁差と記録（docs/m6_5_results.md） / measured (within model)
- 2026-11-02 / claude-code-ci / M6-5R: 局所非接触判定の境界監査。事前計算Φ範囲の適用域は正方形 [-a,a]^2 のみで、外側では実装Φ（最終動径ビン×実θ）を包含しないことを実測で確認（点(0.24,0.05)で16/18自由度が範囲外、板高さで最大1.76e-4m超過）。実際の状態差が出る反例は探索範囲（2,835条件+601点走査）で再現せず、ベンチ条件では外向き区間は省略されていなかったが、保証としては不成立のため安全側に修正: 正方形を出る区間は証明しない、非有限入力は失敗扱い、セル矩形とq範囲に丸め分の微小拡張、着座領域近道にも同じ適用域ガード。Hermite範囲はBernstein形式と相互包含を検証。修正後も速度・省略件数・衝突回数は維持（docs/m6_5r_boundary_audit.md） / measured (within model)
- 2026-11-02 / claude-code-ci / M6-6: 時間スライス非接触証明（2/4/8分割・適応/固定）を試作し物理ビット一致と包含を検証したが、接触探索へ戻る区間が最大21%減っても判定コスト（セル×スライス上限計算）が上回り全条件でM6-5より低速のためデフォルト不採用（localSplitsオプションとして残置、既定1）。代わりに接触サンプル・二分探索の高さ専用評価heightAt（surfAtの高さとビット一致）を採用し全条件で22〜33%短縮。ブラウザでWorker同時実行を実測（5,000粒: 1→8 Workerで物理速度×0.14→×0.45、到着ばらつきはラウンド総時間に比べ小さい）。600,000粒は初期0.02秒で単一10.99秒、リアルタイムまで約2桁差と記録（docs/m6_6_results.md）。作業中にPR #8が閉じられブランチが消える事象があり、PR参照からコミットを復旧して再オープンした / measured (within model)
- 2026-11-02 / claude-code-ci / M6-5R/#8・M6-6/#9 最終監査とマージ: CIワークフローがdisabled_inactivityで停止していたため再有効化。pytest単体実行でreferenceが見えない問題（pythonpath=.）とswift-coreのPEP 668失敗（setup-python化）を修正し、両PRのHEADでfast/webapp成功を確認して順次マージ（#8=7cc8afb、#9=1f52713）。Codexレビューは使用量上限、CodeRabbitは条件によりスキップのため、独立監査記録（docs/m6_5r_final_audit.md, docs/m6_6_final_audit.md）を代替証拠とした。計測フック__physiPoolは開発経路限定に変更 / measured (within model)
- 2026-11-02 / claude-code-ci / M6-7: 初期局面プロファイル（支配項は接触解決の高さ評価、二分探索30回/接触）。着座・起床の高さ+加速度専用評価と二分探索向け融合高さ評価を採用（いずれもビット一致）し、600,000粒初期0.02秒で15,671→11,058ms、100,000粒で9,160→6,039ms。飛行休止の局所上限化は省略量が増えても判定コストが上回り全般に低速化したため既定OFFで残置。600,000粒は物理1秒あたり約550秒相当でリアルタイムまで約2桁差（docs/m6_7_results.md） / measured (within model)
- 2026-11-02 / claude-code-ci / M6-8: ライブ自由度マスクを採用。中心駆動では18自由度中m=0の2つだけが生きるため、Hermite端点が全て0の自由度を高さ・加速度評価から除外（厳密な0項の省略でビット一致、テスト固定）。同一計測内A/Bで初期・形成後とも19〜41%短縮。浅い接触の収束研究では152Hzは既定刻みで収束、649Hzで数%、1241.6Hzで大きく刻み依存（カウントのみでは見逃し分を分離できず、連続衝突検出の参照が必要）と記録。Worker同時実行は5,000粒で1→8 Worker約2.8倍。60万粒リアルタイムには約2桁差（docs/m6_8_results.md） / measured (within model)
