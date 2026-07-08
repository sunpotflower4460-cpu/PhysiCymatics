import React, { useState, useEffect, useRef } from 'react'
import {
  Settings, Play, SlidersHorizontal, Share2, Info, Activity, Database,
  Music, FileJson, Beaker, ArrowRightCircle,
} from 'lucide-react'
import type { Config } from './types'
import * as plate from './physics/plate'
import * as faraday from './physics/faraday'
import { shapeAt } from './physics/shape'
import { ToneEngine } from './physics/audio'

// ============================================================================
// PhysiCymatics Web版 — ported from ui_shell/App.jsx.
// The mock plate provider is REMOVED; plate/fluid now read the real physics
// engines (src/physics/*), which read the committed physics_pack + oracles.
// Rendering-honesty (CLAUDE.md §7): canvases read computed fields only.
// ============================================================================

const FREQ_PLATE_INIT = 152.0 // first center-drive ring (P1)
const FREQ_FLUID_INIT = 60.0  // first hand-experiment frequency (README)
const AMP_MAX = 50            // m/s² slider ceiling (covers 440 Hz threshold 32.8)

const INITIAL_CONFIG: Config = {
  schema: 'physicymatics.config.v0',
  run_id: 'web-init',
  media: 'plate',
  drive: { type: 'point', pos_xy: [0.5, 0.5], freq_hz: FREQ_PLATE_INIT, accel_ms2: 12.0 },
  plate: {
    shape: 'circle', size_m: 0.24, thickness_m: 0.001, material: 'steel',
    E_pa: 2.0e11, rho_kgm3: 7850, nu: 0.3, bc: 'free',
    damping_zeta: plate.zeta().value, // pack PLACEHOLDER 0.002
  },
  fluid: { depth_m: 0.02, viscosity_pas: 1.0e-3, surface_tension_nm: 0.0728 },
  particles: { kind: 'sand', diameter_um: 300, count: 20000 },
  seed: 12345,
  app_version: 'web-0.1',
  physics_core: 'webapp-ts', // real engines connected (not "mock")
}

const LEDGER_DATA = [
  { actor: 'human-app', ts: '2026-07-02T14:31:00+09:00', material: 'steel', size_cm: 24, freq: 152, tier: 'TIER1', claim: 'observed', seed: 12345, media: 'plate' as const, freqInit: FREQ_PLATE_INIT },
  { actor: 'claude-sandbox', ts: '2026-07-02T09:12:00+09:00', material: 'water', size_cm: 24, freq: 60, tier: 'ORACLE', claim: 'measured', seed: 7, media: 'fluid' as const, freqInit: FREQ_FLUID_INIT },
  { actor: 'claude-code-ci', ts: '2026-07-01T22:05:00+09:00', material: 'steel', size_cm: 24, freq: 152, tier: 'TIER1', claim: 'measured', seed: 12345, media: 'plate' as const, freqInit: FREQ_PLATE_INIT },
]

// --- badges / tags ---------------------------------------------------------
const FidelityBadge = ({ config }: { config: Config }) => {
  if (config.physics_core === 'mock') {
    return (
      <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[var(--mock)]/50 bg-[var(--mock)]/20 px-2 py-1 rounded text-[10px] font-mono text-[var(--mock)] backdrop-blur-sm">
        MOCK — 物理未接続
      </div>
    )
  }
  if (config.media === 'fluid') {
    return (
      <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[var(--water)]/50 bg-[var(--water)]/10 px-2 py-1 rounded text-[10px] font-mono text-[var(--water)]">
        <div className="w-2 h-2 rounded-full bg-[var(--water)]" />
        TIER2 水・Faraday閾値（c003）
      </div>
    )
  }
  return (
    <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[var(--ok)]/50 bg-[var(--ok)]/10 px-2 py-1 rounded text-[10px] font-mono text-[var(--ok)]">
      <div className="w-2 h-2 rounded-full bg-[var(--ok)]" />
      TIER1 板・モーダル場（c001）
    </div>
  )
}

const PutInTag = () => <span className="text-[10px] bg-[var(--hair)] px-1 rounded text-[var(--text2)]">PUT-IN</span>
const ComputedTag = () => <span className="text-[10px] border border-[var(--hair)] px-1 rounded text-[var(--ok)]">COMPUTED</span>
const PlaceholderTag = () => <span className="text-[8px] border border-[var(--warn)] text-[var(--warn)] px-1 rounded ml-1 font-mono">PLACEHOLDER</span>

const ActorIcon = ({ actor }: { actor: string }) => {
  switch (actor) {
    case 'human-app': return <span className="text-[var(--signal)]">◐</span>
    case 'claude-sandbox': return <span className="text-[var(--mock)]">◆</span>
    case 'claude-code-ci': return <span className="text-[var(--ok)]">▣</span>
    default: return <span className="text-[var(--text2)]">◇</span>
  }
}

// --- log-frequency slider mapping ------------------------------------------
const LOG_MIN_F = 20
const LOG_MAX_F = 2000
const getLogVal = (f: number) => (Math.log(f / LOG_MIN_F) / Math.log(LOG_MAX_F / LOG_MIN_F)) * 1000
const getFreqFromLog = (v: number) => LOG_MIN_F * Math.pow(LOG_MAX_F / LOG_MIN_F, v / 1000)

// colors (from tokens) for canvas pixel math
const C_PLATE = [22, 29, 38]
const C_SAND = [227, 206, 158]
const C_WATER = [143, 214, 232]
const C_VOID = [11, 15, 20]

// ============================================================================
// Plate canvas — real modal field; node lines (|field|≈0) = where sand collects
// ============================================================================
function drawPlate(canvas: HTMLCanvasElement, config: Config) {
  const zeta = config.plate.damping_zeta
  const field = plate.forcedResponseField(config.drive.freq_hz, config.drive.pos_xy, zeta)
  const ref = plate.referencePeak(config.drive.pos_xy, zeta)
  // display gain: response magnitude relative to this drive point's best
  // reachable resonance. Off-resonance drives (e.g. 90.5 Hz at center) → ~0.
  const gain = ref.peak > 0 ? Math.min(1, field.maxAbs / ref.peak) : 0
  const { data, w, h } = field
  const ctx = canvas.getContext('2d')!
  const img = ctx.createImageData(w, h)
  const nodeEps = 0.12
  for (let i = 0; i < data.length; i++) {
    const idx = i * 4
    const v = data[i]
    if (Number.isNaN(v)) { // outside the disk
      img.data[idx] = C_VOID[0]; img.data[idx + 1] = C_VOID[1]
      img.data[idx + 2] = C_VOID[2]; img.data[idx + 3] = 255
      continue
    }
    // sand accumulates on node lines; brightness scales with response strength
    const nodeStrength = v < nodeEps ? (1 - v / nodeEps) * gain : 0
    for (let c = 0; c < 3; c++) {
      img.data[idx + c] = Math.round(C_PLATE[c] + (C_SAND[c] - C_PLATE[c]) * nodeStrength)
    }
    img.data[idx + 3] = 255
  }
  ctx.putImageData(img, 0, 0)
  return gain
}

// ============================================================================
// Fluid canvas — subharmonic (f/2) standing-wave field above threshold
// ============================================================================
function drawFluid(canvas: HTMLCanvasElement, config: Config, tSec: number) {
  const { freq_hz, accel_ms2 } = config.drive
  const depth = config.fluid.depth_m
  const f = faraday.faradayField(freq_hz, depth, accel_ms2, tSec, 48, 160, 160, config.seed)
  const ctx = canvas.getContext('2d')!
  const img = ctx.createImageData(f.w, f.h)
  for (let i = 0; i < f.data.length; i++) {
    const idx = i * 4
    if (!f.active) { // honest silence below threshold: still, dark water
      img.data[idx] = 12; img.data[idx + 1] = 20; img.data[idx + 2] = 28; img.data[idx + 3] = 255
      continue
    }
    const v = (f.data[i] + 1) / 2 // -1..1 → 0..1
    for (let c = 0; c < 3; c++) {
      img.data[idx + c] = Math.round(12 + (C_WATER[c] - 12) * v)
    }
    img.data[idx + 3] = 255
  }
  ctx.putImageData(img, 0, 0)
}

// ============================================================================
// Simulator screen
// ============================================================================
const SimulatorScreen = ({
  config, setConfig, setRawPanelOpen, tone,
}: {
  config: Config
  setConfig: React.Dispatch<React.SetStateAction<Config>>
  setRawPanelOpen: (b: boolean) => void
  tone: ToneEngine
}) => {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const stageRef = useRef<HTMLDivElement>(null)
  const [showABSheet, setShowABSheet] = useState(false)
  const [plateGain, setPlateGain] = useState(1)

  const detents = plate.eigenfrequencies()

  // Plate: static field, redraw on config change. Fluid: animate f/2.
  useEffect(() => {
    if (!canvasRef.current) return
    if (config.media === 'plate') {
      setPlateGain(drawPlate(canvasRef.current, config))
      return
    }
    let raf = 0
    const t0 = performance.now()
    const loop = () => {
      if (canvasRef.current) drawFluid(canvasRef.current, config, (performance.now() - t0) / 1000)
      raf = requestAnimationFrame(loop)
    }
    loop()
    return () => cancelAnimationFrame(raf)
  }, [config])

  const handlePointerDown = (e: React.PointerEvent) => (e.target as Element).setPointerCapture(e.pointerId)
  const handlePointerMove = (e: React.PointerEvent) => {
    if (!(e.target as Element).hasPointerCapture(e.pointerId) || !stageRef.current) return
    const rect = stageRef.current.getBoundingClientRect()
    const x = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width))
    const y = Math.max(0, Math.min(1, (e.clientY - rect.top) / rect.height))
    setConfig((c) => ({ ...c, drive: { ...c.drive, pos_xy: [x, y] } }))
  }

  const handleFreqChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    let freq = getFreqFromLog(parseFloat(e.target.value))
    if (config.media === 'plate') {
      const closest = detents.find((f) => Math.abs(f - freq) < 3.0)
      if (closest) freq = closest
    }
    setConfig((c) => ({ ...c, drive: { ...c.drive, freq_hz: freq } }))
    if (tone.playing) tone.setFreq(freq)
  }

  const setMedia = (media: 'plate' | 'fluid') =>
    setConfig((c) => ({
      ...c, media,
      drive: { ...c.drive, freq_hz: media === 'plate' ? FREQ_PLATE_INIT : FREQ_FLUID_INIT },
    }))

  const water = faraday.faradayState(config.drive.freq_hz, config.fluid.depth_m, config.drive.accel_ms2)
  const shape = shapeAt(water.lam_mm)

  return (
    <div className="flex flex-col h-full relative">
      {/* Top Bar */}
      <div className="absolute top-0 w-full z-20 p-4 flex justify-between items-center bg-gradient-to-b from-[var(--void)] to-transparent">
        <div className="flex bg-[var(--hair)] rounded p-1 text-sm font-ui">
          <button className={`px-3 py-1 rounded ${config.media === 'plate' ? 'bg-[var(--plate)] text-white' : 'text-[var(--text2)]'}`} onClick={() => setMedia('plate')}>Plate</button>
          <button className={`px-3 py-1 rounded ${config.media === 'fluid' ? 'bg-[var(--plate)] text-white' : 'text-[var(--text2)]'}`} onClick={() => setMedia('fluid')}>Fluid</button>
        </div>
        <div className="text-[12px] text-[var(--text2)] flex items-center gap-1 font-mono bg-[var(--void)]/50 px-2 py-1 rounded backdrop-blur">
          {config.media === 'plate' ? '鋼の円板24cm' : `水 深さ${(config.fluid.depth_m * 1000).toFixed(0)}mm`} <Settings size={12} />
        </div>
      </div>

      <FidelityBadge config={config} />

      {/* Stage */}
      <div ref={stageRef} className="flex-1 relative bg-[var(--void)] overflow-hidden flex items-center justify-center pt-16">
        {config.media === 'plate' ? (
          <>
            <canvas ref={canvasRef} width={160} height={160} className="w-full h-full object-contain opacity-95" />
            <div
              className="absolute w-6 h-6 -ml-3 -mt-3 rounded-full border-2 border-[var(--signal)] flex items-center justify-center cursor-move touch-none bg-[var(--void)]/50"
              style={{ left: `${config.drive.pos_xy[0] * 100}%`, top: `${config.drive.pos_xy[1] * 100}%` }}
              onPointerDown={handlePointerDown} onPointerMove={handlePointerMove}
            >
              <div className="w-1.5 h-1.5 bg-[var(--signal)] rounded-full animate-pulse" />
            </div>
            <div className="absolute bottom-4 w-[90%] left-[5%] text-center p-2 rounded bg-[var(--void)]/70 backdrop-blur border border-[var(--hair)]">
              <div className="text-[10px] text-[var(--text)] font-ui font-bold">
                【節線プレビュー（実モーダル場）{plateGain < 0.05 ? ' — 共鳴せず（ほぼ静止）' : ''}】
              </div>
              <div className="text-[9px] text-[var(--text2)] font-mono mt-1">
                c001 円板固有モードの強制応答。砂は節線に集まる。<br />
                （振幅は形成速度に効くが形は不変＝線形。形成は秒オーダー・c002は未計算）
              </div>
            </div>
          </>
        ) : (
          <>
            <canvas ref={canvasRef} width={160} height={160} className="w-full h-full object-cover opacity-95" />
            <div className="absolute top-24 left-1/2 -translate-x-1/2 bg-[var(--void)]/80 backdrop-blur px-3 py-1.5 rounded-full text-[12px] font-mono flex items-center gap-2 border border-[var(--water)]/30 z-10 whitespace-nowrap">
              <span className="text-[var(--text2)]">駆動 {config.drive.freq_hz.toFixed(1)}Hz</span>
              <ArrowRightCircle size={12} className="text-[var(--water)]" />
              <span className="text-[var(--water)]">応答 {water.respHz.toFixed(1)}Hz (f/2)</span>
            </div>

            <div className="absolute left-6 h-1/2 flex items-center gap-4 max-w-[70%]">
              <div className="h-full w-2 bg-[var(--hair)] rounded-full relative overflow-hidden shrink-0">
                {/* measured error band */}
                <div className="absolute w-full bg-[var(--text2)] opacity-30"
                  style={{ bottom: `${(water.ratio ? (faraday.thresholdAt(config.drive.freq_hz, config.fluid.depth_m).band_ms2[0] / AMP_MAX) * 100 : 0)}%`, height: `${((faraday.thresholdAt(config.drive.freq_hz, config.fluid.depth_m).band_ms2[1] - faraday.thresholdAt(config.drive.freq_hz, config.fluid.depth_m).band_ms2[0]) / AMP_MAX) * 100}%` }} />
                {/* threshold line */}
                <div className="absolute w-full h-[2px] bg-[var(--warn)] shadow-[0_0_8px_var(--warn)]"
                  style={{ bottom: `${Math.min(100, (faraday.thresholdAt(config.drive.freq_hz, config.fluid.depth_m).a_c_ms2 / AMP_MAX) * 100)}%` }} />
                {/* current drive */}
                <div className="absolute bottom-0 w-full bg-[var(--water)] transition-all duration-200"
                  style={{ height: `${Math.min(100, (config.drive.accel_ms2 / AMP_MAX) * 100)}%` }} />
              </div>
              <div className="flex flex-col text-[10px] font-mono gap-1 text-[var(--text2)]">
                <span className="text-[var(--warn)] flex items-center">
                  閾値 {faraday.thresholdAt(config.drive.freq_hz, config.fluid.depth_m).a_c_ms2.toFixed(2)} m/s² <ComputedTag />
                </span>
                <span>現在 {config.drive.accel_ms2.toFixed(1)} m/s²</span>
                <span>λ {water.lam_mm.toFixed(2)} mm</span>
                {water.active ? (
                  water.overLinear ? (
                    <span className="text-[var(--warn)] border border-[var(--warn)] px-1 rounded mt-1">非線形域・計算範囲外（&gt;5×閾値）</span>
                  ) : (
                    <span className="text-[var(--water)] mt-1">波が出現中（振幅は正規化・線形理論は振幅を決めない）</span>
                  )
                ) : (
                  <span className="mt-1 text-[var(--water)] opacity-70 border-l border-[var(--water)] pl-2 leading-tight">閾値下は静寂（本物の静けさ）</span>
                )}
                {/* shape readout (Tier 2.5) */}
                <span className="mt-2 text-[var(--text)] leading-tight">
                  形: {shape.label}
                  <span className="block text-[var(--text2)] opacity-70">{shape.tier}</span>
                </span>
              </div>
            </div>
          </>
        )}
      </div>

      {/* Controls */}
      <div className="bg-[var(--plate)] p-6 flex flex-col gap-6 pb-24 rounded-t-3xl border-t border-[var(--hair)] relative z-30">
        <div className="flex flex-col gap-2">
          <div className="flex justify-between items-end font-mono">
            <span className="text-[12px] text-[var(--text2)]">周波数</span>
            <span className="text-[20px] text-[var(--signal)] font-bold">{config.drive.freq_hz.toFixed(1)} <span className="text-[12px]">Hz</span></span>
          </div>
          <div className="relative h-10 flex items-center">
            <input type="range" min="0" max="1000" step="0.1" value={getLogVal(config.drive.freq_hz)} onChange={handleFreqChange}
              className="w-full absolute z-10 opacity-0 cursor-pointer h-full" />
            <div className="w-full h-1 bg-[var(--hair)] rounded relative">
              {config.media === 'plate' && detents.filter((f) => f >= LOG_MIN_F && f <= LOG_MAX_F).map((f, i) => {
                const pct = (Math.log(f / LOG_MIN_F) / Math.log(LOG_MAX_F / LOG_MIN_F)) * 100
                const isCurrent = Math.abs(f - config.drive.freq_hz) < 0.1
                return <div key={i} className={`absolute top-1/2 -translate-y-1/2 w-[2px] ${isCurrent ? 'bg-[var(--signal)] h-5 shadow-[0_0_8px_var(--signal)]' : 'bg-[var(--text2)]/40 h-3'}`} style={{ left: `${pct}%` }} />
              })}
              <div className="absolute top-1/2 -translate-y-1/2 w-4 h-4 bg-[var(--signal)] rounded-full shadow-lg pointer-events-none transition-all duration-75"
                style={{ left: `calc(${getLogVal(config.drive.freq_hz) / 10}% - 8px)` }} />
            </div>
          </div>
        </div>

        <div className="flex flex-col gap-2">
          <div className="flex justify-between items-end font-mono">
            <span className="text-[12px] text-[var(--text2)]">振幅 {config.media === 'plate' ? '（形成速度／形は不変）' : ''}</span>
            <span className="text-[14px] text-[var(--text)]">{config.drive.accel_ms2.toFixed(1)} <span className="text-[10px]">m/s²</span></span>
          </div>
          <input type="range" min="0" max={AMP_MAX} step="0.5" value={config.drive.accel_ms2}
            onChange={(e) => setConfig((c) => ({ ...c, drive: { ...c.drive, accel_ms2: parseFloat(e.target.value) } }))}
            className="w-full h-1 bg-[var(--hair)] rounded appearance-none cursor-pointer accent-[var(--text)]" />
        </div>

        <div className="flex gap-4">
          <button onClick={() => { if (tone.playing) { tone.stop() } else { tone.start(config.drive.freq_hz) } }}
            className="flex-1 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui active:bg-[var(--hair)]">
            <Play size={16} fill="currentColor" /> {tone.playing ? '停止' : '再生（音）'}
          </button>
          <button className="flex-1 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui active:bg-[var(--hair)]">
            <Activity size={16} /> スイープ
          </button>
          <button onClick={() => setShowABSheet(true)} className="w-16 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center text-sm font-ui font-bold active:bg-[var(--hair)]">A/B</button>
        </div>

        <div className="absolute -top-4 left-1/2 -translate-x-1/2 w-16 h-1.5 bg-[var(--hair)] rounded-full cursor-pointer hover:bg-[var(--text2)]" onClick={() => setRawPanelOpen(true)} />
      </div>

      {showABSheet && (
        <div className="absolute inset-0 z-50 bg-[var(--void)]/90 backdrop-blur-sm flex items-center justify-center p-6">
          <div className="bg-[var(--plate)] border border-[var(--hair)] p-6 rounded-xl text-center w-full shadow-2xl">
            <h3 className="text-lg font-bold font-ui mb-4">A/B 比較</h3>
            <p className="text-sm text-[var(--text)] font-ui mb-6 leading-relaxed">
              どちらが整うかは、板との相対関係で決まります。<br />板サイズを変えて確かめてください。
              <span className="text-[var(--mock)] font-mono text-[10px] mt-4 block">※プロトタイプでは注記のみ</span>
            </p>
            <button onClick={() => setShowABSheet(false)} className="w-full py-3 bg-[var(--hair)] rounded text-sm font-ui font-bold">閉じる</button>
          </div>
        </div>
      )}
    </div>
  )
}

// --- other screens ---------------------------------------------------------
const MusicScreen = () => (
  <div className="flex flex-col h-full bg-[var(--void)] p-6 pt-12 pb-24">
    <h1 className="text-xl font-bold font-ui mb-6">ミュージック</h1>
    <div className="bg-[var(--plate)] p-4 rounded border border-[var(--hair)] mb-6 flex justify-between items-center">
      <div>
        <div className="text-sm">内部トーン（シミュレータの周波数）</div>
        <div className="text-[10px] text-[var(--text2)] font-mono">INTERNAL SYNTH</div>
      </div>
      <button className="p-3 bg-[var(--void)] rounded-full border border-[var(--hair)]"><Play size={16} fill="currentColor" /></button>
    </div>
    <div className="flex-1 bg-[var(--plate)] rounded border border-[var(--hair)] relative overflow-hidden flex flex-col">
      <div className="absolute top-2 right-2 flex items-center gap-2 bg-[var(--void)] px-2 py-1 rounded text-[10px] font-mono border border-[var(--mock)] text-[var(--mock)]">
        MOCK FFT — マイク/ファイル入力は未接続
      </div>
      <div className="flex-1 w-full flex items-end opacity-50 px-2 pb-2 gap-[1px]">
        {Array.from({ length: 40 }).map((_, i) => (
          <div key={i} className="flex-1 bg-[var(--signal)] rounded-t" style={{ height: `${(((i * 37) % 80) + 10)}%`, opacity: 0.6 }} />
        ))}
      </div>
    </div>
    <div className="mt-4 text-[10px] text-[var(--text2)] text-center font-ui leading-relaxed">
      模様は音だけでなく、板と駆動点の性質です。<br />
      <span className="text-[var(--warn)]">模様の形成には秒オーダーかかります（音より遅い）。</span>
    </div>
  </div>
)

const RecordsScreen = ({ setConfig, changeTab }: { setConfig: React.Dispatch<React.SetStateAction<Config>>; changeTab: (t: string) => void }) => {
  const load = (e: (typeof LEDGER_DATA)[number]) => {
    setConfig((c) => ({ ...c, media: e.media, drive: { ...c.drive, freq_hz: e.freqInit }, seed: e.seed }))
    changeTab('simulator')
  }
  return (
    <div className="flex flex-col h-full bg-[var(--void)] p-6 pt-12 pb-24 overflow-y-auto scroll-hide">
      <div className="flex justify-between items-center mb-6">
        <h1 className="text-xl font-bold font-ui">記録ルーム</h1>
        <div className="text-[10px] font-mono bg-[var(--plate)] px-2 py-1 rounded border border-[var(--hair)]">actor: ALL ▾</div>
      </div>
      <div className="flex flex-col gap-4">
        {LEDGER_DATA.map((entry, idx) => (
          <div key={idx} className="bg-[var(--plate)] p-4 rounded border border-[var(--hair)] flex flex-col gap-3">
            <div className="flex justify-between items-start">
              <div className="flex gap-2 items-center text-sm font-mono text-[var(--text2)]"><ActorIcon actor={entry.actor} /> {entry.actor}</div>
              <div className="text-[10px] font-mono text-[var(--text2)]">{entry.ts.split('T')[1].slice(0, 5)}</div>
            </div>
            <div className="text-sm font-ui leading-tight">{entry.material} {entry.size_cm}cm · {entry.freq}Hz</div>
            <div className="flex gap-2 text-[10px] font-mono mt-1">
              <span className={`px-1.5 rounded ${entry.tier === 'MOCK' ? 'bg-[var(--mock)]/20 text-[var(--mock)]' : 'bg-[var(--ok)]/20 text-[var(--ok)]'}`}>{entry.tier}</span>
              <span className="bg-[var(--hair)] px-1.5 rounded text-[var(--text2)]">{entry.claim}</span>
              <span className="text-[var(--text2)]">seed:{entry.seed}</span>
            </div>
            <button onClick={() => load(entry)} className="mt-2 text-[12px] font-ui border border-[var(--hair)] py-1.5 rounded text-[var(--text)] hover:bg-[var(--hair)] flex items-center justify-center gap-1"><Share2 size={12} /> この設定を再現</button>
          </div>
        ))}
      </div>
    </div>
  )
}

const DisciplineScreen = () => (
  <div className="flex flex-col h-full bg-[var(--void)] p-6 pt-12 pb-24 overflow-y-auto scroll-hide font-ui">
    <h1 className="text-xl font-bold mb-6">規律</h1>
    <div className="text-[var(--text2)] leading-loose text-sm">
      <h3 className="text-[var(--text)] font-bold text-lg mb-2">このアプリがやること</h3>
      <p className="mb-4">物理条件に対し、物理方程式の数値解を表示します。すべての模様は方程式の解です（板=c001、水=c003）。</p>
      <h3 className="text-[var(--text)] font-bold text-lg mb-2">やらないこと</h3>
      <ul className="list-disc pl-4 mb-4 space-y-1">
        <li>特定周波数の優遇</li><li>見た目のための非物理的調整</li><li>本物でない物を本物として表示</li>
      </ul>
      <h3 className="text-[var(--text)] font-bold text-lg mb-2">物理階層（TIER）</h3>
      <ul className="mb-4 space-y-2 font-mono text-[11px]">
        <li><span className="text-[var(--ok)]">TIER1 板・モーダル場</span>: c001 固有モードの強制応答（2ルート1.6e-10 / Leissa 0.34%）</li>
        <li><span className="text-[var(--water)]">TIER2 水・Faraday閾値</span>: c003 実測閾値の対数補間、f/2、閾値下は静寂</li>
        <li><span className="text-[var(--mock)]">MOCK</span>: 音声FFT等。この見た目を物理の目標にしない。</li>
      </ul>
      <h3 className="text-[var(--text)] font-bold text-lg mb-2">正直な床</h3>
      <p className="mb-4">板の減衰ζはPLACEHOLDER(0.002・未測定)。水の振幅は線形理論では決まらない（正規化表示）。砂の形成動力学(c002)・光(c006)は本Web版では未計算。</p>
      <div className="text-[10px] mt-8 opacity-50 font-mono border-t border-[var(--hair)] pt-4">
        Kirchhoff–Love Plate Theory<br />Benjamin &amp; Ursell 1954<br />Kumar &amp; Tuckerman 1994 / Chen &amp; Viñals PRE 60,559
      </div>
    </div>
  </div>
)

const RawPanel = ({ open, setOpen, config }: { open: boolean; setOpen: (b: boolean) => void; config: Config }) => {
  if (!open) return null
  const water = faraday.faradayState(config.drive.freq_hz, config.fluid.depth_m, config.drive.accel_ms2)
  const download = () => {
    const dataStr = 'data:text/json;charset=utf-8,' + encodeURIComponent(JSON.stringify(config, null, 2))
    const a = document.createElement('a'); a.setAttribute('href', dataStr); a.setAttribute('download', 'config.json'); a.click()
  }
  const Row = ({ k, v, tag }: { k: string; v: string; tag: React.ReactNode }) => (
    <div className="flex justify-between items-center py-1 border-b border-[var(--hair)]">
      <span className="text-[var(--text2)]">{k}</span>
      <div className="flex items-center gap-2"><span>{v}</span>{tag}</div>
    </div>
  )
  return (
    <div className="absolute inset-0 z-50 bg-[var(--void)]/90 backdrop-blur-sm flex flex-col justify-end">
      <div className="h-full w-full absolute" onClick={() => setOpen(false)} />
      <div className="bg-[var(--plate)] h-[70vh] rounded-t-2xl border-t border-[var(--hair)] p-4 flex flex-col relative">
        <div className="flex justify-between items-center mb-4">
          <div className="font-mono text-[12px] text-[var(--signal)] border-b border-[var(--signal)] pb-2">数値（生の状態）</div>
          <button onClick={download} className="bg-[var(--hair)] p-1.5 rounded text-[var(--text)]"><FileJson size={16} /></button>
        </div>
        <div className="flex-1 overflow-y-auto scroll-hide font-mono text-[11px] text-[var(--text)]">
          <Row k="media" v={config.media} tag={<PutInTag />} />
          <Row k="f_drive" v={`${config.drive.freq_hz.toFixed(2)} Hz`} tag={<PutInTag />} />
          <Row k="accel" v={`${config.drive.accel_ms2.toFixed(2)} m/s²`} tag={<PutInTag />} />
          {config.media === 'plate' ? (
            <Row k="zeta" v={`${config.plate.damping_zeta}`} tag={<PlaceholderTag />} />
          ) : (
            <>
              <Row k="a_c (閾値)" v={`${faraday.thresholdAt(config.drive.freq_hz, config.fluid.depth_m).a_c_ms2.toFixed(3)} m/s²`} tag={<ComputedTag />} />
              <Row k="lambda" v={`${water.lam_mm.toFixed(3)} mm`} tag={<ComputedTag />} />
              <Row k="response" v={`${water.respHz.toFixed(2)} Hz (f/2)`} tag={<ComputedTag />} />
            </>
          )}
          <Row k="physics_core" v={`"${config.physics_core}"`} tag={<PutInTag />} />
          <div className="mt-8 text-[10px] text-[var(--text2)] opacity-50">// 'EMERGED' タグは実物理コア測定に限る。</div>
        </div>
      </div>
    </div>
  )
}

const Onboarding = ({ onDismiss }: { onDismiss: () => void }) => {
  const [step, setStep] = useState(0)
  const steps = [
    'これは物理が先のアプリです。すべての模様は方程式の解です。',
    '模様は周波数だけでは決まりません。板・材質・どこを駆動するかで変わります。',
    '水は駆動の半分の周波数で応答し、閾値以下では静かなままです。それが本物です。',
    'すべての実行は記録され、誰でも（AIでも）再現できます。',
  ]
  return (
    <div className="absolute inset-0 z-[100] bg-[var(--void)] flex items-center justify-center p-6">
      <div className="bg-[var(--plate)] border border-[var(--hair)] p-8 rounded-xl max-w-sm w-full flex flex-col items-center text-center gap-6 shadow-2xl">
        <Beaker size={48} className="text-[var(--signal)] mb-2" />
        <h2 className="font-ui font-bold text-lg">PhysiCymatics</h2>
        <p className="font-ui text-sm text-[var(--text)] leading-relaxed h-16 flex items-center justify-center">{steps[step]}</p>
        <div className="flex gap-2 my-2">{steps.map((_, i) => <div key={i} className={`w-2 h-2 rounded-full ${i === step ? 'bg-[var(--signal)]' : 'bg-[var(--hair)]'}`} />)}</div>
        <button onClick={() => (step < steps.length - 1 ? setStep(step + 1) : onDismiss())}
          className="w-full bg-[var(--text)] text-[var(--void)] font-bold py-3 rounded mt-2 font-ui active:scale-95">
          {step < steps.length - 1 ? '次へ' : '実験室へ'}
        </button>
      </div>
    </div>
  )
}

// --- App root --------------------------------------------------------------
export default function App() {
  const [activeTab, setActiveTab] = useState('simulator')
  const [config, setConfig] = useState<Config>(INITIAL_CONFIG)
  const [rawPanelOpen, setRawPanelOpen] = useState(false)
  const [onboarded, setOnboarded] = useState(false)
  const toneRef = useRef<ToneEngine>(new ToneEngine())

  useEffect(() => () => toneRef.current.stop(), [])

  return (
    <div className="w-full min-h-screen bg-[#000] flex justify-center">
      <div className="w-full max-w-[390px] h-[100dvh] bg-[var(--void)] relative overflow-hidden flex flex-col shadow-2xl ring-1 ring-[var(--hair)]">
        {!onboarded && <Onboarding onDismiss={() => setOnboarded(true)} />}
        <div className="flex-1 overflow-hidden relative">
          {activeTab === 'simulator' && <SimulatorScreen config={config} setConfig={setConfig} setRawPanelOpen={setRawPanelOpen} tone={toneRef.current} />}
          {activeTab === 'music' && <MusicScreen />}
          {activeTab === 'records' && <RecordsScreen setConfig={setConfig} changeTab={setActiveTab} />}
          {activeTab === 'discipline' && <DisciplineScreen />}
        </div>
        <div className="h-20 bg-[var(--void)] border-t border-[var(--hair)] flex justify-between items-center px-6 pb-4 z-40 relative">
          {[
            { id: 'simulator', icon: SlidersHorizontal, label: 'シミュレータ' },
            { id: 'music', icon: Music, label: 'ミュージック' },
            { id: 'records', icon: Database, label: '記録ルーム' },
            { id: 'discipline', icon: Info, label: '規律' },
          ].map(({ id, icon: Icon, label }) => (
            <button key={id} onClick={() => setActiveTab(id)}
              className={`flex flex-col items-center gap-1 transition-colors ${activeTab === id ? 'text-[var(--text)]' : 'text-[var(--text2)]'}`}>
              <Icon size={20} /><span className="text-[10px] font-ui">{label}</span>
            </button>
          ))}
        </div>
        <RawPanel open={rawPanelOpen} setOpen={setRawPanelOpen} config={config} />
      </div>
    </div>
  )
}
