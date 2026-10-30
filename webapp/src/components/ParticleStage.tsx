import { useEffect, useRef, useState } from 'react'
import { Pause, Play, RotateCcw } from 'lucide-react'
import { drawParticles, drawFrameCache, FrameCache } from '../rendering/particleRenderer'
import type { SandSnapshot } from '../simulation/types'
import type { TransientStats } from '../simulation/transientSand'
import { calibratedForceN } from '../simulation/transientSand'
import type { PerformanceScore } from '../simulation/sequence'
import { ScoreStreamAudio } from '../audio/scoreWorklet'

// Sand stage: single-mode (M2) and performance / free-frequency (M3) modes.
// Physics runs in a Worker; this component only sends controls and draws snapshots.
const PLATE_A = 0.12
const COUNTS = [500, 2000, 5000]
const GAMMAS = [80, 40, 0.8]
const DRIVE_CENTER: [number, number] = [0, 0]
const DRIVE_OFF: [number, number] = [0.06, 0]

interface AnySnap { t: number; xy: Float32Array; aliveCount: number; initialCount: number; stats: never }

function performanceScore(forceN: number): PerformanceScore {
  const ev = (start_s: number, frequency_hz: number) => ({
    start_s, end_s: start_s + 12, frequency_hz, forcePeak_N: forceN, attack_s: 0.1, release_s: 0.1,
  })
  return {
    schema: 'physicymatics.score.v1', duration_s: 36, initialPhase_rad: 0, audioGain: 0.15,
    events: [ev(0, 152.0), ev(12, 649.0), ev(24, 368.6)],
  }
}

export default function ParticleStage({ seed0 }: { seed0: number }) {
  const canvasRef = useRef<HTMLCanvasElement>(null)
  const workerRef = useRef<Worker | null>(null)
  const audioRef = useRef<ScoreStreamAudio | null>(null)
  const frameCacheRef = useRef(new FrameCache())
  const [mode, setMode] = useState<'single' | 'perf'>('single')
  const [playing, setPlaying] = useState(false)
  const [count, setCount] = useState(500)
  const [gamma, setGamma] = useState(80)
  const [seed, setSeed] = useState(seed0)
  const [freq, setFreq] = useState(152)
  const [driveOff, setDriveOff] = useState(false)
  const [artwork, setArtwork] = useState(false)
  const [snap, setSnap] = useState<AnySnap | null>(null)
  const [kind, setKind] = useState<'single' | 'transient'>('single')
  const [error, setError] = useState<string | null>(null)

  useEffect(() => {
    const wk = new Worker(new URL('../simulation/worker.ts', import.meta.url), { type: 'module' })
    workerRef.current = wk
    wk.onmessage = (ev: MessageEvent) => {
      if (ev.data.type === 'snapshot') { setKind(ev.data.kind); setSnap(ev.data.snap) }
      else if (ev.data.type === 'frame') {
        // Differential frame (transient): update the display cache, then
        // publish only the small header (time/counts/stats) as state.
        frameCacheRef.current.apply(ev.data.frame)
        const f = ev.data.frame
        setKind('transient')
        setSnap({ t: f.t, xy: frameCacheRef.current.xy, aliveCount: f.aliveCount, initialCount: f.initialCount, stats: f.stats as never })
      }
      else if (ev.data.type === 'error') { setError(String(ev.data.message)); setPlaying(false) }
    }
    wk.onerror = () => { setError('Workerでエラーが起きました'); setPlaying(false) }
    wk.postMessage({ type: 'init', gamma: 80, count: 500, seed: seed0 })
    return () => { wk.terminate(); workerRef.current = null }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [])

  useEffect(() => {
    if (snap && canvasRef.current) {
      if (kind === 'transient') drawFrameCache(canvasRef.current, frameCacheRef.current, PLATE_A)
      else drawParticles(canvasRef.current, snap as unknown as SandSnapshot, PLATE_A)
    }
  }, [snap, kind])

  // Keep audio from running ahead of lagging physics (score time = physics t here).
  useEffect(() => {
    if (!playing || mode !== 'perf') return
    const id = setInterval(() => {
      const au = audioRef.current
      if (!au || !snap) return
      const lag = au.position() - snap.t
      if (lag > 0.75) void au.suspendForLag()
      else if (lag < 0.25) void au.resumeFromLag()
    }, 250)
    return () => clearInterval(id)
  }, [playing, mode, snap])

  const post = (m: unknown) => workerRef.current?.postMessage(m)
  const initSingle = (g: number, c: number, s: number) => {
    setMode('single'); setGamma(g); setCount(c); setSeed(s); setPlaying(false); setError(null)
    void audioRef.current?.pause()
    post({ type: 'init', gamma: g, count: c, seed: s })
  }
  const initPerf = (score: PerformanceScore, c: number, s: number, off: boolean) => {
    setMode('perf'); setCount(c); setSeed(s); setPlaying(false); setError(null)
    if (!audioRef.current) audioRef.current = new ScoreStreamAudio()
    void audioRef.current.load(score)
    post({ type: 'initPerformance', score, count: c, seed: s, driveXY: off ? DRIVE_OFF : DRIVE_CENTER })
  }
  const togglePlay = () => {
    if (playing) {
      post({ type: 'pause' }); void audioRef.current?.pause(); setPlaying(false)
    } else {
      post({ type: 'play' })
      if (mode === 'perf') void audioRef.current?.play()
      setPlaying(true)
    }
  }
  const applyFreq = () => {
    const f = Math.min(2000, Math.max(20, freq))
    setFreq(f)
    post({ type: 'setFrequency', hz: f })
    audioRef.current?.setFrequency(f)
  }
  const applyDrive = (off: boolean) => {
    setDriveOff(off)
    post({ type: 'setDrive', driveXY: off ? DRIVE_OFF : DRIVE_CENTER })
  }

  const ts = kind === 'transient' ? (snap?.stats as unknown as TransientStats) : null
  const ss = kind === 'single' ? snap?.stats as unknown as { t: number; alive: number; stepsPerSec: number } : null
  const speed = (ts?.stepsPerSec ?? ss?.stepsPerSec ?? 0)

  return (
    <div className="flex flex-col h-full relative">
      <div className="flex-1 relative bg-[var(--void)] overflow-hidden flex items-center justify-center pt-16">
        <canvas ref={canvasRef} width={360} height={360} className="w-full h-full object-contain" />
        <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[var(--ok)]/50 bg-[var(--ok)]/10 px-2 py-1 rounded text-[10px] font-mono text-[var(--ok)]">
          {mode === 'single' ? '砂・単一モード（c002）' : '砂・多モード演奏（M3）'}{playing ? ' 実計算中' : ' 一時停止'}
        </div>
        {!artwork && snap && (
          <div className="absolute bottom-4 w-[90%] left-[5%] p-2 rounded bg-[var(--void)]/70 backdrop-blur border border-[var(--hair)] text-[10px] font-mono text-[var(--text2)]">
            物理時刻 {snap.t.toFixed(1)}s ／ 残り {snap.aliveCount}/{snap.initialCount}粒
            {speed > 0 && <> ／ 物理の進み ×{speed.toFixed(2)}{speed < 0.95 ? '（実時間より遅い）' : ''}</>}
            {ts && <> ／ 駆動 {ts.freqHz ? `${ts.freqHz.toFixed(1)}Hz` : '無音区間'} ／ 振幅尺度 {ts.ampUm.toFixed(0)}µm ／ 休止 {Math.round(ts.sleepingFrac * 100)}%</>}
            <br />{mode === 'single' ? `Γ=${gamma}・約152Hz単一モード。` : '同じscoreから音と駆動を生成。'}粒子同士の衝突は未計算。
          </div>
        )}
        {error && (
          <div className="absolute top-28 left-[5%] w-[90%] p-2 rounded border border-[var(--warn)] text-[var(--warn)] text-[11px] font-ui bg-[var(--void)]/80">
            計算が止まりました: {error}（「砂を撒き直す」で最初から再開できます）
          </div>
        )}
      </div>
      <div className="bg-[var(--plate)] p-6 flex flex-col gap-4 pb-24 rounded-t-3xl border-t border-[var(--hair)] relative z-30 max-h-[46%] overflow-y-auto">
        <div className="flex bg-[var(--hair)] rounded p-1 text-[12px] font-ui">
          <button className={`flex-1 px-2 py-1.5 rounded ${mode === 'single' ? 'bg-[var(--plate)] text-white' : 'text-[var(--text2)]'}`}
            onClick={() => initSingle(gamma, count, seed)}>単一モード（152Hz）</button>
          <button className={`flex-1 px-2 py-1.5 rounded ${mode === 'perf' ? 'bg-[var(--plate)] text-white' : 'text-[var(--text2)]'}`}
            onClick={() => initPerf(performanceScore(calibratedForceN(driveOff ? DRIVE_OFF : DRIVE_CENTER)), count, seed, driveOff)}>演奏・自由周波数</button>
        </div>
        <div className="flex gap-4">
          <button onClick={togglePlay} className="flex-1 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui active:bg-[var(--hair)]">
            {playing ? <><Pause size={16} /> 一時停止</> : <><Play size={16} fill="currentColor" /> {snap && snap.t > 0 ? '続きから再開' : mode === 'perf' ? '再生（音と砂）' : '再生'}</>}
          </button>
          <button onClick={() => (mode === 'single'
            ? initSingle(gamma, count, seed + 1)
            : initPerf(performanceScore(calibratedForceN(driveOff ? DRIVE_OFF : DRIVE_CENTER)), count, seed + 1, driveOff))}
            className="flex-1 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui active:bg-[var(--hair)]">
            <RotateCcw size={16} /> 砂を撒き直す
          </button>
        </div>
        {mode === 'perf' && (
          <>
            <div className="flex gap-2 items-end text-[12px] font-mono text-[var(--text2)]">
              <label className="flex-1 flex flex-col gap-1">周波数（20〜2000Hz・吸着なし）
                <input type="number" min={20} max={2000} step={0.1} value={freq}
                  onChange={(e) => setFreq(Number(e.target.value))}
                  className="bg-[var(--void)] border border-[var(--hair)] rounded px-2 py-2 text-[var(--text)]" />
              </label>
              <button onClick={applyFreq} className="border border-[var(--hair)] rounded px-3 py-2 text-[var(--text)] font-ui text-[12px]">
                この周波数へ変更（状態保持）
              </button>
            </div>
            <div className="flex justify-between items-center text-[12px] font-mono text-[var(--text2)]">
              <span>共鳴の目安: 90.5 / 152 / 210 / 346 / 369 / 595 / 649 Hz</span>
              <button onClick={() => applyDrive(!driveOff)} className="border border-[var(--hair)] rounded px-3 py-1.5 text-[var(--text)] font-ui text-[12px]">
                駆動点: {driveOff ? '中心から外す' : '中央'}（押すと切替・状態保持）
              </button>
            </div>
            <div className="text-[11px] text-[var(--text2)] font-ui">
              演奏プリセットは 152→649→368.6Hz（各12秒）。形が出るかは計算結果次第です。
            </div>
          </>
        )}
        <div className="flex gap-4 text-[12px] font-mono text-[var(--text2)]">
          <label className="flex-1 flex flex-col gap-1">粒子数（変更＝撒き直し）
            <select value={count} onChange={(e) => {
              const c = Number(e.target.value)
              if (mode === 'single') initSingle(gamma, c, seed)
              else initPerf(performanceScore(calibratedForceN(driveOff ? DRIVE_OFF : DRIVE_CENTER)), c, seed, driveOff)
            }} className="bg-[var(--void)] border border-[var(--hair)] rounded px-2 py-2 text-[var(--text)]">
              {COUNTS.map((c) => <option key={c} value={c}>{c.toLocaleString()}粒</option>)}
            </select>
          </label>
          {mode === 'single' && (
            <label className="flex-1 flex flex-col gap-1">駆動の強さ Γ（変更＝撒き直し）
              <select value={gamma} onChange={(e) => initSingle(Number(e.target.value), count, seed)}
                className="bg-[var(--void)] border border-[var(--hair)] rounded px-2 py-2 text-[var(--text)]">
                {GAMMAS.map((g) => <option key={g} value={g}>Γ={g}{g === 0.8 ? '（凍結の対照）' : ''}</option>)}
              </select>
            </label>
          )}
        </div>
        <div className="flex justify-between items-center text-[11px] font-mono text-[var(--text2)]">
          <span>seed: {seed}</span>
          <button onClick={() => setArtwork(!artwork)} className="border border-[var(--hair)] rounded px-3 py-1.5 text-[var(--text)] font-ui">
            {artwork ? '数値を表示' : '作品表示（数値を隠す）'}
          </button>
        </div>
      </div>
    </div>
  )
}
