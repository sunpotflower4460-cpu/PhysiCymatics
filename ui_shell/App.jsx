```react
import React, { useState, useEffect, useRef, useMemo } from 'react';
import { Settings, Play, SlidersHorizontal, Share2, Info, Activity, Database, Music, FileJson, Beaker, Download, ArrowRightCircle } from 'lucide-react';

// --- CSS & Design Tokens ---
const STYLES = `
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;700&family=Noto+Sans+JP:wght@400;700&family=Space+Grotesk:wght@400;700&display=swap');

:root {
  --void: #0B0F14;
  --plate: #161D26;
  --hair: #26303B;
  --sand: #E3CE9E;
  --water: #8FD6E8;
  --signal: #F2A33C;
  --ok: #4CC38A;
  --warn: #E5484D;
  --mock: #A78BFA;
  --text: #E9EEF3;
  --text2: #93A1AE;
}

body {
  margin: 0;
  background-color: #000;
  color: var(--text);
  font-family: 'Noto Sans JP', 'Space Grotesk', sans-serif;
  -webkit-font-smoothing: antialiased;
}

.font-ui { font-family: 'Space Grotesk', 'Noto Sans JP', sans-serif; }
.font-mono { font-family: 'IBM Plex Mono', monospace; }

.mock-stripes {
  background: repeating-linear-gradient(
    45deg,
    transparent,
    transparent 5px,
    rgba(167, 139, 250, 0.1) 5px,
    rgba(167, 139, 250, 0.1) 10px
  );
}

.scroll-hide::-webkit-scrollbar { display: none; }
`;

// --- Data Contracts & Mock Provider ---
const INITIAL_CONFIG = {
  schema: "physicymatics.config.v0",
  run_id: "proto-init-123",
  media: "plate",
  drive: { type: "point", pos_xy: [0.5, 0.5], freq_hz: 440.0, accel_ms2: 12.0 },
  plate: {
    shape: "square", size_m: 0.24, thickness_m: 0.001,
    material: "steel", E_pa: 2.0e11, rho_kgm3: 7850, nu: 0.30,
    bc: "free", damping_zeta: 0.005
  },
  fluid: { depth_m: 0.005, viscosity_pas: 1.0e-3, surface_tension_nm: 0.072 },
  particles: { kind: "sand", diameter_um: 300, count: 20000 },
  seed: 12345,
  app_version: "proto-0.1", physics_core: "mock" // ここが"mock"の間はすべての表示がMOCKとなる
};

// NOTE: この見た目を物理の目標にしないこと。これはUIテスト用の単純な解析解です。
class AnalyticMockProvider {
  constructor() {
    this.config = null;
    this.field = null;
    this.dirty = true;
    this.w = 120;
    this.h = 120;
  }
  setConfig(c) {
    if (JSON.stringify(this.config) !== JSON.stringify(c)) {
      this.config = JSON.parse(JSON.stringify(c));
      this.dirty = true;
    }
  }
  eigenfrequencies() {
    // プレースホルダ: c * (m^2 + n^2) に基づく簡易的な固有振動数リスト
    const freqs = [];
    const baseF = 40; 
    for(let m=1; m<=6; m++) {
      for(let n=1; n<=6; n++) {
        freqs.push({ m, n, f_hz: Math.round(baseF * (m*m + n*n)) });
      }
    }
    return freqs.sort((a, b) => a.f_hz - b.f_hz);
  }
  frame(tMs) {
    if (this.config.media === 'fluid') {
      return { kind: "none", reason: "tier2_core_not_connected" };
    }
    
    if (this.dirty) {
      this.field = new Float32Array(this.w * this.h);
      const f = this.config.drive.freq_hz;
      const [x0, y0] = this.config.drive.pos_xy;
      const zeta = this.config.plate.damping_zeta;
      const modes = this.eigenfrequencies();

      let maxAmp = 0;
      for (let y = 0; y < this.h; y++) {
        for (let x = 0; x < this.w; x++) {
          let nx = x / this.w;
          let ny = y / this.h;
          let val = 0;
          for (let mode of modes) {
             let phi0 = Math.sin(mode.m * Math.PI * x0) * Math.sin(mode.n * Math.PI * y0);
             let fn = mode.f_hz;
             let resonance = Math.sqrt(Math.pow(fn*fn - f*f, 2) + Math.pow(2*zeta*fn*f, 2));
             let w_mn = Math.abs(phi0) / (resonance + 0.1);
             let phi = Math.sin(mode.m * Math.PI * nx) * Math.sin(mode.n * Math.PI * ny);
             val += w_mn * phi;
          }
          let absVal = Math.abs(val);
          this.field[y * this.w + x] = absVal;
          if (absVal > maxAmp) maxAmp = absVal;
        }
      }
      if (maxAmp > 0) {
        for(let i=0; i<this.field.length; i++) this.field[i] /= maxAmp;
      }
      this.dirty = false;
    }
    return { kind: "modal_field", source: "analytic", field: this.field, w: this.w, h: this.h };
  }
}

const mockProvider = new AnalyticMockProvider();

// --- Ledger Fixtures ---
const LEDGER_DATA = [
  { schema: "physicymatics.ledger.v0", run_id: "r1", ts: "2026-07-02T14:31:00+09:00", actor: "human-app", config: INITIAL_CONFIG, git_hash: "abc1234", device: "iPhone15,2", tier: "TIER1", claim_tier: "observed", artifacts: ["pattern_0001.png"] },
  { schema: "physicymatics.ledger.v0", run_id: "r2", ts: "2026-07-02T09:12:00+09:00", actor: "claude-sandbox", config: {...INITIAL_CONFIG, media: 'fluid'}, git_hash: "def5678", device: "sandbox-py3.12", tier: "MOCK", claim_tier: "interpretive", artifacts: [] },
  { schema: "physicymatics.ledger.v0", run_id: "r3", ts: "2026-07-01T22:05:00+09:00", actor: "claude-code-ci", config: INITIAL_CONFIG, git_hash: "1a2b3c4", device: "ubuntu-latest", tier: "TIER1", claim_tier: "measured", artifacts: ["test_pass.log"] }
];

// --- Common UI Components ---

// 修正事項1: FidelityBadge は config.physics_core から導出する
const FidelityBadge = ({ config }) => {
  if (config.physics_core === 'mock') {
    return (
      <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[var(--mock)]/50 bg-[var(--mock)]/20 px-2 py-1 rounded text-[10px] font-mono text-[var(--mock)] backdrop-blur-sm">
        MOCK — 物理未接続
      </div>
    );
  }
  
  if (config.media === 'fluid') {
    return (
      <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[#E5484D]/50 bg-[#E5484D]/10 px-2 py-1 rounded text-[10px] font-mono text-[#E5484D]">
        <div className="w-2 h-2 rounded-full bg-[#E5484D] animate-pulse" />
        TIER2 物理コア未接続
      </div>
    );
  }
  return (
    <div className="absolute top-16 right-4 z-10 flex items-center gap-2 border border-[#4CC38A]/50 bg-[#4CC38A]/10 px-2 py-1 rounded text-[10px] font-mono text-[#4CC38A]">
      <div className="w-2 h-2 rounded-full bg-[#4CC38A]" />
      TIER1 解析解・実時間
    </div>
  );
};

const MockWatermark = () => (
  <div className="absolute bottom-2 right-2 text-[10px] font-mono text-[var(--mock)] opacity-50 select-none z-10 pointer-events-none">
    MOCK — 物理未接続
  </div>
);

const PutInTag = () => <span className="text-[10px] bg-[var(--hair)] px-1 rounded text-[var(--text2)]">PUT-IN</span>;
const ComputedTag = () => <span className="text-[10px] border border-[var(--hair)] px-1 rounded text-[var(--ok)]">COMPUTED</span>;
const PlaceholderTag = () => <span className="text-[8px] border border-[var(--warn)] text-[var(--warn)] px-1 rounded ml-1 font-mono">PLACEHOLDER</span>;

const ActorIcon = ({ actor }) => {
  switch(actor) {
    case 'human-app': return <span className="text-[var(--signal)]">◐</span>;
    case 'claude-sandbox': return <span className="text-[var(--mock)]">◆</span>;
    case 'claude-code-ci': return <span className="text-[var(--ok)]">▣</span>;
    default: return <span className="text-[var(--text2)]">◇</span>;
  }
};

// --- Math & Logic helpers ---
const LOG_MIN_F = 20;
const LOG_MAX_F = 2000;
const getLogVal = (freq) => (Math.log(freq / LOG_MIN_F) / Math.log(LOG_MAX_F / LOG_MIN_F)) * 1000;
const getFreqFromLog = (val) => LOG_MIN_F * Math.pow(LOG_MAX_F / LOG_MIN_F, val / 1000);

// --- Screens ---

const SimulatorScreen = ({ config, setConfig, rawPanelOpen, setRawPanelOpen }) => {
  const canvasRef = useRef(null);
  const stageRef = useRef(null);
  const [modes, setModes] = useState([]);
  const [showABSheet, setShowABSheet] = useState(false);
  
  useEffect(() => {
    mockProvider.setConfig(config);
    setModes(mockProvider.eigenfrequencies());
    
    // Draw loop
    let reqId;
    const draw = () => {
      if (config.media === 'fluid') return; 
      
      const frameData = mockProvider.frame(performance.now());
      if (frameData.kind === 'modal_field' && canvasRef.current) {
        const ctx = canvasRef.current.getContext('2d');
        const { field, w, h } = frameData;
        const imgData = ctx.createImageData(w, h);
        
        for (let i = 0; i < field.length; i++) {
          const val = field[i];
          const isNode = val < 0.15; 
          const idx = i * 4;
          if (isNode) {
            imgData.data[idx] = 227;   
            imgData.data[idx+1] = 206; 
            imgData.data[idx+2] = 158; 
            imgData.data[idx+3] = 255; 
          } else {
            imgData.data[idx] = 22;
            imgData.data[idx+1] = 29;
            imgData.data[idx+2] = 38;
            imgData.data[idx+3] = 255;
          }
        }
        ctx.putImageData(imgData, 0, 0);
      }
      reqId = requestAnimationFrame(draw);
    };
    draw();
    return () => cancelAnimationFrame(reqId);
  }, [config]);

  const handlePointerDown = (e) => {
    e.target.setPointerCapture(e.pointerId);
  };
  
  const handlePointerMove = (e) => {
    if (!e.target.hasPointerCapture(e.pointerId)) return;
    const rect = stageRef.current.getBoundingClientRect();
    let x = (e.clientX - rect.left) / rect.width;
    let y = (e.clientY - rect.top) / rect.height;
    x = Math.max(0, Math.min(1, x));
    y = Math.max(0, Math.min(1, y));
    setConfig(c => ({ ...c, drive: { ...c.drive, pos_xy: [x, y] } }));
  };

  const handleFreqChange = (e) => {
    const val = parseFloat(e.target.value); 
    let freq = getFreqFromLog(val);
    
    // 修正事項6: スナップ幅は固定の ±3Hz
    const snapThreshold = 3.0; 
    const closest = modes.find(m => Math.abs(m.f_hz - freq) < snapThreshold);
    if (closest) freq = closest.f_hz;
    
    setConfig(c => ({ ...c, drive: { ...c.drive, freq_hz: freq } }));
  };

  return (
    <div className="flex flex-col h-full relative">
      {/* Top Bar */}
      <div className="absolute top-0 w-full z-20 p-4 flex justify-between items-center bg-gradient-to-b from-[var(--void)] to-transparent">
        <div className="flex bg-[var(--hair)] rounded p-1 text-sm font-ui">
          <button 
            className={`px-3 py-1 rounded ${config.media === 'plate' ? 'bg-[var(--plate)] text-white' : 'text-[var(--text2)]'}`}
            onClick={() => setConfig(c => ({...c, media: 'plate'}))}
          >
            Plate
          </button>
          <button 
            className={`px-3 py-1 rounded ${config.media === 'fluid' ? 'bg-[var(--plate)] text-white' : 'text-[var(--text2)]'}`}
            onClick={() => setConfig(c => ({...c, media: 'fluid'}))}
          >
            Fluid
          </button>
        </div>
        <div className="text-[12px] text-[var(--text2)] flex items-center gap-1 font-mono bg-[var(--void)]/50 px-2 py-1 rounded backdrop-blur">
          鋼板24cm <Settings size={12}/>
        </div>
      </div>
      
      <FidelityBadge config={config} />

      {/* Stage */}
      <div 
        ref={stageRef}
        className="flex-1 relative bg-[var(--void)] overflow-hidden flex items-center justify-center pt-16"
      >
        {config.media === 'plate' ? (
          <>
            <canvas 
              ref={canvasRef} 
              width={120} 
              height={120} 
              className="w-full h-full object-contain opacity-90 mix-blend-screen"
            />
            {/* Drive Point */}
            <div 
              className="absolute w-6 h-6 -ml-3 -mt-3 rounded-full border-2 border-[var(--signal)] flex items-center justify-center cursor-move touch-none bg-[var(--void)]/50"
              style={{ left: `${config.drive.pos_xy[0] * 100}%`, top: `${config.drive.pos_xy[1] * 100}%` }}
              onPointerDown={handlePointerDown}
              onPointerMove={handlePointerMove}
            >
              <div className="w-1.5 h-1.5 bg-[var(--signal)] rounded-full animate-pulse"/>
            </div>
            
            {/* 修正事項3: 即席表示が「節線の場」であることを注記 */}
            <div className="absolute bottom-4 w-[90%] left-[5%] text-center p-2 rounded bg-[var(--void)]/70 backdrop-blur border border-[var(--hair)]">
               <div className="text-[10px] text-[var(--text)] font-ui font-bold">【節線プレビュー（場）】</div>
               <div className="text-[9px] text-[var(--text2)] font-mono mt-1">
                 即時のモーダル場を描画中。実際の砂層の形成には秒オーダーかかります。<br/>
                 (単純支持境界の解析解 — 製品は自由端 FEM)
               </div>
            </div>
          </>
        ) : (
          <div className="w-full h-full bg-gradient-to-t from-[var(--void)] to-[var(--plate)] flex items-center justify-center relative shadow-inner">
            <div className="absolute inset-0 mock-stripes mix-blend-overlay"></div>
            
            {/* 修正事項5: レイアウト重なり解消のためトップ位置調整 */}
            <div className="absolute top-8 left-1/2 -translate-x-1/2 bg-[var(--void)]/80 backdrop-blur px-3 py-1.5 rounded-full text-[12px] font-mono flex items-center gap-2 border border-[var(--water)]/30 z-10">
              <span className="text-[var(--text2)]">駆動 {config.drive.freq_hz.toFixed(1)}Hz</span>
              <ArrowRightCircle size={12} className="text-[var(--water)]" />
              <span className="text-[var(--water)]">応答 {(config.drive.freq_hz / 2).toFixed(1)}Hz (f/2)</span>
            </div>

            <div className="absolute left-6 h-1/2 flex items-center gap-4 max-w-[65%]">
               {/* Threshold Meter (修正事項8: 上限を30に統一) */}
               <div className="h-full w-2 bg-[var(--hair)] rounded-full relative overflow-hidden shrink-0">
                  <div className="absolute bottom-0 w-full bg-[var(--text2)] opacity-30" style={{height: `${(14.0/30)*100}%`}}></div>
                  <div className="absolute bottom-[${(14.0/30)*100}%] w-full h-[2px] bg-[var(--warn)] shadow-[0_0_8px_var(--warn)]" style={{bottom: `${(14.0/30)*100}%`}}></div>
                  <div className="absolute bottom-0 w-full bg-[var(--water)] transition-all duration-200" style={{height: `${Math.min(100, (config.drive.accel_ms2 / 30) * 100)}%`}}></div>
               </div>
               <div className="flex flex-col text-[10px] font-mono gap-1 text-[var(--text2)]">
                  {/* 修正事項4: プレースホルダタグ追加 */}
                  <span className="text-[var(--warn)] flex items-center">
                    閾値 14.0 m/s² <PlaceholderTag />
                  </span>
                  <span>現在 {config.drive.accel_ms2.toFixed(1)} m/s²</span>
                  <span className="mt-2 text-[var(--water)] opacity-70 border-l border-[var(--water)] pl-2 font-ui leading-tight">
                    閾値到達で波が現れます<br/>(物理コア接続後)
                  </span>
               </div>
            </div>
            <MockWatermark />
          </div>
        )}
      </div>

      {/* Controls */}
      <div className="bg-[var(--plate)] p-6 flex flex-col gap-6 pb-24 rounded-t-3xl border-t border-[var(--hair)] relative z-30">
        
        {/* Frequency Dial (修正事項7: 対数スライダー化) */}
        <div className="flex flex-col gap-2">
          <div className="flex justify-between items-end font-mono">
            <span className="text-[12px] text-[var(--text2)]">周波数</span>
            <span className="text-[20px] text-[var(--signal)] font-bold">{config.drive.freq_hz.toFixed(1)} <span className="text-[12px]">Hz</span></span>
          </div>
          <div className="relative h-10 flex items-center">
             <input 
              type="range" 
              min="0" max="1000" step="0.1" 
              value={getLogVal(config.drive.freq_hz)}
              onChange={handleFreqChange}
              className="w-full absolute z-10 opacity-0 cursor-pointer h-full"
            />
            {/* Custom Track rendering (Log Scale) */}
            <div className="w-full h-1 bg-[var(--hair)] rounded relative">
              {modes.filter(m => m.f_hz >= LOG_MIN_F && m.f_hz <= LOG_MAX_F).map((m, i) => {
                 const pct = (Math.log(m.f_hz / LOG_MIN_F) / Math.log(LOG_MAX_F / LOG_MIN_F)) * 100;
                 const isCurrent = Math.abs(m.f_hz - config.drive.freq_hz) < 0.1;
                 return (
                   <div key={i} className={`absolute top-1/2 -translate-y-1/2 w-[2px] h-3 transition-all ${isCurrent ? 'bg-[var(--signal)] h-5 shadow-[0_0_8px_var(--signal)]' : 'bg-[var(--text2)]/40'}`} style={{left: `${pct}%`}} />
                 );
              })}
              <div className="absolute top-1/2 -translate-y-1/2 w-4 h-4 bg-[var(--signal)] rounded-full shadow-lg pointer-events-none transition-all duration-75" style={{left: `calc(${(Math.log(config.drive.freq_hz / LOG_MIN_F) / Math.log(LOG_MAX_F / LOG_MIN_F)) * 100}% - 8px)`}} />
            </div>
          </div>
        </div>

        {/* Amplitude */}
        <div className="flex flex-col gap-2">
           <div className="flex justify-between items-end font-mono">
            <span className="text-[12px] text-[var(--text2)]">振幅</span>
            <span className="text-[14px] text-[var(--text)]">{config.drive.accel_ms2.toFixed(1)} <span className="text-[10px]">m/s²</span></span>
          </div>
          <input 
              type="range" 
              min="0" max="30" step="0.5" 
              value={config.drive.accel_ms2}
              onChange={(e) => setConfig(c => ({...c, drive: {...c.drive, accel_ms2: parseFloat(e.target.value)}}))}
              className="w-full h-1 bg-[var(--hair)] rounded appearance-none cursor-pointer accent-[var(--text)]"
            />
        </div>

        {/* Buttons */}
        <div className="flex gap-4">
           <button className="flex-1 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui active:bg-[var(--hair)]">
             <Play size={16} fill="currentColor" /> 再生
           </button>
           <button className="flex-1 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui active:bg-[var(--hair)]">
             <Activity size={16} /> スイープ
           </button>
           <button onClick={() => setShowABSheet(true)} className="w-16 bg-[var(--void)] border border-[var(--hair)] rounded py-3 flex items-center justify-center gap-2 text-sm font-ui font-bold active:bg-[var(--hair)]">
             A/B
           </button>
        </div>
        
        {/* Raw Panel Handle */}
        <div 
          className="absolute -top-4 left-1/2 -translate-x-1/2 w-16 h-1.5 bg-[var(--hair)] rounded-full cursor-pointer hover:bg-[var(--text2)]"
          onClick={() => setRawPanelOpen(true)}
        />
      </div>

      {/* 修正事項9: A/B比較シートのインライン化 */}
      {showABSheet && (
        <div className="absolute inset-0 z-50 bg-[var(--void)]/90 backdrop-blur-sm flex items-center justify-center p-6 animate-in fade-in">
           <div className="bg-[var(--plate)] border border-[var(--hair)] p-6 rounded-xl text-center w-full shadow-2xl relative">
              <h3 className="text-lg font-bold font-ui mb-4">A/B 比較</h3>
              <p className="text-sm text-[var(--text)] font-ui mb-6 leading-relaxed">
                 どちらが整うかは、<br/>
                 板との相対関係で決まります。<br/>
                 板サイズを変えて確かめてください。<br/>
                 <span className="text-[var(--mock)] font-mono text-[10px] mt-4 block">※プロトタイプでは注記のみです</span>
              </p>
              <button onClick={() => setShowABSheet(false)} className="w-full py-3 bg-[var(--hair)] rounded text-sm font-ui font-bold hover:bg-[var(--text2)]/50 transition-colors">
                閉じる
              </button>
           </div>
        </div>
      )}
    </div>
  );
};

const MusicScreen = () => {
  return (
    <div className="flex flex-col h-full bg-[var(--void)] p-6 pt-12 pb-24">
       <h1 className="text-xl font-bold font-ui mb-6">ミュージック</h1>
       
       <div className="bg-[var(--plate)] p-4 rounded border border-[var(--hair)] mb-6 flex justify-between items-center">
         <div>
           <div className="text-sm">Sine Sweep (Demo)</div>
           <div className="text-[10px] text-[var(--text2)] font-mono">INTERNAL SYNTH</div>
         </div>
         <button className="p-3 bg-[var(--void)] rounded-full border border-[var(--hair)]"><Play size={16} fill="currentColor" /></button>
       </div>

       {/* Spectrogram Mock */}
       <div className="flex-1 bg-[var(--plate)] rounded border border-[var(--hair)] relative overflow-hidden flex flex-col">
          {/* 修正事項2: MOCK FFTの表示 */}
          <div className="absolute top-2 right-2 flex items-center gap-2 bg-[var(--void)] px-2 py-1 rounded text-[10px] font-mono border border-[var(--mock)] text-[var(--mock)]">
             MOCK FFT - 音声未接続
          </div>
          
          <div className="flex-1 w-full flex items-end opacity-50 px-2 pb-2 gap-[1px]">
             {/* Mock FFT Bars */}
             {Array.from({length: 40}).map((_, i) => (
                <div key={i} className="flex-1 bg-[var(--signal)] rounded-t transition-all duration-300" style={{height: `${Math.random() * 80 + 10}%`, opacity: Math.random() * 0.5 + 0.5}}></div>
             ))}
          </div>
       </div>

       <div className="mt-4 text-[10px] text-[var(--text2)] text-center font-ui leading-relaxed">
         模様は音だけでなく、板と駆動点の性質です。<br/>
         <span className="text-[var(--warn)]">模様の形成には秒オーダーかかります（音楽より遅い）。</span>
       </div>
    </div>
  );
};

const RecordsScreen = ({ setConfig, changeTab }) => {
  const loadConfig = (conf) => {
    setConfig(conf);
    changeTab('simulator');
  };

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
                 <div className="flex gap-2 items-center text-sm font-mono text-[var(--text2)]">
                   <ActorIcon actor={entry.actor} /> {entry.actor}
                 </div>
                 <div className="text-[10px] font-mono text-[var(--text2)]">{entry.ts.split('T')[1].slice(0,5)}</div>
              </div>
              
              <div className="text-sm font-ui leading-tight">
                 {entry.config.plate.material} {entry.config.plate.size_m * 100}cm · {entry.config.drive.freq_hz}Hz
              </div>
              
              <div className="flex gap-2 text-[10px] font-mono mt-1">
                <span className={`px-1.5 rounded ${entry.tier === 'TIER1' ? 'bg-[var(--ok)]/20 text-[var(--ok)]' : 'bg-[var(--mock)]/20 text-[var(--mock)]'}`}>
                  {entry.tier}
                </span>
                <span className="bg-[var(--hair)] px-1.5 rounded text-[var(--text2)]">{entry.claim_tier}</span>
                <span className="text-[var(--text2)]">seed:{entry.config.seed}</span>
              </div>

              <button 
                onClick={() => loadConfig(entry.config)}
                className="mt-2 text-[12px] font-ui border border-[var(--hair)] py-1.5 rounded text-[var(--text)] hover:bg-[var(--hair)] flex items-center justify-center gap-1"
              >
                <Share2 size={12} /> この設定を再現
              </button>
           </div>
         ))}
       </div>
    </div>
  );
};

const DisciplineScreen = () => {
  return (
    <div className="flex flex-col h-full bg-[var(--void)] p-6 pt-12 pb-24 overflow-y-auto scroll-hide font-ui">
       <h1 className="text-xl font-bold mb-6">規律</h1>
       
       <div className="prose prose-invert prose-sm text-[var(--text2)] leading-loose">
         <h3 className="text-[var(--text)] font-bold text-lg mb-2">このアプリがやること</h3>
         <p className="mb-4">ユーザーの設定した物理条件に対し、物理方程式の数値解を表示します。すべての模様は方程式の解から導かれます。</p>
         
         <h3 className="text-[var(--text)] font-bold text-lg mb-2">このアプリがやらないこと</h3>
         <ul className="list-disc pl-4 mb-4 space-y-1">
           <li>特定周波数の優遇</li>
           <li>見た目のための非物理的調整</li>
           <li>特定の効果の示唆や確約</li>
         </ul>

         <h3 className="text-[var(--text)] font-bold text-lg mb-2">物理階層（TIER）</h3>
         <ul className="mb-4 space-y-2 font-mono text-[11px]">
           <li><span className="text-[var(--ok)]">TIER1 解析解</span>: 閉形式で解ける厳密解のプレビュー</li>
           <li><span className="text-[var(--warn)]">TIER2 物理コア未接続</span>: モデルは存在するが計算リソース未接続</li>
           <li><span className="text-[var(--mock)]">MOCK</span>: プレースホルダ。この見た目を物理の目標にしてはならない。</li>
         </ul>

         <h3 className="text-[var(--text)] font-bold text-lg mb-2">再現性の約束と正直な床</h3>
         <p className="mb-4">参照実装（倍精度）はビット厳密です。GPU実装は機種内で決定的、機種間は許容誤差内とします。</p>

         <div className="text-[10px] mt-8 opacity-50 font-mono border-t border-[var(--hair)] pt-4">
           Ref: Kirchhoff–Love Plate Theory<br/>
           Benjamin & Ursell 1954<br/>
           Kumar & Tuckerman 1994
         </div>
       </div>
    </div>
  );
};

const RawPanel = ({ open, setOpen, config }) => {
  const [activeTab, setActiveTab] = useState('数値');

  const downloadConfig = () => {
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(config, null, 2));
    const dlAnchorElem = document.createElement('a');
    dlAnchorElem.setAttribute("href",     dataStr     );
    dlAnchorElem.setAttribute("download", "config.json");
    dlAnchorElem.click();
  };

  if (!open) return null;

  return (
    <div className="absolute inset-0 z-50 bg-[var(--void)]/90 backdrop-blur-sm flex flex-col justify-end animate-in fade-in duration-200">
       <div className="h-full w-full absolute" onClick={() => setOpen(false)}></div>
       <div className="bg-[var(--plate)] h-[70vh] rounded-t-2xl border-t border-[var(--hair)] p-4 flex flex-col relative">
          
          <div className="flex justify-between items-center mb-4">
            <div className="flex gap-4 border-b border-[var(--hair)] font-mono text-[12px]">
              {['モード', 'スペクトル', '数値', '閾値'].map(tab => (
                <button 
                  key={tab}
                  className={`pb-2 ${activeTab === tab ? 'text-[var(--signal)] border-b border-[var(--signal)]' : 'text-[var(--text2)]'}`}
                  onClick={() => setActiveTab(tab)}
                >{tab}</button>
              ))}
            </div>
            <button onClick={downloadConfig} className="bg-[var(--hair)] p-1.5 rounded hover:bg-[var(--text2)] text-[var(--text)] transition-colors">
              <FileJson size={16} />
            </button>
          </div>

          <div className="flex-1 overflow-y-auto scroll-hide font-mono text-[11px] text-[var(--text)]">
            <div className="flex flex-col gap-2">
              <div className="flex justify-between items-center py-1 border-b border-[var(--hair)]">
                 <span className="text-[var(--text2)]">f_drive</span>
                 <div className="flex items-center gap-2">
                   <span>{config.drive.freq_hz.toFixed(2)} <span className="text-[var(--text2)]">Hz</span></span>
                   <PutInTag />
                 </div>
              </div>
              <div className="flex justify-between items-center py-1 border-b border-[var(--hair)]">
                 <span className="text-[var(--text2)]">accel</span>
                 <div className="flex items-center gap-2">
                   <span>{config.drive.accel_ms2.toFixed(2)} <span className="text-[var(--text2)]">m/s²</span></span>
                   <PutInTag />
                 </div>
              </div>
              <div className="flex justify-between items-center py-1 border-b border-[var(--hair)]">
                 <span className="text-[var(--text2)]">resonance_dist</span>
                 <div className="flex items-center gap-2">
                   <span>0.15 <span className="text-[var(--text2)]">Hz</span></span>
                   <ComputedTag />
                 </div>
              </div>
              <div className="flex justify-between items-center py-1 border-b border-[var(--hair)]">
                 <span className="text-[var(--text2)]">physics_core</span>
                 <div className="flex items-center gap-2">
                   <span className="text-[var(--mock)]">"{config.physics_core}"</span>
                   <PutInTag />
                 </div>
              </div>
            </div>
            <div className="mt-8 text-[10px] text-[var(--text2)] opacity-50">
               // Raw JSON state dumped above.<br/>
               // 'EMERGED' tag is strictly reserved for physical core measurements.
            </div>
          </div>
       </div>
    </div>
  );
};

const Onboarding = ({ onDismiss }) => {
  const [step, setStep] = useState(0);
  const steps = [
    "これは物理が先のアプリです。すべての模様は方程式の解です。",
    "模様は周波数だけでは決まりません。板・材質・どこを駆動するかで変わります。",
    "水は駆動の半分の周波数で応答し、閾値以下では静かなままです。それが本物です。",
    "すべての実行は記録され、誰でも（AIでも）再現できます。"
  ];

  return (
    <div className="absolute inset-0 z-[100] bg-[var(--void)] flex items-center justify-center p-6">
      <div className="bg-[var(--plate)] border border-[var(--hair)] p-8 rounded-xl max-w-sm w-full flex flex-col items-center text-center gap-6 shadow-2xl">
        <Beaker size={48} className="text-[var(--signal)] mb-2" />
        <h2 className="font-ui font-bold text-lg">PhysiCymatics</h2>
        <p className="font-ui text-sm text-[var(--text)] leading-relaxed h-16 flex items-center justify-center">
          {steps[step]}
        </p>
        <div className="flex gap-2 my-2">
          {steps.map((_, i) => (
             <div key={i} className={`w-2 h-2 rounded-full ${i === step ? 'bg-[var(--signal)]' : 'bg-[var(--hair)]'}`} />
          ))}
        </div>
        <button 
          onClick={() => {
            if (step < steps.length - 1) setStep(step + 1);
            else onDismiss();
          }}
          className="w-full bg-[var(--text)] text-[var(--void)] font-bold py-3 rounded mt-2 font-ui transition-transform active:scale-95"
        >
          {step < steps.length - 1 ? '次へ' : '実験室へ'}
        </button>
      </div>
    </div>
  );
};

// --- App Root ---

export default function App() {
  const [activeTab, setActiveTab] = useState('simulator');
  const [config, setConfig] = useState(INITIAL_CONFIG);
  const [rawPanelOpen, setRawPanelOpen] = useState(false);
  const [onboarded, setOnboarded] = useState(false);

  return (
    <div className="w-full min-h-screen bg-[#000] flex justify-center">
      <style dangerouslySetInnerHTML={{__html: STYLES}} />
      
      {/* Mobile Container */}
      <div className="w-full max-w-[390px] h-[100dvh] bg-[var(--void)] relative overflow-hidden flex flex-col shadow-2xl ring-1 ring-[var(--hair)]">
        
        {!onboarded && <Onboarding onDismiss={() => setOnboarded(true)} />}

        {/* Content Area */}
        <div className="flex-1 overflow-hidden relative">
           {activeTab === 'simulator' && <SimulatorScreen config={config} setConfig={setConfig} rawPanelOpen={rawPanelOpen} setRawPanelOpen={setRawPanelOpen} />}
           {activeTab === 'music' && <MusicScreen />}
           {activeTab === 'records' && <RecordsScreen setConfig={setConfig} changeTab={setActiveTab} />}
           {activeTab === 'discipline' && <DisciplineScreen />}
        </div>

        {/* Bottom Navigation */}
        <div className="h-20 bg-[var(--void)] border-t border-[var(--hair)] flex justify-between items-center px-6 pb-4 z-40 relative">
           <button onClick={() => setActiveTab('simulator')} className={`flex flex-col items-center gap-1 transition-colors ${activeTab === 'simulator' ? 'text-[var(--text)]' : 'text-[var(--text2)]'}`}>
             <SlidersHorizontal size={20} />
             <span className="text-[10px] font-ui">シミュレータ</span>
           </button>
           <button onClick={() => setActiveTab('music')} className={`flex flex-col items-center gap-1 transition-colors ${activeTab === 'music' ? 'text-[var(--text)]' : 'text-[var(--text2)]'}`}>
             <Music size={20} />
             <span className="text-[10px] font-ui">ミュージック</span>
           </button>
           <button onClick={() => setActiveTab('records')} className={`flex flex-col items-center gap-1 transition-colors ${activeTab === 'records' ? 'text-[var(--text)]' : 'text-[var(--text2)]'}`}>
             <Database size={20} />
             <span className="text-[10px] font-ui">記録ルーム</span>
           </button>
           <button onClick={() => setActiveTab('discipline')} className={`flex flex-col items-center gap-1 transition-colors ${activeTab === 'discipline' ? 'text-[var(--text)]' : 'text-[var(--text2)]'}`}>
             <Info size={20} />
             <span className="text-[10px] font-ui">規律</span>
           </button>
        </div>

        {/* Drawers */}
        <RawPanel open={rawPanelOpen} setOpen={setRawPanelOpen} config={config} />
      </div>
    </div>
  );
}


```
