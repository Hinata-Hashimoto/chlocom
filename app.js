'use strict';
/* Chloroplast genome synteny viewer.
   Top/bottom genomes, each with its own center; shared zoom span.
   Homology links: "gene name" (same normalized name) or "BLAST" (precomputed
   per-gene blastn, thresholded live by % identity). Click a gene to align the
   other genome to its homolog. Zoom in to read DNA. Invert SSC, notes, etc. */

const NS = 'http://www.w3.org/2000/svg';
const CAT = {
  PSI:            ['#1b5e20', 'photosystem I (psa)'],
  PSII:           ['#4caf50', 'photosystem II (psb)'],
  cytb6f:         ['#00acc1', 'cytochrome b6f (pet)'],
  photosynthesis: ['#8bc34a', 'ATP synthase / RuBisCO'],
  ndh:            ['#c0ca33', 'NADH dehydrog. (ndh)'],
  ribosomal:      ['#42A5F5', 'ribosomal protein'],
  RNApol:         ['#AB47BC', 'RNA polymerase'],
  rRNA:           ['#EF5350', 'rRNA'],
  tRNA:           ['#FFA726', 'tRNA'],
  chl:            ['#26A69A', 'Chl biosynthesis'],
  other:          ['#BDBDBD', 'other'],
  intergenic:     ['#a1887f', 'intergenic (IGR)'],
};
const REGION_COLOR = { LSC: '#cfe2f3', SSC: '#fce5cd', IRa: '#d9ead3', IRb: '#d9ead3' };
const BASE_COLOR = { A: '#2e7d32', C: '#1565c0', G: '#ef6c00', T: '#c62828', U: '#c62828', N: '#999' };
const COMP = { A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', N: 'N' };

const L = { ml: 196, mr: 24, topAxisY: 28, topSeqY: 48, topY: 64, geneH: 28,
            botY: 300, botSeqY: 342, botAxisY: 360 };
const SVG_H = 380;             // SVG の高さ（上下の余白を詰めた）
const SEQ_MIN_PX = 7;
const shortName = d => d.replace(/\s*\(.*/, '').replace(/^(\w)\w+\s+/, '$1. ');

const svg = document.getElementById('svg');
const tooltip = document.getElementById('tooltip');
const SEQ = () => window.GENOME_SEQ || {};
const state = {
  top: null, bot: null,
  span: 1, centerTop: 0.5, centerBot: 0.5,
  hover: null, pin: null,                  // {which:'top'|'bot', key}
  onlyShared: false, showLabels: false,
  hiddenCats: new Set(['intergenic']),      // IGR は既定で非表示（凡例でON）
  sel: null,
  revcomp: false, invSSC: false,
  mode: 'name', minPid: 50,                 // homology mode + BLAST identity threshold
  notes: loadNotes(),
  _tb: new Map(), _bt: new Map(),           // link maps top->bot, bot->top (each {key,pid})
};

function loadNotes() { try { return JSON.parse(localStorage.getItem('cpv_notes') || '{}'); } catch (e) { return {}; } }
function saveNotes() { try { localStorage.setItem('cpv_notes', JSON.stringify(state.notes)); } catch (e) {} }
function noteKey(label, key) { return label + '|' + key; }
function hasNote(label, key) { const n = state.notes[noteKey(label, key)]; return n && n.trim(); }
function exportNotes() {
  if (!Object.keys(state.notes).some(k => (state.notes[k] || '').trim())) { alert('No notes to export yet.'); return; }
  download('chlocom_notes.json', JSON.stringify({ format: 'cpview-notes', version: 1, notes: state.notes }, null, 2));
}
function importNotes(file) {
  const r = new FileReader();
  r.onload = () => {
    let incoming; try { const j = JSON.parse(r.result); incoming = (j && j.notes) ? j.notes : j; } catch (e) {}
    if (!incoming || typeof incoming !== 'object') { alert('Could not read notes JSON.'); return; }
    let n = 0; for (const k in incoming) if (typeof incoming[k] === 'string') { state.notes[k] = incoming[k]; n++; }
    saveNotes(); updateSeqPanel(); render(); alert(`Imported ${n} note(s) (merged).`);
  };
  r.readAsText(file);
}

/* ---------- geometry helpers ---------- */
function el(name, attrs = {}, parent = null) {
  const e = document.createElementNS(NS, name);
  for (const k in attrs) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function plotW() { return svg.clientWidth - L.ml - L.mr; }
function centerOf(which) { return which === 'top' ? state.centerTop : state.centerBot; }
function viewOf(which) { const c = centerOf(which); return { start: c - state.span / 2, end: c + state.span / 2 }; }
function fracToX(f, which) { const v = viewOf(which); return L.ml + (f - v.start) / (v.end - v.start) * plotW(); }
function xToFrac(px, which) { const v = viewOf(which); return v.start + (px - L.ml) / plotW() * (v.end - v.start); }
function geneXs(g, len, which) {
  let x0 = fracToX(g.start / len, which), x1 = fracToX(g.end / len, which);
  if (x1 - x0 < 2.2) { const c = (x0 + x1) / 2; x0 = c - 1.1; x1 = c + 1.1; }
  return [x0, x1];
}
function clampCenters() {
  const h = state.span / 2;
  if (state.span >= 1) { state.centerTop = state.centerBot = 0.5; return; }
  state.centerTop = Math.min(1 - h, Math.max(h, state.centerTop));
  state.centerBot = Math.min(1 - h, Math.max(h, state.centerBot));
}
function whichTrack(clientY) { const y = clientY - svg.getBoundingClientRect().top; return y < (L.topY + L.botY) / 2 ? 'top' : 'bot'; }
function genomeOf(which) { return GENOME_DATA[which === 'top' ? state.top : state.bot]; }

/* ---- SSC inversion (flip-flop isomer): reverse-complement the SSC of the bottom genome ---- */
function invRegions(which) {
  if (!(state.invSSC && which === 'bot')) return [];
  return genomeOf(which).regions.filter(r => r.label === 'SSC').map(r => [r.start, r.end]);
}
function invContains(regs, p) { for (const m of regs) if (p >= m[0] && p < m[1]) return m; return null; }
function dispGene(which, gene) {
  const m = invContains(invRegions(which), (gene.start + gene.end) / 2);
  if (!m) return { start: gene.start, end: gene.end, strand: gene.strand };
  return { start: m[0] + m[1] - gene.end, end: m[0] + m[1] - gene.start, strand: -gene.strand };
}

/* ---------- homology links (name or BLAST) ---------- */
function computeLinks() {
  const top = GENOME_DATA[state.top], bot = GENOME_DATA[state.bot];
  const tb = new Map(), bt = new Map();
  if (state.mode === 'name') {
    const botSet = new Set(bot.genes.map(g => g.key));
    for (const g of top.genes) if (botSet.has(g.key)) { tb.set(g.key, { key: g.key, pid: null }); bt.set(g.key, { key: g.key, pid: null }); }
  } else {
    const hits = (window.GENOME_BLAST || {})[`${state.top}::${state.bot}`] || [];
    for (const h of hits) {
      if (h.pid < state.minPid) continue;
      const cur = tb.get(h.q);
      if (!cur || cur.pid < h.pid) tb.set(h.q, { key: h.s, pid: h.pid });
    }
    for (const [tk, lk] of tb) {           // reverse map: best per bottom gene
      const cur = bt.get(lk.key);
      if (!cur || cur.pid < lk.pid) bt.set(lk.key, { key: tk, pid: lk.pid });
    }
  }
  state._tb = tb; state._bt = bt;
}
function activeKeys() {
  const a = state.hover || state.pin;
  if (!a) return { t: null, b: null };
  if (a.which === 'top') return { t: a.key, b: (state._tb.get(a.key) || {}).key || null };
  return { b: a.key, t: (state._bt.get(a.key) || {}).key || null };
}
const visibleCat = c => !state.hiddenCats.has(c);
function geneClass(linked, isActive) { return 'gene' + (state.onlyShared && !linked ? ' dim' : '') + (isActive ? ' hl' : ''); }
function ribbonClass(isActive, anyActive) { return 'ribbon' + (anyActive ? (isActive ? ' hl' : ' dim') : ''); }
function applyHighlight() {
  const { t: aT, b: aB } = activeKeys(), anyActive = !!(state.hover || state.pin);
  for (const r of (state._genes || [])) r.el.setAttribute('class', geneClass(r.linked, r.which === 'top' ? r.gene.key === aT : r.gene.key === aB));
  for (const r of (state._ribbons || [])) r.el.setAttribute('class', ribbonClass(r.tk === aT, anyActive));
}

/* ---------- shape builders ---------- */
function arrowPoints(x0, x1, y, h, strand) {
  const w = x1 - x0, tip = Math.min(w * 0.42, 8);
  return strand >= 0
    ? `${x0},${y} ${x1 - tip},${y} ${x1},${y + h / 2} ${x1 - tip},${y + h} ${x0},${y + h}`
    : `${x1},${y} ${x0 + tip},${y} ${x0},${y + h / 2} ${x0 + tip},${y + h} ${x1},${y + h}`;
}
function ribbonPath(xt0, xt1, xb0, xb1, ytop, ybot) {
  const m = (ytop + ybot) / 2;
  return `M${xt0},${ytop} L${xt1},${ytop} C${xt1},${m} ${xb1},${m} ${xb1},${ybot} L${xb0},${ybot} C${xb0},${m} ${xt0},${m} ${xt0},${ytop} Z`;
}

/* ---------- render ---------- */
function render() {
  computeLinks();
  while (svg.firstChild) svg.removeChild(svg.firstChild);
  const W = svg.clientWidth, H = SVG_H;
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const top = GENOME_DATA[state.top], bot = GENOME_DATA[state.bot];
  if (!top || !bot) return;
  const topKeys = new Map(top.genes.map(g => [g.key, g]));
  const botKeys = new Map(bot.genes.map(g => [g.key, g]));
  const { t: aT, b: aB } = activeKeys(), anyActive = !!(state.hover || state.pin);
  state._genes = []; state._ribbons = [];

  const defs = el('defs', {}, svg);
  const clip = el('clipPath', { id: 'plot' }, defs);
  el('rect', { x: L.ml, y: 0, width: plotW(), height: H }, clip);
  const plot = el('g', { 'clip-path': 'url(#plot)' }, svg);

  // region bands
  for (const [g, y, w] of [[top, L.topY, 'top'], [bot, L.botY, 'bot']]) {
    for (const r of g.regions) {
      const x0 = fracToX(r.start / g.length, w), x1 = fracToX(r.end / g.length, w);
      el('rect', { class: 'region', x: x0, y: y - 4, width: Math.max(0, x1 - x0),
        height: L.geneH + 8, fill: REGION_COLOR[r.label] || '#eee', 'fill-opacity': .35 }, plot);
      if (x1 - x0 > 34) el('text', { class: 'region-label', x: (x0 + x1) / 2, y: y - 7, 'text-anchor': 'middle' }, plot).textContent = r.label;
    }
  }

  // ribbons from link map
  const ytop = L.topY + L.geneH, ybot = L.botY, vT = viewOf('top'), vB = viewOf('bot');
  for (const [tk, lk] of state._tb) {
    const a = topKeys.get(tk), b = botKeys.get(lk.key);
    if (!a || !b || !visibleCat(a.cat)) continue;
    const da = dispGene('top', a), db = dispGene('bot', b);
    if ((da.end / top.length < vT.start || da.start / top.length > vT.end) &&
        (db.end / bot.length < vB.start || db.start / bot.length > vB.end)) continue;
    const [at0, at1] = geneXs(da, top.length, 'top'), [bt0, bt1] = geneXs(db, bot.length, 'bot');
    const rb = el('path', { class: ribbonClass(tk === aT, anyActive), d: ribbonPath(at0, at1, bt0, bt1, ytop, ybot), fill: CAT[a.cat][0] }, plot);
    state._ribbons.push({ el: rb, tk });
  }

  // gene arrows
  for (const [g, y, which, isTop] of [[top, L.topY, 'top', true], [bot, L.botY, 'bot', false]]) {
    const len = g.length, v = viewOf(which);
    for (const gene of g.genes) {
      if (!visibleCat(gene.cat)) continue;
      const dg = dispGene(which, gene);
      if (dg.end / len < v.start || dg.start / len > v.end) continue;
      const [x0, x1] = geneXs(dg, len, which);
      const linked = which === 'top' ? state._tb.has(gene.key) : state._bt.has(gene.key);
      const isActive = which === 'top' ? gene.key === aT : gene.key === aB;
      let poly;
      if (gene.cat === 'intergenic')                     // IGR は矢印でなく低めの矩形で描く
        poly = el('rect', { class: geneClass(linked, isActive), x: Math.min(x0, x1), y: y + L.geneH * 0.28,
          width: Math.max(1, Math.abs(x1 - x0)), height: L.geneH * 0.44, fill: CAT[gene.cat][0] }, plot);
      else
        poly = el('polygon', { class: geneClass(linked, isActive),
          points: arrowPoints(x0, x1, y, L.geneH, dg.strand), fill: CAT[gene.cat][0] }, plot);
      poly.__gene = gene; poly.__which = which;
      state._genes.push({ el: poly, gene, which, linked });
      poly.addEventListener('mouseenter', ev => onGeneHover(gene, which, ev));
      poly.addEventListener('mousemove', moveTooltip);
      poly.addEventListener('mouseleave', onGeneLeave);
      if (isActive || (state.showLabels && x1 - x0 > 9))
        el('text', { class: 'gene-name' + (isActive ? ' hl' : ''), x: (x0 + x1) / 2, y: isTop ? y + L.geneH + 13 : y - 7, 'text-anchor': 'middle' }, plot).textContent = gene.name;
      if (hasNote(g.label, gene.key))
        el('text', { class: 'notemark', x: (x0 + x1) / 2, y: isTop ? y - 4 : y + L.geneH + 16, 'text-anchor': 'middle' }, plot).textContent = '✎';
    }
    drawSequence(g, which, isTop ? L.topSeqY : L.botSeqY, plot);
  }

  if (state.sel) {
    const g = genomeOf(state.sel.which), y = state.sel.which === 'top' ? L.topY : L.botY;
    const xa = fracToX(state.sel.fA, state.sel.which), xb = fracToX(state.sel.fB, state.sel.which);
    el('rect', { class: 'selrect', x: Math.min(xa, xb), y: y - 8, width: Math.abs(xb - xa), height: L.geneH + 16 }, plot);
  }

  for (const [g, y, which, isTop] of [[top, L.topY, 'top', true], [bot, L.botY, 'bot', false]]) {
    const ly = y + L.geneH / 2;
    el('text', { class: 'genome-label', x: L.ml - 12, y: ly - 2, 'text-anchor': 'end' }, svg).textContent = shortName(g.display);
    el('text', { class: 'genome-sub', x: L.ml - 12, y: ly + 12, 'text-anchor': 'end' }, svg).textContent = `${g.length.toLocaleString()} bp`;
    drawAxis(g, which, isTop ? L.topAxisY : L.botAxisY, isTop);
  }
  updateStatus();
}

function drawSequence(g, which, y, plot) {
  const seqs = SEQ()[g.label]; if (!seqs) return;
  const v = viewOf(which), len = g.length;
  if (plotW() / (state.span * len) < SEQ_MIN_PX) return;
  const regs = invRegions(which);
  const d0 = Math.max(0, Math.floor(v.start * len)), d1 = Math.min(len, Math.ceil(v.end * len));
  el('rect', { x: L.ml, y: y - 11, width: plotW(), height: 15, fill: '#fff', 'fill-opacity': .95 }, plot);
  for (let d = d0; d < d1; d++) {
    const m = invContains(regs, d), gi = m ? (m[0] + m[1] - 1 - d) : d;
    let ch = seqs[gi] || 'N';
    if (m) ch = COMP[ch] || 'N';
    if (state.revcomp) ch = COMP[ch] || 'N';
    el('text', { class: 'base', x: fracToX((d + 0.5) / len, which), y, 'text-anchor': 'middle', fill: BASE_COLOR[ch] || '#555' }, plot).textContent = ch;
  }
}
function drawAxis(g, which, y, above) {
  const v = viewOf(which), len = g.length;
  const startBp = Math.max(0, v.start) * len, endBp = Math.min(1, v.end) * len, span = endBp - startBp;
  const step = niceStep(span / 6);
  el('line', { class: 'axis', x1: L.ml, y1: y, x2: L.ml + plotW(), y2: y }, svg);
  for (let bp = Math.ceil(startBp / step) * step; bp <= endBp; bp += step) {
    const x = fracToX(bp / len, which);
    el('line', { class: 'axis', x1: x, y1: y - 3, x2: x, y2: y + 3 }, svg);
    el('text', { class: 'axis-tick', x, y: above ? y - 6 : y + 14, 'text-anchor': 'middle' }, svg)
      .textContent = span > 20000 ? (bp / 1000).toFixed(0) + 'k' : span > 2000 ? (bp / 1000).toFixed(1) + 'k' : Math.round(bp).toLocaleString();
  }
}
function niceStep(x) { const p = Math.pow(10, Math.floor(Math.log10(x))); const n = x / p; return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * p; }

/* ---------- sequence utilities ---------- */
function revComp(s) { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'N'; return o; }
function geneSeq(label, gene) { const s = (SEQ()[label] || '').slice(gene.start, gene.end); return gene.strand < 0 ? revComp(s) : s; }
function fasta(header, seq) { let o = '>' + header + '\n'; for (let i = 0; i < seq.length; i += 70) o += seq.slice(i, i + 70) + '\n'; return o; }
function download(name, text) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain' }));
  a.download = name; document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
function copyText(text, btn) {
  const done = () => { const t = btn.textContent; btn.textContent = 'Copied'; setTimeout(() => btn.textContent = t, 1200); };
  if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(done, () => fallbackCopy(text, done));
  else fallbackCopy(text, done);
}
function fallbackCopy(text, done) {
  const ta = document.createElement('textarea'); ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
  document.body.appendChild(ta); ta.select(); try { document.execCommand('copy'); done(); } catch (e) {} ta.remove();
}

/* ---------- clicked-gene info panel + cross-species presence ---------- */
// 系統順（おおまかな分類群でグループ化、ツリー風表示用）
const PHYLO = [
  ['Red lineage', ['Cyanidioschyzon_merolae', 'Phaeodactylum_tricornutum']],
  ['Green alga', ['Chlamydomonas_reinhardtii']],
  ['Liverwort', ['Marchantia_polymorpha']],
  ['Fern', ['Ceratopteris_thalictroides']],
  ['Gymnosperm', ['Cryptomeria_japonica']],
  ['Monocots', ['Oryza_sativa', 'Zea_mays', 'Sorghum_bicolor']],
  ['Eudicots', ['Arabidopsis_thaliana', 'Manihot_esculenta', 'Nicotiana_benthamiana']],
];
// 各種に「この遺伝子があるか」を BLASTN（事前計算）+ 現在のしきい値で判定
function presenceHTML(gene) {
  const src = state.pin.which === 'top' ? state.top : state.bot;   // クリック元ゲノムをクエリに
  const base = gene.key.replace(/__ir2$/, '');
  const cands = new Set([gene.key, base, base + '__ir2']);          // IR コピーも同一遺伝子とみなす
  const B = window.GENOME_BLAST || {};
  const info = lab => {
    if (lab === src) return { ok: true, self: true };
    let best = -1;
    for (const h of (B[`${src}::${lab}`] || []))
      if (cands.has(h.q) && h.pid >= state.minPid && h.pid > best) best = h.pid;
    return { ok: best >= 0, pid: best };
  };
  let n = 0, total = 0, rows = '';
  for (const [clade, labels] of PHYLO) {
    const present = labels.filter(l => GENOME_DATA[l]);
    if (!present.length) continue;
    let items = '';
    for (const lab of present) {
      total++; const r = info(lab); if (r.ok) n++;
      const cur = (lab === state.top || lab === state.bot) ? ' cur' : '';
      const pidtxt = (r.ok && !r.self) ? ` <span class="meta">${r.pid}%</span>` : '';
      items += `<span class="sp${r.ok ? ' yes' : ' no'}${cur}">${r.ok ? '●' : '○'} ${shortName(GENOME_DATA[lab].display)}${pidtxt}</span>`;
    }
    rows += `<div class="clade"><span class="cl">${clade}</span>${items}</div>`;
  }
  return `<div class="presence"><div class="pres-head">Present in ${n}/${total} species ` +
         `<span class="meta">(BLASTN ≥ ${state.minPid}% identity, query = ${shortName(GENOME_DATA[src].display)})</span></div>${rows}</div>`;
}
function updateGeneInfo() {
  const panel = document.getElementById('geneInfo');
  if (!state.pin) { panel.hidden = true; panel.innerHTML = ''; return; }
  const me = genomeOf(state.pin.which), gene = me.genes.find(g => g.key === state.pin.key);
  if (!gene) { panel.hidden = true; return; }
  const link = state.pin.which === 'top' ? state._tb.get(state.pin.key) : state._bt.get(state.pin.key);
  const hom = link ? genomeOf(state.pin.which === 'top' ? 'bot' : 'top').genes.find(g => g.key === link.key) : null;
  const T = GENOME_DATA[state.top], B = GENOME_DATA[state.bot];
  const fmt = x => x ? `${x.start.toLocaleString()}–${x.end.toLocaleString()} bp (${x.strand > 0 ? '+' : '−'})` : '— none';
  const tG = state.pin.which === 'top' ? gene : hom, bG = state.pin.which === 'top' ? hom : gene;
  const idtxt = (state.mode === 'blast' && link && link.pid != null) ? ` · BLASTN ${link.pid}% id` : '';
  panel.innerHTML =
    `<div class="gi-head"><b>${gene.name}</b> <span class="meta">[${CAT[gene.cat][1]}]${idtxt}</span></div>` +
    `<div class="gi-pos">${shortName(T.display)}: ${fmt(tG)}　／　${shortName(B.display)}: ${fmt(bG)}` +
    (hom ? '' : ' <span class="meta">— no homolog at current setting</span>') + `</div>` +
    presenceHTML(gene);
  panel.hidden = false;
}

/* ---------- sequence panel (clicked gene + homolog) ---------- */
function updateSeqPanel() {
  updateGeneInfo();
  const panel = document.getElementById('seqPanel');
  if (!state.pin) { panel.hidden = true; panel.innerHTML = ''; return; }
  const clickedG = genomeOf(state.pin.which), clicked = clickedG.genes.find(g => g.key === state.pin.key);
  if (!clicked) { panel.hidden = true; return; }
  const link = state.pin.which === 'top' ? state._tb.get(state.pin.key) : state._bt.get(state.pin.key);
  const otherWhich = state.pin.which === 'top' ? 'bot' : 'top';
  const hom = link ? genomeOf(otherWhich).genes.find(g => g.key === link.key) : null;
  const topGene = state.pin.which === 'top' ? clicked : hom;
  const botGene = state.pin.which === 'top' ? hom : clicked;
  const idtxt = (state.mode === 'blast' && link && link.pid != null) ? ` <span class="meta">— BLAST ${link.pid}% identity</span>` : '';
  const homtxt = (hom && hom.name !== clicked.name) ? ` <span class="meta">↔ ${hom.name}</span>` : '';
  panel.innerHTML = `<button class="close" id="seqClose">✕ Close</button><h3>Sequence: <i>${clicked.name}</i>${homtxt}${idtxt}</h3>`;
  for (const [g, gene] of [[GENOME_DATA[state.top], topGene], [GENOME_DATA[state.bot], botGene]]) {
    const row = document.createElement('div'); row.className = 'row';
    if (!gene) { row.innerHTML = `<span class="gname">${shortName(g.display)}</span><span class="meta">${link === undefined ? '' : 'no homolog'}</span>`; panel.appendChild(row); continue; }
    const seq = geneSeq(g.label, gene);
    const head = `${gene.name}_${shortName(g.display).replace(/[^A-Za-z0-9]/g, '')} ${gene.start}-${gene.end}(${gene.strand > 0 ? '+' : '-'}) ${seq.length}bp`;
    row.innerHTML = `<span class="gname">${shortName(g.display)} — ${gene.name}</span>` +
      `<span class="meta">${gene.start.toLocaleString()}–${gene.end.toLocaleString()} bp (${gene.strand > 0 ? '+' : '−'}) / ${seq.length} bp</span>`;
    const bCopy = document.createElement('button'); bCopy.textContent = 'Copy'; bCopy.onclick = () => copyText(seq, bCopy);
    const bFa = document.createElement('button'); bFa.textContent = 'FASTA'; bFa.onclick = () => download(`${gene.name}_${g.label}.fasta`, fasta(head, seq));
    row.append(bCopy, bFa);
    const pre = document.createElement('pre'); pre.textContent = seq;
    const nk = noteKey(g.label, gene.key);
    const note = document.createElement('textarea'); note.className = 'note';
    note.placeholder = `Note for ${gene.name} in ${shortName(g.display)} …`;
    note.value = state.notes[nk] || '';
    note.addEventListener('input', () => { state.notes[nk] = note.value; saveNotes(); });
    note.addEventListener('change', render);
    panel.append(row, pre, note);
  }
  panel.hidden = false;
  document.getElementById('seqClose').onclick = () => { state.pin = null; updateSeqPanel(); render(); };
}

/* ---------- range extraction ---------- */
function buildRangePanel() {
  const sel = document.getElementById('rangeGenome'); sel.innerHTML = '';
  for (const which of ['top', 'bot']) sel.add(new Option(`${which === 'top' ? 'Top' : 'Bottom'}: ${shortName(genomeOf(which).display)}`, which));
  document.getElementById('rangeGet').onclick = () => extractRange();
  for (const id of ['rangeGenome', 'rangeFrom', 'rangeLen']) document.getElementById(id).addEventListener('input', updateRangeTo);
  updateRangeTo();
}
function updateRangeTo() {
  const g = genomeOf(document.getElementById('rangeGenome').value || 'top');
  const from = Math.max(0, parseInt(document.getElementById('rangeFrom').value || '0', 10));
  const len = Math.max(1, parseInt(document.getElementById('rangeLen').value || '1', 10));
  document.getElementById('rangeTo').textContent = `→ ends at ${Math.min(g.length, from + len).toLocaleString()} bp  (From increases downstream, + direction)`;
}
let lastRange = null;
function setRange(which, from, len) {
  document.getElementById('rangeGenome').value = which;
  document.getElementById('rangeFrom').value = Math.round(from);
  document.getElementById('rangeLen').value = Math.max(1, Math.round(len));
  updateRangeTo(); extractRange();
}
function extractRange() {
  const which = document.getElementById('rangeGenome').value || 'top', g = genomeOf(which);
  let from = Math.max(0, Math.min(g.length - 1, parseInt(document.getElementById('rangeFrom').value || '0', 10)));
  let len = Math.max(1, parseInt(document.getElementById('rangeLen').value || '1', 10));
  const to = Math.min(g.length, from + len);
  let seq = (SEQ()[g.label] || '').slice(from, to);
  if (state.revcomp) seq = revComp(seq);
  const rc = state.revcomp ? ' rev-comp' : '';
  lastRange = { label: g.label, from, to, seq, header: `${g.label} ${from}-${to}${rc} (rotated LSC-first coords) ${seq.length}bp` };
  const out = document.getElementById('rangeOut'); out.hidden = false;
  out.innerHTML = `<div class="meta">${shortName(g.display)} : ${from.toLocaleString()}–${to.toLocaleString()} bp${rc} / ${seq.length} bp</div>`;
  const pre = document.createElement('pre'); pre.textContent = seq; out.appendChild(pre);
}
function copyRange() { if (!lastRange) extractRange(); if (lastRange) copyText(lastRange.seq, document.getElementById('rangeCopy')); }
function dlRange() { if (!lastRange) extractRange(); if (lastRange) download(`range_${lastRange.label}_${lastRange.from}-${lastRange.to}.fasta`, fasta(lastRange.header, lastRange.seq)); }
function downloadVisibleRange() {
  let out = '';
  for (const which of ['top', 'bot']) {
    const g = genomeOf(which), v = viewOf(which);
    const a = Math.max(0, Math.round(v.start * g.length)), b = Math.min(g.length, Math.round(v.end * g.length));
    out += fasta(`${g.label} ${a}-${b} (rotated LSC-first coords) ${(b - a)}bp`, (SEQ()[g.label] || '').slice(a, b));
  }
  download(`view_${state.top}_vs_${state.bot}.fasta`, out);
}

/* ---------- screenshot ---------- */
const SVG_CSS = `
.gene{stroke:rgba(0,0,0,.35);stroke-width:.5}.gene.dim{opacity:.18}.gene.hl{stroke:#111;stroke-width:1.6}
.ribbon{fill-opacity:.22;stroke:none}.ribbon.dim{fill-opacity:.05}.ribbon.hl{fill-opacity:.72}
.region{stroke:#fff;stroke-width:1}.region-label{font:700 11px sans-serif;fill:#555}
.genome-label{font:italic 13px sans-serif;fill:#333}.genome-sub{font:11px sans-serif;fill:#777}
.axis{stroke:#bbb;stroke-width:1}.axis-tick{font:10px sans-serif;fill:#999}
.gene-name{font:10px sans-serif;fill:#333}.gene-name.hl{font-weight:700;fill:#000}
.base{font:600 12px ui-monospace,monospace}.selrect{fill:rgba(21,101,192,.16);stroke:#1565c0;stroke-width:1}
.notemark{font:11px sans-serif;fill:#e65100}text{font-family:sans-serif}`;
function saveScreenshot() {
  const W = svg.clientWidth, H = SVG_H, scale = 2;
  const clone = svg.cloneNode(true);
  clone.setAttribute('xmlns', NS); clone.setAttribute('width', W); clone.setAttribute('height', H);
  const st = document.createElementNS(NS, 'style'); st.textContent = SVG_CSS; clone.insertBefore(st, clone.firstChild);
  const img = new Image();
  img.onload = () => {
    const c = document.createElement('canvas'); c.width = W * scale; c.height = H * scale;
    const ctx = c.getContext('2d'); ctx.scale(scale, scale); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, W, H); ctx.drawImage(img, 0, 0);
    c.toBlob(b => { const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = `chlocom_${state.top}_vs_${state.bot}.png`; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 1000); }, 'image/png');
  };
  img.onerror = () => alert('Screenshot failed.');
  img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(new XMLSerializer().serializeToString(clone));
}

/* ---------- interactions ---------- */
function onGeneHover(gene, which, ev) {
  state.hover = { which, key: gene.key }; applyHighlight();
  const link = which === 'top' ? state._tb.get(gene.key) : state._bt.get(gene.key);
  const hom = link ? genomeOf(which === 'top' ? 'bot' : 'top').genes.find(x => x.key === link.key) : null;
  const T = GENOME_DATA[state.top], B = GENOME_DATA[state.bot];
  const fmt = x => x ? `${x.start.toLocaleString()}–${x.end.toLocaleString()} bp (${x.strand > 0 ? '+' : '−'})` : '— none';
  const tG = which === 'top' ? gene : hom, bG = which === 'top' ? hom : gene;
  const idtxt = (state.mode === 'blast' && link && link.pid != null) ? ` <span class="mono">${link.pid}% id</span>` : '';
  tooltip.innerHTML = `<b>${gene.name}</b> <span class="mono">[${CAT[gene.cat][1]}]</span>${idtxt}<br>` +
    `${shortName(T.display)}: ${fmt(tG)}<br>${shortName(B.display)}: ${fmt(bG)}` +
    (hom ? (hom.name !== gene.name ? `<br><span class="mono">homolog: ${hom.name} — click to align</span>` : '<br><span class="mono">click to align &amp; show sequence</span>')
         : '<br><span class="mono">no homolog at current setting</span>');
  tooltip.hidden = false; moveTooltip(ev);
}
function onGeneLeave() { state.hover = null; tooltip.hidden = true; applyHighlight(); }
function moveTooltip(ev) {
  const pad = 14; let x = ev.clientX + pad, y = ev.clientY + pad;
  if (x + tooltip.offsetWidth > window.innerWidth) x = ev.clientX - tooltip.offsetWidth - pad;
  tooltip.style.left = x + 'px'; tooltip.style.top = y + 'px';
}
function onGeneClick(gene, which) {
  state.pin = (state.pin && state.pin.which === which && state.pin.key === gene.key) ? null : { which, key: gene.key };
  if (state.pin) alignOther(gene, which);
  updateSeqPanel(); render();
}
function alignOther(gene, which) {
  const link = which === 'top' ? state._tb.get(gene.key) : state._bt.get(gene.key);
  if (!link) return;
  const me = genomeOf(which), otherWhich = which === 'top' ? 'bot' : 'top', other = genomeOf(otherWhich);
  const hom = other.genes.find(x => x.key === link.key); if (!hom) return;
  const dc = dispGene(which, gene), dh = dispGene(otherWhich, hom);
  const newCenter = (dh.start + dh.end) / 2 / other.length - (dc.start + dc.end) / 2 / me.length + centerOf(which);
  if (otherWhich === 'top') state.centerTop = newCenter; else state.centerBot = newCenter;
  clampCenters();
}

svg.addEventListener('wheel', ev => {
  ev.preventDefault();
  const u = (ev.clientX - svg.getBoundingClientRect().left - L.ml) / plotW();
  if (u < 0 || u > 1) return;
  const spanNew = Math.min(1, Math.max(0.0008, state.span * (ev.deltaY < 0 ? 0.82 : 1.22)));
  for (const w of ['top', 'bot']) {
    const fUnder = centerOf(w) - state.span / 2 + u * state.span, c = fUnder + spanNew * (0.5 - u);
    if (w === 'top') state.centerTop = c; else state.centerBot = c;
  }
  state.span = spanNew; clampCenters(); render();
}, { passive: false });

let drag = null; const DRAG_PX = 4;
svg.addEventListener('pointerdown', ev => {
  const which = whichTrack(ev.clientY);
  if (ev.shiftKey) { const f = xToFrac(ev.clientX - svg.getBoundingClientRect().left, which); drag = { mode: 'sel', which, moved: false }; state.sel = { which, fA: f, fB: f }; }
  else drag = { mode: 'pan', x: ev.clientX, cT: state.centerTop, cB: state.centerBot, moved: false };
  svg.setPointerCapture(ev.pointerId);
});
svg.addEventListener('pointermove', ev => {
  if (!drag) return;
  if (drag.mode === 'sel') { drag.moved = true; state.sel.fB = xToFrac(ev.clientX - svg.getBoundingClientRect().left, drag.which); render(); }
  else {
    if (!drag.moved && Math.abs(ev.clientX - drag.x) < DRAG_PX) return;
    drag.moved = true; svg.classList.add('dragging');
    const dFrac = (ev.clientX - drag.x) / plotW() * state.span;
    state.centerTop = drag.cT - dFrac; state.centerBot = drag.cB - dFrac; clampCenters(); render();
  }
});
svg.addEventListener('pointerup', ev => {
  const d = drag; drag = null; svg.classList.remove('dragging');
  if (!d) return;
  if (d.mode === 'sel') {
    if (d.moved) { const g = genomeOf(d.which), fA = Math.min(state.sel.fA, state.sel.fB), fB = Math.max(state.sel.fA, state.sel.fB); setRange(d.which, fA * g.length, (fB - fA) * g.length); }
    else { state.sel = null; render(); }
  } else if (!d.moved) clickAt(ev.clientX, ev.clientY);
});
svg.addEventListener('pointercancel', () => { drag = null; svg.classList.remove('dragging'); });
function clickAt(cx, cy) {
  let t = document.elementFromPoint(cx, cy);
  while (t && !t.__gene && t !== document.body) t = t.parentNode;
  if (t && t.__gene) { onGeneClick(t.__gene, t.__which); return; }
  state.sel = null; state.pin = null; document.getElementById('rangeOut').hidden = true; updateSeqPanel(); render();
}

/* ---------- controls ---------- */
function buildSelects() {
  const top = document.getElementById('topSelect'), bot = document.getElementById('botSelect');
  for (const g of GENOME_INDEX) { top.add(new Option(`${g.display}  (${g.length.toLocaleString()} bp)`, g.label)); bot.add(new Option(`${g.display}  (${g.length.toLocaleString()} bp)`, g.label)); }
  top.value = state.top; bot.value = state.bot;
  top.onchange = () => { state.top = top.value; afterGenomeChange(); };
  bot.onchange = () => { state.bot = bot.value; afterGenomeChange(); };
}
function afterGenomeChange() { resetView(); buildGeneList(); buildRangePanel(); render(); }
function buildLegend() {
  const box = document.getElementById('legend'); box.innerHTML = '';
  const head = document.createElement('span'); head.className = 'leghead'; head.textContent = 'Gene categories:'; box.appendChild(head);
  const all = document.createElement('button'); all.className = 'legbtn'; all.textContent = 'Select all'; all.onclick = () => { state.hiddenCats.clear(); buildLegend(); render(); };
  const none = document.createElement('button'); none.className = 'legbtn'; none.textContent = 'Clear all'; none.onclick = () => { state.hiddenCats = new Set(Object.keys(CAT)); buildLegend(); render(); };
  box.append(all, none);
  for (const k in CAT) {
    const it = document.createElement('span'); it.className = 'item' + (state.hiddenCats.has(k) ? ' off' : '');
    it.innerHTML = `<span class="sw" style="background:${CAT[k][0]}"></span>${CAT[k][1]}`; it.title = 'Click to toggle';
    it.onclick = () => { state.hiddenCats.has(k) ? state.hiddenCats.delete(k) : state.hiddenCats.add(k); buildLegend(); render(); };
    box.appendChild(it);
  }
}
function buildGeneList() {
  const dl = document.getElementById('geneList'); dl.innerHTML = '';
  const names = new Set(); for (const lab of [state.top, state.bot]) for (const g of GENOME_DATA[lab].genes) names.add(g.name);
  [...names].sort().forEach(n => dl.appendChild(new Option(n)));
}
function doSearch(q) {
  q = q.trim().toLowerCase(); if (!q) return;
  for (const which of ['top', 'bot']) {
    const g = genomeOf(which);
    const hit = g.genes.find(x => x.name.toLowerCase() === q) || g.genes.find(x => x.name.toLowerCase().startsWith(q)) || g.genes.find(x => x.key === q.replace(/[^a-z0-9]/g, ''));
    if (hit) {
      state.span = 0.08; const dh = dispGene(which, hit), cF = (dh.start + dh.end) / 2 / g.length;
      if (which === 'top') state.centerTop = cF; else state.centerBot = cF;
      clampCenters(); state.pin = { which, key: hit.key }; alignOther(hit, which); updateSeqPanel(); render(); return;
    }
  }
  document.getElementById('status').textContent = `"${q}" not found`;
}
function resetView() { state.span = 1; state.centerTop = state.centerBot = 0.5; state.pin = null; state.sel = null; updateSeqPanel(); }
function updateThrUI() {
  // スライダーは BLASTN モードのときだけ操作可（presence は最後に設定した値を使う）
  const on = state.mode === 'blast';
  document.getElementById('thrWrap').classList.toggle('off', !on);
  document.getElementById('minPid').disabled = !on;
  document.getElementById('minPidVal').textContent = state.minPid + '%';
}
const nGenes = g => g.genes.reduce((n, x) => n + (x.cat !== 'intergenic' && !x.key.endsWith('__ir2')), 0);
function updateStatus() {
  const t = GENOME_DATA[state.top], b = GENOME_DATA[state.bot];
  const modetxt = state.mode === 'blast' ? `BLASTN ≥${state.minPid}%` : 'gene name';
  document.getElementById('status').textContent =
    `${modetxt}: ${state._tb.size} homolog link(s)  |  ${shortName(t.display)} ${nGenes(t)} · ${shortName(b.display)} ${nGenes(b)} genes  |  ` +
    `zoom ${(1 / state.span).toFixed(1)}×  |  view top ${Math.round(state.span * t.length).toLocaleString()} bp / bottom ${Math.round(state.span * b.length).toLocaleString()} bp`;
}

/* ---------- init ---------- */
function init() {
  if (!window.GENOME_INDEX || !GENOME_INDEX.length) { document.getElementById('status').textContent = 'No data. Run generate_data.py first.'; return; }
  const byLabel = Object.fromEntries(GENOME_INDEX.map(g => [g.label, g]));
  state.top = byLabel['Nicotiana_benthamiana'] ? 'Nicotiana_benthamiana' : GENOME_INDEX[0].label;
  state.bot = byLabel['Chlamydomonas_reinhardtii'] ? 'Chlamydomonas_reinhardtii' : (GENOME_INDEX[1] || GENOME_INDEX[0]).label;

  buildSelects(); buildLegend(); buildGeneList(); buildRangePanel();
  document.getElementById('onlyShared').onchange = e => { state.onlyShared = e.target.checked; render(); };
  document.getElementById('showLabels').onchange = e => { state.showLabels = e.target.checked; render(); };
  document.getElementById('revcomp').onchange = e => { state.revcomp = e.target.checked; render(); if (!document.getElementById('rangeOut').hidden) extractRange(); };
  document.getElementById('invSSC').onchange = e => { state.invSSC = e.target.checked; render(); };
  const modeSel = document.getElementById('mode'); modeSel.value = state.mode;
  modeSel.onchange = () => { state.mode = modeSel.value; updateThrUI(); render(); updateGeneInfo(); };
  const pidR = document.getElementById('minPid');
  const minFloor = (window.GENOME_BLAST_META || {}).minPid || 0;   // 内部でシンテニーとする最小%
  pidR.min = minFloor;                                             // それ未満はデータが無いので選べない
  if (state.minPid < minFloor) state.minPid = minFloor;
  pidR.value = state.minPid;
  document.getElementById('minPidVal').textContent = state.minPid + '%';
  pidR.oninput = () => { state.minPid = +pidR.value; document.getElementById('minPidVal').textContent = pidR.value + '%'; render(); updateGeneInfo(); };
  document.getElementById('resetBtn').onclick = () => { resetView(); render(); };
  document.getElementById('rangeBtn').onclick = downloadVisibleRange;
  document.getElementById('shotBtn').onclick = saveScreenshot;
  document.getElementById('notesExport').onclick = exportNotes;
  document.getElementById('notesImport').onclick = () => document.getElementById('notesFile').click();
  document.getElementById('notesFile').onchange = e => { if (e.target.files[0]) importNotes(e.target.files[0]); e.target.value = ''; };
  document.getElementById('swapBtn').onclick = () => {
    [state.top, state.bot] = [state.bot, state.top]; [state.centerTop, state.centerBot] = [state.centerBot, state.centerTop];
    document.getElementById('topSelect').value = state.top; document.getElementById('botSelect').value = state.bot;
    buildGeneList(); buildRangePanel(); updateSeqPanel(); render();
  };
  const search = document.getElementById('search');
  search.addEventListener('change', () => doSearch(search.value));
  search.addEventListener('keydown', e => { if (e.key === 'Enter') doSearch(search.value); });
  window.addEventListener('resize', render);
  updateThrUI(); render();
}
window.addEventListener('DOMContentLoaded', init);
