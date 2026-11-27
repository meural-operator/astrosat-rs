/* AstroSat Studio: Minimalist Offline Scientific Workbench Controller */
(() => {
  'use strict';

  const C = window.AstroCharts;
  const $ = (id) => document.getElementById(id);

  const STUDIO_CHART_COLORS = {
    light: {
      bg: '#ffffff',
      plotBg: '#f5f7fa',
      grid: 'rgba(14, 23, 38, 0.08)',
      axis: '#b9c5d4',
      text: '#0e1726',
      muted: '#4a596e',
      crosshair: 'rgba(18, 72, 179, 0.65)',
      gti: 'rgba(15, 118, 110, 0.10)',
      gtiBorder: 'rgba(15, 118, 110, 0.35)',
      refLine: 'rgba(185, 28, 28, 0.76)',
      badPix: [235, 240, 245],
    },
    dark: {
      bg: '#0e192b',
      plotBg: '#15243c',
      grid: 'rgba(148, 166, 192, 0.13)',
      axis: '#31486d',
      text: '#eef3f9',
      muted: '#94a6c0',
      crosshair: 'rgba(90, 156, 255, 0.70)',
      gti: 'rgba(45, 212, 191, 0.12)',
      gtiBorder: 'rgba(45, 212, 191, 0.35)',
      refLine: 'rgba(248, 113, 113, 0.80)',
      badPix: [21, 36, 60],
    },
  };

  const RUNNER_CONTEXT_NOTES = {
    czti: 'End-to-end CZTI Level-1 to Level-2 reduction: event cleaning, GTIs, DPH/DPI, sky images, light curves, and spectra.',
    laxpc: 'LAXPC Level-1 to Level-2 reduction with deadtime correction, MKF/SAA GTIs, Format A/B events, light curves, and spectra.',
    bary: 'Solar System barycentric time correction (TDB) for CZTI or LAXPC event and light-curve FITS files.',
    'czti-addspec': 'Co-add multiple CZTI OGIP .pha spectra and optional .rsp matrices across quadrants or orbits.',
    'laxpc-backshift': 'Compute LAXPC secular and HV-jump gain shifts and shift a background .pha spectrum to the observation epoch.',
  };

  const state = {
    theme: localStorage.getItem('astrosat-studio-theme') || 'light',
    health: null,
    activeTab: 'runner',
    runnerPipe: 'czti',
    jobPollTimer: null,
    jobPollFailures: 0,
    jobsCache: [],
    expandedJobOutputs: new Set(),
    collapsedJobOutputs: new Set(),
    lastSubmittedJobId: null,
    addspecRows: [{ pha: '', rsp: '' }, { pha: '', rsp: '' }],
    backshiftEpochMode: 'fits',
    cztiQuadsExplicit: false,
    lcMode: 'lightcurve',
    lcItems: [],
    lcPlot: null,
    imgPlot: null,
    imgData: null,
    pulsePlot: null,
    pgramPlot: null,
    pulsePreset: 'laxpc4',
    pulseInputs: [],
    pulseLastResult: null,
    inspDir: '',
    inspParent: '',
    inspEntries: [],
    inspData: null,
    inspHdu: 0,
    inspOffset: 0,
    inspLimit: 50,
    modalTarget: null,
    modalMode: 'file',
    modalDir: '',
    modalParent: '',
    modalEntries: [],
    modalCallback: null,
    modalReturnFocus: null,
    tomlReturnFocus: null,
    toastTimer: null,
  };

  // =========================================================================
  // Utilities
  // =========================================================================
  function esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  function basename(p) {
    if (!p) return '';
    return String(p).replace(/\\/g, '/').split('/').filter(Boolean).pop() || p;
  }

  function fmtBytes(n) {
    if (!Number.isFinite(n) || n < 0) return 'N/A';
    if (n < 1024) return `${n} B`;
    if (n < 1048576) return `${(n / 1024).toFixed(1)} KB`;
    if (n < 1073741824) return `${(n / 1048576).toFixed(1)} MB`;
    return `${(n / 1073741824).toFixed(2)} GB`;
  }

  function fmtFixed(v, d = 2, fallback = 'N/A') {
    return Number.isFinite(v) ? v.toFixed(d) : fallback;
  }

  function numOrNull(v) {
    if (v === '' || v == null) return null;
    const n = Number(v);
    return Number.isFinite(n) ? n : null;
  }

  function showBanner(prefix, msg) {
    const el = $(`${prefix}-banner`);
    const msgEl = $(`${prefix}-banner-msg`);
    if (!el || !msgEl) return;
    msgEl.textContent = msg;
    el.classList.remove('hidden');
  }

  function hideBanner(prefix) {
    const el = $(`${prefix}-banner`);
    if (el) el.classList.add('hidden');
  }

  function showToast(msg, { isError = false, actionLabel = '', onAction = null, durationMs = 5200 } = {}) {
    const toast = $('studio-toast');
    const msgEl = $('toast-msg');
    const actBtn = $('toast-action');
    if (!toast || !msgEl) return;
    clearTimeout(state.toastTimer);
    msgEl.textContent = msg;
    toast.classList.toggle('toast-err', isError);
    if (actBtn) {
      if (actionLabel && onAction) {
        actBtn.textContent = actionLabel;
        actBtn.classList.remove('hidden');
        actBtn.onclick = () => {
          onAction();
          hideToast();
        };
      } else {
        actBtn.classList.add('hidden');
        actBtn.onclick = null;
      }
    }
    toast.classList.remove('hidden');
    if (durationMs > 0) {
      state.toastTimer = setTimeout(hideToast, durationMs);
    }
  }

  function hideToast() {
    clearTimeout(state.toastTimer);
    const toast = $('studio-toast');
    if (toast) toast.classList.add('hidden');
  }

  async function apiGet(url) {
    const resp = await fetch(url);
    const data = await resp.json().catch(() => ({ error: `HTTP ${resp.status}` }));
    if (!resp.ok) throw new Error(data.error || `Request failed (${resp.status})`);
    return data;
  }

  async function apiPost(url, body) {
    const resp = await fetch(url, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    const data = await resp.json().catch(() => ({ error: `HTTP ${resp.status}` }));
    if (!resp.ok) throw new Error(data.error || `Request failed (${resp.status})`);
    return data;
  }

  // =========================================================================
  // Theme & Navigation
  // =========================================================================
  function applyTheme(theme) {
    state.theme = theme === 'dark' ? 'dark' : 'light';
    document.documentElement.setAttribute('data-theme', state.theme);
    localStorage.setItem('astrosat-studio-theme', state.theme);

    const pal = C.getPalette(state.theme);
    state.lcItems.forEach((it, idx) => { it.color = pal[idx % pal.length]; });
    if (state.lcPlot) {
      state.lcPlot.setTheme(state.theme);
      if (state.lcItems.length) {
        renderLcSeriesList();
        redrawLcPlot();
      }
    }
    if (state.imgPlot) state.imgPlot.setTheme(state.theme);
    if (state.pulsePlot) {
      state.pulsePlot.setTheme(state.theme);
      if (state.pulseLastResult) renderPulseResult(state.pulseLastResult);
    }
    if (state.pgramPlot) state.pgramPlot.setTheme(state.theme);

    const toggle = $('theme-toggle');
    if (toggle) toggle.textContent = state.theme === 'dark' ? 'Light' : 'Dark';
  }

  function safeResizeVisiblePlots(tab) {
    requestAnimationFrame(() => {
      requestAnimationFrame(() => {
        if (tab === 'lc-spec' && state.lcPlot) state.lcPlot.resize();
        if (tab === 'imaging' && state.imgPlot) state.imgPlot.resize();
        if (tab === 'pulse') {
          if (state.pulsePlot) state.pulsePlot.resize();
          if (state.pgramPlot && !$('periodogram-card').classList.contains('hidden')) {
            state.pgramPlot.resize();
          }
        }
      });
    });
  }

  function switchTab(tab) {
    state.activeTab = tab;
    document.querySelectorAll('.nav-link').forEach((b) => {
      const act = b.dataset.tab === tab;
      b.classList.toggle('active', act);
      b.setAttribute('aria-selected', act ? 'true' : 'false');
    });
    document.querySelectorAll('.view').forEach((p) => {
      p.classList.toggle('active', p.id === `tab-${tab}`);
    });
    safeResizeVisiblePlots(tab);
  }

  // =========================================================================
  // Shared File Browser Modal
  // =========================================================================
  function openFileModal({ targetId = null, mode = 'file', initialDir = '', onSelect = null }) {
    state.modalReturnFocus = document.activeElement instanceof HTMLElement ? document.activeElement : null;
    state.modalTarget = targetId;
    state.modalMode = mode;
    state.modalCallback = onSelect;
    $('file-modal-title').textContent = mode === 'dir' ? 'Select directory' : 'Select FITS or data file';
    $('modal-select-dir').classList.toggle('hidden', mode !== 'dir');
    $('modal-hint').textContent = mode === 'dir'
      ? 'Open a folder and click Use this directory'
      : 'Click a file to select it, or click a folder to open';
    $('modal-filter').value = '';
    $('file-modal').classList.remove('hidden');
    requestAnimationFrame(() => $('modal-filter').focus());

    const jumps = $('modal-quick-jumps');
    jumps.innerHTML = '';
    const h = state.health;
    if (h) {
      const candidates = [
        { label: 'Data root', path: h.data_dir },
        { label: 'CZTI L1', path: h.presets && h.presets.czti_l1_dir },
        { label: 'CZTI L2', path: h.presets && h.presets.czti_l2_dir },
        { label: 'LAXPC L1', path: h.presets && h.presets.laxpc_l1_dir },
        { label: 'LAXPC L2', path: h.presets && h.presets.laxpc_l2_dir },
      ].filter((c) => c.path);
      for (const c of candidates) {
        const b = document.createElement('button');
        b.type = 'button';
        b.className = 'text-action';
        b.textContent = c.label;
        b.onclick = () => loadModalDir(c.path);
        jumps.appendChild(b);
      }
    }
    loadModalDir(initialDir || (h && h.data_dir) || '');
  }

  function closeFileModal() {
    $('file-modal').classList.add('hidden');
    if (state.modalReturnFocus && typeof state.modalReturnFocus.focus === 'function') {
      state.modalReturnFocus.focus();
    }
    state.modalReturnFocus = null;
  }

  async function loadModalDir(dirPath) {
    try {
      const q = dirPath ? `?path=${encodeURIComponent(dirPath)}` : '';
      const data = await apiGet(`/api/browse${q}`);
      state.modalDir = data.path;
      state.modalParent = data.parent || '';
      state.modalEntries = data.entries || [];
      $('modal-path-input').value = data.path;
      renderModalEntries();
    } catch (err) {
      $('modal-entries-body').innerHTML = `<tr><td colspan="3" class="status-err">${esc(err.message)}</td></tr>`;
    }
  }

  function renderModalEntries() {
    const q = $('modal-filter').value.trim().toLowerCase();
    const body = $('modal-entries-body');
    const items = state.modalEntries.filter((e) => !q || e.name.toLowerCase().includes(q));
    if (!items.length) {
      body.innerHTML = '<tr><td colspan="3" class="table-empty">Empty directory</td></tr>';
      return;
    }
    body.innerHTML = items.map((e, i) => {
      const kind = e.is_dir ? 'Folder' : (e.product && e.product !== 'other' ? e.product : e.kind);
      const prefix = e.is_dir ? '+ ' : '';
      return `<tr class="clickable" data-idx="${i}">
        <td class="mono-cell">${esc(prefix + e.name)}</td>
        <td class="meta-muted">${esc(kind)}</td>
        <td class="num">${e.is_dir ? '' : fmtBytes(e.size_bytes)}</td>
      </tr>`;
    }).join('');
    body.querySelectorAll('tr.clickable').forEach((tr) => {
      tr.onclick = () => {
        const entry = items[Number(tr.dataset.idx)];
        if (entry.is_dir) {
          loadModalDir(entry.path);
        } else if (state.modalMode === 'file') {
          if (state.modalTarget && $(state.modalTarget)) {
            $(state.modalTarget).value = entry.path;
            $(state.modalTarget).dispatchEvent(new Event('input', { bubbles: true }));
          }
          if (state.modalCallback) state.modalCallback(entry.path);
          closeFileModal();
        }
      };
    });
  }

  // =========================================================================
  // VIEW 1: PIPELINE RUNNER & STANDALONE TOOLS
  // =========================================================================
  function setRunnerPipe(pipe) {
    state.runnerPipe = pipe;
    document.querySelectorAll('#runner-pipeline-seg button, #runner-tools-seg button').forEach((x) => {
      x.classList.toggle('active', x.dataset.pipe === pipe);
    });
    const forms = ['czti', 'laxpc', 'bary', 'czti-addspec', 'laxpc-backshift'];
    for (const f of forms) {
      const el = $(`form-${f}`);
      if (el) el.classList.toggle('hidden', pipe !== f);
    }
    const note = $('runner-context-note');
    if (note) note.textContent = RUNNER_CONTEXT_NOTES[pipe] || '';

    const submit = $('runner-submit');
    if (submit) {
      if (pipe === 'czti-addspec') submit.textContent = 'Co-add spectra';
      else if (pipe === 'laxpc-backshift') submit.textContent = 'Compute backshift';
      else submit.textContent = 'Run pipeline';
    }

    const isCliTool = pipe === 'czti-addspec' || pipe === 'laxpc-backshift';
    $('copy-toml-btn').disabled = isCliTool;
    $('download-toml-btn').disabled = isCliTool;
    $('preview-toml-caption').textContent = isCliTool
      ? 'TOML configuration (not used by standalone CLI subcommands)'
      : 'TOML configuration';

    updateCliPreview();
  }

  function setBackshiftEpochMode(mode) {
    state.backshiftEpochMode = mode === 'day' ? 'day' : 'fits';
    document.querySelectorAll('#backshift-epoch-seg button').forEach((b) => {
      b.classList.toggle('active', b.dataset.epoch === state.backshiftEpochMode);
    });
    $('backshift-fits-wrap').classList.toggle('hidden', state.backshiftEpochMode !== 'fits');
    $('backshift-day-wrap').classList.toggle('hidden', state.backshiftEpochMode !== 'day');
    scheduleCliPreview();
  }

  function stemWithoutExt(p) {
    const base = basename(p);
    return base.replace(/\.gz$/i, '').replace(/\.[^.]+$/, '');
  }

  function findPairedRspFromCache(phaPath) {
    const targetStem = stemWithoutExt(phaPath);
    const parentDir = phaPath.replace(/\\/g, '/').replace(/\/[^/]+$/, '');
    for (const j of state.jobsCache) {
      const outs = j.outputs || [];
      for (const o of outs) {
        if ((o.product === 'response' || /\.rsp$/i.test(o.path)) && stemWithoutExt(o.path) === targetStem) {
          return o.path;
        }
      }
    }
    if (targetStem) {
      return `${parentDir}/${targetStem}.rsp`;
    }
    return '';
  }

  async function resolveVerifiedRsp(phaPath) {
    const targetStem = stemWithoutExt(phaPath);
    for (const j of state.jobsCache) {
      for (const o of (j.outputs || [])) {
        if ((o.product === 'response' || /\.rsp$/i.test(o.path)) && stemWithoutExt(o.path) === targetStem) {
          return o.path;
        }
      }
    }
    const norm = phaPath.replace(/\\/g, '/');
    const slash = norm.lastIndexOf('/');
    if (slash > 0) {
      const dir = norm.slice(0, slash);
      try {
        const res = await apiGet(`/api/browse?path=${encodeURIComponent(dir)}`);
        const match = (res.entries || []).find((e) => !e.is_dir && /\.rsp(\.gz)?$/i.test(e.name) && stemWithoutExt(e.name) === targetStem);
        if (match) return match.path;
      } catch {
        // Ignore directory lookup failures
      }
    }
    return '';
  }

  function populateAddspecQuickPick() {
    const sel = $('addspec-quick-pick');
    if (!sel) return;
    sel.innerHTML = '<option value="">Add a detected .pha (and paired .rsp)...</option>';
    const seen = new Set();
    const addOpt = (path, label) => {
      if (!path || seen.has(path)) return;
      seen.add(path);
      const opt = document.createElement('option');
      opt.value = path;
      opt.textContent = label;
      opt.title = path;
      sel.appendChild(opt);
    };
    for (const j of state.jobsCache) {
      if (j.pipeline === 'czti' && j.status === 'succeeded') {
        for (const o of (j.outputs || [])) {
          if (o.product === 'spectrum' || /\.pha$/i.test(o.path)) {
            addOpt(o.path, `${j.id}: ${basename(o.path)}`);
          }
        }
      }
    }
    const samples = (state.health && state.health.presets && state.health.presets.sample_spectra) || [];
    for (const s of samples) {
      if (/czt/i.test(s) || /\.pha$/i.test(s)) {
        addOpt(s, `Sample: ${basename(s)}`);
      }
    }
  }

  function readAddspecRowsFromDom() {
    const domRows = document.querySelectorAll('#addspec-rows .pair-row');
    if (domRows.length > 0) {
      const collected = Array.from(domRows).map((el) => {
        const phaEl = el.querySelector('.addspec-pha');
        const rspEl = el.querySelector('.addspec-rsp');
        return {
          pha: phaEl ? phaEl.value : '',
          rsp: rspEl ? rspEl.value : '',
        };
      });
      state.addspecRows = collected;
      return collected;
    }
    return state.addspecRows;
  }

  function evaluateAddspecRows() {
    const rawRows = readAddspecRowsFromDom();
    const indexedNonEmpty = rawRows
      .map((r, i) => ({
        rowNum: i + 1,
        pha: (r.pha || '').trim(),
        rsp: (r.rsp || '').trim(),
      }))
      .filter((r) => r.pha || r.rsp);

    const orphanRsp = indexedNonEmpty.find((r) => !r.pha && r.rsp);
    const hasAnyRsp = indexedNonEmpty.some((r) => Boolean(r.rsp));
    const missingRsp = hasAnyRsp ? indexedNonEmpty.find((r) => r.pha && !r.rsp) : null;

    let error = '';
    if (orphanRsp) {
      error = `Row #${orphanRsp.rowNum} specifies a .rsp without a matching .pha spectrum.`;
    } else if (missingRsp) {
      error = `Row #${missingRsp.rowNum} is missing a .rsp (provide a matching .rsp for all ${indexedNonEmpty.length} rows, or leave all .rsp fields blank).`;
    }

    return {
      rows: indexedNonEmpty,
      pha: indexedNonEmpty.map((r) => r.pha),
      rsp: hasAnyRsp ? indexedNonEmpty.map((r) => r.rsp) : [],
      hasAnyRsp,
      error,
    };
  }

  function updateAddspecHint() {
    const hint = $('addspec-row-hint');
    if (!hint) return;
    const evalRes = evaluateAddspecRows();
    if (evalRes.error) {
      hint.textContent = evalRes.error;
      hint.classList.add('status-warn');
    } else {
      hint.textContent = 'Provide a matching .rsp for every row, or leave all .rsp fields blank.';
      hint.classList.remove('status-warn');
    }
  }

  function renderAddspecRows() {
    const box = $('addspec-rows');
    if (!box) return;
    box.innerHTML = state.addspecRows.map((row, idx) => `
      <div class="pair-row" data-idx="${idx}">
        <span class="pair-index">#${idx + 1}</span>
        <div class="pair-fields">
          <div class="path-field">
            <input type="text" class="input-path addspec-pha" value="${esc(row.pha)}" placeholder="Input spectrum .pha (required)" />
            <button type="button" class="subtle-btn addspec-browse-pha">Browse</button>
          </div>
          <div class="path-field">
            <input type="text" class="input-path addspec-rsp" value="${esc(row.rsp)}" placeholder="Optional response .rsp (all rows or none)" />
            <button type="button" class="subtle-btn addspec-browse-rsp">Browse</button>
          </div>
        </div>
        <div class="pair-actions">
          <button type="button" class="text-action addspec-up" ${idx === 0 ? 'disabled' : ''} title="Move row up">Up</button>
          <button type="button" class="text-action addspec-down" ${idx === state.addspecRows.length - 1 ? 'disabled' : ''} title="Move row down">Down</button>
          ${state.addspecRows.length > 1 ? '<button type="button" class="text-action addspec-del">Remove</button>' : ''}
        </div>
      </div>
    `).join('');

    box.querySelectorAll('.pair-row').forEach((el) => {
      const idx = Number(el.dataset.idx);
      const row = state.addspecRows[idx];
      const syncPha = (e) => {
        row.pha = e.target.value;
        updateAddspecHint();
        scheduleCliPreview();
      };
      const syncRsp = (e) => {
        row.rsp = e.target.value;
        updateAddspecHint();
        scheduleCliPreview();
      };
      el.querySelector('.addspec-pha').oninput = syncPha;
      el.querySelector('.addspec-pha').onchange = syncPha;
      el.querySelector('.addspec-rsp').oninput = syncRsp;
      el.querySelector('.addspec-rsp').onchange = syncRsp;
      el.querySelector('.addspec-browse-pha').onclick = () => openFileModal({
        mode: 'file',
        onSelect: async (p) => {
          row.pha = p;
          if (!row.rsp) {
            const rsp = await resolveVerifiedRsp(p);
            if (rsp) row.rsp = rsp;
          }
          renderAddspecRows();
          scheduleCliPreview();
        },
      });
      el.querySelector('.addspec-browse-rsp').onclick = () => openFileModal({
        mode: 'file',
        onSelect: (p) => {
          row.rsp = p;
          renderAddspecRows();
          scheduleCliPreview();
        },
      });
      el.querySelector('.addspec-up').onclick = () => {
        if (idx > 0) {
          const tmp = state.addspecRows[idx - 1];
          state.addspecRows[idx - 1] = state.addspecRows[idx];
          state.addspecRows[idx] = tmp;
          renderAddspecRows();
          scheduleCliPreview();
        }
      };
      el.querySelector('.addspec-down').onclick = () => {
        if (idx < state.addspecRows.length - 1) {
          const tmp = state.addspecRows[idx + 1];
          state.addspecRows[idx + 1] = state.addspecRows[idx];
          state.addspecRows[idx] = tmp;
          renderAddspecRows();
          scheduleCliPreview();
        }
      };
      const delBtn = el.querySelector('.addspec-del');
      if (delBtn) {
        delBtn.onclick = () => {
          state.addspecRows.splice(idx, 1);
          renderAddspecRows();
          scheduleCliPreview();
        };
      }
    });
    updateAddspecHint();
  }

  async function discoverWorkspaceCztiSpectra() {
    const dataDir = state.health && state.health.data_dir;
    if (!dataDir) return [];
    const tryDir = async (dir) => {
      try {
        const res = await apiGet(`/api/browse?path=${encodeURIComponent(dir)}`);
        const entries = res.entries || [];
        const phas = entries.filter((e) => !e.is_dir && /\.pha$/i.test(e.name));
        if (phas.length) {
          return phas.slice(0, 4).map((p) => {
            const stem = stemWithoutExt(p.name);
            const rsp = entries.find((e) => !e.is_dir && /\.rsp$/i.test(e.name) && stemWithoutExt(e.name) === stem);
            return { pha: p.path, rsp: rsp ? rsp.path : '' };
          });
        }
        return entries.filter((e) => e.is_dir);
      } catch {
        return [];
      }
    };
    const root = dataDir.replace(/[\\/]+$/, '');
    for (const base of [`${root}/sample/czti_2018`, `${root}/czti_2018`]) {
      const top = await tryDir(base);
      if (top.length && top[0].pha) return top;
      for (const sub of top) {
        if (!sub.path) continue;
        const m0 = await tryDir(`${sub.path}/czti/modeM0`);
        if (m0.length && m0[0].pha) return m0;
      }
    }
    return [];
  }

  async function useLatestCztiOutputsForAddspec() {
    const baseOut = state.health && state.health.cwd
      ? `${state.health.cwd.replace(/[\\/]+$/, '')}/.tmp`
      : './.tmp';
    const cztiJob = state.jobsCache.find((j) => j.pipeline === 'czti' && j.status === 'succeeded'
      && (j.outputs || []).some((o) => o.product === 'spectrum' || /\.pha$/i.test(o.path)));

    let rows = [];
    if (cztiJob) {
      const phas = (cztiJob.outputs || []).filter((o) => o.product === 'spectrum' || /\.pha$/i.test(o.path));
      const rsps = (cztiJob.outputs || []).filter((o) => o.product === 'response' || /\.rsp$/i.test(o.path));
      rows = phas.map((p) => {
        const stem = stemWithoutExt(p.path);
        const match = rsps.find((r) => stemWithoutExt(r.path) === stem);
        return { pha: p.path, rsp: match ? match.path : '' };
      });
    } else {
      const samples = ((state.health && state.health.presets && state.health.presets.sample_spectra) || [])
        .filter((s) => /czt/i.test(s) && /\.pha$/i.test(s))
        .slice(0, 4);
      for (const s of samples) {
        const rsp = await resolveVerifiedRsp(s);
        rows.push({ pha: s, rsp });
      }
      if (!rows.length) {
        rows = await discoverWorkspaceCztiSpectra();
      }
    }
    if (!rows.length) {
      showToast('No CZTI .pha spectra found in recent jobs or workspace samples.', { isError: true });
      return;
    }
    state.addspecRows = rows;
    const allHaveRsp = rows.every((r) => Boolean(r.rsp));
    $('addspec-out-pha').value = `${baseOut}/studio_coadd.pha`;
    $('addspec-out-rsp').value = allHaveRsp ? `${baseOut}/studio_coadd.rsp` : '';
    renderAddspecRows();
    updateCliPreview();
  }

  function collectRunnerPayload() {
    const pipe = state.runnerPipe;
    const threads = numOrNull($('runner-threads').value);
    const payload = { pipeline: pipe };
    if (threads != null && threads >= 1) payload.threads = Math.floor(threads);

    if (pipe === 'czti') {
      const czti = {};
      const setStr = (k, id) => { const v = $(id).value.trim(); if (v) czti[k] = v; };
      const setNum = (k, id) => { const v = numOrNull($(id).value); if (v != null) czti[k] = v; };
      setStr('input_dir', 'czti-input-dir');
      setStr('m0', 'czti-m0');
      setStr('ss', 'czti-ss');
      setStr('tct', 'czti-tct');
      setStr('mkf', 'czti-mkf');
      setStr('caldb_dir', 'czti-caldb-dir');
      setStr('output_dir', 'czti-output-dir');
      setStr('output_prefix', 'czti-prefix');
      setStr('mkf_threshold_file', 'czti-mkf-thresh');
      setStr('user_gti_file', 'czti-user-gti');
      czti.entry_stage = Number($('czti-entry-stage').value);
      czti.exit_stage = Number($('czti-exit-stage').value);
      czti.algo = $('czti-algo').value;
      czti.gti_type = $('czti-gti-type').value;
      setNum('ra_deg', 'czti-ra');
      setNum('dec_deg', 'czti-dec');
      setNum('time_binsize', 'czti-tbin');
      setNum('emin_kev', 'czti-emin');
      setNum('emax_kev', 'czti-emax');
      const outType = $('czti-output-type').value;
      if (outType && outType !== 'both') czti.output_type = outType;
      czti.image_algo = $('czti-img-algo').value;
      const ov = numOrNull($('czti-oversampling').value);
      if (ov != null) czti.oversampling = Math.floor(ov);
      const selectedQuads = Array.from(document.querySelectorAll('.czti-quad-cb:checked')).map((cb) => Number(cb.value));
      if (state.cztiQuadsExplicit || selectedQuads.length < 4) {
        czti.quads = selectedQuads;
      }
      czti.mask_weight = $('czti-mask-weight').checked;
      czti.gen_resp = $('czti-gen-resp').checked;
      if ($('czti-gen-evt').checked) czti.generate_eventfile = true;
      payload.czti = czti;
    } else if (pipe === 'laxpc') {
      const laxpc = {};
      const setStr = (k, id) => { const v = $(id).value.trim(); if (v) laxpc[k] = v; };
      const setNum = (k, id) => { const v = numOrNull($(id).value); if (v != null) laxpc[k] = v; };
      setStr('input_dir', 'laxpc-input-dir');
      setStr('caldb_dir', 'laxpc-caldb-dir');
      setStr('output_dir', 'laxpc-output-dir');
      const unitsStr = $('laxpc-units').value.trim();
      if (unitsStr) {
        const units = unitsStr.split(/[,\s]+/).map((s) => s.trim()).filter(Boolean);
        laxpc.units = units;
      }
      setNum('tbin', 'laxpc-tbin');
      setNum('deadtime_sec', 'laxpc-deadtime');
      laxpc.spectrum_channels = Number($('laxpc-channels').value);
      setNum('min_elv_deg', 'laxpc-min-elv');
      setStr('anodes', 'laxpc-anodes');
      setStr('energy_range', 'laxpc-energy-range');
      const pMin = numOrNull($('laxpc-pha-min').value);
      if (pMin != null) laxpc.pha_min = Math.floor(pMin);
      const pMax = numOrNull($('laxpc-pha-max').value);
      if (pMax != null) laxpc.pha_max = Math.floor(pMax);
      setNum('saa_entry_lon_deg', 'laxpc-saa-entry');
      setNum('saa_exit_lon_deg', 'laxpc-saa-exit');
      setNum('merge_gap_sec', 'laxpc-merge-gap');
      setNum('min_duration_sec', 'laxpc-min-dur');
      setStr('user_gti_file', 'laxpc-user-gti');
      if ($('laxpc-user-gti-only').checked) laxpc.user_gti_only = true;
      if ($('laxpc-ascii').checked) laxpc.write_ascii_products = true;
      setStr('background_model_path', 'laxpc-bkg-model');
      setStr('background_pha_path', 'laxpc-bkg-pha');
      setNum('background_peak_ch', 'laxpc-bkg-peak');
      laxpc.write_format_a_combined = $('laxpc-fmt-a').checked;
      laxpc.write_format_b_events = $('laxpc-fmt-b').checked;
      laxpc.write_orbit_file = $('laxpc-orbit').checked;
      laxpc.single_only = $('laxpc-single').checked;
      laxpc.layer1_only = $('laxpc-layer1').checked;
      payload.laxpc = laxpc;
    } else if (pipe === 'bary') {
      const bary = {};
      const setStr = (k, id) => { const v = $(id).value.trim(); if (v) bary[k] = v; };
      const setNum = (k, id) => { const v = numOrNull($(id).value); if (v != null) bary[k] = v; };
      setStr('input', 'bary-input');
      setStr('output', 'bary-output');
      const orb = $('bary-orbit').value.trim();
      const mkf = $('bary-mkf').value.trim();
      if (orb) bary.orbit = orb;
      else if (mkf) bary.mkf = mkf;
      setStr('ephemeris', 'bary-ephem');
      setNum('ra_deg', 'bary-ra');
      setNum('dec_deg', 'bary-dec');
      bary.frame = $('bary-frame').value;
      bary.write_checksums = $('bary-checksums').checked;
      payload.bary = bary;
    } else if (pipe === 'czti-addspec') {
      const evalRes = evaluateAddspecRows();
      const addspec = { pha: evalRes.pha };
      if (evalRes.hasAnyRsp) addspec.rsp = evalRes.rsp;
      const outPha = $('addspec-out-pha').value.trim();
      if (outPha) addspec.out_pha = outPha;
      const outRsp = $('addspec-out-rsp').value.trim();
      if (outRsp) addspec.out_rsp = outRsp;
      addspec.add_exposure = $('addspec-add-exposure').checked;
      payload.czti_addspec = addspec;
    } else if (pipe === 'laxpc-backshift') {
      const bs = {
        unit: $('backshift-unit').value,
        layer1_only: $('backshift-layer1').checked,
        single_only: $('backshift-single').checked,
      };
      const caldb = $('backshift-caldb').value.trim();
      if (caldb) bs.caldb_dir = caldb;
      if (state.backshiftEpochMode === 'fits') {
        const f = $('backshift-fits').value.trim();
        if (f) bs.fits = f;
      } else {
        const md = numOrNull($('backshift-mission-day').value);
        if (md != null) bs.mission_day = md;
      }
      const bkg = $('backshift-bkg-pha').value.trim();
      if (bkg) bs.bkg_pha = bkg;
      const bpk = numOrNull($('backshift-bkg-peak').value);
      if (bpk != null) bs.bkg_peak_ch = bpk;
      const out = $('backshift-output').value.trim();
      if (out) bs.output = out;
      const pMin = numOrNull($('backshift-pha-min').value);
      if (pMin != null) bs.pha_min = Math.floor(pMin);
      const pMax = numOrNull($('backshift-pha-max').value);
      if (pMax != null) bs.pha_max = Math.floor(pMax);
      payload.laxpc_backshift = bs;
    }
    return payload;
  }

  let previewTimer = null;
  function scheduleCliPreview() {
    clearTimeout(previewTimer);
    previewTimer = setTimeout(updateCliPreview, 160);
  }

  async function updateCliPreview() {
    if (state.runnerPipe === 'czti') {
      const selectedQuads = document.querySelectorAll('.czti-quad-cb:checked');
      if (selectedQuads.length === 0) {
        const msg = 'Select at least one CZTI quadrant (Q0 to Q3).';
        $('preview-cli').textContent = `# Preview error: ${msg}`;
        $('preview-toml').textContent = `# Preview error: ${msg}`;
        return;
      }
    }
    if (state.runnerPipe === 'czti-addspec') {
      const evalRes = evaluateAddspecRows();
      if (evalRes.error) {
        $('preview-cli').textContent = `# Preview error: ${evalRes.error}`;
        $('preview-toml').textContent = '# Standalone CLI tool: astrosat czti addspec uses CLI flags (no AppConfig TOML table).\n';
        return;
      }
    }
    try {
      const data = await apiPost('/api/cli-preview', collectRunnerPayload());
      $('preview-cli').textContent = data.cli_command || '';
      $('preview-toml').textContent = data.toml_config || '';
    } catch (err) {
      $('preview-cli').textContent = `# Preview error: ${err.message}`;
      $('preview-toml').textContent = `# Preview error: ${err.message}`;
    }
  }

  async function autofillRunnerSample() {
    const p = state.health && state.health.presets;
    if (!p) return;
    const baseOut = state.health && state.health.cwd
      ? `${state.health.cwd.replace(/[\\/]+$/, '')}/.tmp`
      : './.tmp';
    if (state.runnerPipe === 'czti') {
      if (p.czti_l1_dir) $('czti-input-dir').value = p.czti_l1_dir;
      if (p.czti_caldb) $('czti-caldb-dir').value = p.czti_caldb;
      $('czti-output-dir').value = `${baseOut}/studio_czti_out`;
      $('czti-prefix').value = 'czti_studio';
    } else if (state.runnerPipe === 'laxpc') {
      if (p.laxpc_l1_dir) $('laxpc-input-dir').value = p.laxpc_l1_dir;
      if (p.laxpc_caldb) $('laxpc-caldb-dir').value = p.laxpc_caldb;
      $('laxpc-output-dir').value = `${baseOut}/studio_laxpc_out`;
    } else if (state.runnerPipe === 'bary') {
      const ev = (p.sample_events && p.sample_events[0]) || '';
      if (ev) $('bary-input').value = ev;
      if (p.laxpc_orb || p.czti_orb) {
        $('bary-orbit').value = p.laxpc_orb || p.czti_orb;
        $('bary-mkf').value = '';
      }
      if (p.jpleph_405) $('bary-ephem').value = p.jpleph_405;
      $('bary-output').value = `${baseOut}/studio_bary.evt`;
      $('bary-ra').value = '83.633083';
      $('bary-dec').value = '22.014500';
    } else if (state.runnerPipe === 'czti-addspec') {
      await useLatestCztiOutputsForAddspec();
      return;
    } else if (state.runnerPipe === 'laxpc-backshift') {
      $('backshift-unit').value = '10';
      if (p.laxpc_caldb) $('backshift-caldb').value = p.laxpc_caldb;
      const laxEv = findSampleEvent('laxpc');
      if (laxEv) {
        $('backshift-fits').value = laxEv;
        setBackshiftEpochMode('fits');
      } else {
        $('backshift-mission-day').value = '250';
        setBackshiftEpochMode('day');
      }
      const dataRoot = (state.health && state.health.data_dir || '').replace(/[\\/]+$/, '');
      const bkgDirs = [p.laxpc_caldb, dataRoot ? `${dataRoot}/ref_src/laxpc` : ''].filter(Boolean);
      const browsed = await Promise.all(bkgDirs.map((bd) => apiGet(`/api/browse?path=${encodeURIComponent(bd)}`).catch(() => null)));
      for (const dir of browsed) {
        if (!dir) continue;
        const bkg = (dir.entries || []).find((e) => !e.is_dir && /^backlxp10.*\.pha$/i.test(e.name));
        if (bkg) {
          $('backshift-bkg-pha').value = bkg.path;
          $('backshift-output').value = `${baseOut}/studio_shifted_bkg.pha`;
          break;
        }
      }
    }
    updateCliPreview();
  }

  async function submitJob() {
    hideBanner('runner');
    if (state.runnerPipe === 'czti') {
      const selectedQuads = document.querySelectorAll('.czti-quad-cb:checked');
      if (selectedQuads.length === 0) {
        const msg = 'Select at least one CZTI quadrant (Q0 to Q3).';
        showBanner('runner', msg);
        const statusEl = $('runner-inline-status');
        if (statusEl) statusEl.textContent = `Validation error: ${msg}`;
        showToast(msg, { isError: true });
        return;
      }
    }
    if (state.runnerPipe === 'czti-addspec') {
      const evalRes = evaluateAddspecRows();
      if (evalRes.error) {
        showBanner('runner', evalRes.error);
        const statusEl = $('runner-inline-status');
        if (statusEl) statusEl.textContent = `Validation error: ${evalRes.error}`;
        showToast(evalRes.error, { isError: true });
        return;
      }
    }
    const btn = $('runner-submit');
    btn.disabled = true;
    try {
      const rec = await apiPost('/api/jobs/run', collectRunnerPayload());
      state.lastSubmittedJobId = rec.id;
      const statusEl = $('runner-inline-status');
      if (statusEl) statusEl.textContent = `Started ${rec.id} (${rec.pipeline}).`;
      showToast(`Started ${rec.id} (${rec.pipeline}).`, {
        actionLabel: 'View job',
        onAction: () => {
          const stage = $('jobs-stage');
          if (stage) stage.scrollIntoView({ behavior: 'smooth', block: 'start' });
        },
      });
      await refreshJobs();
      startJobPolling();
      if (window.innerWidth <= 1024) {
        const stage = $('jobs-stage');
        if (stage) stage.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
      }
    } catch (err) {
      showBanner('runner', err.message);
      const statusEl = $('runner-inline-status');
      if (statusEl) statusEl.textContent = `Submission error: ${err.message}`;
      showToast(err.message, { isError: true });
    } finally {
      btn.disabled = false;
    }
  }

  function startJobPolling() {
    if (state.jobPollTimer) return;
    state.jobPollTimer = setInterval(async () => {
      const jobs = await refreshJobs();
      if (jobs === null) {
        if ((state.jobPollFailures || 0) >= 5) {
          clearInterval(state.jobPollTimer);
          state.jobPollTimer = null;
        }
        return;
      }
      if (!jobs.some((j) => j.status === 'running' || j.status === 'queued')) {
        clearInterval(state.jobPollTimer);
        state.jobPollTimer = null;
      }
    }, 1200);
  }

  async function refreshJobs() {
    try {
      const data = await apiGet('/api/jobs');
      state.jobPollFailures = 0;
      const jobs = data.jobs || [];
      state.jobsCache = jobs;
      populateAddspecQuickPick();
      if (state.lastSubmittedJobId) {
        const tracked = jobs.find((j) => j.id === state.lastSubmittedJobId);
        const statusEl = $('runner-inline-status');
        if (tracked && statusEl) {
          if (tracked.status === 'succeeded') {
            const sec = Number.isFinite(tracked.elapsed_ms) ? `${(tracked.elapsed_ms / 1000).toFixed(2)} s` : '';
            statusEl.textContent = `${tracked.id} succeeded${sec ? ` in ${sec}` : ''}.`;
          } else if (tracked.status === 'failed') {
            statusEl.textContent = `${tracked.id} failed: ${tracked.error || 'see job details'}`;
          }
        }
      }
      renderJobs(jobs);
      if (jobs.some((j) => j.status === 'running' || j.status === 'queued')) {
        startJobPolling();
      }
      return jobs;
    } catch {
      state.jobPollFailures = (state.jobPollFailures || 0) + 1;
      return null;
    }
  }

  function renderJobSummaryHtml(j) {
    const s = j.summary;
    if (s == null) return '';
    if (typeof s === 'string') return `<div class="job-summary">${esc(s)}</div>`;
    if (typeof s !== 'object') return `<div class="job-summary">${esc(String(s))}</div>`;

    const kvGrid = (pairs) => {
      const items = pairs.filter(([, val]) => val != null && val !== '');
      if (!items.length) return '';
      return `<dl class="job-kv-grid">${items.map(([label, val]) => `
        <div class="job-kv-item"><dt>${esc(label)}</dt><dd>${esc(val)}</dd></div>
      `).join('')}</dl>`;
    };

    if (j.pipeline === 'laxpc-backshift') {
      return kvGrid([
        ['Unit', s.unit],
        ['Source / bkg day', `${fmtFixed(s.source_day, 2)} / ${fmtFixed(s.background_day, 2)}`],
        ['Delta ch (30 keV)', Number.isFinite(s.delta_ch_30kev) ? `${s.delta_ch_30kev >= 0 ? '+' : ''}${s.delta_ch_30kev.toFixed(3)} ch` : null],
        ['Gain ratio', fmtFixed(s.gain_ratio, 5, null)],
        ['Secular scale', fmtFixed(s.secular_scale, 5, null)],
        ['Crossed HV jumps', s.crossed_hv_jumps != null ? String(s.crossed_hv_jumps) : null],
        ['Recommended RMF', s.recommended_rmf],
        ['LC bkg shift frac', Number.isFinite(s.lcbk_shift_fraction) ? C.formatSci(s.lcbk_shift_fraction, 4) : null],
      ]);
    }

    if (j.pipeline === 'czti-addspec') {
      return kvGrid([
        ['Input spectra', s.inputs != null ? String(s.inputs) : null],
        ['Input responses', s.rsp_inputs != null ? String(s.rsp_inputs) : null],
        ['Exposure mode', s.add_exposure ? 'Summed (multi-orbit)' : 'Averaged (multi-quadrant)'],
        ['Exposure', Number.isFinite(s.exposure_sec) ? `${s.exposure_sec.toFixed(2)} s` : null],
        ['Channels', s.channels != null ? String(s.channels) : null],
        ['Total counts', Number.isFinite(s.total_counts) ? C.formatSci(s.total_counts, 5) : null],
      ]);
    }

    if (j.pipeline === 'laxpc' && Array.isArray(s.units)) {
      const unitSummaries = s.units.map((u) => {
        const exp = Number.isFinite(u.exposure_sec) ? `${u.exposure_sec.toFixed(1)}s` : 'N/A';
        const rate = Number.isFinite(u.lightcurve_mean_rate_cps) ? `${C.formatSci(u.lightcurve_mean_rate_cps, 3)} c/s` : 'N/A';
        return `${u.unit}: ${u.gti_filtered_events || 0} evts, ${exp}, ${rate}`;
      }).join(' | ');
      return `<div class="job-summary">${esc(`Units: ${s.units_processed || s.units.length}, Format A events: ${s.format_a_total_events || 0}${unitSummaries ? ' | ' + unitSummaries : ''}`)}</div>`;
    }

    const parts = Object.entries(s)
      .filter(([, v]) => v != null && v !== '' && typeof v !== 'object')
      .slice(0, 8)
      .map(([k, v]) => `${k}: ${typeof v === 'number' ? C.formatSci(v, 4) : String(v)}`);
    return parts.length ? `<div class="job-summary">${esc(parts.join(', '))}</div>` : '';
  }

  function renderJobOutputButtons(outputs) {
    const groups = [
      { key: 'lightcurve', title: 'Light curves', items: [] },
      { key: 'spectrum', title: 'Spectra', items: [] },
      { key: 'image', title: '2D images', items: [] },
      { key: 'events', title: 'Event files', items: [] },
      { key: 'aux', title: 'Responses & auxiliary', items: [] },
    ];
    for (const o of outputs) {
      const prod = o.product || 'other';
      const bucket = groups.find((g) => g.key === prod) || groups[4];
      bucket.items.push(o);
    }
    const nonEmpty = groups.filter((g) => g.items.length > 0);
    return nonEmpty.map((g) => {
      const links = g.items.map((o) => {
        const prod = o.product || 'other';
        const name = o.label || basename(o.path);
        const primary = (prod === 'lightcurve' || prod === 'spectrum')
          ? `<button type="button" class="text-action job-open" data-act="${prod}" data-path="${esc(o.path)}">Plot ${esc(basename(o.path))}</button>`
          : (prod === 'image'
            ? `<button type="button" class="text-action job-open" data-act="image" data-path="${esc(o.path)}">Image ${esc(basename(o.path))}</button>`
            : (prod === 'events'
              ? `<button type="button" class="text-action job-open" data-act="pulse" data-path="${esc(o.path)}">Fold ${esc(basename(o.path))}</button>`
              : ''));
        const inspLabel = primary ? 'Inspect' : `Inspect ${esc(name)} (${esc(basename(o.path))})`;
        return `${primary}<button type="button" class="text-action job-open" data-act="inspect" data-path="${esc(o.path)}">${inspLabel}</button>`;
      }).join('');
      return `<div class="job-output-group">
        <div class="job-output-group-title">${esc(g.title)} (${g.items.length})</div>
        <div class="job-links">${links}</div>
      </div>`;
    }).join('');
  }

  function renderJobs(jobs) {
    const box = $('jobs-list');
    if (!jobs.length) {
      box.innerHTML = '<p class="empty-note">No jobs started yet. Choose a pipeline or tool on the left and run it to inspect generated FITS products.</p>';
      return;
    }
    box.innerHTML = jobs.map((j, idx) => {
      const stCls = j.status === 'succeeded' ? 'status-ok' : (j.status === 'failed' ? 'status-err' : 'status-warn');
      const elapsed = Number.isFinite(j.elapsed_ms) ? `${(j.elapsed_ms / 1000).toFixed(2)} s` : 'Running';
      const summaryHtml = renderJobSummaryHtml(j);
      const outputs = j.outputs || [];
      let outputsBlock = '';
      if (outputs.length > 0) {
        const groupedHtml = renderJobOutputButtons(outputs);
        if (outputs.length > 8) {
          const jobId = String(j.id);
          const isOpen = state.expandedJobOutputs.has(jobId)
            || (idx === 0 && !state.collapsedJobOutputs.has(jobId));
          outputsBlock = `<details class="job-outputs-details" data-job-id="${esc(jobId)}" ${isOpen ? 'open' : ''}>
            <summary class="job-outputs-summary">${outputs.length} generated FITS products</summary>
            <div class="job-outputs-wrap">${groupedHtml}</div>
          </details>`;
        } else {
          outputsBlock = `<div class="job-outputs-wrap">${groupedHtml}</div>`;
        }
      }
      return `<div class="job-row">
        <div class="job-top">
          <div class="job-title">
            <span class="${stCls}">${esc(j.status)}</span>
            <strong>${esc(j.pipeline.toUpperCase())}</strong>
            <span class="meta-muted">${esc(j.id)}</span>
          </div>
          <span class="meta-muted">${elapsed}</span>
        </div>
        ${summaryHtml}
        ${j.error ? `<div class="status-err">${esc(j.error)}</div>` : ''}
        ${outputsBlock}
      </div>`;
    }).join('');

    box.querySelectorAll('.job-open').forEach((btn) => {
      btn.onclick = () => openProductInTab(btn.dataset.act, btn.dataset.path);
    });
    box.querySelectorAll('.job-outputs-details').forEach((det) => {
      det.addEventListener('toggle', () => {
        const id = String(det.dataset.jobId || '');
        if (!id) return;
        if (det.open) {
          state.expandedJobOutputs.add(id);
          state.collapsedJobOutputs.delete(id);
        } else {
          state.expandedJobOutputs.delete(id);
          state.collapsedJobOutputs.add(id);
        }
      });
    });
  }

  // =========================================================================
  // Two-Way TOML Config Export & Load/Import
  // =========================================================================
  function downloadTomlConfig() {
    const text = $('preview-toml').textContent || '';
    if (!text.trim()) return;
    const blob = new Blob([text], { type: 'application/toml;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `astrosat_${state.runnerPipe}.toml`;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => {
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }, 150);
  }

  function openTomlModal() {
    state.tomlReturnFocus = document.activeElement instanceof HTMLElement ? document.activeElement : null;
    $('toml-modal-err').classList.add('hidden');
    $('toml-modal-err').textContent = '';
    $('toml-modal').classList.remove('hidden');
    requestAnimationFrame(() => $('toml-paste-input').focus());
  }

  function closeTomlModal() {
    $('toml-modal').classList.add('hidden');
    if (state.tomlReturnFocus && typeof state.tomlReturnFocus.focus === 'function') {
      state.tomlReturnFocus.focus();
    }
    state.tomlReturnFocus = null;
  }

  function setSelectChecked(selectEl, rawValue, fieldName, warnings) {
    if (!selectEl || rawValue == null) return false;
    const target = String(rawValue).trim();
    const match = Array.from(selectEl.options).find(
      (o) => o.value === target || o.value.toLowerCase() === target.toLowerCase()
    );
    if (match) {
      selectEl.value = match.value;
      return true;
    }
    if (Array.isArray(warnings)) {
      warnings.push(`Unsupported ${fieldName} value "${target}"`);
    }
    return false;
  }

  function setSelectValueWithDynamicOption(selectEl, value) {
    const strVal = String(value);
    const exists = Array.from(selectEl.options).some((o) => o.value === strVal);
    if (!exists) {
      const opt = document.createElement('option');
      opt.value = strVal;
      opt.textContent = strVal;
      selectEl.appendChild(opt);
    }
    selectEl.value = strVal;
  }

  function resetRunnerAndTimingFields() {
    $('runner-threads').value = '';

    // CZTI defaults
    $('czti-input-dir').value = '';
    $('czti-m0').value = '';
    $('czti-ss').value = '';
    $('czti-tct').value = '';
    $('czti-mkf').value = '';
    $('czti-caldb-dir').value = '';
    $('czti-output-dir').value = '';
    $('czti-prefix').value = 'czti_l2';
    $('czti-mkf-thresh').value = '';
    $('czti-user-gti').value = '';
    $('czti-entry-stage').value = '1';
    $('czti-exit-stage').value = '3';
    $('czti-gti-type').value = 'quad';
    $('czti-algo').value = 'default';
    $('czti-ra').value = '';
    $('czti-dec').value = '';
    $('czti-tbin').value = '1.0';
    $('czti-emin').value = '20';
    $('czti-emax').value = '200';
    $('czti-output-type').value = 'both';
    $('czti-mask-weight').checked = true;
    $('czti-gen-resp').checked = true;
    $('czti-gen-evt').checked = false;
    $('czti-img-algo').value = 'balcc';
    $('czti-oversampling').value = '';
    state.cztiQuadsExplicit = false;
    document.querySelectorAll('.czti-quad-cb').forEach((cb) => {
      cb.checked = true;
    });
    $('czti-advanced').open = false;

    // LAXPC defaults
    $('laxpc-input-dir').value = '';
    $('laxpc-caldb-dir').value = '';
    $('laxpc-output-dir').value = '';
    $('laxpc-units').value = '1,2,3';
    $('laxpc-tbin').value = '1.0';
    $('laxpc-deadtime').value = '';
    $('laxpc-channels').value = '512';
    $('laxpc-min-elv').value = '0';
    $('laxpc-fmt-a').checked = true;
    $('laxpc-fmt-b').checked = true;
    $('laxpc-orbit').checked = true;
    $('laxpc-single').checked = false;
    $('laxpc-layer1').checked = false;
    $('laxpc-anodes').value = '';
    $('laxpc-energy-range').value = '';
    $('laxpc-pha-min').value = '';
    $('laxpc-pha-max').value = '';
    $('laxpc-saa-entry').value = '';
    $('laxpc-saa-exit').value = '';
    $('laxpc-merge-gap').value = '';
    $('laxpc-min-dur').value = '';
    $('laxpc-user-gti').value = '';
    $('laxpc-user-gti-only').checked = false;
    $('laxpc-ascii').checked = false;
    $('laxpc-bkg-model').value = '';
    $('laxpc-bkg-pha').value = '';
    $('laxpc-bkg-peak').value = '';
    $('laxpc-advanced').open = false;

    // Barycenter defaults
    $('bary-input').value = '';
    $('bary-output').value = '';
    $('bary-orbit').value = '';
    $('bary-mkf').value = '';
    $('bary-ephem').value = '';
    $('bary-ra').value = '83.633083';
    $('bary-dec').value = '22.014500';
    $('bary-frame').value = 'ICRS';
    $('bary-checksums').checked = true;
    $('bary-advanced').open = false;

    // Timing defaults
    $('pulse-freq').value = '29.658051471';
    $('pulse-period').value = (1 / 29.658051471).toFixed(12).replace(/\.?0+$/, '');
    $('pulse-fdot').value = '0';
    $('pulse-fddot').value = '0';
    $('pulse-epoch').value = '';
    $('pulse-nbins').value = '32';
    $('pulse-harmonics').value = '2';
    $('pulse-ephem-advanced').open = false;
  }

  function hydrateFromParsedConfig(res) {
    const cfg = res.config || {};
    const sections = res.sections || [];
    const warnings = [];
    resetRunnerAndTimingFields();
    if (cfg.threads != null) $('runner-threads').value = String(cfg.threads);

    if (sections.includes('czti') && cfg.czti) {
      const c = cfg.czti;
      const setIf = (id, v) => { if (v != null) $(id).value = String(v); };
      setIf('czti-input-dir', c.input_dir);
      setIf('czti-m0', c.l1_m0_file);
      setIf('czti-ss', c.l1_ss_file);
      setIf('czti-tct', c.tct_file);
      setIf('czti-mkf', c.mkf_file);
      setIf('czti-caldb-dir', c.caldb_dir);
      setIf('czti-output-dir', c.output_dir);
      setIf('czti-prefix', c.output_prefix);
      setIf('czti-mkf-thresh', c.mkf_threshold_file);
      setIf('czti-user-gti', c.user_gti_file);
      if (c.entry_stage != null) setSelectChecked($('czti-entry-stage'), c.entry_stage, 'czti.entry_stage', warnings);
      if (c.exit_stage != null) setSelectChecked($('czti-exit-stage'), c.exit_stage, 'czti.exit_stage', warnings);
      if (c.gti_type != null) setSelectChecked($('czti-gti-type'), c.gti_type, 'czti.gti_type', warnings);
      if (c.algo != null) setSelectChecked($('czti-algo'), c.algo, 'czti.algo', warnings);
      setIf('czti-ra', c.ra_deg);
      setIf('czti-dec', c.dec_deg);
      setIf('czti-tbin', c.time_binsize);
      setIf('czti-emin', c.emin_kev);
      setIf('czti-emax', c.emax_kev);
      if (c.output_type != null) {
        const normOut = {
          spectrum: 'spec',
          spec: 'spec',
          lightcurve: 'lc',
          lc: 'lc',
          all: 'both',
          both: 'both',
        }[String(c.output_type).trim().toLowerCase()] || c.output_type;
        setSelectChecked($('czti-output-type'), normOut, 'czti.output_type', warnings);
      }
      if (c.mask_weight != null) $('czti-mask-weight').checked = Boolean(c.mask_weight);
      if (c.gen_resp != null) $('czti-gen-resp').checked = Boolean(c.gen_resp);
      if (c.generate_eventfile != null) $('czti-gen-evt').checked = Boolean(c.generate_eventfile);
      if (c.image_algo != null) setSelectChecked($('czti-img-algo'), c.image_algo, 'czti.image_algo', warnings);
      setIf('czti-oversampling', c.oversampling);
      if (Array.isArray(c.quads)) {
        state.cztiQuadsExplicit = true;
        document.querySelectorAll('.czti-quad-cb').forEach((cb) => {
          cb.checked = c.quads.includes(Number(cb.value));
        });
      }
      const hasAdv = [
        c.l1_m0_file, c.l1_ss_file, c.tct_file, c.mkf_file, c.mkf_threshold_file,
        c.user_gti_file, c.output_type, c.generate_eventfile, c.oversampling, c.quads,
      ].some((v) => v != null);
      $('czti-advanced').open = hasAdv;
    }

    if (sections.includes('laxpc') && cfg.laxpc) {
      const l = cfg.laxpc;
      const setIf = (id, v) => { if (v != null) $(id).value = String(v); };
      setIf('laxpc-input-dir', l.input_dir);
      setIf('laxpc-caldb-dir', l.caldb_dir);
      setIf('laxpc-output-dir', l.output_dir);
      if (Array.isArray(l.units)) $('laxpc-units').value = l.units.join(',');
      setIf('laxpc-tbin', l.tbin);
      setIf('laxpc-deadtime', l.deadtime_sec);
      if (l.spectrum_channels != null) setSelectValueWithDynamicOption($('laxpc-channels'), l.spectrum_channels);
      setIf('laxpc-min-elv', l.min_elv_deg);
      if (l.write_format_a_combined != null) $('laxpc-fmt-a').checked = Boolean(l.write_format_a_combined);
      if (l.write_format_b_events != null) $('laxpc-fmt-b').checked = Boolean(l.write_format_b_events);
      if (l.write_orbit_file != null) $('laxpc-orbit').checked = Boolean(l.write_orbit_file);
      if (l.single_only != null) $('laxpc-single').checked = Boolean(l.single_only);
      if (l.layer1_only != null) $('laxpc-layer1').checked = Boolean(l.layer1_only);
      setIf('laxpc-anodes', l.anodes);
      setIf('laxpc-energy-range', l.energy_range);
      setIf('laxpc-pha-min', l.pha_min);
      setIf('laxpc-pha-max', l.pha_max);
      setIf('laxpc-saa-entry', l.saa_entry_lon_deg);
      setIf('laxpc-saa-exit', l.saa_exit_lon_deg);
      setIf('laxpc-merge-gap', l.merge_gap_sec);
      setIf('laxpc-min-dur', l.min_duration_sec);
      setIf('laxpc-user-gti', l.user_gti_file);
      if (l.user_gti_only != null) $('laxpc-user-gti-only').checked = Boolean(l.user_gti_only);
      if (l.write_ascii_products != null) $('laxpc-ascii').checked = Boolean(l.write_ascii_products);
      setIf('laxpc-bkg-model', l.background_model_path);
      setIf('laxpc-bkg-pha', l.background_pha_path);
      setIf('laxpc-bkg-peak', l.background_peak_ch);
      const hasAdv = [
        l.deadtime_sec, l.energy_range, l.pha_min, l.pha_max, l.anodes,
        l.saa_entry_lon_deg, l.saa_exit_lon_deg, l.merge_gap_sec, l.min_duration_sec,
        l.user_gti_file, l.user_gti_only, l.write_ascii_products,
        l.background_model_path, l.background_pha_path, l.background_peak_ch,
      ].some((v) => v != null);
      $('laxpc-advanced').open = hasAdv;
    }

    if (sections.includes('bary') && cfg.bary) {
      const b = cfg.bary;
      if (b.orbit_file != null) {
        $('bary-orbit').value = b.orbit_file;
        $('bary-mkf').value = '';
      }
      if (b.ra_deg != null) $('bary-ra').value = String(b.ra_deg);
      if (b.dec_deg != null) $('bary-dec').value = String(b.dec_deg);
      if (b.frame != null) setSelectChecked($('bary-frame'), b.frame, 'bary.frame', warnings);
      if (b.ephemeris != null) $('bary-ephem').value = b.ephemeris;
      if (b.write_checksums != null) $('bary-checksums').checked = Boolean(b.write_checksums);
      $('bary-advanced').open = false;
    }

    if (sections.includes('timing') && cfg.timing) {
      const t = cfg.timing;
      if (t.freq_hz != null) {
        $('pulse-freq').value = String(t.freq_hz);
        if (t.freq_hz > 0) $('pulse-period').value = (1 / t.freq_hz).toFixed(12).replace(/\.?0+$/, '');
      }
      if (t.fdot != null) $('pulse-fdot').value = String(t.fdot);
      if (t.fddot != null) {
        $('pulse-fddot').value = String(t.fddot);
        $('pulse-ephem-advanced').open = true;
      }
      if (t.harmonics != null) setSelectValueWithDynamicOption($('pulse-harmonics'), t.harmonics);
      if (t.nbins != null) setSelectValueWithDynamicOption($('pulse-nbins'), t.nbins);
    }

    const firstPipe = sections.find((s) => ['czti', 'laxpc', 'bary'].includes(s));
    if (firstPipe) setRunnerPipe(firstPipe);
    else updateCliPreview();

    const isTimingOnly = !firstPipe && sections.includes('timing');
    if (isTimingOnly) {
      switchTab('pulse');
      state.tomlReturnFocus = $('pulse-freq');
    }
    const baseSummary = isTimingOnly
      ? 'Loaded timing configuration into Pulse alignment'
      : (sections.length ? `Loaded TOML sections: ${sections.join(', ')}` : 'Loaded TOML configuration');
    const summary = warnings.length ? `${baseSummary} (warnings: ${warnings.join('; ')})` : baseSummary;
    const statusEl = $('runner-inline-status');
    if (statusEl) statusEl.textContent = summary;
    if (warnings.length) {
      showBanner('runner', warnings.join('; '));
      showToast(summary, { isError: true });
    } else {
      hideBanner('runner');
      showToast(summary);
    }
  }

  async function applyTomlImport() {
    const raw = $('toml-paste-input').value;
    const errBox = $('toml-modal-err');
    errBox.classList.add('hidden');
    if (!raw.trim()) {
      errBox.textContent = 'Paste TOML configuration or choose a .toml file first.';
      errBox.classList.remove('hidden');
      return;
    }
    const btn = $('toml-apply-btn');
    btn.disabled = true;
    try {
      const res = await apiPost('/api/config/parse', { content: raw, format: 'auto' });
      hydrateFromParsedConfig(res);
      closeTomlModal();
    } catch (err) {
      errBox.textContent = err.message;
      errBox.classList.remove('hidden');
    } finally {
      btn.disabled = false;
    }
  }

  function openProductInTab(action, path) {
    if (action === 'lightcurve' || action === 'spectrum') {
      switchTab('lc-spec');
      setLcMode(action);
      $('lc-path-input').value = path;
      addLcOrSpecFile(path);
    } else if (action === 'image') {
      switchTab('imaging');
      $('img-path-input').value = path;
      loadImageFile(path);
    } else if (action === 'pulse') {
      switchTab('pulse');
      state.pulseInputs = [{
        label: basename(path),
        path,
        caldb: (state.health && state.health.presets && state.health.presets.laxpc_caldb) || '',
        use_weights: false,
        bands: [{ label: 'Full band', emin_kev: '', emax_kev: '', unit: '' }],
      }];
      renderPulseInputs();
    } else {
      switchTab('inspector');
      $('insp-path-input').value = path;
      inspectFitsFile(path, 0, 0);
    }
  }

  // =========================================================================
  // VIEW 2: LIGHT CURVES & SPECTRA
  // =========================================================================
  function setLcMode(mode) {
    if (state.lcMode !== mode) {
      state.lcMode = mode;
      state.lcItems = [];
      state.lcPlot.clear();
      $('lc-empty').classList.remove('hidden');
      renderLcSeriesList();
      redrawLcPlot();
    }
    document.querySelectorAll('#lc-mode-seg button').forEach((b) => {
      b.classList.toggle('active', b.dataset.mode === mode);
    });
    $('lc-opts-box').classList.toggle('hidden', mode !== 'lightcurve');
    $('spec-opts-box').classList.toggle('hidden', mode !== 'spectrum');
    $('lc-gti-wrap').classList.toggle('hidden', mode !== 'lightcurve');
    if (mode === 'spectrum') {
      $('lc-log-y').checked = true;
      $('lc-step').checked = true;
      state.lcPlot.setOptions({ logY: true, stepMode: true });
    } else {
      $('lc-log-y').checked = false;
      $('lc-step').checked = false;
      state.lcPlot.setOptions({ logY: false, stepMode: false });
    }
    populateLcSampleDropdown();
  }

  function populateLcSampleDropdown() {
    const sel = $('lc-sample-select');
    sel.innerHTML = '<option value="">Select detected product...</option>';
    const p = state.health && state.health.presets;
    if (!p) return;
    const list = state.lcMode === 'lightcurve' ? (p.sample_lightcurves || []) : (p.sample_spectra || []);
    for (const f of list) {
      const opt = document.createElement('option');
      opt.value = f;
      opt.textContent = basename(f);
      opt.title = f;
      sel.appendChild(opt);
    }
  }

  function cleanTitle(t, fallback) {
    return String(t || fallback || '').replace(/\s*\u2014\s*/g, ' / ');
  }

  async function addLcOrSpecFile(pathOverride = null) {
    const path = (pathOverride != null ? pathOverride : $('lc-path-input').value).trim();
    if (!path) return;
    hideBanner('lc');
    $('lc-loading').classList.remove('hidden');
    try {
      const existingIdx = state.lcItems.findIndex((it) => it.path === path);
      const pal = C.getPalette(state.theme);
      const color = existingIdx >= 0 ? state.lcItems[existingIdx].color : pal[state.lcItems.length % pal.length];
      let raw;
      if (state.lcMode === 'lightcurve') {
        const rebin = Number($('lc-rebin').value) || 1;
        raw = await apiGet(`/api/lightcurve?path=${encodeURIComponent(path)}&rebin=${rebin}`);
      } else {
        const goodOnly = $('spec-good-only').checked;
        const caldb = $('spec-caldb-input').value.trim();
        let url = `/api/spectrum?path=${encodeURIComponent(path)}&good_only=${goodOnly}`;
        if (caldb) url += `&caldb=${encodeURIComponent(caldb)}`;
        raw = await apiGet(url);
      }
      const entry = {
        id: existingIdx >= 0 ? state.lcItems[existingIdx].id : `${Date.now()}_${Math.random().toString(36).slice(2, 6)}`,
        path,
        label: cleanTitle(raw.title, basename(path)),
        color,
        hidden: false,
        raw,
      };
      if (existingIdx >= 0) state.lcItems[existingIdx] = entry;
      else state.lcItems.push(entry);
      renderLcSeriesList();
      redrawLcPlot();
    } catch (err) {
      showBanner('lc', err.message);
    } finally {
      $('lc-loading').classList.add('hidden');
    }
  }

  async function reloadAllLcItems() {
    if (!state.lcItems.length) return;
    const prev = state.lcItems.map((it) => ({ path: it.path, color: it.color, hidden: it.hidden }));
    state.lcItems = [];
    $('lc-loading').classList.remove('hidden');
    hideBanner('lc');
    try {
      for (const item of prev) {
        let raw;
        if (state.lcMode === 'lightcurve') {
          const rebin = Number($('lc-rebin').value) || 1;
          raw = await apiGet(`/api/lightcurve?path=${encodeURIComponent(item.path)}&rebin=${rebin}`);
        } else {
          const goodOnly = $('spec-good-only').checked;
          const caldb = $('spec-caldb-input').value.trim();
          let url = `/api/spectrum?path=${encodeURIComponent(item.path)}&good_only=${goodOnly}`;
          if (caldb) url += `&caldb=${encodeURIComponent(caldb)}`;
          raw = await apiGet(url);
        }
        state.lcItems.push({
          id: `${Date.now()}_${Math.random().toString(36).slice(2, 6)}`,
          path: item.path,
          label: cleanTitle(raw.title, basename(item.path)),
          color: item.color,
          hidden: item.hidden,
          raw,
        });
      }
      renderLcSeriesList();
      redrawLcPlot();
    } catch (err) {
      showBanner('lc', err.message);
    } finally {
      $('lc-loading').classList.add('hidden');
    }
  }

  function renderLcSeriesList() {
    const box = $('lc-series-list');
    if (!state.lcItems.length) {
      box.innerHTML = '<p class="empty-note">No files loaded yet.</p>';
      return;
    }
    box.innerHTML = state.lcItems.map((it, idx) => {
      const r = it.raw;
      const meta = state.lcMode === 'lightcurve'
        ? `${(r.time || []).length} pts, ${C.formatSci(r.mean_rate, 3)} c/s`
        : `${r.good_channels || (r.channel || []).length} ch, ${C.formatSci(r.total_rate, 3)} c/s`;
      return `<div class="series-item">
        <div class="series-left">
          <span class="swatch" style="background:${it.color};opacity:${it.hidden ? 0.3 : 1};margin-right:0;"></span>
          <div style="min-width:0;">
            <div class="series-name" title="${esc(it.path)}">${esc(basename(it.path))}</div>
            <div class="series-meta">${esc(meta)}</div>
          </div>
        </div>
        <div class="inline-actions">
          <button type="button" class="text-action lc-toggle-vis" data-idx="${idx}">${it.hidden ? 'Show' : 'Hide'}</button>
          <button type="button" class="text-action lc-remove" data-idx="${idx}">Remove</button>
        </div>
      </div>`;
    }).join('');
    box.querySelectorAll('.lc-toggle-vis').forEach((b) => {
      b.onclick = () => {
        const it = state.lcItems[Number(b.dataset.idx)];
        it.hidden = !it.hidden;
        renderLcSeriesList();
        redrawLcPlot();
      };
    });
    box.querySelectorAll('.lc-remove').forEach((b) => {
      b.onclick = () => {
        state.lcItems.splice(Number(b.dataset.idx), 1);
        renderLcSeriesList();
        redrawLcPlot();
      };
    });
  }

  function redrawLcPlot() {
    const items = state.lcItems;
    if (!items.length) {
      state.lcPlot.clear();
      $('lc-empty').classList.remove('hidden');
      $('lc-stat-count').textContent = '0';
      $('lc-stat-exp').textContent = 'N/A';
      $('lc-stat-rate').textContent = 'N/A';
      $('lc-stat-meta').textContent = 'N/A';
      return;
    }
    $('lc-empty').classList.add('hidden');
    const first = items[0].raw;

    if (state.lcMode === 'lightcurve') {
      const rel = $('lc-rel-time').checked;
      const t0 = Number.isFinite(first.tstart) ? first.tstart : ((first.time && first.time[0]) || 0);
      const series = items.map((it) => ({
        id: it.id,
        label: basename(it.path),
        color: it.color,
        hidden: it.hidden,
        x: (it.raw.time || []).map((t) => (rel ? t - t0 : t)),
        y: it.raw.rate || [],
        yerr: it.raw.error || [],
      }));
      const gtis = (first.gtis || []).map((g) => [rel ? g.start - t0 : g.start, rel ? g.stop - t0 : g.stop]);
      state.lcPlot.setData({
        series,
        gtis,
        title: items.length === 1 ? cleanTitle(first.title, basename(items[0].path)) : `Overlay (${items.length} light curves)`,
        xLabel: rel ? `Time - ${t0.toFixed(1)}` : `Time (${first.timesys || 'MET'})`,
        xUnit: 's',
        yLabel: 'Count rate',
        yUnit: first.units || 'counts/s',
      });
      $('lc-stat-count').textContent = String(items.length);
      $('lc-stat-exp').textContent = Number.isFinite(first.exposure) ? `${first.exposure.toFixed(2)} s` : 'N/A';
      $('lc-stat-rate').textContent = Number.isFinite(first.mean_rate) ? `${C.formatSci(first.mean_rate, 4)} c/s` : 'N/A';
      $('lc-stat-meta').textContent = Number.isFinite(first.effective_timedel)
        ? `${C.formatSci(first.effective_timedel, 3)} s (${first.effective_rebin || 1}x)`
        : 'N/A';
    } else {
      const useEnergy = $('spec-xaxis').value === 'energy' && first.energy_mapped;
      const series = items.map((it) => {
        const r = it.raw;
        const xArr = (useEnergy && r.energy_kev && r.energy_kev.length) ? r.energy_kev : r.channel;
        return {
          id: it.id,
          label: basename(it.path),
          color: it.color,
          hidden: it.hidden,
          x: xArr || [],
          y: r.rate_or_counts || [],
          yerr: r.stat_err || [],
        };
      });
      state.lcPlot.setData({
        series,
        gtis: [],
        title: items.length === 1 ? cleanTitle(first.title, basename(items[0].path)) : `Overlay (${items.length} spectra)`,
        xLabel: useEnergy ? 'Energy' : (first.chantype || 'Channel'),
        xUnit: useEnergy ? 'keV' : '',
        yLabel: first.is_rate ? 'Count rate' : 'Counts',
        yUnit: first.units || (first.is_rate ? 'counts/s/ch' : 'counts/ch'),
      });
      $('lc-stat-count').textContent = String(items.length);
      $('lc-stat-exp').textContent = Number.isFinite(first.exposure) ? `${first.exposure.toFixed(2)} s` : 'N/A';
      $('lc-stat-rate').textContent = Number.isFinite(first.total_rate) ? `${C.formatSci(first.total_rate, 4)} c/s` : 'N/A';
      $('lc-stat-meta').textContent = `${first.energy_source || 'channel'}, ${first.good_channels}/${first.total_channels} good`;
    }
  }

  // =========================================================================
  // VIEW 3: 2D IMAGING
  // =========================================================================
  async function loadImageFile(pathOverride = null, hduOverride = null) {
    const path = (pathOverride != null ? pathOverride : $('img-path-input').value).trim();
    if (!path) return;
    hideBanner('img');
    $('img-loading').classList.remove('hidden');
    try {
      let url = `/api/image?path=${encodeURIComponent(path)}`;
      if (hduOverride !== null && hduOverride !== '') url += `&hdu=${encodeURIComponent(hduOverride)}`;
      const data = await apiGet(url);
      state.imgData = data;
      $('img-empty').classList.add('hidden');

      const hduSel = $('img-hdu-select');
      const avail = data.available_image_hdus || [];
      hduSel.innerHTML = avail.map((h) =>
        `<option value="${h.index}" ${h.index === data.hdu_index ? 'selected' : ''}>#${h.index} ${esc(h.extname || 'PRIMARY')} (${h.width}x${h.height})</option>`
      ).join('') || `<option value="${data.hdu_index}">#${data.hdu_index} ${esc(data.hdu || 'IMAGE')}</option>`;

      const quadDefault = Boolean(data.quad_overlay || (data.width === 128 && data.height === 128));
      $('img-quads').checked = quadDefault;
      state.imgPlot.setOptions({ showQuads: quadDefault });

      state.imgPlot.setImage({
        width: data.width,
        height: data.height,
        data: data.pixels,
        stats: data.stats,
        wcs: data.wcs && data.wcs.has_wcs ? data.wcs : null,
        peak: data.stats && data.stats.peak ? data.stats.peak : null,
        title: `${basename(data.path)} [${data.hdu || ('HDU ' + data.hdu_index)}]`,
      });

      const st = data.stats || {};
      const pk = st.peak || {};
      $('img-stat-dim').textContent = `${data.width} x ${data.height}`;
      $('img-stat-snr').textContent = Number.isFinite(pk.snr) ? `${pk.snr.toFixed(2)} sigma` : 'N/A';
      $('img-stat-peak').textContent = Number.isFinite(pk.x) ? `(${pk.x}, ${pk.y})` : 'N/A';
      $('img-stat-max').textContent = C.formatSci(pk.val != null ? pk.val : st.max, 4);
      $('img-stat-minmed').textContent = `${C.formatSci(st.min, 3)} / ${C.formatSci(st.p50, 3)}`;
      $('img-stat-meanstd').textContent = `${C.formatSci(st.mean, 3)} +/- ${C.formatSci(st.std, 3)}`;

      $('img-wcs-label').textContent = data.wcs && data.wcs.has_wcs
        ? `WCS: ${data.wcs.ctype1 || 'X'} / ${data.wcs.ctype2 || 'Y'}`
        : 'Detector pixel coordinates';
    } catch (err) {
      showBanner('img', err.message);
    } finally {
      $('img-loading').classList.add('hidden');
    }
  }

  function handleImageHover(hit) {
    if (!hit) {
      $('img-cur-pix').textContent = 'N/A';
      $('img-cur-val').textContent = 'N/A';
      $('img-cur-wcs1').textContent = 'N/A';
      $('img-cur-wcs2').textContent = 'N/A';
      return;
    }
    $('img-cur-pix').textContent = `(${hit.fx}, ${hit.fy})`;
    $('img-cur-val').textContent = C.formatSci(hit.val, 5);
    if (hit.world) {
      $('img-cur-wcs1-lbl').textContent = hit.world.labelX;
      $('img-cur-wcs2-lbl').textContent = hit.world.labelY;
      $('img-cur-wcs1').textContent = `${hit.world.x.toFixed(4)} deg`;
      $('img-cur-wcs2').textContent = `${hit.world.y.toFixed(4)} deg`;
    } else {
      $('img-cur-wcs1-lbl').textContent = 'FITS 1-based X';
      $('img-cur-wcs2-lbl').textContent = 'FITS 1-based Y';
      $('img-cur-wcs1').textContent = String(hit.fx + 1);
      $('img-cur-wcs2').textContent = String(hit.fy + 1);
    }
  }

  // =========================================================================
  // VIEW 4: PHASE-FOLDED & MULTI-BAND PULSE ALIGNMENT
  // =========================================================================
  function findSampleEvent(kind) {
    const evs = (state.health && state.health.presets && state.health.presets.sample_events) || [];
    if (kind === 'laxpc') {
      return evs.find((p) => /laxpc|lxp|level2\.event/i.test(p) && /bary/i.test(p))
        || evs.find((p) => /laxpc|lxp|level2\.event/i.test(p))
        || '';
    }
    if (kind === 'czti') {
      return evs.find((p) => /czt/i.test(p) && /bary/i.test(p))
        || evs.find((p) => /czt/i.test(p))
        || '';
    }
    return evs[0] || '';
  }

  function applyPulsePreset(kind) {
    state.pulsePreset = kind;
    document.querySelectorAll('#pulse-preset-seg button').forEach((b) => {
      b.classList.toggle('active', b.dataset.preset === kind);
    });
    const p = (state.health && state.health.presets) || {};
    const laxpcPath = findSampleEvent('laxpc');
    const cztiPath = findSampleEvent('czti');
    const laxpcCaldb = p.laxpc_caldb || '';

    if (kind === 'laxpc4') {
      state.pulseInputs = [{
        label: 'LAXPC',
        path: laxpcPath,
        caldb: laxpcCaldb,
        use_weights: false,
        bands: [
          { label: '3-10 keV', emin_kev: 3, emax_kev: 10, unit: '' },
          { label: '10-20 keV', emin_kev: 10, emax_kev: 20, unit: '' },
          { label: '20-40 keV', emin_kev: 20, emax_kev: 40, unit: '' },
          { label: '40-80 keV', emin_kev: 40, emax_kev: 80, unit: '' },
        ],
      }];
    } else if (kind === 'czti3') {
      state.pulseInputs = [{
        label: 'CZTI',
        path: cztiPath,
        caldb: '',
        use_weights: false,
        bands: [
          { label: '20-60 keV', emin_kev: 20, emax_kev: 60, unit: '' },
          { label: '60-120 keV', emin_kev: 60, emax_kev: 120, unit: '' },
          { label: '120-200 keV', emin_kev: 120, emax_kev: 200, unit: '' },
        ],
      }];
    } else if (kind === 'cross') {
      state.pulseInputs = [
        {
          label: 'LAXPC',
          path: laxpcPath,
          caldb: laxpcCaldb,
          use_weights: false,
          bands: [
            { label: '3-10 keV', emin_kev: 3, emax_kev: 10, unit: '' },
            { label: '10-20 keV', emin_kev: 10, emax_kev: 20, unit: '' },
            { label: '20-40 keV', emin_kev: 20, emax_kev: 40, unit: '' },
            { label: '40-80 keV', emin_kev: 40, emax_kev: 80, unit: '' },
          ],
        },
        {
          label: 'CZTI',
          path: cztiPath,
          caldb: '',
          use_weights: false,
          bands: [
            { label: '20-60 keV', emin_kev: 20, emax_kev: 60, unit: '' },
            { label: '60-120 keV', emin_kev: 60, emax_kev: 120, unit: '' },
            { label: '120-200 keV', emin_kev: 120, emax_kev: 200, unit: '' },
          ],
        },
      ];
    }
    renderPulseInputs();
  }

  function renderPulseInputs() {
    const box = $('pulse-inputs-container');
    const pal = C.getPalette(state.theme);
    let colorIdx = 0;

    if (!state.pulseInputs.length) {
      box.innerHTML = '<p class="empty-note">Select a preset above or add an event file.</p>';
      return;
    }

    box.innerHTML = state.pulseInputs.map((inp, ii) => {
      const bandsHtml = inp.bands.map((b, bi) => {
        const col = pal[colorIdx++ % pal.length];
        return `<div class="band-line" data-ii="${ii}" data-bi="${bi}">
          <span class="swatch" style="background:${col};margin-right:0;"></span>
          <input type="text" class="b-label" value="${esc(b.label)}" placeholder="Band" />
          <input type="number" step="any" class="b-emin" value="${esc(b.emin_kev)}" placeholder="Emin" title="Min energy (keV)" />
          <input type="number" step="any" class="b-emax" value="${esc(b.emax_kev)}" placeholder="Emax" title="Max energy (keV)" />
          <input type="number" min="1" max="3" class="b-unit" value="${esc(b.unit)}" placeholder="Unit" title="LAXPC unit (1-3)" />
          <button type="button" class="text-action b-del" title="Remove band">x</button>
        </div>`;
      }).join('');

      return `<div class="pulse-input-group" data-ii="${ii}">
        <div class="section-head">
          <input type="text" class="inp-label" value="${esc(inp.label)}" placeholder="Instrument" style="width:128px;font-weight:600;height:32px;" />
          <div class="inline-actions">
            <label class="check-label"><input type="checkbox" class="inp-weights" ${inp.use_weights ? 'checked' : ''} /> Weights</label>
            <button type="button" class="text-action inp-add-band">Add band</button>
            <button type="button" class="text-action inp-del">Remove</button>
          </div>
        </div>
        <div class="path-field">
          <input type="text" class="input-path inp-path" value="${esc(inp.path)}" placeholder="Barycentered event FITS path" />
          <button type="button" class="subtle-btn inp-browse">Browse</button>
        </div>
        <div class="path-field">
          <input type="text" class="input-path inp-caldb" value="${esc(inp.caldb || '')}" placeholder="Optional CALDB dir (LAXPC keV)" />
          <button type="button" class="subtle-btn inp-browse-caldb">CALDB</button>
        </div>
        <div class="field-stack" style="gap:6px;">${bandsHtml}</div>
      </div>`;
    }).join('');

    box.querySelectorAll('.pulse-input-group').forEach((grp) => {
      const ii = Number(grp.dataset.ii);
      const inp = state.pulseInputs[ii];
      grp.querySelector('.inp-label').oninput = (e) => { inp.label = e.target.value; };
      grp.querySelector('.inp-path').oninput = (e) => { inp.path = e.target.value; };
      grp.querySelector('.inp-caldb').oninput = (e) => { inp.caldb = e.target.value; };
      grp.querySelector('.inp-weights').onchange = (e) => { inp.use_weights = e.target.checked; };
      grp.querySelector('.inp-del').onclick = () => { state.pulseInputs.splice(ii, 1); renderPulseInputs(); };
      grp.querySelector('.inp-add-band').onclick = () => {
        inp.bands.push({ label: `Band ${inp.bands.length + 1}`, emin_kev: '', emax_kev: '', unit: '' });
        renderPulseInputs();
      };
      grp.querySelector('.inp-browse').onclick = () => openFileModal({
        mode: 'file',
        onSelect: (p) => { inp.path = p; renderPulseInputs(); },
      });
      grp.querySelector('.inp-browse-caldb').onclick = () => openFileModal({
        mode: 'dir',
        onSelect: (p) => { inp.caldb = p; renderPulseInputs(); },
      });
      grp.querySelectorAll('.band-line').forEach((row) => {
        const bi = Number(row.dataset.bi);
        const b = inp.bands[bi];
        row.querySelector('.b-label').oninput = (e) => { b.label = e.target.value; };
        row.querySelector('.b-emin').oninput = (e) => { b.emin_kev = e.target.value; };
        row.querySelector('.b-emax').oninput = (e) => { b.emax_kev = e.target.value; };
        row.querySelector('.b-unit').oninput = (e) => { b.unit = e.target.value; };
        row.querySelector('.b-del').onclick = () => { inp.bands.splice(bi, 1); renderPulseInputs(); };
      });
    });
  }

  async function runPeriodSearch() {
    hideBanner('pulse');
    const path = $('search-path').value.trim() || (state.pulseInputs[0] && state.pulseInputs[0].path) || '';
    if (!path) {
      showBanner('pulse', 'Select a barycentered event file before running period search.');
      return;
    }
    const freqMinRaw = $('search-freq-min').value.trim();
    const freqMaxRaw = $('search-freq-max').value.trim();
    const freqMin = numOrNull(freqMinRaw);
    const freqMax = numOrNull(freqMaxRaw);
    const freqHz = numOrNull($('pulse-freq').value);
    const periodS = numOrNull($('pulse-period').value);
    const trialsElVal = $('search-trials').value.trim();
    const trialsRaw = numOrNull(trialsElVal);
    const dfElVal = $('search-df').value.trim();
    const dfHz = numOrNull(dfElVal);

    const hasFreqMin = freqMinRaw !== '' || freqMin != null;
    const hasFreqMax = freqMaxRaw !== '' || freqMax != null;
    if (hasFreqMin !== hasFreqMax) {
      showBanner('pulse', 'Provide both Min frequency and Max frequency (Hz) for a range search, or leave both blank.');
      return;
    }
    if (hasFreqMin && hasFreqMax) {
      if (freqMin == null || freqMax == null || freqMin <= 0 || freqMin >= freqMax) {
        showBanner('pulse', 'Frequency range must satisfy 0 < Min frequency < Max frequency (Hz).');
        return;
      }
    }

    const body = {
      path,
      search_span_hz: numOrNull($('search-span').value) ?? 0.002,
      fdot: numOrNull($('pulse-fdot').value) ?? 0,
      fddot: numOrNull($('pulse-fddot').value) ?? 0,
      harmonics: Number($('pulse-harmonics').value) || 2,
      nbins: Number($('pulse-nbins').value) || 32,
    };
    if (freqMin != null && freqMax != null) {
      body.freq_min = freqMin;
      body.freq_max = freqMax;
    } else if (freqHz != null && freqHz > 0) {
      body.freq_hz = freqHz;
    } else if (periodS != null && periodS > 0) {
      body.period_s = periodS;
    } else {
      showBanner('pulse', 'Enter a valid positive spin frequency (Hz), spin period (s), or frequency range.');
      return;
    }
    if (trialsElVal !== '') {
      if (trialsRaw == null || trialsRaw < 1) {
        showBanner('pulse', 'Number of trials must be >= 1 (or leave blank to use df_hz).');
        return;
      }
      body.n_trials = Math.floor(trialsRaw);
    } else if (dfElVal !== '') {
      if (dfHz == null || dfHz <= 0) {
        showBanner('pulse', 'Frequency step df (Hz) must be finite and > 0 (or leave blank to use default trials).');
        return;
      }
      body.df_hz = dfHz;
    } else {
      body.n_trials = 300;
    }

    const ep = numOrNull($('pulse-epoch').value);
    if (ep != null) body.epoch = ep;
    const gti = $('search-gti').value.trim();
    if (gti) body.gti = gti;
    const hdu = $('search-hdu').value.trim();
    if (hdu) body.hdu = hdu;
    const piMin = numOrNull($('search-pi-min').value);
    if (piMin != null) body.pi_min = Math.floor(piMin);
    const piMax = numOrNull($('search-pi-max').value);
    if (piMax != null) body.pi_max = Math.floor(piMax);
    const unit = numOrNull($('search-unit').value);
    if (unit != null) body.unit = Math.floor(unit);
    const layer = numOrNull($('search-layer').value);
    if (layer != null) body.layer = Math.floor(layer);

    const btn = $('pulse-run-search');
    btn.disabled = true;
    $('search-result-meta').textContent = 'Searching...';
    try {
      const res = await apiPost('/api/pulse/search', body);
      if (Number.isFinite(res.best_freq_hz) && res.best_freq_hz > 0) {
        $('pulse-freq').value = res.best_freq_hz.toFixed(9);
        $('pulse-period').value = (1 / res.best_freq_hz).toFixed(12).replace(/\.?0+$/, '');
      }
      $('search-result-meta').textContent = `Peak ${fmtFixed(res.best_freq_hz, 7)} Hz, Z2 ${C.formatSci(res.max_z2_n, 4)} (${res.elapsed_ms} ms)`;
      $('periodogram-card').classList.remove('hidden');
      $('pgram-summary').textContent = `${fmtFixed(res.best_freq_hz, 7)} Hz, Z2 = ${C.formatSci(res.max_z2_n, 4)} (${fmtFixed(res.detection_sigma, 1)} sigma)`;

      requestAnimationFrame(() => {
        state.pgramPlot.resize();
        const pal = C.getPalette(state.theme);
        state.pgramPlot.setData({
          series: [{
            id: 'z2',
            label: `Z2_${res.harmonics}`,
            color: pal[0],
            x: res.periodogram.freq_hz,
            y: res.periodogram.z2_n,
            yerr: null,
          }],
          vlines: [{ x: res.best_freq_hz, label: `Peak ${res.best_freq_hz.toFixed(6)} Hz` }],
          title: `Z2_${res.harmonics} periodogram (${basename(res.path)}, ${res.num_events.toLocaleString()} events)`,
          xLabel: 'Frequency',
          xUnit: 'Hz',
          yLabel: `Z2_${res.harmonics} power`,
          yUnit: '',
        });
      });
    } catch (err) {
      showBanner('pulse', err.message);
      $('search-result-meta').textContent = 'Search failed';
    } finally {
      btn.disabled = false;
    }
  }

  async function runMultiBandFold() {
    hideBanner('pulse');
    const freq = numOrNull($('pulse-freq').value);
    if (freq == null || freq <= 0) {
      showBanner('pulse', 'Enter a valid positive spin frequency in Hz.');
      return;
    }
    const validInputs = state.pulseInputs
      .filter((inp) => inp.path && inp.path.trim() && inp.bands.length)
      .map((inp) => ({
        label: inp.label.trim() || basename(inp.path),
        path: inp.path.trim(),
        caldb: inp.caldb && inp.caldb.trim() ? inp.caldb.trim() : undefined,
        use_weights: Boolean(inp.use_weights),
        bands: inp.bands.map((b) => {
          const obj = { label: b.label.trim() || 'Band' };
          const emin = numOrNull(b.emin_kev);
          const emax = numOrNull(b.emax_kev);
          const unit = numOrNull(b.unit);
          if (emin != null) obj.emin_kev = emin;
          if (emax != null) obj.emax_kev = emax;
          if (unit != null) obj.unit = Math.floor(unit);
          return obj;
        }),
      }));

    if (!validInputs.length) {
      showBanner('pulse', 'Provide at least one event file with one energy band.');
      return;
    }

    const body = {
      inputs: validInputs,
      freq_hz: freq,
      fdot: numOrNull($('pulse-fdot').value) || 0,
      fddot: numOrNull($('pulse-fddot').value) || 0,
      nbins: Number($('pulse-nbins').value) || 32,
      harmonics: Number($('pulse-harmonics').value) || 2,
    };
    const ep = numOrNull($('pulse-epoch').value);
    if (ep != null) body.epoch = ep;

    $('pulse-loading').classList.remove('hidden');
    try {
      const res = await apiPost('/api/pulse/multiband', body);
      state.pulseLastResult = res;
      renderPulseResult(res);
    } catch (err) {
      showBanner('pulse', err.message);
    } finally {
      $('pulse-loading').classList.add('hidden');
    }
  }

  function renderPulseResult(res) {
    const pal = C.getPalette(state.theme);
    const bands = (res.bands || []).map((b, i) => ({
      label: b.full_label || `${b.input_label}: ${b.band_label}`,
      color: pal[i % pal.length],
      values: b.counts,
      errors: b.errors,
      peakPhase: b.peak_phase,
      z2: b.z2_n,
      chi2Red: b.reduced_chi2,
      pf: b.pulsed_fraction_minmax,
      pfRms: b.pulsed_fraction_rms,
      lagCycles: b.phase_lag_cycles,
      lagMs: b.phase_lag_ms,
      numEvents: b.num_events,
    }));

    $('pulse-empty').classList.toggle('hidden', bands.length > 0);
    let refIdx = -1;
    if (res.reference_band) {
      refIdx = bands.findIndex((b) => b.label === res.reference_band && (b.numEvents || 0) > 0);
    }
    if (refIdx < 0) {
      refIdx = bands.findIndex((b) => (b.numEvents || 0) > 0);
    }
    if (refIdx < 0 && res.reference_band) {
      refIdx = bands.findIndex((b) => b.label === res.reference_band);
    }
    if (refIdx < 0 && bands.length > 0) {
      refIdx = 0;
    }
    const refBand = refIdx >= 0 ? bands[refIdx] : null;
    const refPeak = refBand ? refBand.peakPhase : null;
    const periodMs = Number.isFinite(res.period_sec) ? (res.period_sec * 1000).toFixed(4) : 'N/A';
    const fddotVal = res.fddot ?? numOrNull($('pulse-fddot').value) ?? 0;
    const fddotPart = fddotVal !== 0 ? `, fddot = ${fddotVal} Hz/s^2` : '';
    const sub = `f = ${fmtFixed(res.freq_hz, 9)} Hz (P = ${periodMs} ms), fdot = ${res.fdot ?? 0} Hz/s${fddotPart}, T0 = ${fmtFixed(res.epoch, 2)} s, ${res.nbins} bins`;

    state.pulsePlot.setData({
      bands,
      refPeakPhase: refPeak,
      title: 'Multi-band phase-folded pulse alignment',
      subtitle: sub,
    });

    $('pulse-elapsed').textContent = `Folded in ${res.elapsed_ms} ms, reference: ${res.reference_band || (refBand && refBand.label) || 'N/A'}`;

    $('pulse-table-body').innerHTML = bands.map((b, idx) => {
      const isRef = idx === refIdx;
      const sign = (b.lagCycles ?? 0) >= 0 ? '+' : '';
      const lagCycStr = isRef ? '0.0000 (ref)' : (Number.isFinite(b.lagCycles) ? `${sign}${b.lagCycles.toFixed(4)}` : 'N/A');
      const lagMsStr = isRef ? '0.000' : (Number.isFinite(b.lagMs) ? `${sign}${b.lagMs.toFixed(3)}` : 'N/A');
      return `<tr>
        <td><span class="swatch" style="background:${b.color}"></span><strong>${esc(b.label)}</strong></td>
        <td class="num">${(b.numEvents || 0).toLocaleString()}</td>
        <td class="num">${fmtFixed(b.peakPhase, 4)}</td>
        <td class="num">${lagCycStr}</td>
        <td class="num">${lagMsStr}</td>
        <td class="num">${C.formatSci(b.z2, 5)}</td>
        <td class="num">${C.formatSci(b.chi2Red, 4)}</td>
        <td class="num">${Number.isFinite(b.pf) ? (b.pf * 100).toFixed(2) + '%' : 'N/A'}</td>
        <td class="num">${Number.isFinite(b.pfRms) ? (b.pfRms * 100).toFixed(2) + '%' : 'N/A'}</td>
      </tr>`;
    }).join('');
  }

  // =========================================================================
  // VIEW 5: FITS INSPECTOR
  // =========================================================================
  async function loadInspectorBrowser(dirPath) {
    try {
      const q = dirPath ? `?path=${encodeURIComponent(dirPath)}` : '';
      const data = await apiGet(`/api/browse${q}`);
      state.inspDir = data.path;
      state.inspParent = data.parent || '';
      state.inspEntries = data.entries || [];
      $('insp-dir-input').value = data.path;
      renderInspectorBrowser();
    } catch (err) {
      showBanner('insp', err.message);
    }
  }

  function renderInspectorBrowser() {
    const q = $('insp-file-filter').value.trim().toLowerCase();
    const body = $('insp-browser-body');
    const items = state.inspEntries.filter((e) => !q || e.name.toLowerCase().includes(q));
    if (!items.length) {
      body.innerHTML = '<tr><td colspan="2" class="table-empty">No matching entries</td></tr>';
      return;
    }
    body.innerHTML = items.map((e, i) => `
      <tr class="clickable" data-idx="${i}">
        <td class="mono-cell" title="${esc(e.name)}">${esc((e.is_dir ? '+ ' : '') + e.name)}</td>
        <td class="num meta-muted">${e.is_dir ? 'Folder' : fmtBytes(e.size_bytes)}</td>
      </tr>`).join('');
    body.querySelectorAll('tr.clickable').forEach((tr) => {
      tr.onclick = () => {
        const e = items[Number(tr.dataset.idx)];
        if (e.is_dir) loadInspectorBrowser(e.path);
        else {
          $('insp-path-input').value = e.path;
          inspectFitsFile(e.path, 0, 0);
        }
      };
    });
  }

  async function inspectFitsFile(pathOverride = null, hdu = 0, offset = 0) {
    const path = (pathOverride != null ? pathOverride : $('insp-path-input').value).trim();
    if (!path) return;
    hideBanner('insp');
    state.inspHdu = hdu;
    state.inspOffset = offset;
    try {
      const data = await apiGet(`/api/inspect?path=${encodeURIComponent(path)}&full=true&hdu=${hdu}&offset=${offset}&limit=${state.inspLimit}`);
      state.inspData = data;
      $('insp-val-box').classList.add('hidden');
      renderInspectorData(data);
    } catch (err) {
      showBanner('insp', err.message);
    }
  }

  async function validateFitsFile() {
    const path = $('insp-path-input').value.trim();
    if (!path) return;
    hideBanner('insp');
    try {
      const res = await apiGet(`/api/validate?path=${encodeURIComponent(path)}&strict=true`);
      $('insp-val-box').classList.remove('hidden');
      $('insp-val-summary').textContent = `${res.valid ? 'PASS' : 'FAIL'}: ${res.errors} errors, ${res.warnings} warnings`;
      const lines = (res.findings || []).map((f) => `[HDU ${f.hdu_index}] ${f.severity.toUpperCase()}: ${f.message}`);
      $('insp-val-findings').textContent = lines.length ? lines.join('\n') : '0 errors, 0 warnings. Strict OGIP and AstroSat validation passed.';
    } catch (err) {
      showBanner('insp', err.message);
    }
  }

  function renderInspectorData(d) {
    $('insp-m-inst').textContent = [d.instrume, d.detnam].filter(Boolean).join(' / ') || 'N/A';
    $('insp-m-obj').textContent = [d.object, d.obs_id].filter(Boolean).join(', ') || 'N/A';
    $('insp-m-exp').textContent = Number.isFinite(d.exposure) ? `${d.exposure.toFixed(2)} s` : 'N/A';
    $('insp-m-time').textContent = [d.timesys, d.timeref].filter(Boolean).join(' / ') || 'N/A';
    $('insp-m-radec').textContent = (Number.isFinite(d.ra_obj) && Number.isFinite(d.dec_obj))
      ? `${d.ra_obj.toFixed(4)} deg, ${d.dec_obj.toFixed(4)} deg` : 'N/A';

    const h = d.health || {};
    const csHtml = h.checksummed_hdus > 0
      ? (h.all_checksums_ok
        ? `<span class="status-ok">Checksums verified (${h.checksummed_hdus}/${d.num_hdus} HDUs)</span>`
        : '<span class="status-err">Checksum mismatch detected</span>')
      : '<span>No checksums recorded</span>';
    $('insp-status-line').innerHTML = `${csHtml}<span>File size: ${fmtBytes(d.file_size_bytes)}</span>`;

    const gtis = d.gtis || [];
    const gtiWrap = $('insp-gti-wrap');
    if (gtis.length && Number.isFinite(d.tstart) && Number.isFinite(d.tstop) && d.tstop > d.tstart) {
      gtiWrap.classList.remove('hidden');
      const span = d.tstop - d.tstart;
      const gtiExp = Number.isFinite(h.gti_exposure_sec) ? h.gti_exposure_sec : gtis.reduce((a, g) => a + (g.stop - g.start), 0);
      const duty = (gtiExp / span) * 100;
      $('insp-gti-summary').textContent = `${gtis.length} intervals, ${gtiExp.toFixed(1)} s / ${span.toFixed(1)} s (${duty.toFixed(1)}% duty)`;
      $('insp-gti-bar').innerHTML = gtis.slice(0, 400).map((g) => {
        const left = Math.max(0, Math.min(100, ((g.start - d.tstart) / span) * 100));
        const width = Math.max(0.3, Math.min(100 - left, ((g.stop - g.start) / span) * 100));
        return `<div class="gti-seg" style="left:${left.toFixed(3)}%;width:${width.toFixed(3)}%;" title="[${g.start.toFixed(1)}, ${g.stop.toFixed(1)}] (${(g.stop - g.start).toFixed(1)} s)"></div>`;
      }).join('');
    } else {
      gtiWrap.classList.add('hidden');
    }

    const qbox = $('insp-quick-open');
    qbox.innerHTML = '';
    const hdus = d.hdus || [];
    const isImageHdu = (x) => /image/i.test(x.kind || '');
    const addQ = (label, act) => {
      const b = document.createElement('button');
      b.type = 'button';
      b.className = 'text-action';
      b.textContent = label;
      b.onclick = () => openProductInTab(act, d.path);
      qbox.appendChild(b);
    };
    if (hdus.some((x) => /RATE|LIGHTCURVE/i.test(x.extname || '') || /LIGHTCURVE/i.test(x.hduclas1 || ''))) addQ('Open light curve', 'lightcurve');
    if (hdus.some((x) => /SPECTRUM/i.test(x.extname || '') || /SPECTRUM/i.test(x.hduclas1 || ''))) addQ('Open spectrum', 'spectrum');
    if (hdus.some((x) => isImageHdu(x) && x.shape && x.shape.length === 2)) addQ('Open in 2D imaging', 'image');
    if (hdus.some((x) => /EVENT|Q0|Q1/i.test(x.extname || '') || /EVENTS/i.test(x.hduclas1 || ''))) addQ('Fold pulse profile', 'pulse');

    $('insp-hdu-count').textContent = `${d.num_hdus} HDUs`;
    const hduBody = $('insp-hdu-body');
    hduBody.innerHTML = hdus.map((x) => {
      const sel = x.index === state.inspHdu ? 'selected' : '';
      const shapeStr = isImageHdu(x)
        ? (x.shape && x.shape.length ? x.shape.join('x') : 'Empty')
        : `${(x.nrows || 0).toLocaleString()} rows`;
      const cs = x.has_checksum
        ? (x.checksum_ok && x.datasum_ok ? '<span class="status-ok">OK</span>' : '<span class="status-err">Fail</span>')
        : '<span class="meta-muted">None</span>';
      return `<tr class="clickable ${sel}" data-hdu="${x.index}">
        <td class="num">${x.index}</td>
        <td class="mono-cell">${esc(x.extname || 'PRIMARY')}</td>
        <td class="meta-muted">${esc(x.kind)}</td>
        <td class="num">${esc(shapeStr)}</td>
        <td class="num">${Number.isFinite(x.exposure) ? x.exposure.toFixed(2) : 'N/A'}</td>
        <td>${cs}</td>
      </tr>`;
    }).join('');
    hduBody.querySelectorAll('tr').forEach((tr) => {
      tr.onclick = () => inspectFitsFile(d.path, Number(tr.dataset.hdu), 0);
    });

    const sh = d.selected_hdu;
    if (sh) {
      $('insp-selected-hdu-label').textContent = `HDU #${sh.index}: ${sh.extname || 'PRIMARY'} (${sh.kind})`;
      renderInspectorCards();
      renderInspectorRows(sh.table);
    }
  }

  function renderInspectorCards() {
    const sh = state.inspData && state.inspData.selected_hdu;
    const body = $('insp-cards-body');
    if (!sh || !sh.cards) {
      body.innerHTML = '<tr><td colspan="3" class="table-empty">No header cards</td></tr>';
      return;
    }
    const q = $('insp-card-search').value.trim().toLowerCase();
    const cards = sh.cards.filter((c) => !q
      || (c.keyword || '').toLowerCase().includes(q)
      || (c.value || '').toLowerCase().includes(q)
      || (c.comment || '').toLowerCase().includes(q));
    body.innerHTML = cards.map((c) => `<tr>
      <td class="kw-cell">${esc(c.keyword)}</td>
      <td class="mono-cell">${esc(c.value)}</td>
      <td class="meta-muted">${esc(c.comment)}</td>
    </tr>`).join('') || '<tr><td colspan="3" class="table-empty">No matching cards</td></tr>';
  }

  function renderInspectorRows(tbl) {
    const head = $('insp-rows-head');
    const body = $('insp-rows-body');
    if (!tbl || !tbl.columns || !tbl.columns.length) {
      $('insp-rows-info').textContent = 'Selected HDU is not a binary table.';
      head.innerHTML = '';
      body.innerHTML = '';
      $('insp-rows-prev').disabled = true;
      $('insp-rows-next').disabled = true;
      return;
    }
    const end = Math.min(tbl.total_rows, tbl.offset + tbl.rows.length);
    $('insp-rows-info').textContent = `Rows ${tbl.offset + 1} to ${end} of ${tbl.total_rows.toLocaleString()}`;
    $('insp-rows-prev').disabled = tbl.offset <= 0;
    $('insp-rows-next').disabled = end >= tbl.total_rows;

    head.innerHTML = `<tr><th class="num">#</th>${tbl.columns.map((c) => `<th>${esc(c)}</th>`).join('')}</tr>`;
    body.innerHTML = tbl.rows.map((r, ri) =>
      `<tr><td class="num meta-muted">${tbl.offset + ri + 1}</td>${r.map((cell) => `<td class="mono-cell row-cell" title="${esc(cell)}">${esc(cell)}</td>`).join('')}</tr>`
    ).join('');
  }

  // =========================================================================
  // Initialization & Event Wiring
  // =========================================================================
  async function initHealth() {
    try {
      const h = await apiGet('/api/health');
      state.health = h;
      $('workspace-meta').textContent = `${basename(h.data_dir) || h.data_dir}, ${h.threads || 1} threads, v${h.version || '0.1'}`;
      $('workspace-meta').title = `Data directory: ${h.data_dir}`;

      populateLcSampleDropdown();
      const imgSel = $('img-sample-select');
      for (const f of (h.presets && h.presets.sample_images) || []) {
        const opt = document.createElement('option');
        opt.value = f;
        opt.textContent = basename(f);
        opt.title = f;
        imgSel.appendChild(opt);
      }
      if (h.presets && h.presets.laxpc_caldb) $('spec-caldb-input').value = h.presets.laxpc_caldb;

      applyPulsePreset(findSampleEvent('laxpc') ? 'laxpc4' : (findSampleEvent('czti') ? 'czti3' : 'laxpc4'));
      autofillRunnerSample();
      loadInspectorBrowser(h.data_dir || '');
      refreshJobs();
    } catch {
      $('workspace-meta').textContent = 'Backend offline';
    }
  }

  function wireEvents() {
    $('theme-toggle').onclick = () => applyTheme(state.theme === 'dark' ? 'light' : 'dark');

    document.querySelectorAll('[data-dismiss]').forEach((btn) => {
      btn.onclick = () => hideBanner(btn.dataset.dismiss);
    });

    document.querySelectorAll('.nav-link').forEach((btn) => {
      btn.onclick = () => switchTab(btn.dataset.tab);
    });

    window.addEventListener('keydown', (ev) => {
      const fileOpen = !$('file-modal').classList.contains('hidden');
      const tomlOpen = !$('toml-modal').classList.contains('hidden');
      if (ev.key === 'Escape') {
        if (fileOpen) { closeFileModal(); return; }
        if (tomlOpen) { closeTomlModal(); return; }
      }
      if (fileOpen || tomlOpen) return;
      if (ev.target && /INPUT|TEXTAREA|SELECT/.test(ev.target.tagName)) return;
      const map = { '1': 'runner', '2': 'lc-spec', '3': 'imaging', '4': 'pulse', '5': 'inspector' };
      if (map[ev.key]) switchTab(map[ev.key]);
    });

    document.querySelectorAll('.browse-btn').forEach((btn) => {
      btn.onclick = () => openFileModal({ targetId: btn.dataset.target, mode: btn.dataset.mode || 'file' });
    });

    $('file-modal').onclick = (e) => { if (e.target === $('file-modal')) closeFileModal(); };
    $('file-modal-close').onclick = closeFileModal;
    $('modal-cancel').onclick = closeFileModal;
    $('modal-up-btn').onclick = () => { if (state.modalParent) loadModalDir(state.modalParent); };
    $('modal-go-btn').onclick = () => loadModalDir($('modal-path-input').value.trim());
    $('modal-path-input').onkeydown = (e) => { if (e.key === 'Enter') loadModalDir($('modal-path-input').value.trim()); };
    $('modal-filter').oninput = renderModalEntries;
    $('modal-select-dir').onclick = () => {
      if (state.modalTarget && $(state.modalTarget)) {
        $(state.modalTarget).value = state.modalDir;
        $(state.modalTarget).dispatchEvent(new Event('input', { bubbles: true }));
      }
      if (state.modalCallback) state.modalCallback(state.modalDir);
      closeFileModal();
    };

    // TOML modal & Toast
    $('runner-load-toml').onclick = openTomlModal;
    $('toml-modal').onclick = (e) => { if (e.target === $('toml-modal')) closeTomlModal(); };
    $('toml-modal-close').onclick = closeTomlModal;
    $('toml-modal-cancel').onclick = closeTomlModal;
    $('toml-clear-btn').onclick = () => {
      $('toml-paste-input').value = '';
      $('toml-file-input').value = '';
      $('toml-file-name').textContent = 'Or drag and drop a .toml file here';
      $('toml-modal-err').classList.add('hidden');
    };
    const loadTomlFile = async (file) => {
      if (!file) return;
      $('toml-file-name').textContent = file.name;
      $('toml-paste-input').value = await file.text();
    };
    $('toml-file-input').onchange = (e) => {
      if (e.target.files && e.target.files[0]) loadTomlFile(e.target.files[0]);
    };
    const dropzone = $('toml-dropzone');
    dropzone.addEventListener('dragover', (e) => {
      e.preventDefault();
      dropzone.classList.add('drag-over');
    });
    dropzone.addEventListener('dragleave', () => dropzone.classList.remove('drag-over'));
    dropzone.addEventListener('drop', (e) => {
      e.preventDefault();
      dropzone.classList.remove('drag-over');
      const file = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
      if (file) loadTomlFile(file);
    });
    $('toml-apply-btn').onclick = applyTomlImport;
    $('toast-close').onclick = hideToast;

    // View 1: Runner & Tools
    document.querySelectorAll('#runner-pipeline-seg button, #runner-tools-seg button').forEach((b) => {
      b.onclick = () => setRunnerPipe(b.dataset.pipe);
    });
    document.querySelectorAll('.czti-quad-cb').forEach((cb) => {
      cb.onchange = () => {
        state.cztiQuadsExplicit = true;
        scheduleCliPreview();
      };
    });
    $('bary-orbit').addEventListener('input', () => {
      if ($('bary-orbit').value.trim()) $('bary-mkf').value = '';
    });
    $('bary-mkf').addEventListener('input', () => {
      if ($('bary-mkf').value.trim()) $('bary-orbit').value = '';
    });
    document.querySelectorAll('#backshift-epoch-seg button').forEach((b) => {
      b.onclick = () => setBackshiftEpochMode(b.dataset.epoch);
    });
    $('addspec-add-row').onclick = () => {
      state.addspecRows.push({ pha: '', rsp: '' });
      renderAddspecRows();
      scheduleCliPreview();
    };
    $('addspec-from-job-btn').onclick = useLatestCztiOutputsForAddspec;
    $('addspec-quick-pick').onchange = async (e) => {
      const p = e.target.value;
      if (!p) return;
      e.target.value = '';
      const rsp = await resolveVerifiedRsp(p);
      const emptyRow = state.addspecRows.find((r) => !r.pha.trim());
      if (emptyRow) {
        emptyRow.pha = p;
        if (!emptyRow.rsp && rsp) emptyRow.rsp = rsp;
      } else {
        state.addspecRows.push({ pha: p, rsp });
      }
      renderAddspecRows();
      scheduleCliPreview();
    };
    $('tab-runner').addEventListener('input', scheduleCliPreview);
    $('tab-runner').addEventListener('change', scheduleCliPreview);
    $('runner-autofill').onclick = autofillRunnerSample;
    $('runner-submit').onclick = submitJob;
    $('jobs-refresh-btn').onclick = refreshJobs;
    const copyWithFeedback = (btn, text) => {
      if (navigator.clipboard) navigator.clipboard.writeText(text).catch(() => {});
      const orig = btn.textContent;
      btn.textContent = 'Copied';
      setTimeout(() => { btn.textContent = orig; }, 1100);
    };
    $('copy-cli-btn').onclick = (e) => copyWithFeedback(e.currentTarget, $('preview-cli').textContent);
    $('copy-toml-btn').onclick = (e) => copyWithFeedback(e.currentTarget, $('preview-toml').textContent);
    $('download-toml-btn').onclick = downloadTomlConfig;

    // View 2: Light Curves & Spectra
    document.querySelectorAll('#lc-mode-seg button').forEach((b) => {
      b.onclick = () => setLcMode(b.dataset.mode);
    });
    $('lc-add-btn').onclick = () => addLcOrSpecFile();
    $('lc-clear-btn').onclick = () => { state.lcItems = []; renderLcSeriesList(); redrawLcPlot(); };
    $('lc-sample-select').onchange = (e) => {
      if (e.target.value) {
        $('lc-path-input').value = e.target.value;
        addLcOrSpecFile(e.target.value);
        e.target.value = '';
      }
    };
    $('lc-rebin').onchange = reloadAllLcItems;
    $('spec-good-only').onchange = reloadAllLcItems;
    $('lc-rel-time').onchange = redrawLcPlot;
    $('spec-xaxis').onchange = redrawLcPlot;
    $('lc-show-err').onchange = (e) => state.lcPlot.setOptions({ showErrors: e.target.checked });
    $('lc-log-x').onchange = (e) => state.lcPlot.setOptions({ logX: e.target.checked });
    $('lc-log-y').onchange = (e) => state.lcPlot.setOptions({ logY: e.target.checked });
    $('lc-step').onchange = (e) => state.lcPlot.setOptions({ stepMode: e.target.checked });
    $('lc-show-gti').onchange = (e) => state.lcPlot.setOptions({ showGti: e.target.checked });
    $('lc-reset-btn').onclick = () => state.lcPlot.resetView();
    $('lc-svg-btn').onclick = () => state.lcPlot.exportSVG({ filename: `astrosat_${state.lcMode}.svg` });
    $('lc-png-btn').onclick = () => state.lcPlot.exportPNG({ filename: `astrosat_${state.lcMode}.png` });

    // View 3: 2D Imaging
    $('img-load-btn').onclick = () => loadImageFile();
    $('img-sample-select').onchange = (e) => {
      if (e.target.value) {
        $('img-path-input').value = e.target.value;
        loadImageFile(e.target.value);
        e.target.value = '';
      }
    };
    $('img-hdu-select').onchange = (e) => { if (state.imgData) loadImageFile(state.imgData.path, e.target.value); };
    $('img-cmap').onchange = (e) => state.imgPlot.setOptions({ colormap: e.target.value });
    $('img-stretch').onchange = (e) => state.imgPlot.setOptions({ stretch: e.target.value });
    $('img-quads').onchange = (e) => state.imgPlot.setOptions({ showQuads: e.target.checked });
    $('img-peak').onchange = (e) => state.imgPlot.setOptions({ showPeak: e.target.checked });
    $('img-reset-btn').onclick = () => state.imgPlot.resetView();
    $('img-png-btn').onclick = () => state.imgPlot.exportPNG({ filename: 'astrosat_2d_image.png' });

    // View 4: Pulse Alignment
    const syncPeriodFromFreq = () => {
      const f = numOrNull($('pulse-freq').value);
      $('pulse-period').value = (f != null && f > 0) ? (1 / f).toFixed(12).replace(/\.?0+$/, '') : '';
    };
    const syncFreqFromPeriod = () => {
      const p = numOrNull($('pulse-period').value);
      if (p != null && p > 0) $('pulse-freq').value = (1 / p).toFixed(10).replace(/\.?0+$/, '');
    };
    syncPeriodFromFreq();
    $('pulse-freq').addEventListener('input', syncPeriodFromFreq);
    $('pulse-period').addEventListener('input', syncFreqFromPeriod);
    $('pulse-search-toggle').onclick = () => $('pulse-search-box').classList.toggle('hidden');
    $('pulse-run-search').onclick = runPeriodSearch;
    $('pgram-close').onclick = () => $('periodogram-card').classList.add('hidden');
    document.querySelectorAll('#pulse-preset-seg button').forEach((b) => {
      b.onclick = () => applyPulsePreset(b.dataset.preset);
    });
    $('pulse-add-input').onclick = () => {
      state.pulseInputs.push({
        label: `Input ${state.pulseInputs.length + 1}`,
        path: '',
        caldb: '',
        use_weights: false,
        bands: [{ label: 'Full band', emin_kev: '', emax_kev: '', unit: '' }],
      });
      renderPulseInputs();
    };
    $('pulse-fold-btn').onclick = runMultiBandFold;
    document.querySelectorAll('#pulse-mode-seg button').forEach((b) => {
      b.onclick = () => {
        document.querySelectorAll('#pulse-mode-seg button').forEach((x) => x.classList.toggle('active', x === b));
        state.pulsePlot.setOptions({ mode: b.dataset.mode });
      };
    });
    $('pulse-show-err').onchange = (e) => state.pulsePlot.setOptions({ showErrors: e.target.checked });
    $('pulse-show-guide').onchange = (e) => state.pulsePlot.setOptions({ showPeakGuide: e.target.checked });
    $('pulse-svg-btn').onclick = () => state.pulsePlot.exportSVG({ filename: 'astrosat_pulse_alignment.svg' });
    $('pulse-png-btn').onclick = () => state.pulsePlot.exportPNG({ filename: 'astrosat_pulse_alignment.png' });

    // View 5: FITS Inspector
    $('insp-up-btn').onclick = () => { if (state.inspParent) loadInspectorBrowser(state.inspParent); };
    $('insp-dir-go').onclick = () => loadInspectorBrowser($('insp-dir-input').value.trim());
    $('insp-dir-input').onkeydown = (e) => { if (e.key === 'Enter') loadInspectorBrowser($('insp-dir-input').value.trim()); };
    $('insp-file-filter').oninput = renderInspectorBrowser;
    $('insp-load-btn').onclick = () => inspectFitsFile(null, 0, 0);
    $('insp-path-input').onkeydown = (e) => { if (e.key === 'Enter') inspectFitsFile(null, 0, 0); };
    $('insp-validate-btn').onclick = validateFitsFile;
    $('insp-card-search').oninput = renderInspectorCards;
    document.querySelectorAll('#insp-detail-seg button').forEach((b) => {
      b.onclick = () => {
        document.querySelectorAll('#insp-detail-seg button').forEach((x) => x.classList.toggle('active', x === b));
        const sub = b.dataset.sub;
        $('insp-cards-panel').classList.toggle('hidden', sub !== 'cards');
        $('insp-rows-panel').classList.toggle('hidden', sub !== 'rows');
      };
    });
    $('insp-rows-prev').onclick = () => {
      if (state.inspData && state.inspOffset > 0) {
        inspectFitsFile(state.inspData.path, state.inspHdu, Math.max(0, state.inspOffset - state.inspLimit));
      }
    };
    $('insp-rows-next').onclick = () => {
      if (state.inspData) {
        inspectFitsFile(state.inspData.path, state.inspHdu, state.inspOffset + state.inspLimit);
      }
    };
  }

  function boot() {
    state.lcPlot = new C.SeriesPlot($('lc-canvas'), $('lc-tooltip'), {
      theme: state.theme,
      themeColors: STUDIO_CHART_COLORS,
    });
    state.imgPlot = new C.ImagePlot2D($('img-canvas'), {
      theme: state.theme,
      themeColors: STUDIO_CHART_COLORS,
      onHover: handleImageHover,
    });
    state.pulsePlot = new C.PulseAlignmentPlot($('pulse-canvas'), $('pulse-tooltip'), {
      theme: state.theme,
      themeColors: STUDIO_CHART_COLORS,
    });
    state.pgramPlot = new C.SeriesPlot($('pgram-canvas'), $('pgram-tooltip'), {
      theme: state.theme,
      themeColors: STUDIO_CHART_COLORS,
      showErrors: false,
    });

    applyTheme(state.theme);
    renderAddspecRows();
    setRunnerPipe(state.runnerPipe);
    setBackshiftEpochMode(state.backshiftEpochMode);
    wireEvents();
    initHealth();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else {
    boot();
  }
})();
