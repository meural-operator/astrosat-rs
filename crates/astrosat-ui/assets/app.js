/* astrosat-ui — Zero-dependency offline Scientific Web UI Application */
(() => {
  'use strict';

  const C = window.AstroCharts;
  const $ = (id) => document.getElementById(id);

  // =========================================================================
  // Global App State
  // =========================================================================
  const state = {
    theme: localStorage.getItem('astrosat-theme') || 'dark',
    health: null,
    activeTab: 'runner',
    runnerPipe: 'czti',
    jobPollTimer: null,
    // Light Curves & Spectra
    lcMode: 'lightcurve',
    lcItems: [], // [{ id, path, label, color, hidden, raw }]
    lcPlot: null,
    // 2D Imaging
    imgPlot: null,
    imgData: null,
    // Pulse Alignment
    pulsePlot: null,
    pgramPlot: null,
    pulseInputs: [], // [{ label, path, caldb, use_weights, bands: [{ label, emin_kev, emax_kev, unit }] }]
    pulseLastResult: null,
    // Inspector
    inspDir: '',
    inspEntries: [],
    inspData: null,
    inspHdu: 0,
    inspOffset: 0,
    inspLimit: 50,
    // Modal file browser
    modalTarget: null,
    modalMode: 'file', // 'file' | 'dir'
    modalDir: '',
    modalEntries: [],
    modalCallback: null,
  };

  // =========================================================================
  // Helpers
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
    if (!Number.isFinite(n) || n < 0) return '—';
    if (n < 1024) return `${n} B`;
    if (n < 1048576) return `${(n / 1024).toFixed(1)} KB`;
    if (n < 1073741824) return `${(n / 1048576).toFixed(1)} MB`;
    return `${(n / 1073741824).toFixed(2)} GB`;
  }

  function fmtFixed(v, d = 2, fallback = '—') {
    return Number.isFinite(v) ? v.toFixed(d) : fallback;
  }

  function numOrNull(v) {
    if (v === '' || v == null) return null;
    const n = Number(v);
    return Number.isFinite(n) ? n : null;
  }

  function showBanner(prefix, msg, kind = 'err') {
    const el = $(`${prefix}-banner`);
    const msgEl = $(`${prefix}-banner-msg`);
    if (!el || !msgEl) return;
    el.className = `banner banner-${kind}`;
    msgEl.textContent = msg;
  }

  function hideBanner(prefix) {
    const el = $(`${prefix}-banner`);
    if (el) el.classList.add('hidden');
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
    state.theme = theme === 'light' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', state.theme);
    localStorage.setItem('astrosat-theme', state.theme);
    const pal = C.getPalette(state.theme);
    state.lcItems.forEach((it, idx) => { it.color = pal[idx % pal.length]; });
    if (state.lcPlot) {
      state.lcPlot.setTheme(state.theme);
      if (state.lcItems.length) { renderLcSeriesList(); redrawLcPlot(); }
    }
    if (state.imgPlot) state.imgPlot.setTheme(state.theme);
    if (state.pulsePlot) {
      state.pulsePlot.setTheme(state.theme);
      if (state.pulseLastResult) renderPulseResult(state.pulseLastResult);
    }
    if (state.pgramPlot) state.pgramPlot.setTheme(state.theme);
    const icon = $('theme-icon-sun');
    if (icon) {
      icon.innerHTML = state.theme === 'dark'
        ? '<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.93 4.93l1.41 1.41M17.66 17.66l1.41 1.41M2 12h2M20 12h2M4.93 19.07l1.41-1.41M17.66 6.34l1.41-1.41"/>'
        : '<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/>';
    }
  }

  function switchTab(tab) {
    state.activeTab = tab;
    document.querySelectorAll('.nav-tab').forEach((b) => {
      const act = b.dataset.tab === tab;
      b.classList.toggle('active', act);
      b.setAttribute('aria-selected', act ? 'true' : 'false');
    });
    document.querySelectorAll('.tab-panel').forEach((p) => {
      p.classList.toggle('active', p.id === `tab-${tab}`);
    });
    // Resize visible canvases cleanly after becoming visible
    requestAnimationFrame(() => {
      if (tab === 'lc-spec' && state.lcPlot) state.lcPlot.resize();
      if (tab === 'imaging' && state.imgPlot) state.imgPlot.resize();
      if (tab === 'pulse') {
        if (state.pulsePlot) state.pulsePlot.resize();
        if (state.pgramPlot) state.pgramPlot.resize();
      }
    });
  }

  // =========================================================================
  // Shared Workspace File Browser Modal
  // =========================================================================
  function openFileModal({ targetId = null, mode = 'file', initialDir = '', onSelect = null }) {
    state.modalTarget = targetId;
    state.modalMode = mode;
    state.modalCallback = onSelect;
    $('file-modal-title').textContent = mode === 'dir' ? 'Select Directory' : 'Select FITS / Data File';
    $('modal-select-dir').style.display = mode === 'dir' ? '' : 'none';
    $('modal-hint').textContent = mode === 'dir' ? 'Navigate into a folder and click Select This Directory' : 'Click a file to select, or click a folder to browse';
    $('modal-filter').value = '';
    $('file-modal').classList.remove('hidden');

    // Populate quick jumps from health presets
    const jumps = $('modal-quick-jumps');
    jumps.innerHTML = '';
    const h = state.health;
    if (h) {
      const candidates = [
        { label: 'Data Dir', path: h.data_dir },
        { label: 'CZTI L1', path: h.presets && h.presets.czti_l1_dir },
        { label: 'CZTI L2', path: h.presets && h.presets.czti_l2_dir },
        { label: 'LAXPC L1', path: h.presets && h.presets.laxpc_l1_dir },
        { label: 'LAXPC L2', path: h.presets && h.presets.laxpc_l2_dir },
      ].filter((c) => c.path);
      for (const c of candidates) {
        const b = document.createElement('button');
        b.className = 'btn btn-sm';
        b.textContent = c.label;
        b.onclick = () => loadModalDir(c.path);
        jumps.appendChild(b);
      }
    }

    const startDir = initialDir || (h && h.data_dir) || '';
    loadModalDir(startDir);
  }

  function closeFileModal() {
    $('file-modal').classList.add('hidden');
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
      $('modal-entries-body').innerHTML = `<tr><td colspan="3" style="color:var(--err);">${esc(err.message)}</td></tr>`;
    }
  }

  function renderModalEntries() {
    const q = $('modal-filter').value.trim().toLowerCase();
    const body = $('modal-entries-body');
    const items = state.modalEntries.filter((e) => !q || e.name.toLowerCase().includes(q));
    if (!items.length) {
      body.innerHTML = `<tr><td colspan="3" class="muted">Empty directory</td></tr>`;
      return;
    }
    body.innerHTML = items.map((e, i) => {
      const icon = e.is_dir ? '▸' : '·';
      const badge = e.is_dir ? 'dir' : (e.product && e.product !== 'other' ? e.product : e.kind);
      return `<tr class="file-row" data-idx="${i}">
        <td><span class="file-name">${icon} ${esc(e.name)}</span></td>
        <td><span class="badge">${esc(badge)}</span></td>
        <td class="num">${e.is_dir ? '—' : fmtBytes(e.size_bytes)}</td>
      </tr>`;
    }).join('');
    body.querySelectorAll('tr.file-row').forEach((tr) => {
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
  // TAB 1: PIPELINE RUNNER
  // =========================================================================
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
      setStr('caldb_dir', 'czti-caldb-dir');
      setStr('output_dir', 'czti-output-dir');
      setStr('output_prefix', 'czti-prefix');
      czti.entry_stage = Number($('czti-entry-stage').value);
      czti.exit_stage = Number($('czti-exit-stage').value);
      czti.algo = $('czti-algo').value;
      czti.gti_type = $('czti-gti-type').value;
      setNum('ra_deg', 'czti-ra');
      setNum('dec_deg', 'czti-dec');
      setNum('time_binsize', 'czti-tbin');
      setNum('emin_kev', 'czti-emin');
      setNum('emax_kev', 'czti-emax');
      czti.image_algo = $('czti-img-algo').value;
      czti.mask_weight = $('czti-mask-weight').checked;
      czti.gen_resp = $('czti-gen-resp').checked;
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
        if (units.length) laxpc.units = units;
      }
      setNum('tbin', 'laxpc-tbin');
      laxpc.spectrum_channels = Number($('laxpc-channels').value);
      setNum('min_elv_deg', 'laxpc-min-elv');
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
      setStr('orbit', 'bary-orbit');
      setStr('ephemeris', 'bary-ephem');
      setNum('ra_deg', 'bary-ra');
      setNum('dec_deg', 'bary-dec');
      bary.frame = $('bary-frame').value;
      bary.write_checksums = $('bary-checksums').checked;
      payload.bary = bary;
    }
    return payload;
  }

  let previewTimer = null;
  function scheduleCliPreview() {
    clearTimeout(previewTimer);
    previewTimer = setTimeout(updateCliPreview, 180);
  }

  async function updateCliPreview() {
    try {
      const data = await apiPost('/api/cli-preview', collectRunnerPayload());
      $('preview-cli').textContent = data.cli_command || '';
      $('preview-toml').textContent = data.toml_config || '';
    } catch (err) {
      $('preview-cli').textContent = `# Preview error: ${err.message}`;
    }
  }

  function autofillRunnerSample() {
    const p = state.health && state.health.presets;
    if (!p) return;
    const baseOut = state.health && state.health.cwd
      ? `${state.health.cwd.replace(/[\\/]+$/, '')}/.tmp`
      : './.tmp';
    if (state.runnerPipe === 'czti') {
      if (p.czti_l1_dir) $('czti-input-dir').value = p.czti_l1_dir;
      if (p.czti_caldb) $('czti-caldb-dir').value = p.czti_caldb;
      $('czti-output-dir').value = `${baseOut}/ui_czti_out`;
      $('czti-prefix').value = 'czti_ui';
    } else if (state.runnerPipe === 'laxpc') {
      if (p.laxpc_l1_dir) $('laxpc-input-dir').value = p.laxpc_l1_dir;
      if (p.laxpc_caldb) $('laxpc-caldb-dir').value = p.laxpc_caldb;
      $('laxpc-output-dir').value = `${baseOut}/ui_laxpc_out`;
    } else if (state.runnerPipe === 'bary') {
      const ev = (p.sample_events && p.sample_events[0]) || '';
      if (ev) $('bary-input').value = ev;
      if (p.laxpc_orb || p.czti_orb) $('bary-orbit').value = p.laxpc_orb || p.czti_orb;
      if (p.jpleph_405) $('bary-ephem').value = p.jpleph_405;
      $('bary-output').value = `${baseOut}/ui_bary.evt`;
      $('bary-ra').value = '83.633083';
      $('bary-dec').value = '22.014500';
    }
    updateCliPreview();
  }

  async function submitJob() {
    hideBanner('runner');
    const btn = $('runner-submit');
    btn.disabled = true;
    try {
      await apiPost('/api/jobs/run', collectRunnerPayload());
      await refreshJobs();
      startJobPolling();
    } catch (err) {
      showBanner('runner', err.message, 'err');
    } finally {
      btn.disabled = false;
    }
  }

  function startJobPolling() {
    if (state.jobPollTimer) return;
    state.jobPollTimer = setInterval(async () => {
      const jobs = await refreshJobs();
      if (!jobs.some((j) => j.status === 'running')) {
        clearInterval(state.jobPollTimer);
        state.jobPollTimer = null;
      }
    }, 1200);
  }

  async function refreshJobs() {
    try {
      const data = await apiGet('/api/jobs');
      const jobs = data.jobs || [];
      renderJobs(jobs);
      return jobs;
    } catch {
      return [];
    }
  }

  function formatJobSummary(summary) {
    if (summary == null) return '';
    if (typeof summary === 'string') return summary;
    if (typeof summary === 'object') {
      return Object.entries(summary)
        .filter(([, v]) => v != null && v !== '')
        .map(([k, v]) => {
          const valStr = typeof v === 'number' ? C.formatSci(v, 4)
            : (typeof v === 'object' ? JSON.stringify(v) : String(v));
          return `${k}: ${valStr}`;
        })
        .join(' · ');
    }
    return String(summary);
  }

  function renderJobs(jobs) {
    const box = $('jobs-list');
    if (!jobs.length) {
      box.innerHTML = `<div class="muted text-sm">No pipeline jobs run yet. Configure a pipeline and click <strong>Run Pipeline</strong>.</div>`;
      return;
    }
    // Backend returns jobs newest-first
    box.innerHTML = jobs.map((j) => {
      const badgeCls = j.status === 'succeeded' ? 'badge-ok' : (j.status === 'failed' ? 'badge-err' : 'badge-warn');
      const elapsed = Number.isFinite(j.elapsed_ms) ? `${(j.elapsed_ms / 1000).toFixed(2)}s` : '…';
      const summaryText = formatJobSummary(j.summary);
      const outputsHtml = (j.outputs || []).slice(0, 12).map((o) => {
        const prod = o.product || 'other';
        let openAction = '';
        if (prod === 'lightcurve' || prod === 'spectrum') {
          openAction = `<button class="btn btn-sm job-open" data-act="${prod}" data-path="${esc(o.path)}">Plot ${esc(o.label || basename(o.path))}</button>`;
        } else if (prod === 'image') {
          openAction = `<button class="btn btn-sm job-open" data-act="image" data-path="${esc(o.path)}">2D View ${esc(o.label || basename(o.path))}</button>`;
        } else if (prod === 'events') {
          openAction = `<button class="btn btn-sm job-open" data-act="pulse" data-path="${esc(o.path)}">Fold ${esc(o.label || basename(o.path))}</button>`;
        }
        return `${openAction}<button class="btn btn-sm btn-ghost job-open" data-act="inspect" data-path="${esc(o.path)}" title="Inspect ${esc(o.path)}">Inspect ${esc(basename(o.path))}</button>`;
      }).join('');
      return `<div class="job-item">
        <div class="row-between">
          <div class="row">
            <span class="badge ${badgeCls}">${esc(j.status)}</span>
            <strong class="mono">${esc(j.pipeline.toUpperCase())}</strong>
            <span class="text-sm mono muted">#${esc(j.id)}</span>
          </div>
          <span class="text-sm mono muted">${elapsed}</span>
        </div>
        ${summaryText ? `<div class="text-sm mono">${esc(summaryText)}</div>` : ''}
        ${j.error ? `<div class="text-sm" style="color:var(--err);">${esc(j.error)}</div>` : ''}
        ${outputsHtml ? `<div class="job-outputs">${outputsHtml}</div>` : ''}
      </div>`;
    }).join('');

    box.querySelectorAll('.job-open').forEach((btn) => {
      btn.onclick = () => openProductInTab(btn.dataset.act, btn.dataset.path);
    });
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
        bands: [{ label: 'Full Band', emin_kev: '', emax_kev: '', unit: '' }],
      }];
      renderPulseInputs();
    } else {
      switchTab('inspector');
      $('insp-path-input').value = path;
      inspectFitsFile(path, 0, 0);
    }
  }

  // =========================================================================
  // TAB 2: LIGHT CURVES & SPECTRA
  // =========================================================================
  function setLcMode(mode) {
    if (state.lcMode !== mode) {
      state.lcMode = mode;
      state.lcItems = [];
      state.lcPlot.clear();
      $('lc-empty').classList.remove('hidden');
      renderLcSeriesList();
    }
    document.querySelectorAll('#lc-mode-seg button').forEach((b) => {
      b.classList.toggle('active', b.dataset.mode === mode);
    });
    $('lc-opts-box').style.display = mode === 'lightcurve' ? '' : 'none';
    $('spec-opts-box').style.display = mode === 'spectrum' ? '' : 'none';
    $('lc-gti-wrap').style.display = mode === 'lightcurve' ? '' : 'none';
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
    sel.innerHTML = `<option value="">— Select detected sample —</option>`;
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

  async function addLcOrSpecFile(pathOverride = null) {
    const path = (pathOverride != null ? pathOverride : $('lc-path-input').value).trim();
    if (!path) return;
    hideBanner('lc');
    $('lc-loading').classList.remove('hidden');
    try {
      const pal = C.getPalette(state.theme);
      const color = pal[state.lcItems.length % pal.length];
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
      state.lcItems.push({
        id: `${Date.now()}_${Math.random().toString(36).slice(2, 6)}`,
        path,
        label: raw.title || basename(path),
        color,
        hidden: false,
        raw,
      });
      renderLcSeriesList();
      redrawLcPlot();
    } catch (err) {
      showBanner('lc', err.message, 'err');
    } finally {
      $('lc-loading').classList.add('hidden');
    }
  }

  async function reloadAllLcItems() {
    if (!state.lcItems.length) return;
    const paths = state.lcItems.map((it) => ({ path: it.path, color: it.color, hidden: it.hidden }));
    state.lcItems = [];
    $('lc-loading').classList.remove('hidden');
    hideBanner('lc');
    try {
      for (const item of paths) {
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
          label: raw.title || basename(item.path),
          color: item.color,
          hidden: item.hidden,
          raw,
        });
      }
      renderLcSeriesList();
      redrawLcPlot();
    } catch (err) {
      showBanner('lc', err.message, 'err');
    } finally {
      $('lc-loading').classList.add('hidden');
    }
  }

  function renderLcSeriesList() {
    const box = $('lc-series-list');
    if (!state.lcItems.length) {
      box.innerHTML = `<div class="muted text-sm">No files loaded yet.</div>`;
      return;
    }
    box.innerHTML = state.lcItems.map((it, idx) => {
      const r = it.raw;
      const meta = state.lcMode === 'lightcurve'
        ? `${(r.time || []).length} pts · ${C.formatSci(r.mean_rate, 3)} c/s`
        : `${r.good_channels || (r.channel || []).length} ch · ${C.formatSci(r.total_rate, 3)} c/s`;
      return `<div class="series-chip">
        <div class="series-chip-left">
          <span class="color-dot" style="background:${it.color}; opacity:${it.hidden ? 0.3 : 1};"></span>
          <div style="min-width:0;">
            <div class="series-chip-name" title="${esc(it.path)}">${esc(basename(it.path))}</div>
            <div class="series-chip-meta">${esc(meta)}</div>
          </div>
        </div>
        <div class="row" style="flex-wrap:nowrap;">
          <button class="btn btn-ghost btn-sm lc-toggle-vis" data-idx="${idx}" title="Toggle visibility">${it.hidden ? 'Show' : 'Hide'}</button>
          <button class="btn btn-ghost btn-sm lc-remove" data-idx="${idx}" title="Remove">×</button>
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
      $('lc-stat-exp').textContent = '—';
      $('lc-stat-rate').textContent = '—';
      $('lc-stat-meta').textContent = '—';
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
        title: items.length === 1 ? (first.title || basename(items[0].path)) : `Overlay (${items.length} Light Curves)`,
        xLabel: rel ? `Time − ${t0.toFixed(1)}` : `Time (${first.timesys || 'MET'})`,
        xUnit: 's',
        yLabel: 'Count Rate',
        yUnit: first.units || 'counts/s',
      });
      $('lc-stat-count').textContent = String(items.length);
      $('lc-stat-exp').textContent = Number.isFinite(first.exposure) ? `${first.exposure.toFixed(2)} s` : '—';
      $('lc-stat-rate').textContent = Number.isFinite(first.mean_rate) ? `${C.formatSci(first.mean_rate, 4)} c/s` : '—';
      $('lc-stat-meta').textContent = Number.isFinite(first.effective_timedel)
        ? `Δt = ${C.formatSci(first.effective_timedel, 3)} s (${first.effective_rebin || 1}×)`
        : '—';
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
        title: items.length === 1 ? (first.title || basename(items[0].path)) : `Overlay (${items.length} Spectra)`,
        xLabel: useEnergy ? 'Energy' : (first.chantype || 'Channel'),
        xUnit: useEnergy ? 'keV' : '',
        yLabel: first.is_rate ? 'Count Rate' : 'Counts',
        yUnit: first.units || (first.is_rate ? 'counts/s/ch' : 'counts/ch'),
      });
      $('lc-stat-count').textContent = String(items.length);
      $('lc-stat-exp').textContent = Number.isFinite(first.exposure) ? `${first.exposure.toFixed(2)} s` : '—';
      $('lc-stat-rate').textContent = Number.isFinite(first.total_rate) ? `${C.formatSci(first.total_rate, 4)} c/s` : '—';
      $('lc-stat-meta').textContent = `${first.energy_source || 'channel'} · ${first.good_channels}/${first.total_channels} good`;
    }
  }

  // =========================================================================
  // TAB 3: 2D IMAGING (DPH / DPI / SKY MAPS)
  // =========================================================================
  async function loadImageFile(pathOverride = null, hduOverride = null) {
    const path = (pathOverride != null ? pathOverride : $('img-path-input').value).trim();
    if (!path) return;
    hideBanner('img');
    $('img-loading').classList.remove('hidden');
    try {
      let url = `/api/image?path=${encodeURIComponent(path)}`;
      if (hduOverride !== null && hduOverride !== '') {
        url += `&hdu=${encodeURIComponent(hduOverride)}`;
      }
      const data = await apiGet(url);
      state.imgData = data;
      $('img-empty').classList.add('hidden');

      // Populate HDU dropdown
      const hduSel = $('img-hdu-select');
      const avail = data.available_image_hdus || [];
      hduSel.innerHTML = avail.map((h) =>
        `<option value="${h.index}" ${h.index === data.hdu_index ? 'selected' : ''}>#${h.index} ${esc(h.extname || 'PRIMARY')} (${h.width}×${h.height})</option>`
      ).join('') || `<option value="${data.hdu_index}">#${data.hdu_index} ${esc(data.hdu || 'IMAGE')}</option>`;

      // Auto-toggle Q0-Q3 overlay if 128x128 CZTI detector image
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

      // Stats cards
      const st = data.stats || {};
      const pk = st.peak || {};
      $('img-stat-dim').textContent = `${data.width} × ${data.height}`;
      $('img-stat-snr').textContent = Number.isFinite(pk.snr) ? `${pk.snr.toFixed(2)} σ` : '—';
      $('img-stat-peak').textContent = Number.isFinite(pk.x) ? `(${pk.x}, ${pk.y})` : '—';
      $('img-stat-max').textContent = C.formatSci(pk.val != null ? pk.val : st.max, 4);
      $('img-stat-minmed').textContent = `${C.formatSci(st.min, 3)} / ${C.formatSci(st.p50, 3)}`;
      $('img-stat-meanstd').textContent = `${C.formatSci(st.mean, 3)} ± ${C.formatSci(st.std, 3)}`;

      if (data.wcs && data.wcs.has_wcs) {
        $('img-wcs-badge').textContent = `WCS: ${data.wcs.ctype1 || 'X'} / ${data.wcs.ctype2 || 'Y'}`;
        $('img-wcs-badge').className = 'badge badge-accent';
      } else {
        $('img-wcs-badge').textContent = 'Detector Pixel Coords';
        $('img-wcs-badge').className = 'badge';
      }
    } catch (err) {
      showBanner('img', err.message, 'err');
    } finally {
      $('img-loading').classList.add('hidden');
    }
  }

  function handleImageHover(hit) {
    if (!hit) {
      $('img-cur-pix').textContent = '—';
      $('img-cur-val').textContent = '—';
      $('img-cur-wcs1').textContent = '—';
      $('img-cur-wcs2').textContent = '—';
      return;
    }
    $('img-cur-pix').textContent = `(${hit.fx}, ${hit.fy})`;
    $('img-cur-val').textContent = C.formatSci(hit.val, 5);
    if (hit.world) {
      $('img-cur-wcs1-lbl').textContent = hit.world.labelX;
      $('img-cur-wcs2-lbl').textContent = hit.world.labelY;
      $('img-cur-wcs1').textContent = `${hit.world.x.toFixed(4)}°`;
      $('img-cur-wcs2').textContent = `${hit.world.y.toFixed(4)}°`;
    } else {
      $('img-cur-wcs1-lbl').textContent = 'FITS 1-based X';
      $('img-cur-wcs2-lbl').textContent = 'FITS 1-based Y';
      $('img-cur-wcs1').textContent = String(hit.fx + 1);
      $('img-cur-wcs2').textContent = String(hit.fy + 1);
    }
  }

  // =========================================================================
  // TAB 4: PHASE-FOLDED & MULTI-BAND PULSE ALIGNMENT
  // =========================================================================
  function findSampleEvent(kind) {
    const evs = (state.health && state.health.presets && state.health.presets.sample_events) || [];
    if (kind === 'laxpc') {
      return evs.find((p) => /laxpc|lxp|level2\.event/i.test(p) && /bary/i.test(p))
        || evs.find((p) => /laxpc|lxp|level2\.event/i.test(p))
        || evs[0] || '';
    }
    if (kind === 'czti') {
      return evs.find((p) => /czt/i.test(p) && /bary/i.test(p))
        || evs.find((p) => /czt/i.test(p))
        || '';
    }
    return evs[0] || '';
  }

  function applyPulsePreset(kind) {
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
          { label: '3–10 keV', emin_kev: 3, emax_kev: 10, unit: '' },
          { label: '10–20 keV', emin_kev: 10, emax_kev: 20, unit: '' },
          { label: '20–40 keV', emin_kev: 20, emax_kev: 40, unit: '' },
          { label: '40–80 keV', emin_kev: 40, emax_kev: 80, unit: '' },
        ],
      }];
    } else if (kind === 'czti3') {
      state.pulseInputs = [{
        label: 'CZTI',
        path: cztiPath,
        caldb: '',
        use_weights: false,
        bands: [
          { label: '20–60 keV', emin_kev: 20, emax_kev: 60, unit: '' },
          { label: '60–120 keV', emin_kev: 60, emax_kev: 120, unit: '' },
          { label: '120–200 keV', emin_kev: 120, emax_kev: 200, unit: '' },
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
            { label: '3–10 keV', emin_kev: 3, emax_kev: 10, unit: '' },
            { label: '10–20 keV', emin_kev: 10, emax_kev: 20, unit: '' },
            { label: '20–40 keV', emin_kev: 20, emax_kev: 40, unit: '' },
            { label: '40–80 keV', emin_kev: 40, emax_kev: 80, unit: '' },
          ],
        },
        {
          label: 'CZTI',
          path: cztiPath,
          caldb: '',
          use_weights: false,
          bands: [
            { label: '20–60 keV', emin_kev: 20, emax_kev: 60, unit: '' },
            { label: '60–120 keV', emin_kev: 60, emax_kev: 120, unit: '' },
            { label: '120–200 keV', emin_kev: 120, emax_kev: 200, unit: '' },
          ],
        },
      ];
    }
    renderPulseInputs();
  }

  function renderPulseInputs() {
    const box = $('pulse-inputs-container');
    const pal = C.getPalette(state.theme);
    let bandColorIdx = 0;

    if (!state.pulseInputs.length) {
      box.innerHTML = `<div class="muted text-sm">Click a preset above or <strong>+ Event File</strong> to configure pulse bands.</div>`;
      return;
    }

    box.innerHTML = state.pulseInputs.map((inp, ii) => {
      const bandsHtml = inp.bands.map((b, bi) => {
        const col = pal[bandColorIdx++ % pal.length];
        return `<div class="band-row" data-ii="${ii}" data-bi="${bi}">
          <span class="color-dot" style="background:${col}"></span>
          <input type="text" class="b-label" value="${esc(b.label)}" placeholder="Label" />
          <input type="number" step="any" class="b-emin" value="${esc(b.emin_kev)}" placeholder="Emin" title="Min Energy (keV)" />
          <input type="number" step="any" class="b-emax" value="${esc(b.emax_kev)}" placeholder="Emax" title="Max Energy (keV)" />
          <input type="number" min="1" max="3" class="b-unit" value="${esc(b.unit)}" placeholder="Unit" title="Optional LAXPC unit (1-3)" />
          <button class="btn btn-ghost btn-sm b-del" title="Remove band" style="padding:2px 5px;">×</button>
        </div>`;
      }).join('');

      return `<div class="job-item" data-ii="${ii}">
        <div class="row-between">
          <input type="text" class="inp-label" value="${esc(inp.label)}" placeholder="Instrument label" style="width:120px; font-weight:600;" />
          <div class="row">
            <label class="check-pill" title="Use mask WEIGHT column (CZTI)"><input type="checkbox" class="inp-weights" ${inp.use_weights ? 'checked' : ''} /> Weights</label>
            <button class="btn btn-sm inp-add-band">+ Band</button>
            <button class="btn btn-sm btn-ghost inp-del" title="Remove input">×</button>
          </div>
        </div>
        <div class="input-with-btn">
          <input type="text" class="mono inp-path" value="${esc(inp.path)}" placeholder="Barycentered .evn / .evt / .fits path" />
          <button class="btn btn-sm inp-browse">Browse</button>
        </div>
        <div class="input-with-btn">
          <input type="text" class="mono inp-caldb" value="${esc(inp.caldb || '')}" placeholder="Optional CALDB dir (for LAXPC keV)" />
          <button class="btn btn-sm inp-browse-caldb">CALDB</button>
        </div>
        <div class="stack" style="gap:4px;">${bandsHtml}</div>
      </div>`;
    }).join('');

    // Wire input listeners
    box.querySelectorAll('.job-item').forEach((card) => {
      const ii = Number(card.dataset.ii);
      const inp = state.pulseInputs[ii];
      card.querySelector('.inp-label').oninput = (e) => { inp.label = e.target.value; };
      card.querySelector('.inp-path').oninput = (e) => { inp.path = e.target.value; };
      card.querySelector('.inp-caldb').oninput = (e) => { inp.caldb = e.target.value; };
      card.querySelector('.inp-weights').onchange = (e) => { inp.use_weights = e.target.checked; };
      card.querySelector('.inp-del').onclick = () => { state.pulseInputs.splice(ii, 1); renderPulseInputs(); };
      card.querySelector('.inp-add-band').onclick = () => {
        inp.bands.push({ label: `Band ${inp.bands.length + 1}`, emin_kev: '', emax_kev: '', unit: '' });
        renderPulseInputs();
      };
      card.querySelector('.inp-browse').onclick = () => {
        openFileModal({
          mode: 'file',
          onSelect: (p) => { inp.path = p; renderPulseInputs(); },
        });
      };
      card.querySelector('.inp-browse-caldb').onclick = () => {
        openFileModal({
          mode: 'dir',
          onSelect: (p) => { inp.caldb = p; renderPulseInputs(); },
        });
      };
      card.querySelectorAll('.band-row').forEach((row) => {
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
      showBanner('pulse', 'Please select a barycentered event file first.', 'err');
      return;
    }
    const freq = numOrNull($('pulse-freq').value);
    const body = {
      path,
      freq_hz: freq,
      search_span_hz: numOrNull($('search-span').value) || 0.002,
      n_trials: Math.floor(numOrNull($('search-trials').value) || 300),
      fdot: numOrNull($('pulse-fdot').value) || 0,
      harmonics: Number($('pulse-harmonics').value) || 2,
      nbins: Number($('pulse-nbins').value) || 32,
    };
    const ep = numOrNull($('pulse-epoch').value);
    if (ep != null) body.epoch = ep;

    const btn = $('pulse-run-search');
    btn.disabled = true;
    $('search-result-meta').textContent = 'Searching…';
    try {
      const res = await apiPost('/api/pulse/search', body);
      if (Number.isFinite(res.best_freq_hz)) {
        $('pulse-freq').value = res.best_freq_hz.toFixed(9);
      }
      $('search-result-meta').textContent = `Peak ν = ${fmtFixed(res.best_freq_hz, 7)} Hz · Z² = ${C.formatSci(res.max_z2_n, 4)} (${res.elapsed_ms} ms)`;
      $('periodogram-card').style.display = '';
      $('pgram-badge').textContent = `ν = ${fmtFixed(res.best_freq_hz, 7)} Hz · Z² = ${C.formatSci(res.max_z2_n, 4)} · ${fmtFixed(res.detection_sigma, 1)}σ`;

      requestAnimationFrame(() => {
        state.pgramPlot.resize();
        const pal = C.getPalette(state.theme);
        state.pgramPlot.setData({
          series: [{
            id: 'z2',
            label: `Z²_${res.harmonics}`,
            color: pal[0],
            x: res.periodogram.freq_hz,
            y: res.periodogram.z2_n,
            yerr: null,
          }],
          vlines: [{ x: res.best_freq_hz, label: `Peak ${res.best_freq_hz.toFixed(6)} Hz` }],
          title: `Z²_${res.harmonics} Periodogram — ${basename(res.path)} (${res.num_events.toLocaleString()} events)`,
          xLabel: 'Frequency',
          xUnit: 'Hz',
          yLabel: `Z²_${res.harmonics} Power`,
          yUnit: '',
        });
      });
    } catch (err) {
      showBanner('pulse', err.message, 'err');
      $('search-result-meta').textContent = 'Search failed';
    } finally {
      btn.disabled = false;
    }
  }

  async function runMultiBandFold() {
    hideBanner('pulse');
    const freq = numOrNull($('pulse-freq').value);
    if (freq == null || freq <= 0) {
      showBanner('pulse', 'Please enter a valid positive spin frequency ν (Hz).', 'err');
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
      showBanner('pulse', 'Add at least one event file with one energy band (or click a preset button).', 'err');
      return;
    }

    const body = {
      inputs: validInputs,
      freq_hz: freq,
      fdot: numOrNull($('pulse-fdot').value) || 0,
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
      showBanner('pulse', err.message, 'err');
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
    const refPeak = bands.length ? bands[0].peakPhase : null;
    const periodMs = Number.isFinite(res.period_sec) ? (res.period_sec * 1000).toFixed(4) : '—';
    const sub = `ν = ${fmtFixed(res.freq_hz, 9)} Hz (P = ${periodMs} ms)  ·  ν̇ = ${res.fdot ?? 0} Hz/s  ·  T₀ = ${fmtFixed(res.epoch, 2)} s  ·  ${res.nbins} bins`;

    state.pulsePlot.setData({
      bands,
      refPeakPhase: refPeak,
      title: 'Multi-Band Phase-Folded Pulse Profile Alignment',
      subtitle: sub,
    });

    $('pulse-elapsed').textContent = `Folded in ${res.elapsed_ms} ms · Ref: ${res.reference_band || (bands[0] && bands[0].label) || '—'}`;

    const tbody = $('pulse-table-body');
    tbody.innerHTML = bands.map((b, idx) => {
      const sign = (b.lagCycles ?? 0) >= 0 ? '+' : '';
      const lagCycStr = idx === 0 ? '0.0000 (ref)' : (Number.isFinite(b.lagCycles) ? `${sign}${b.lagCycles.toFixed(4)}` : '—');
      const lagMsStr = idx === 0 ? '0.000' : (Number.isFinite(b.lagMs) ? `${sign}${b.lagMs.toFixed(3)}` : '—');
      return `<tr>
        <td><span class="tt-swatch" style="background:${b.color}"></span><strong>${esc(b.label)}</strong></td>
        <td class="num">${(b.numEvents || 0).toLocaleString()}</td>
        <td class="num">${fmtFixed(b.peakPhase, 4)}</td>
        <td class="num">${lagCycStr}</td>
        <td class="num">${lagMsStr}</td>
        <td class="num">${C.formatSci(b.z2, 5)}</td>
        <td class="num">${C.formatSci(b.chi2Red, 4)}</td>
        <td class="num">${Number.isFinite(b.pf) ? (b.pf * 100).toFixed(2) + '%' : '—'}</td>
        <td class="num">${Number.isFinite(b.pfRms) ? (b.pfRms * 100).toFixed(2) + '%' : '—'}</td>
      </tr>`;
    }).join('');
  }

  // =========================================================================
  // TAB 5: FITS INSPECTOR
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
      showBanner('insp', err.message, 'err');
    }
  }

  function renderInspectorBrowser() {
    const q = $('insp-file-filter').value.trim().toLowerCase();
    const body = $('insp-browser-body');
    const items = state.inspEntries.filter((e) => !q || e.name.toLowerCase().includes(q));
    if (!items.length) {
      body.innerHTML = `<tr><td colspan="2" class="muted">No matching entries</td></tr>`;
      return;
    }
    body.innerHTML = items.map((e, i) => `
      <tr class="file-row" data-idx="${i}">
        <td><span class="file-name">${e.is_dir ? '▸' : '·'} ${esc(e.name)}</span></td>
        <td class="num">${e.is_dir ? '<span class="badge">dir</span>' : fmtBytes(e.size_bytes)}</td>
      </tr>`).join('');
    body.querySelectorAll('tr.file-row').forEach((tr) => {
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
      showBanner('insp', err.message, 'err');
    }
  }

  async function validateFitsFile() {
    const path = $('insp-path-input').value.trim();
    if (!path) return;
    hideBanner('insp');
    try {
      const res = await apiGet(`/api/validate?path=${encodeURIComponent(path)}&strict=true`);
      $('insp-val-box').classList.remove('hidden');
      $('insp-val-summary').textContent = `${res.valid ? 'PASS' : 'FAIL'} · ${res.errors} errors · ${res.warnings} warnings`;
      const lines = (res.findings || []).map((f) => `[HDU ${f.hdu_index}] ${f.severity.toUpperCase()}: ${f.message}`);
      $('insp-val-findings').textContent = lines.length ? lines.join('\n') : '0 errors, 0 warnings — strict OGIP & AstroSat validation passed.';
      const badgeCls = res.valid && res.warnings === 0 ? 'badge-ok' : (res.valid ? 'badge-warn' : 'badge-err');
      const existing = $('insp-val-badge');
      if (existing) existing.remove();
      const span = document.createElement('span');
      span.id = 'insp-val-badge';
      span.className = `badge ${badgeCls}`;
      span.textContent = res.valid ? `validate --strict: PASS (${res.warnings}W)` : `validate --strict: FAIL (${res.errors}E)`;
      $('insp-badges').appendChild(span);
    } catch (err) {
      showBanner('insp', err.message, 'err');
    }
  }

  function renderInspectorData(d) {
    $('insp-m-inst').textContent = [d.instrume, d.detnam].filter(Boolean).join(' / ') || '—';
    $('insp-m-obj').textContent = [d.object, d.obs_id].filter(Boolean).join(' · ') || '—';
    $('insp-m-exp').textContent = Number.isFinite(d.exposure) ? `${d.exposure.toFixed(2)} s` : '—';
    $('insp-m-time').textContent = [d.timesys, d.timeref].filter(Boolean).join(' / ') || '—';
    $('insp-m-radec').textContent = (Number.isFinite(d.ra_obj) && Number.isFinite(d.dec_obj))
      ? `${d.ra_obj.toFixed(4)}°, ${d.dec_obj.toFixed(4)}°` : '—';

    // Checksum badges
    const badges = $('insp-badges');
    const h = d.health || {};
    const csBadge = h.checksummed_hdus > 0
      ? (h.all_checksums_ok
        ? `<span class="badge badge-ok">CHECKSUM/DATASUM OK (${h.checksummed_hdus}/${d.num_hdus})</span>`
        : `<span class="badge badge-err">CHECKSUM MISMATCH</span>`)
      : `<span class="badge">No Checksums</span>`;
    badges.innerHTML = `${csBadge}<span class="badge">${fmtBytes(d.file_size_bytes)}</span>`;

    // GTI summary bar
    const gtis = d.gtis || [];
    const gtiWrap = $('insp-gti-wrap');
    if (gtis.length && Number.isFinite(d.tstart) && Number.isFinite(d.tstop) && d.tstop > d.tstart) {
      gtiWrap.style.display = '';
      const span = d.tstop - d.tstart;
      const gtiExp = Number.isFinite(h.gti_exposure_sec) ? h.gti_exposure_sec : gtis.reduce((a, g) => a + (g.stop - g.start), 0);
      const duty = (gtiExp / span) * 100;
      $('insp-gti-summary').textContent = `${gtis.length} intervals · ${gtiExp.toFixed(1)} s / ${span.toFixed(1)} s (${duty.toFixed(1)}% duty)`;
      $('insp-gti-bar').innerHTML = gtis.slice(0, 400).map((g) => {
        const left = Math.max(0, Math.min(100, ((g.start - d.tstart) / span) * 100));
        const width = Math.max(0.3, Math.min(100 - left, ((g.stop - g.start) / span) * 100));
        return `<div class="gti-seg" style="left:${left.toFixed(3)}%;width:${width.toFixed(3)}%;" title="[${g.start.toFixed(1)}, ${g.stop.toFixed(1)}] (${(g.stop - g.start).toFixed(1)}s)"></div>`;
      }).join('');
    } else {
      gtiWrap.style.display = 'none';
    }

    // Quick-open buttons based on HDU contents
    const qbox = $('insp-quick-open');
    qbox.innerHTML = '';
    const hdus = d.hdus || [];
    const isImageHdu = (x) => /image/i.test(x.kind || '');
    const hasImg = hdus.some((x) => isImageHdu(x) && x.shape && x.shape.length === 2);
    const hasSpec = hdus.some((x) => /SPECTRUM/i.test(x.extname || '') || /SPECTRUM/i.test(x.hduclas1 || ''));
    const hasLc = hdus.some((x) => /RATE|LIGHTCURVE/i.test(x.extname || '') || /LIGHTCURVE/i.test(x.hduclas1 || ''));
    const hasEvt = hdus.some((x) => /EVENT|Q0|Q1/i.test(x.extname || '') || /EVENTS/i.test(x.hduclas1 || ''));
    const addQ = (label, act) => {
      const b = document.createElement('button');
      b.className = 'btn btn-sm';
      b.textContent = label;
      b.onclick = () => openProductInTab(act, d.path);
      qbox.appendChild(b);
    };
    if (hasLc) addQ('Open Light Curve', 'lightcurve');
    if (hasSpec) addQ('Open Spectrum', 'spectrum');
    if (hasImg) addQ('Open in 2D Imaging', 'image');
    if (hasEvt) addQ('Fold in Pulse Alignment', 'pulse');

    // HDU table
    $('insp-hdu-count').textContent = `${d.num_hdus} HDUs`;
    const hduBody = $('insp-hdu-body');
    hduBody.innerHTML = hdus.map((x) => {
      const sel = x.index === state.inspHdu ? 'selected' : '';
      const shapeStr = isImageHdu(x) ? (x.shape && x.shape.length ? x.shape.join('×') : 'empty') : `${(x.nrows || 0).toLocaleString()} rows`;
      const cs = x.has_checksum
        ? (x.checksum_ok && x.datasum_ok ? '<span class="badge badge-ok">OK</span>' : '<span class="badge badge-err">FAIL</span>')
        : '<span class="muted">—</span>';
      return `<tr class="file-row ${sel}" data-hdu="${x.index}">
        <td class="num">${x.index}</td>
        <td class="mono">${esc(x.extname || 'PRIMARY')}</td>
        <td><span class="badge">${esc(x.kind)}</span></td>
        <td class="num">${esc(shapeStr)}</td>
        <td class="num">${Number.isFinite(x.exposure) ? x.exposure.toFixed(2) : '—'}</td>
        <td>${cs}</td>
      </tr>`;
    }).join('');
    hduBody.querySelectorAll('tr').forEach((tr) => {
      tr.onclick = () => inspectFitsFile(d.path, Number(tr.dataset.hdu), 0);
    });

    // Selected HDU cards & sample rows
    const sh = d.selected_hdu;
    if (sh) {
      $('insp-selected-hdu-label').textContent = `HDU #${sh.index} · ${sh.extname || 'PRIMARY'} (${sh.kind})`;
      renderInspectorCards();
      renderInspectorRows(sh.table);
    }
  }

  function renderInspectorCards() {
    const sh = state.inspData && state.inspData.selected_hdu;
    const body = $('insp-cards-body');
    if (!sh || !sh.cards) {
      body.innerHTML = `<tr><td colspan="3" class="muted">No header cards</td></tr>`;
      return;
    }
    const q = $('insp-card-search').value.trim().toLowerCase();
    const cards = sh.cards.filter((c) => !q ||
      (c.keyword || '').toLowerCase().includes(q) ||
      (c.value || '').toLowerCase().includes(q) ||
      (c.comment || '').toLowerCase().includes(q));
    body.innerHTML = cards.map((c) => `<tr>
      <td class="mono" style="font-weight:600;color:var(--accent);">${esc(c.keyword)}</td>
      <td class="mono">${esc(c.value)}</td>
      <td class="muted">${esc(c.comment)}</td>
    </tr>`).join('') || `<tr><td colspan="3" class="muted">No matching cards</td></tr>`;
  }

  function renderInspectorRows(tbl) {
    const head = $('insp-rows-head');
    const body = $('insp-rows-body');
    if (!tbl || !tbl.columns || !tbl.columns.length) {
      $('insp-rows-info').textContent = 'Selected HDU is not a BINTABLE (no tabular rows).';
      head.innerHTML = '';
      body.innerHTML = '';
      $('insp-rows-prev').disabled = true;
      $('insp-rows-next').disabled = true;
      return;
    }
    const end = Math.min(tbl.total_rows, tbl.offset + tbl.rows.length);
    $('insp-rows-info').textContent = `Rows ${tbl.offset + 1}–${end} of ${tbl.total_rows.toLocaleString()}`;
    $('insp-rows-prev').disabled = tbl.offset <= 0;
    $('insp-rows-next').disabled = end >= tbl.total_rows;

    head.innerHTML = `<tr><th>#</th>${tbl.columns.map((c) => `<th>${esc(c)}</th>`).join('')}</tr>`;
    body.innerHTML = tbl.rows.map((r, ri) =>
      `<tr><td class="num muted">${tbl.offset + ri + 1}</td>${r.map((cell) => `<td class="mono">${esc(cell)}</td>`).join('')}</tr>`
    ).join('');
  }

  // =========================================================================
  // Initialization & Event Wiring
  // =========================================================================
  async function initHealth() {
    try {
      const h = await apiGet('/api/health');
      state.health = h;
      $('version-pill').textContent = `v${h.version || '0.1'}`;
      $('workspace-text').textContent = `${basename(h.data_dir) || h.data_dir} · ${h.threads || 1}T`;
      $('workspace-pill').title = `Data dir: ${h.data_dir} (${h.threads} threads)`;
      $('status-dot').classList.remove('offline');

      // Populate samples
      populateLcSampleDropdown();
      const imgSel = $('img-sample-select');
      for (const f of (h.presets && h.presets.sample_images) || []) {
        const opt = document.createElement('option');
        opt.value = f;
        opt.textContent = basename(f);
        opt.title = f;
        imgSel.appendChild(opt);
      }
      if (h.presets && h.presets.laxpc_caldb) {
        $('spec-caldb-input').value = h.presets.laxpc_caldb;
      }

      // Initialize default Pulse preset & Inspector browser & CLI preview
      applyPulsePreset('laxpc4');
      autofillRunnerSample();
      loadInspectorBrowser(h.data_dir || '');
      refreshJobs();
    } catch (err) {
      $('status-dot').classList.add('offline');
      $('workspace-text').textContent = 'Offline';
    }
  }

  function wireEvents() {
    // Theme toggle
    $('theme-toggle').onclick = () => applyTheme(state.theme === 'dark' ? 'light' : 'dark');

    // Main 5-tab nav + keyboard shortcuts 1..5
    document.querySelectorAll('.nav-tab').forEach((btn) => {
      btn.onclick = () => switchTab(btn.dataset.tab);
    });
    window.addEventListener('keydown', (ev) => {
      if (ev.target && /INPUT|TEXTAREA|SELECT/.test(ev.target.tagName)) return;
      const map = { '1': 'runner', '2': 'lc-spec', '3': 'imaging', '4': 'pulse', '5': 'inspector' };
      if (map[ev.key]) switchTab(map[ev.key]);
    });

    // Shared browse buttons
    document.querySelectorAll('.browse-btn').forEach((btn) => {
      btn.onclick = () => openFileModal({
        targetId: btn.dataset.target,
        mode: btn.dataset.mode || 'file',
      });
    });

    // Modal events
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

    // Tab 1: Runner
    document.querySelectorAll('#runner-pipeline-seg button').forEach((b) => {
      b.onclick = () => {
        state.runnerPipe = b.dataset.pipe;
        document.querySelectorAll('#runner-pipeline-seg button').forEach((x) => x.classList.toggle('active', x === b));
        $('form-czti').style.display = state.runnerPipe === 'czti' ? '' : 'none';
        $('form-laxpc').style.display = state.runnerPipe === 'laxpc' ? '' : 'none';
        $('form-bary').style.display = state.runnerPipe === 'bary' ? '' : 'none';
        updateCliPreview();
      };
    });
    $('tab-runner').addEventListener('input', scheduleCliPreview);
    $('tab-runner').addEventListener('change', scheduleCliPreview);
    $('runner-autofill').onclick = autofillRunnerSample;
    $('runner-submit').onclick = submitJob;
    $('jobs-refresh-btn').onclick = refreshJobs;
    const copyWithFeedback = (btn, text) => {
      if (navigator.clipboard) navigator.clipboard.writeText(text).catch(() => {});
      const orig = btn.textContent;
      btn.textContent = 'Copied!';
      setTimeout(() => { btn.textContent = orig; }, 1100);
    };
    $('copy-cli-btn').onclick = (e) => copyWithFeedback(e.currentTarget, $('preview-cli').textContent);
    $('copy-toml-btn').onclick = (e) => copyWithFeedback(e.currentTarget, $('preview-toml').textContent);

    // Tab 2: Light Curves & Spectra
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

    // Tab 3: 2D Imaging
    $('img-load-btn').onclick = () => loadImageFile();
    $('img-sample-select').onchange = (e) => {
      if (e.target.value) {
        $('img-path-input').value = e.target.value;
        loadImageFile(e.target.value);
        e.target.value = '';
      }
    };
    $('img-hdu-select').onchange = (e) => {
      if (state.imgData) loadImageFile(state.imgData.path, e.target.value);
    };
    $('img-cmap').onchange = (e) => state.imgPlot.setOptions({ colormap: e.target.value });
    $('img-stretch').onchange = (e) => state.imgPlot.setOptions({ stretch: e.target.value });
    $('img-quads').onchange = (e) => state.imgPlot.setOptions({ showQuads: e.target.checked });
    $('img-peak').onchange = (e) => state.imgPlot.setOptions({ showPeak: e.target.checked });
    $('img-reset-btn').onclick = () => state.imgPlot.resetView();
    $('img-png-btn').onclick = () => state.imgPlot.exportPNG({ filename: 'astrosat_2d_image.png' });

    // Tab 4: Pulse Alignment
    $('pulse-search-toggle').onclick = () => {
      const box = $('pulse-search-box');
      box.style.display = box.style.display === 'none' ? '' : 'none';
    };
    $('pulse-run-search').onclick = runPeriodSearch;
    $('pgram-close').onclick = () => { $('periodogram-card').style.display = 'none'; };
    $('preset-laxpc4').onclick = () => applyPulsePreset('laxpc4');
    $('preset-czti3').onclick = () => applyPulsePreset('czti3');
    $('preset-cross').onclick = () => applyPulsePreset('cross');
    $('pulse-add-input').onclick = () => {
      state.pulseInputs.push({
        label: `Input ${state.pulseInputs.length + 1}`,
        path: '',
        caldb: '',
        use_weights: false,
        bands: [{ label: 'Full Band', emin_kev: '', emax_kev: '', unit: '' }],
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

    // Tab 5: FITS Inspector
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
        $('insp-cards-panel').style.display = sub === 'cards' ? '' : 'none';
        $('insp-rows-panel').style.display = sub === 'rows' ? '' : 'none';
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
    state.lcPlot = new C.SeriesPlot($('lc-canvas'), $('lc-tooltip'), { theme: state.theme });
    state.imgPlot = new C.ImagePlot2D($('img-canvas'), { theme: state.theme, onHover: handleImageHover });
    state.pulsePlot = new C.PulseAlignmentPlot($('pulse-canvas'), $('pulse-tooltip'), { theme: state.theme });
    state.pgramPlot = new C.SeriesPlot($('pgram-canvas'), $('pgram-tooltip'), { theme: state.theme, showErrors: false });

    applyTheme(state.theme);
    wireEvents();
    initHealth();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else {
    boot();
  }
})();
