/* astrosat-ui — Zero-dependency offline HiDPI Canvas + SVG Scientific Charting Engine */
(() => {
  'use strict';

  // Colorblind-safe Okabe-Ito-inspired palette tuned for dark & light backgrounds
  const PALETTE_DARK = [
    '#56b4e9', // sky blue
    '#e69f00', // amber
    '#009e73', // teal green
    '#f07178', // coral red
    '#cc79a7', // magenta
    '#f5e062', // warm yellow
    '#8be9fd', // cyan
    '#a78bfa', // soft violet
  ];

  const PALETTE_LIGHT = [
    '#0072b2', // deep blue
    '#d55e00', // vermilion
    '#009e73', // teal green
    '#cc79a7', // reddish purple
    '#e69f00', // orange
    '#56b4e9', // sky blue
    '#6d28d9', // violet
    '#334155', // slate
  ];

  function getPalette(theme) {
    return theme === 'light' ? PALETTE_LIGHT : PALETTE_DARK;
  }

  function getThemeColors(theme) {
    if (theme === 'light') {
      return {
        bg: '#ffffff',
        plotBg: '#fafbfd',
        grid: 'rgba(15, 23, 42, 0.08)',
        axis: '#94a3b8',
        text: '#1e293b',
        muted: '#64748b',
        crosshair: 'rgba(37, 99, 235, 0.55)',
        gti: 'rgba(5, 150, 105, 0.11)',
        gtiBorder: 'rgba(5, 150, 105, 0.35)',
        refLine: 'rgba(220, 38, 38, 0.72)',
        quadLine: 'rgba(255, 255, 255, 0.80)',
        quadShadow: 'rgba(0, 0, 0, 0.65)',
        badPix: [226, 232, 240],
      };
    }
    return {
      bg: '#0b0d12',
      plotBg: '#0e1118',
      grid: 'rgba(148, 163, 184, 0.11)',
      axis: '#475569',
      text: '#e2e8f0',
      muted: '#8b96ae',
      crosshair: 'rgba(96, 165, 250, 0.65)',
      gti: 'rgba(52, 211, 153, 0.10)',
      gtiBorder: 'rgba(52, 211, 153, 0.32)',
      refLine: 'rgba(248, 113, 113, 0.78)',
      quadLine: 'rgba(255, 255, 255, 0.72)',
      quadShadow: 'rgba(0, 0, 0, 0.75)',
      badPix: [24, 29, 42],
    };
  }

  function resolveThemeColors(activeTheme, optsTheme, themeColors) {
    const key = activeTheme === 'light' ? 'light' : 'dark';
    const base = getThemeColors(key);
    if (!themeColors || typeof themeColors !== 'object') return base;
    if ('light' in themeColors || 'dark' in themeColors) {
      const scoped = themeColors[key];
      return scoped && typeof scoped === 'object' ? Object.assign({}, base, scoped) : base;
    }
    if (!activeTheme || activeTheme === optsTheme) {
      return Object.assign({}, base, themeColors);
    }
    return base;
  }

  // Exact perceptual colormap control points (matplotlib Viridis, Inferno, Cividis, Grayscale)
  const CMAP_STOPS = {
    viridis: [
      [68, 1, 84], [72, 24, 106], [71, 45, 123], [63, 66, 134],
      [53, 86, 140], [43, 105, 142], [35, 123, 142], [29, 140, 139],
      [31, 158, 137], [53, 175, 127], [90, 192, 108], [134, 207, 82],
      [181, 222, 43], [223, 227, 24], [253, 231, 37],
    ],
    inferno: [
      [0, 0, 4], [14, 11, 38], [40, 11, 84], [73, 11, 106],
      [106, 23, 110], [138, 34, 106], [170, 47, 94], [198, 65, 76],
      [221, 89, 55], [238, 119, 34], [249, 153, 16], [252, 189, 31],
      [246, 223, 82], [244, 244, 150], [252, 255, 164],
    ],
    cividis: [
      [0, 32, 77], [0, 44, 106], [18, 56, 111], [47, 68, 110],
      [70, 80, 108], [91, 93, 107], [112, 107, 108], [133, 122, 109],
      [154, 138, 109], [175, 155, 106], [197, 173, 98], [219, 191, 86],
      [239, 211, 69], [253, 234, 69],
    ],
    grayscale: [
      [0, 0, 0], [255, 255, 255],
    ],
  };

  function buildLut(stops) {
    const lut = new Uint8Array(256 * 3);
    const n = stops.length - 1;
    for (let i = 0; i < 256; i++) {
      const pos = (i / 255) * n;
      const idx = Math.min(n - 1, Math.floor(pos));
      const f = pos - idx;
      const a = stops[idx];
      const b = stops[idx + 1];
      lut[i * 3] = Math.round(a[0] + (b[0] - a[0]) * f);
      lut[i * 3 + 1] = Math.round(a[1] + (b[1] - a[1]) * f);
      lut[i * 3 + 2] = Math.round(a[2] + (b[2] - a[2]) * f);
    }
    return lut;
  }

  const COLORMAPS = {
    viridis: buildLut(CMAP_STOPS.viridis),
    inferno: buildLut(CMAP_STOPS.inferno),
    cividis: buildLut(CMAP_STOPS.cividis),
    grayscale: buildLut(CMAP_STOPS.grayscale),
  };

  function stretchValue(v, vmin, vmax, mode) {
    if (!Number.isFinite(v)) return NaN;
    const span = vmax - vmin;
    let t = span > 0 ? (v - vmin) / span : 0.5;
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    if (mode === 'log') {
      return Math.log1p(999 * t) / Math.log(1000);
    }
    if (mode === 'asinh') {
      return Math.asinh(10 * t) / Math.asinh(10);
    }
    return t; // linear or p1_p99 (vmin/vmax already set to p1/p99)
  }

  function formatSci(v, sig = 4) {
    if (!Number.isFinite(v)) return '—';
    if (v === 0) return '0';
    const av = Math.abs(v);
    if (av >= 1e5 || (av > 0 && av < 1e-3)) {
      return v.toExponential(Math.max(1, sig - 1)).replace(/e\+?/, 'e');
    }
    // Strip trailing zeros cleanly
    const s = Number(v.toPrecision(sig)).toString();
    return s;
  }

  function formatTick(v, step) {
    if (!Number.isFinite(v)) return '';
    if (Math.abs(v) < 1e-12 * Math.max(1, Math.abs(step || 1))) return '0';
    const av = Math.abs(v);
    if (av >= 1e5 || (av > 0 && av < 1e-3)) {
      return v.toExponential(1).replace(/\.0e/, 'e').replace(/e\+/, 'e');
    }
    if (step && step > 0) {
      const dec = Math.max(0, Math.min(6, -Math.floor(Math.log10(step)) + 1));
      const fixed = v.toFixed(dec);
      return fixed.includes('.') ? fixed.replace(/\.?0+$/, '') : fixed;
    }
    return formatSci(v, 4);
  }

  function niceStep(range, targetCount) {
    if (!Number.isFinite(range) || range <= 0) return 1;
    const rough = range / Math.max(2, targetCount);
    const pow10 = Math.pow(10, Math.floor(Math.log10(rough)));
    const frac = rough / pow10;
    let nice;
    if (frac <= 1.2) nice = 1;
    else if (frac <= 2.2) nice = 2;
    else if (frac <= 3.5) nice = 2.5;
    else if (frac <= 7) nice = 5;
    else nice = 10;
    return nice * pow10;
  }

  function niceTicks(min, max, targetCount = 6) {
    if (!Number.isFinite(min) || !Number.isFinite(max)) return { ticks: [], step: 1 };
    if (max <= min) return { ticks: [min], step: 1 };
    const step = niceStep(max - min, targetCount);
    const start = Math.ceil(min / step - 1e-9) * step;
    const ticks = [];
    for (let v = start; v <= max + step * 1e-9 && ticks.length < 32; v += step) {
      const clean = Math.abs(v) < step * 1e-9 ? 0 : v;
      ticks.push(clean);
    }
    return { ticks, step };
  }

  function logTicks(minVal, maxVal) {
    const lo = Math.max(minVal, 1e-30);
    const hi = Math.max(maxVal, lo * 1.1);
    const pMin = Math.floor(Math.log10(lo));
    const pMax = Math.ceil(Math.log10(hi));
    const decades = pMax - pMin;
    const mults = decades <= 2 ? [1, 2, 5] : [1];
    const ticks = [];
    for (let p = pMin; p <= pMax && ticks.length < 40; p++) {
      const base = Math.pow(10, p);
      for (const m of mults) {
        const v = m * base;
        if (v >= lo * 0.9999 && v <= hi * 1.0001) {
          ticks.push(v);
        }
      }
    }
    return { ticks, step: 0 };
  }

  function escapeXml(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  function triggerDownload(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    setTimeout(() => {
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }, 150);
  }

  // Exact FITS 0-based pixel <-> Canvas pixel conversions (FITS y=0 is BOTTOM row!)
  function canvasToFitsPixel(mx, my, imgRect, width, height) {
    if (!imgRect || imgRect.w <= 0 || imgRect.h <= 0 || width <= 0 || height <= 0) {
      return null;
    }
    const u = (mx - imgRect.x) / imgRect.w;
    const v = (my - imgRect.y) / imgRect.h;
    if (u < 0 || u >= 1 || v < 0 || v >= 1) return null;
    const fx = Math.min(width - 1, Math.max(0, Math.floor(u * width)));
    const fy = Math.min(height - 1, Math.max(0, (height - 1) - Math.floor(v * height)));
    return { fx, fy };
  }

  function fitsPixelToCanvas(fx, fy, imgRect, width, height) {
    const cx = imgRect.x + ((fx + 0.5) / width) * imgRect.w;
    const cy = imgRect.y + (((height - 1 - fy) + 0.5) / height) * imgRect.h;
    return { cx, cy };
  }

  // Exact FITS WCS conversion from 0-based (fx, fy) to world coordinates
  function fitsPixelToWorld(fx, fy, wcs) {
    if (!wcs) return null;
    const crpix1 = Number(wcs.crpix1);
    const crpix2 = Number(wcs.crpix2);
    const crval1 = Number(wcs.crval1);
    const crval2 = Number(wcs.crval2);
    if (![crpix1, crpix2, crval1, crval2].every(Number.isFinite)) return null;

    // FITS 1-based pixel coordinates
    const dx = (fx + 1) - crpix1;
    const dy = (fy + 1) - crpix2;

    const unpackMat2x2 = (m) => {
      if (!Array.isArray(m)) return null;
      if (m.length === 2 && Array.isArray(m[0]) && Array.isArray(m[1]) && m[0].length === 2 && m[1].length === 2) {
        const flat = [Number(m[0][0]), Number(m[0][1]), Number(m[1][0]), Number(m[1][1])];
        return flat.every(Number.isFinite) ? flat : null;
      }
      if (m.length === 4) {
        const flat = m.map(Number);
        return flat.every(Number.isFinite) ? flat : null;
      }
      return null;
    };

    let m11, m12, m21, m22;
    const cdMat = unpackMat2x2(wcs.cd);
    if (cdMat) {
      [m11, m12, m21, m22] = cdMat;
    } else if ([wcs.cd1_1, wcs.cd1_2, wcs.cd2_1, wcs.cd2_2].every(Number.isFinite)) {
      m11 = wcs.cd1_1; m12 = wcs.cd1_2; m21 = wcs.cd2_1; m22 = wcs.cd2_2;
    } else {
      const cdelt1 = Number.isFinite(wcs.cdelt1) ? wcs.cdelt1 : 1;
      const cdelt2 = Number.isFinite(wcs.cdelt2) ? wcs.cdelt2 : 1;
      const pcMat = unpackMat2x2(wcs.pc);
      if (pcMat) {
        m11 = cdelt1 * pcMat[0];
        m12 = cdelt1 * pcMat[1];
        m21 = cdelt2 * pcMat[2];
        m22 = cdelt2 * pcMat[3];
      } else if ([wcs.pc1_1, wcs.pc1_2, wcs.pc2_1, wcs.pc2_2].every(Number.isFinite)) {
        m11 = cdelt1 * wcs.pc1_1;
        m12 = cdelt1 * wcs.pc1_2;
        m21 = cdelt2 * wcs.pc2_1;
        m22 = cdelt2 * wcs.pc2_2;
      } else {
        const rot = Number.isFinite(wcs.crota2) ? (wcs.crota2 * Math.PI) / 180 : 0;
        const c = Math.cos(rot);
        const s = Math.sin(rot);
        m11 = cdelt1 * c;
        m12 = -cdelt2 * s;
        m21 = cdelt1 * s;
        m22 = cdelt2 * c;
      }
    }

    const xi = m11 * dx + m12 * dy;
    const eta = m21 * dx + m22 * dy;

    const ctype1 = String(wcs.ctype1 || '').toUpperCase();
    const ctype2 = String(wcs.ctype2 || '').toUpperCase();
    const isTan = ctype1.includes('TAN') || ctype2.includes('TAN');

    if (isTan) {
      const deg2rad = Math.PI / 180;
      const ra0 = crval1 * deg2rad;
      const dec0 = crval2 * deg2rad;
      const xRad = xi * deg2rad;
      const yRad = eta * deg2rad;
      const r = Math.hypot(xRad, yRad);
      let raDeg = crval1;
      let decDeg = crval2;
      if (r > 0) {
        const c = Math.atan(r);
        const sinC = Math.sin(c);
        const cosC = Math.cos(c);
        const sinD0 = Math.sin(dec0);
        const cosD0 = Math.cos(dec0);
        const dec = Math.asin(Math.max(-1, Math.min(1, cosC * sinD0 + (yRad * sinC * cosD0) / r)));
        const ra = ra0 + Math.atan2(xRad * sinC, r * cosD0 * cosC - yRad * sinD0 * sinC);
        raDeg = ((ra / deg2rad) % 360 + 360) % 360;
        decDeg = dec / deg2rad;
      }
      return {
        x: raDeg,
        y: decDeg,
        labelX: 'RA',
        labelY: 'Dec',
        unit: 'deg',
        isRaDec: true,
        offsetX: xi,
        offsetY: eta,
      };
    }

    return {
      x: crval1 + xi,
      y: crval2 + eta,
      labelX: ctype1 || 'θx',
      labelY: ctype2 || 'θy',
      unit: wcs.cunit1 || 'deg',
      isRaDec: false,
      offsetX: xi,
      offsetY: eta,
    };
  }

  // =========================================================================
  // 1. Interactive 1D Scientific Plotter (Light Curves, Spectra, Periodograms)
  // =========================================================================
  class SeriesPlot {
    constructor(canvas, tooltipEl, opts = {}) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.tooltipEl = tooltipEl || null;
      this.opts = Object.assign({
        theme: 'dark',
        themeColors: null,
        title: '',
        xLabel: 'X',
        yLabel: 'Y',
        xUnit: '',
        yUnit: '',
        logX: false,
        logY: false,
        showErrors: true,
        showGti: true,
        stepMode: false,
        margins: { top: 28, right: 20, bottom: 44, left: 68 },
      }, opts);

      this.series = []; // { id, label, color, x[], y[], yerr[], step?, lineWidth?, markers? }
      this.gtis = [];   // [[start, stop], ...]
      this.vlines = []; // [{ x, label, color, dash }]
      this.view = null; // { xMin, xMax, yMin, yMax }
      this.dataBounds = null;
      this.hover = null;
      this.dragging = false;
      this.dragStart = null;
      this.dragView = null;
      this.width = 600;
      this.height = 400;

      this._bindEvents();
      this._ro = new ResizeObserver(() => this.resize());
      if (canvas.parentElement) this._ro.observe(canvas.parentElement);
      this.resize();
    }

    setTheme(theme) {
      this.opts.theme = theme;
      this.render();
    }

    setOptions(patch) {
      const logChanged = ('logX' in patch && patch.logX !== this.opts.logX) ||
                         ('logY' in patch && patch.logY !== this.opts.logY);
      Object.assign(this.opts, patch);
      if (logChanged) {
        this._computeDataBounds();
        this.resetView();
      } else {
        this.render();
      }
    }

    setData({ series = [], gtis = [], vlines = [], title, xLabel, yLabel, xUnit, yUnit }) {
      this.series = series;
      this.gtis = gtis;
      this.vlines = vlines;
      if (title !== undefined) this.opts.title = title;
      if (xLabel !== undefined) this.opts.xLabel = xLabel;
      if (yLabel !== undefined) this.opts.yLabel = yLabel;
      if (xUnit !== undefined) this.opts.xUnit = xUnit;
      if (yUnit !== undefined) this.opts.yUnit = yUnit;
      this._computeDataBounds();
      this.resetView();
    }

    clear() {
      this.series = [];
      this.gtis = [];
      this.vlines = [];
      this.dataBounds = null;
      this.view = null;
      this.hover = null;
      if (this.tooltipEl) this.tooltipEl.classList.add('hidden');
      this.render();
    }

    _computeDataBounds() {
      let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
      const logX = this.opts.logX;
      const logY = this.opts.logY;
      for (const s of this.series) {
        if (s.hidden) continue;
        const n = Math.min(s.x.length, s.y.length);
        for (let i = 0; i < n; i++) {
          const x = s.x[i];
          const y = s.y[i];
          if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
          if (logX && x <= 0) continue;
          if (logY && y <= 0) continue;
          if (x < xMin) xMin = x;
          if (x > xMax) xMax = x;
          const e = (this.opts.showErrors && s.yerr && Number.isFinite(s.yerr[i])) ? s.yerr[i] : 0;
          const yLo = logY ? Math.max(y * 0.5, y - e) : (y - e);
          const yHi = y + e;
          if (yLo < yMin) yMin = yLo;
          if (yHi > yMax) yMax = yHi;
        }
      }
      if (!Number.isFinite(xMin) || !Number.isFinite(yMin)) {
        this.dataBounds = null;
        return;
      }
      if (logX) {
        xMin = Math.max(1e-30, xMin);
        xMax = Math.max(xMin * 1.1, xMax);
        const pad = Math.pow(xMax / xMin, 0.04);
        xMin /= pad;
        xMax *= pad;
      } else {
        const dx = (xMax - xMin) || Math.max(1, Math.abs(xMin) * 0.1);
        xMin -= dx * 0.025;
        xMax += dx * 0.025;
      }
      if (logY) {
        yMin = Math.max(1e-30, yMin);
        yMax = Math.max(yMin * 1.2, yMax);
        const pad = Math.pow(yMax / yMin, 0.08);
        yMin /= pad;
        yMax *= pad;
      } else {
        const dy = (yMax - yMin) || Math.max(1, Math.abs(yMin) * 0.1);
        yMin -= dy * 0.07;
        yMax += dy * 0.09;
      }
      this.dataBounds = { xMin, xMax, yMin, yMax };
    }

    resetView() {
      this.view = this.dataBounds ? { ...this.dataBounds } : null;
      this.render();
    }

    resize() {
      const parent = this.canvas.parentElement;
      if (!parent) return;
      const rect = parent.getBoundingClientRect();
      if (rect.width < 10 || rect.height < 10) return;
      const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
      this.width = Math.round(rect.width);
      this.height = Math.round(rect.height);
      this.canvas.width = Math.round(this.width * dpr);
      this.canvas.height = Math.round(this.height * dpr);
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.render();
    }

    _plotRect() {
      const m = this.opts.margins;
      return {
        x: m.left,
        y: m.top,
        w: Math.max(20, this.width - m.left - m.right),
        h: Math.max(20, this.height - m.top - m.bottom),
      };
    }

    _xToPx(x, pr, v) {
      if (this.opts.logX) {
        const a = Math.log10(Math.max(1e-30, v.xMin));
        const b = Math.log10(Math.max(1e-30, v.xMax));
        return pr.x + ((Math.log10(Math.max(1e-30, x)) - a) / (b - a || 1)) * pr.w;
      }
      return pr.x + ((x - v.xMin) / (v.xMax - v.xMin || 1)) * pr.w;
    }

    _yToPx(y, pr, v) {
      if (this.opts.logY) {
        const a = Math.log10(Math.max(1e-30, v.yMin));
        const b = Math.log10(Math.max(1e-30, v.yMax));
        return pr.y + pr.h - ((Math.log10(Math.max(1e-30, y)) - a) / (b - a || 1)) * pr.h;
      }
      return pr.y + pr.h - ((y - v.yMin) / (v.yMax - v.yMin || 1)) * pr.h;
    }

    _pxToX(px, pr, v) {
      const f = (px - pr.x) / (pr.w || 1);
      if (this.opts.logX) {
        const a = Math.log10(Math.max(1e-30, v.xMin));
        const b = Math.log10(Math.max(1e-30, v.xMax));
        return Math.pow(10, a + f * (b - a));
      }
      return v.xMin + f * (v.xMax - v.xMin);
    }

    _pxToY(py, pr, v) {
      const f = (pr.y + pr.h - py) / (pr.h || 1);
      if (this.opts.logY) {
        const a = Math.log10(Math.max(1e-30, v.yMin));
        const b = Math.log10(Math.max(1e-30, v.yMax));
        return Math.pow(10, a + f * (b - a));
      }
      return v.yMin + f * (v.yMax - v.yMin);
    }

    _bindEvents() {
      const c = this.canvas;
      c.addEventListener('wheel', (ev) => {
        if (!this.view) return;
        ev.preventDefault();
        const rect = c.getBoundingClientRect();
        const mx = ev.clientX - rect.left;
        const my = ev.clientY - rect.top;
        const pr = this._plotRect();
        const fx = Math.max(0, Math.min(1, (mx - pr.x) / pr.w));
        const fy = Math.max(0, Math.min(1, (pr.y + pr.h - my) / pr.h));
        const factor = ev.deltaY > 0 ? 1.18 : 1 / 1.18;
        const v = this.view;

        if (this.opts.logX) {
          const a = Math.log10(v.xMin), b = Math.log10(v.xMax);
          const pivot = a + fx * (b - a);
          const span = (b - a) * factor;
          v.xMin = Math.pow(10, pivot - fx * span);
          v.xMax = Math.pow(10, pivot + (1 - fx) * span);
        } else {
          const pivot = v.xMin + fx * (v.xMax - v.xMin);
          const span = (v.xMax - v.xMin) * factor;
          v.xMin = pivot - fx * span;
          v.xMax = pivot + (1 - fx) * span;
        }
        if (this.opts.logY) {
          const a = Math.log10(v.yMin), b = Math.log10(v.yMax);
          const pivot = a + fy * (b - a);
          const span = (b - a) * factor;
          v.yMin = Math.pow(10, pivot - fy * span);
          v.yMax = Math.pow(10, pivot + (1 - fy) * span);
        } else {
          const pivot = v.yMin + fy * (v.yMax - v.yMin);
          const span = (v.yMax - v.yMin) * factor;
          v.yMin = pivot - fy * span;
          v.yMax = pivot + (1 - fy) * span;
        }
        this.render();
      }, { passive: false });

      c.addEventListener('mousedown', (ev) => {
        if (!this.view || ev.button !== 0) return;
        this.dragging = true;
        this.dragStart = { x: ev.clientX, y: ev.clientY };
        this.dragView = { ...this.view };
      });

      window.addEventListener('mousemove', (ev) => {
        if (this.dragging && this.dragView) {
          const pr = this._plotRect();
          const dxPx = ev.clientX - this.dragStart.x;
          const dyPx = ev.clientY - this.dragStart.y;
          const dv = this.dragView;
          if (this.opts.logX) {
            const a = Math.log10(dv.xMin), b = Math.log10(dv.xMax);
            const shift = -(dxPx / pr.w) * (b - a);
            this.view.xMin = Math.pow(10, a + shift);
            this.view.xMax = Math.pow(10, b + shift);
          } else {
            const shift = -(dxPx / pr.w) * (dv.xMax - dv.xMin);
            this.view.xMin = dv.xMin + shift;
            this.view.xMax = dv.xMax + shift;
          }
          if (this.opts.logY) {
            const a = Math.log10(dv.yMin), b = Math.log10(dv.yMax);
            const shift = (dyPx / pr.h) * (b - a);
            this.view.yMin = Math.pow(10, a + shift);
            this.view.yMax = Math.pow(10, b + shift);
          } else {
            const shift = (dyPx / pr.h) * (dv.yMax - dv.yMin);
            this.view.yMin = dv.yMin + shift;
            this.view.yMax = dv.yMax + shift;
          }
          this.render();
          return;
        }
        const rect = c.getBoundingClientRect();
        if (ev.clientX < rect.left || ev.clientX > rect.right || ev.clientY < rect.top || ev.clientY > rect.bottom) {
          return;
        }
        this._handleHover(ev.clientX - rect.left, ev.clientY - rect.top);
      });

      window.addEventListener('mouseup', () => {
        this.dragging = false;
        this.dragView = null;
      });

      c.addEventListener('mouseleave', () => {
        this.hover = null;
        if (this.tooltipEl) this.tooltipEl.classList.add('hidden');
        this.render();
      });

      c.addEventListener('dblclick', () => {
        this.resetView();
      });
    }

    _handleHover(mx, my) {
      if (!this.view || this.series.length === 0) return;
      const pr = this._plotRect();
      if (mx < pr.x || mx > pr.x + pr.w || my < pr.y || my > pr.y + pr.h) {
        this.hover = null;
        if (this.tooltipEl) this.tooltipEl.classList.add('hidden');
        this.render();
        return;
      }
      const xVal = this._pxToX(mx, pr, this.view);
      const yVal = this._pxToY(my, pr, this.view);
      const hits = [];
      for (const s of this.series) {
        if (s.hidden || !s.x.length) continue;
        // Binary search assuming sorted X, fallback to linear if unsorted
        let bestIdx = -1;
        let bestDist = Infinity;
        const n = s.x.length;
        if (n > 2 && s.x[0] <= s.x[n - 1]) {
          let lo = 0, hi = n - 1;
          while (lo < hi) {
            const mid = (lo + hi) >> 1;
            if (s.x[mid] < xVal) lo = mid + 1;
            else hi = mid;
          }
          for (const cand of [lo - 1, lo, lo + 1]) {
            if (cand >= 0 && cand < n && Number.isFinite(s.x[cand]) && Number.isFinite(s.y[cand])) {
              const d = Math.abs(s.x[cand] - xVal);
              if (d < bestDist) { bestDist = d; bestIdx = cand; }
            }
          }
        } else {
          for (let i = 0; i < n; i++) {
            if (!Number.isFinite(s.x[i]) || !Number.isFinite(s.y[i])) continue;
            const d = Math.abs(s.x[i] - xVal);
            if (d < bestDist) { bestDist = d; bestIdx = i; }
          }
        }
        if (bestIdx >= 0) {
          hits.push({
            label: s.label,
            color: s.color,
            x: s.x[bestIdx],
            y: s.y[bestIdx],
            yerr: s.yerr ? s.yerr[bestIdx] : null,
          });
        }
      }
      this.hover = { mx, my, xVal, yVal, hits };
      this._updateTooltip();
      this.render();
    }

    _updateTooltip() {
      const el = this.tooltipEl;
      if (!el) return;
      if (!this.hover || !this.hover.hits.length) {
        el.classList.add('hidden');
        return;
      }
      const { mx, my, hits } = this.hover;
      const xUnit = this.opts.xUnit ? ` ${this.opts.xUnit}` : '';
      const yUnit = this.opts.yUnit ? ` ${this.opts.yUnit}` : '';
      const xRef = hits[0].x;
      let html = `<div class="tt-head">${escapeXml(this.opts.xLabel)}: ${formatSci(xRef, 6)}${escapeXml(xUnit)}</div>`;
      for (const h of hits.slice(0, 6)) {
        const errStr = (h.yerr != null && Number.isFinite(h.yerr)) ? ` ± ${formatSci(h.yerr, 3)}` : '';
        html += `<div class="tt-row"><span><span class="tt-swatch" style="background:${h.color}"></span>${escapeXml(h.label)}</span><span>${formatSci(h.y, 4)}${errStr}${escapeXml(yUnit)}</span></div>`;
      }
      el.innerHTML = html;
      el.classList.remove('hidden');
      const pad = 12;
      const tw = el.offsetWidth || 210;
      const th = el.offsetHeight || 70;
      let left = mx + pad;
      let top = my + pad;
      if (left + tw > this.width - 8) left = Math.max(8, mx - tw - pad);
      if (top + th > this.height - 8) top = Math.max(8, my - th - pad);
      el.style.left = `${left}px`;
      el.style.top = `${top}px`;
    }

    render(targetCtx = null, width = this.width, height = this.height, overrideTheme = null) {
      const ctx = targetCtx || this.ctx;
      const tc = resolveThemeColors(overrideTheme || this.opts.theme, this.opts.theme, this.opts.themeColors);
      ctx.save();
      ctx.fillStyle = tc.bg;
      ctx.fillRect(0, 0, width, height);

      if (!this.view || !this.series.length) {
        ctx.restore();
        return;
      }

      const m = this.opts.margins;
      const pr = {
        x: m.left,
        y: m.top,
        w: Math.max(20, width - m.left - m.right),
        h: Math.max(20, height - m.top - m.bottom),
      };
      const v = this.view;

      // Plot background
      ctx.fillStyle = tc.plotBg;
      ctx.fillRect(pr.x, pr.y, pr.w, pr.h);

      // Ticks & grid
      const xt = this.opts.logX ? logTicks(v.xMin, v.xMax) : niceTicks(v.xMin, v.xMax, Math.max(4, Math.floor(pr.w / 100)));
      const yt = this.opts.logY ? logTicks(v.yMin, v.yMax) : niceTicks(v.yMin, v.yMax, Math.max(4, Math.floor(pr.h / 65)));

      ctx.strokeStyle = tc.grid;
      ctx.lineWidth = 1;
      for (const tx of xt.ticks) {
        const px = Math.round(this._xToPx(tx, pr, v)) + 0.5;
        if (px < pr.x || px > pr.x + pr.w) continue;
        ctx.beginPath(); ctx.moveTo(px, pr.y); ctx.lineTo(px, pr.y + pr.h); ctx.stroke();
      }
      for (const ty of yt.ticks) {
        const py = Math.round(this._yToPx(ty, pr, v)) + 0.5;
        if (py < pr.y || py > pr.y + pr.h) continue;
        ctx.beginPath(); ctx.moveTo(pr.x, py); ctx.lineTo(pr.x + pr.w, py); ctx.stroke();
      }

      // Clipped plot interior
      ctx.save();
      ctx.beginPath();
      ctx.rect(pr.x, pr.y, pr.w, pr.h);
      ctx.clip();

      // GTI shaded regions
      if (this.opts.showGti && this.gtis && this.gtis.length) {
        ctx.fillStyle = tc.gti;
        for (const [g0, g1] of this.gtis) {
          if (!Number.isFinite(g0) || !Number.isFinite(g1)) continue;
          const x0 = this._xToPx(g0, pr, v);
          const x1 = this._xToPx(g1, pr, v);
          if (x1 < pr.x || x0 > pr.x + pr.w) continue;
          ctx.fillRect(x0, pr.y, Math.max(1, x1 - x0), pr.h);
        }
      }

      // Vertical reference lines
      for (const vl of this.vlines) {
        if (!Number.isFinite(vl.x)) continue;
        const px = this._xToPx(vl.x, pr, v);
        if (px < pr.x || px > pr.x + pr.w) continue;
        ctx.save();
        ctx.strokeStyle = vl.color || tc.refLine;
        ctx.lineWidth = 1.25;
        ctx.setLineDash(vl.dash || [5, 4]);
        ctx.beginPath(); ctx.moveTo(px, pr.y); ctx.lineTo(px, pr.y + pr.h); ctx.stroke();
        if (vl.label) {
          ctx.fillStyle = vl.color || tc.refLine;
          ctx.font = '11px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
          ctx.fillText(vl.label, px + 5, pr.y + 14);
        }
        ctx.restore();
      }

      // Data series
      for (const s of this.series) {
        if (s.hidden) continue;
        const n = Math.min(s.x.length, s.y.length);
        if (!n) continue;
        const stepMode = s.step !== undefined ? s.step : this.opts.stepMode;

        // Error bars first
        if (this.opts.showErrors && s.yerr && s.yerr.length === n) {
          ctx.strokeStyle = s.color;
          ctx.globalAlpha = 0.45;
          ctx.lineWidth = 1;
          ctx.beginPath();
          for (let i = 0; i < n; i++) {
            const x = s.x[i], y = s.y[i], e = s.yerr[i];
            if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(e) || e <= 0) continue;
            if (this.opts.logX && x <= 0) continue;
            if (this.opts.logY && y <= 0) continue;
            const px = this._xToPx(x, pr, v);
            if (px < pr.x - 4 || px > pr.x + pr.w + 4) continue;
            const yLo = this.opts.logY ? Math.max(v.yMin * 0.5, y - e) : (y - e);
            const yHi = y + e;
            const pyLo = this._yToPx(yLo, pr, v);
            const pyHi = this._yToPx(yHi, pr, v);
            ctx.moveTo(px, pyLo);
            ctx.lineTo(px, pyHi);
          }
          ctx.stroke();
          ctx.globalAlpha = 1;
        }

        // Line / step path
        ctx.strokeStyle = s.color;
        ctx.lineWidth = s.lineWidth || 1.6;
        ctx.beginPath();
        let started = false;
        let prevPx = 0, prevPy = 0;
        for (let i = 0; i < n; i++) {
          const x = s.x[i], y = s.y[i];
          if (!Number.isFinite(x) || !Number.isFinite(y) || (this.opts.logX && x <= 0) || (this.opts.logY && y <= 0)) {
            started = false;
            continue;
          }
          const px = this._xToPx(x, pr, v);
          const py = this._yToPx(y, pr, v);
          if (!started) {
            ctx.moveTo(px, py);
            started = true;
          } else if (stepMode) {
            const midX = (prevPx + px) * 0.5;
            ctx.lineTo(midX, prevPy);
            ctx.lineTo(midX, py);
            ctx.lineTo(px, py);
          } else {
            ctx.lineTo(px, py);
          }
          prevPx = px; prevPy = py;
        }
        ctx.stroke();

        // Point markers when sparse
        if (s.markers || n <= 120) {
          ctx.fillStyle = s.color;
          const r = n <= 60 ? 2.5 : 1.8;
          for (let i = 0; i < n; i++) {
            const x = s.x[i], y = s.y[i];
            if (!Number.isFinite(x) || !Number.isFinite(y) || (this.opts.logX && x <= 0) || (this.opts.logY && y <= 0)) continue;
            const px = this._xToPx(x, pr, v);
            if (px < pr.x - 4 || px > pr.x + pr.w + 4) continue;
            const py = this._yToPx(y, pr, v);
            ctx.beginPath();
            ctx.arc(px, py, r, 0, Math.PI * 2);
            ctx.fill();
          }
        }
      }

      // Crosshair on live canvas only
      if (!targetCtx && this.hover) {
        ctx.save();
        ctx.strokeStyle = tc.crosshair;
        ctx.lineWidth = 1;
        ctx.setLineDash([4, 3]);
        ctx.beginPath();
        ctx.moveTo(this.hover.mx, pr.y); ctx.lineTo(this.hover.mx, pr.y + pr.h);
        ctx.moveTo(pr.x, this.hover.my); ctx.lineTo(pr.x + pr.w, this.hover.my);
        ctx.stroke();
        ctx.restore();
      }

      ctx.restore(); // end clip

      // Plot border
      ctx.strokeStyle = tc.axis;
      ctx.lineWidth = 1;
      ctx.strokeRect(pr.x + 0.5, pr.y + 0.5, pr.w - 1, pr.h - 1);

      // Tick labels
      ctx.fillStyle = tc.muted;
      ctx.font = '11px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      for (const tx of xt.ticks) {
        const px = this._xToPx(tx, pr, v);
        if (px < pr.x + 8 || px > pr.x + pr.w - 8) continue;
        ctx.fillText(formatTick(tx, xt.step), px, pr.y + pr.h + 6);
      }

      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      for (const ty of yt.ticks) {
        const py = this._yToPx(ty, pr, v);
        if (py < pr.y + 6 || py > pr.y + pr.h - 6) continue;
        ctx.fillText(formatTick(ty, yt.step), pr.x - 7, py);
      }

      // Axis labels & Title
      ctx.fillStyle = tc.text;
      ctx.font = '500 12px system-ui, -apple-system, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      const xTitle = this.opts.xUnit ? `${this.opts.xLabel} (${this.opts.xUnit})` : this.opts.xLabel;
      ctx.fillText(xTitle, pr.x + pr.w / 2, height - 4);

      ctx.save();
      ctx.translate(14, pr.y + pr.h / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      const yTitle = this.opts.yUnit ? `${this.opts.yLabel} (${this.opts.yUnit})` : this.opts.yLabel;
      ctx.fillText(yTitle, 0, 0);
      ctx.restore();

      if (this.opts.title) {
        ctx.fillStyle = tc.text;
        ctx.font = '600 12.5px system-ui, -apple-system, sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        ctx.fillText(this.opts.title, pr.x, 7);
      }

      // Compact top-right legend when 1+ series
      const visSeries = this.series.filter((s) => !s.hidden);
      if (visSeries.length > 0) {
        ctx.font = '11px system-ui, -apple-system, sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'middle';
        let lx = pr.x + pr.w - 10;
        const ly = 14;
        for (let i = visSeries.length - 1; i >= 0; i--) {
          const s = visSeries[i];
          const tw = ctx.measureText(s.label).width;
          lx -= tw;
          if (lx < pr.x + 180) break;
          ctx.fillStyle = tc.muted;
          ctx.fillText(s.label, lx, ly);
          lx -= 14;
          ctx.fillStyle = s.color;
          ctx.fillRect(lx, ly - 3, 9, 6);
          lx -= 12;
        }
      }

      ctx.restore();
    }

    exportSVG({ filename = 'plot.svg', theme = null, width = 960, height = 540 } = {}) {
      if (!this.view || !this.series.length) return;
      const th = theme || this.opts.theme;
      const tc = resolveThemeColors(th, this.opts.theme, this.opts.themeColors);
      const m = this.opts.margins;
      const pr = { x: m.left, y: m.top + 6, w: width - m.left - m.right, h: height - m.top - m.bottom - 6 };
      const v = this.view;
      const xt = this.opts.logX ? logTicks(v.xMin, v.xMax) : niceTicks(v.xMin, v.xMax, 7);
      const yt = this.opts.logY ? logTicks(v.yMin, v.yMax) : niceTicks(v.yMin, v.yMax, 6);

      const parts = [];
      parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`);
      parts.push(`<defs><clipPath id="plotClip"><rect x="${pr.x}" y="${pr.y}" width="${pr.w}" height="${pr.h}"/></clipPath></defs>`);
      parts.push(`<rect width="${width}" height="${height}" fill="${tc.bg}"/>`);
      parts.push(`<rect x="${pr.x}" y="${pr.y}" width="${pr.w}" height="${pr.h}" fill="${tc.plotBg}" stroke="${tc.axis}" stroke-width="1"/>`);

      // Grid & ticks
      for (const tx of xt.ticks) {
        const px = this._xToPx(tx, pr, v);
        if (px < pr.x || px > pr.x + pr.w) continue;
        parts.push(`<line x1="${px.toFixed(2)}" y1="${pr.y}" x2="${px.toFixed(2)}" y2="${pr.y + pr.h}" stroke="${tc.grid}"/>`);
        parts.push(`<text x="${px.toFixed(2)}" y="${pr.y + pr.h + 16}" fill="${tc.muted}" font-family="monospace" font-size="11" text-anchor="middle">${escapeXml(formatTick(tx, xt.step))}</text>`);
      }
      for (const ty of yt.ticks) {
        const py = this._yToPx(ty, pr, v);
        if (py < pr.y || py > pr.y + pr.h) continue;
        parts.push(`<line x1="${pr.x}" y1="${py.toFixed(2)}" x2="${pr.x + pr.w}" y2="${py.toFixed(2)}" stroke="${tc.grid}"/>`);
        parts.push(`<text x="${pr.x - 8}" y="${(py + 4).toFixed(2)}" fill="${tc.muted}" font-family="monospace" font-size="11" text-anchor="end">${escapeXml(formatTick(ty, yt.step))}</text>`);
      }

      parts.push(`<g clip-path="url(#plotClip)">`);
      if (this.opts.showGti && this.gtis.length) {
        for (const [g0, g1] of this.gtis) {
          const x0 = this._xToPx(g0, pr, v);
          const x1 = this._xToPx(g1, pr, v);
          if (x1 < pr.x || x0 > pr.x + pr.w) continue;
          parts.push(`<rect x="${x0.toFixed(2)}" y="${pr.y}" width="${Math.max(1, x1 - x0).toFixed(2)}" height="${pr.h}" fill="${tc.gti}"/>`);
        }
      }

      for (const vl of this.vlines) {
        const px = this._xToPx(vl.x, pr, v);
        parts.push(`<line x1="${px.toFixed(2)}" y1="${pr.y}" x2="${px.toFixed(2)}" y2="${pr.y + pr.h}" stroke="${vl.color || tc.refLine}" stroke-dasharray="5,4" stroke-width="1.2"/>`);
      }

      for (const s of this.series) {
        if (s.hidden) continue;
        const n = Math.min(s.x.length, s.y.length);
        const stepMode = s.step !== undefined ? s.step : this.opts.stepMode;
        if (this.opts.showErrors && s.yerr && s.yerr.length === n) {
          let errD = '';
          for (let i = 0; i < n; i++) {
            const x = s.x[i], y = s.y[i], e = s.yerr[i];
            if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(e) || e <= 0) continue;
            if ((this.opts.logX && x <= 0) || (this.opts.logY && y <= 0)) continue;
            const px = this._xToPx(x, pr, v);
            const yLo = this.opts.logY ? Math.max(v.yMin * 0.5, y - e) : (y - e);
            const pyLo = this._yToPx(yLo, pr, v);
            const pyHi = this._yToPx(y + e, pr, v);
            errD += `M${px.toFixed(1)},${pyLo.toFixed(1)}V${pyHi.toFixed(1)}`;
          }
          if (errD) parts.push(`<path d="${errD}" stroke="${s.color}" stroke-opacity="0.45" stroke-width="1" fill="none"/>`);
        }
        let d = '';
        let started = false;
        let prevPx = 0, prevPy = 0;
        for (let i = 0; i < n; i++) {
          const x = s.x[i], y = s.y[i];
          if (!Number.isFinite(x) || !Number.isFinite(y) || (this.opts.logX && x <= 0) || (this.opts.logY && y <= 0)) {
            started = false; continue;
          }
          const px = this._xToPx(x, pr, v);
          const py = this._yToPx(y, pr, v);
          if (!started) {
            d += `M${px.toFixed(2)},${py.toFixed(2)}`;
            started = true;
          } else if (stepMode) {
            const mx = ((prevPx + px) * 0.5).toFixed(2);
            d += `H${mx}V${py.toFixed(2)}H${px.toFixed(2)}`;
          } else {
            d += `L${px.toFixed(2)},${py.toFixed(2)}`;
          }
          prevPx = px; prevPy = py;
        }
        if (d) parts.push(`<path d="${d}" stroke="${s.color}" stroke-width="${s.lineWidth || 1.6}" fill="none"/>`);
      }
      parts.push(`</g>`);

      const xTitle = this.opts.xUnit ? `${this.opts.xLabel} (${this.opts.xUnit})` : this.opts.xLabel;
      const yTitle = this.opts.yUnit ? `${this.opts.yLabel} (${this.opts.yUnit})` : this.opts.yLabel;
      parts.push(`<text x="${(pr.x + pr.w / 2).toFixed(1)}" y="${height - 8}" fill="${tc.text}" font-family="sans-serif" font-size="12" font-weight="500" text-anchor="middle">${escapeXml(xTitle)}</text>`);
      parts.push(`<text x="16" y="${(pr.y + pr.h / 2).toFixed(1)}" fill="${tc.text}" font-family="sans-serif" font-size="12" font-weight="500" text-anchor="middle" transform="rotate(-90 16 ${(pr.y + pr.h / 2).toFixed(1)})">${escapeXml(yTitle)}</text>`);
      if (this.opts.title) {
        parts.push(`<text x="${pr.x}" y="20" fill="${tc.text}" font-family="sans-serif" font-size="13" font-weight="600">${escapeXml(this.opts.title)}</text>`);
      }

      // Legend
      const vis = this.series.filter((s) => !s.hidden);
      let lx = pr.x + pr.w;
      for (let i = vis.length - 1; i >= 0; i--) {
        const s = vis[i];
        const estW = s.label.length * 6.8 + 20;
        lx -= estW;
        if (lx < pr.x + 220) break;
        parts.push(`<rect x="${lx}" y="11" width="10" height="7" fill="${s.color}"/>`);
        parts.push(`<text x="${lx + 14}" y="18" fill="${tc.muted}" font-family="sans-serif" font-size="11">${escapeXml(s.label)}</text>`);
      }

      parts.push(`</svg>`);
      triggerDownload(new Blob([parts.join('\n')], { type: 'image/svg+xml;charset=utf-8' }), filename);
    }

    exportPNG({ filename = 'plot.png', theme = null, scale = 2 } = {}) {
      if (!this.view || !this.series.length) return;
      const w = Math.max(800, this.width);
      const h = Math.max(460, this.height);
      const off = document.createElement('canvas');
      off.width = w * scale;
      off.height = h * scale;
      const octx = off.getContext('2d');
      octx.setTransform(scale, 0, 0, scale, 0, 0);
      this.render(octx, w, h, theme || this.opts.theme);
      off.toBlob((blob) => {
        if (blob) triggerDownload(blob, filename);
      }, 'image/png');
    }
  }

  // =========================================================================
  // 2. Interactive 2D Scientific FITS Image Viewer (DPH / DPI / Sky Maps)
  // =========================================================================
  class ImagePlot2D {
    constructor(canvas, opts = {}) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.opts = Object.assign({
        theme: 'dark',
        themeColors: null,
        colormap: 'viridis',
        stretch: 'linear', // linear | log | asinh | p1_p99
        showQuads: true,
        showPeak: true,
        onHover: null,
      }, opts);

      this.img = null;       // { width, height, data, stats, wcs, title, peak }
      this.offscreen = null; // offscreen canvas holding rendered FITS pixels (already Y-inverted!)
      this.zoom = 1;
      this.panX = 0;
      this.panY = 0;
      this.lastImgRect = null;
      this.hoverPix = null;
      this.dragging = false;
      this.dragStart = null;
      this.width = 600;
      this.height = 500;

      this._bindEvents();
      this._ro = new ResizeObserver(() => this.resize());
      if (canvas.parentElement) this._ro.observe(canvas.parentElement);
      this.resize();
    }

    setTheme(theme) {
      this.opts.theme = theme;
      if (this.img) this._rebuildBitmap();
      this.render();
    }

    setOptions(patch) {
      const needBitmap = ('colormap' in patch) || ('stretch' in patch) || ('themeColors' in patch);
      Object.assign(this.opts, patch);
      if (needBitmap && this.img) this._rebuildBitmap();
      this.render();
    }

    setImage(img) {
      this.img = img;
      this.zoom = 1;
      this.panX = 0;
      this.panY = 0;
      this.hoverPix = null;
      if (img) this._rebuildBitmap();
      this.render();
    }

    resetView() {
      this.zoom = 1;
      this.panX = 0;
      this.panY = 0;
      this.render();
    }

    resize() {
      const parent = this.canvas.parentElement;
      if (!parent) return;
      const rect = parent.getBoundingClientRect();
      if (rect.width < 10 || rect.height < 10) return;
      const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
      this.width = Math.round(rect.width);
      this.height = Math.round(rect.height);
      this.canvas.width = Math.round(this.width * dpr);
      this.canvas.height = Math.round(this.height * dpr);
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.render();
    }

    _valueRange() {
      if (!this.img || !this.img.stats) return [0, 1];
      const st = this.img.stats;
      if (this.opts.stretch === 'p1_p99' && Number.isFinite(st.p1) && Number.isFinite(st.p99) && st.p99 > st.p1) {
        return [st.p1, st.p99];
      }
      const lo = Number.isFinite(st.min) ? st.min : 0;
      const hi = Number.isFinite(st.max) && st.max > lo ? st.max : lo + 1;
      return [lo, hi];
    }

    _rebuildBitmap() {
      const { width, height, data } = this.img;
      const off = document.createElement('canvas');
      off.width = width;
      off.height = height;
      const octx = off.getContext('2d');
      const idata = octx.createImageData(width, height);
      const buf = idata.data;
      const lut = COLORMAPS[this.opts.colormap] || COLORMAPS.viridis;
      const tc = resolveThemeColors(this.opts.theme, this.opts.theme, this.opts.themeColors);
      const [vmin, vmax] = this._valueRange();
      const mode = this.opts.stretch === 'p1_p99' ? 'linear' : this.opts.stretch;

      // CRITICAL FITS CONVENTION: FITS row y=0 is BOTTOM; canvas imgY = (height - 1) - y
      for (let fy = 0; fy < height; fy++) {
        const imgY = (height - 1) - fy;
        const srcRow = fy * width;
        const dstRow = imgY * width * 4;
        for (let fx = 0; fx < width; fx++) {
          const val = data[srcRow + fx];
          const off4 = dstRow + fx * 4;
          if (!Number.isFinite(val)) {
            buf[off4] = tc.badPix[0];
            buf[off4 + 1] = tc.badPix[1];
            buf[off4 + 2] = tc.badPix[2];
            buf[off4 + 3] = 255;
          } else {
            const t = stretchValue(val, vmin, vmax, mode);
            const idx = Math.min(255, Math.max(0, Math.round(t * 255))) * 3;
            buf[off4] = lut[idx];
            buf[off4 + 1] = lut[idx + 1];
            buf[off4 + 2] = lut[idx + 2];
            buf[off4 + 3] = 255;
          }
        }
      }
      octx.putImageData(idata, 0, 0);
      this.offscreen = off;
    }

    _computeLayout(width, height) {
      const m = { top: 30, right: 86, bottom: 42, left: 58 };
      const availW = Math.max(40, width - m.left - m.right);
      const availH = Math.max(40, height - m.top - m.bottom);
      const aspect = this.img ? (this.img.width / this.img.height) : 1;
      let baseW = availW;
      let baseH = availW / aspect;
      if (baseH > availH) {
        baseH = availH;
        baseW = availH * aspect;
      }
      const box = {
        x: m.left + (availW - baseW) / 2,
        y: m.top + (availH - baseH) / 2,
        w: baseW,
        h: baseH,
      };
      const imgRect = {
        x: box.x + (box.w * (1 - this.zoom)) / 2 + this.panX,
        y: box.y + (box.h * (1 - this.zoom)) / 2 + this.panY,
        w: box.w * this.zoom,
        h: box.h * this.zoom,
      };
      const cbar = {
        x: box.x + box.w + 18,
        y: box.y,
        w: 14,
        h: box.h,
      };
      return { box, imgRect, cbar };
    }

    _bindEvents() {
      const c = this.canvas;
      c.addEventListener('wheel', (ev) => {
        if (!this.img) return;
        ev.preventDefault();
        const rect = c.getBoundingClientRect();
        const mx = ev.clientX - rect.left;
        const my = ev.clientY - rect.top;
        const { imgRect } = this._computeLayout(this.width, this.height);
        const factor = ev.deltaY > 0 ? 1 / 1.2 : 1.2;
        const nextZoom = Math.max(1, Math.min(24, this.zoom * factor));
        const ratio = nextZoom / this.zoom;
        // Zoom toward cursor
        const cx = imgRect.x + imgRect.w / 2;
        const cy = imgRect.y + imgRect.h / 2;
        this.panX = this.panX * ratio + (mx - cx) * (1 - ratio);
        this.panY = this.panY * ratio + (my - cy) * (1 - ratio);
        this.zoom = nextZoom;
        if (this.zoom <= 1.001) { this.zoom = 1; this.panX = 0; this.panY = 0; }
        this.render();
      }, { passive: false });

      c.addEventListener('mousedown', (ev) => {
        if (!this.img || ev.button !== 0) return;
        this.dragging = true;
        this.dragStart = { x: ev.clientX, y: ev.clientY, panX: this.panX, panY: this.panY };
      });

      window.addEventListener('mousemove', (ev) => {
        if (this.dragging && this.dragStart) {
          this.panX = this.dragStart.panX + (ev.clientX - this.dragStart.x);
          this.panY = this.dragStart.panY + (ev.clientY - this.dragStart.y);
          this.render();
          return;
        }
        if (!this.img || !this.lastImgRect) return;
        const rect = c.getBoundingClientRect();
        const mx = ev.clientX - rect.left;
        const my = ev.clientY - rect.top;
        const { box } = this._computeLayout(this.width, this.height);
        if (mx < box.x || mx > box.x + box.w || my < box.y || my > box.y + box.h) {
          if (this.hoverPix) {
            this.hoverPix = null;
            if (this.opts.onHover) this.opts.onHover(null);
            this.render();
          }
          return;
        }
        const pix = canvasToFitsPixel(mx, my, this.lastImgRect, this.img.width, this.img.height);
        if (!pix) return;
        const val = this.img.data[pix.fy * this.img.width + pix.fx];
        const world = fitsPixelToWorld(pix.fx, pix.fy, this.img.wcs);
        this.hoverPix = { ...pix, val, world };
        if (this.opts.onHover) this.opts.onHover(this.hoverPix);
        this.render();
      });

      window.addEventListener('mouseup', () => {
        this.dragging = false;
        this.dragStart = null;
      });

      c.addEventListener('mouseleave', () => {
        this.hoverPix = null;
        if (this.opts.onHover) this.opts.onHover(null);
        this.render();
      });

      c.addEventListener('dblclick', () => this.resetView());
    }

    render(targetCtx = null, width = this.width, height = this.height, overrideTheme = null) {
      const ctx = targetCtx || this.ctx;
      const tc = resolveThemeColors(overrideTheme || this.opts.theme, this.opts.theme, this.opts.themeColors);
      ctx.save();
      ctx.fillStyle = tc.bg;
      ctx.fillRect(0, 0, width, height);

      if (!this.img || !this.offscreen) {
        ctx.restore();
        return;
      }

      const { box, imgRect, cbar } = this._computeLayout(width, height);
      if (!targetCtx) this.lastImgRect = imgRect;

      // Clip to display viewport box
      ctx.save();
      ctx.beginPath();
      ctx.rect(box.x, box.y, box.w, box.h);
      ctx.clip();

      ctx.fillStyle = tc.plotBg;
      ctx.fillRect(box.x, box.y, box.w, box.h);
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(this.offscreen, imgRect.x, imgRect.y, imgRect.w, imgRect.h);

      // CZTI 128x128 Q0..Q3 quadrant overlay (exact detxy_to_locxy geometry!)
      // Q0: x=[0,64), y=[64,128) (top-left)
      // Q1: x=[64,128), y=[64,128) (top-right)
      // Q2: x=[64,128), y=[0,64) (bottom-right)
      // Q3: x=[0,64), y=[0,64) (bottom-left)
      if (this.opts.showQuads && this.img.width === 128 && this.img.height === 128) {
        const midX = imgRect.x + imgRect.w * 0.5;
        const midY = imgRect.y + imgRect.h * 0.5;
        ctx.save();
        ctx.strokeStyle = tc.quadLine;
        ctx.lineWidth = 1.25;
        ctx.setLineDash([5, 4]);
        ctx.beginPath();
        ctx.moveTo(midX, imgRect.y); ctx.lineTo(midX, imgRect.y + imgRect.h);
        ctx.moveTo(imgRect.x, midY); ctx.lineTo(imgRect.x + imgRect.w, midY);
        ctx.stroke();
        ctx.setLineDash([]);

        const labels = [
          { name: 'Q0', x: imgRect.x + 8, y: imgRect.y + 8, align: 'left', base: 'top' },
          { name: 'Q1', x: imgRect.x + imgRect.w - 8, y: imgRect.y + 8, align: 'right', base: 'top' },
          { name: 'Q2', x: imgRect.x + imgRect.w - 8, y: imgRect.y + imgRect.h - 8, align: 'right', base: 'bottom' },
          { name: 'Q3', x: imgRect.x + 8, y: imgRect.y + imgRect.h - 8, align: 'left', base: 'bottom' },
        ];
        ctx.font = '700 12px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
        for (const lb of labels) {
          ctx.textAlign = lb.align;
          ctx.textBaseline = lb.base;
          ctx.fillStyle = tc.quadShadow;
          ctx.fillText(lb.name, lb.x + 1, lb.y + 1);
          ctx.fillStyle = '#ffffff';
          ctx.fillText(lb.name, lb.x, lb.y);
        }
        ctx.restore();
      }

      // Peak SNR crosshair
      if (this.opts.showPeak && this.img.peak && Number.isFinite(this.img.peak.x) && Number.isFinite(this.img.peak.y)) {
        const { cx, cy } = fitsPixelToCanvas(this.img.peak.x, this.img.peak.y, imgRect, this.img.width, this.img.height);
        ctx.save();
        ctx.strokeStyle = tc.peak || '#fbbf24';
        ctx.lineWidth = 1.5;
        const r = 8;
        ctx.beginPath();
        ctx.arc(cx, cy, r, 0, Math.PI * 2);
        ctx.moveTo(cx - r - 6, cy); ctx.lineTo(cx - r + 2, cy);
        ctx.moveTo(cx + r - 2, cy); ctx.lineTo(cx + r + 6, cy);
        ctx.moveTo(cx, cy - r - 6); ctx.lineTo(cx, cy - r + 2);
        ctx.moveTo(cx, cy + r - 2); ctx.lineTo(cx, cy + r + 6);
        ctx.stroke();
        ctx.restore();
      }

      // Live hover pixel highlight
      if (!targetCtx && this.hoverPix) {
        const { cx, cy } = fitsPixelToCanvas(this.hoverPix.fx, this.hoverPix.fy, imgRect, this.img.width, this.img.height);
        ctx.save();
        ctx.strokeStyle = tc.hoverLine || 'rgba(255,255,255,0.7)';
        ctx.lineWidth = 1;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(cx, box.y); ctx.lineTo(cx, box.y + box.h);
        ctx.moveTo(box.x, cy); ctx.lineTo(box.x + box.w, cy);
        ctx.stroke();
        ctx.restore();
      }

      ctx.restore(); // end image clip

      // Box frame
      ctx.strokeStyle = tc.axis;
      ctx.lineWidth = 1;
      ctx.strokeRect(box.x + 0.5, box.y + 0.5, box.w - 1, box.h - 1);

      // Pixel / WCS axis ticks along box edges
      const visibleFx0 = Math.max(0, ((box.x - imgRect.x) / imgRect.w) * this.img.width);
      const visibleFx1 = Math.min(this.img.width, ((box.x + box.w - imgRect.x) / imgRect.w) * this.img.width);
      const visibleFy1 = Math.min(this.img.height, (1 - (box.y - imgRect.y) / imgRect.h) * this.img.height);
      const visibleFy0 = Math.max(0, (1 - (box.y + box.h - imgRect.y) / imgRect.h) * this.img.height);

      const xTicks = niceTicks(visibleFx0, visibleFx1, 5);
      const yTicks = niceTicks(visibleFy0, visibleFy1, 5);
      ctx.fillStyle = tc.muted;
      ctx.font = '11px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';

      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      for (const fx of xTicks.ticks) {
        const px = imgRect.x + (fx / this.img.width) * imgRect.w;
        if (px < box.x + 4 || px > box.x + box.w - 4) continue;
        ctx.beginPath(); ctx.moveTo(px, box.y + box.h); ctx.lineTo(px, box.y + box.h + 4); ctx.stroke();
        ctx.fillText(String(Math.round(fx)), px, box.y + box.h + 6);
      }

      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      for (const fy of yTicks.ticks) {
        const py = imgRect.y + (1 - fy / this.img.height) * imgRect.h;
        if (py < box.y + 4 || py > box.y + box.h - 4) continue;
        ctx.beginPath(); ctx.moveTo(box.x - 4, py); ctx.lineTo(box.x, py); ctx.stroke();
        ctx.fillText(String(Math.round(fy)), box.x - 7, py);
      }

      // Axis labels & Title
      ctx.fillStyle = tc.text;
      ctx.font = '500 11.5px system-ui, -apple-system, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText('FITS X (pixel)', box.x + box.w / 2, height - 4);

      ctx.save();
      ctx.translate(14, box.y + box.h / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      ctx.fillText('FITS Y (pixel, bottom = 0)', 0, 0);
      ctx.restore();

      if (this.img.title) {
        ctx.fillStyle = tc.text;
        ctx.font = '600 12.5px system-ui, -apple-system, sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'top';
        ctx.fillText(this.img.title, box.x, 8);
      }

      // Vertical Colorbar
      const lut = COLORMAPS[this.opts.colormap] || COLORMAPS.viridis;
      const steps = Math.max(32, Math.round(cbar.h));
      for (let i = 0; i < steps; i++) {
        const t = 1 - i / (steps - 1);
        const idx = Math.min(255, Math.max(0, Math.round(t * 255))) * 3;
        ctx.fillStyle = `rgb(${lut[idx]},${lut[idx + 1]},${lut[idx + 2]})`;
        ctx.fillRect(cbar.x, cbar.y + (i / steps) * cbar.h, cbar.w, Math.ceil(cbar.h / steps) + 0.5);
      }
      ctx.strokeStyle = tc.axis;
      ctx.strokeRect(cbar.x + 0.5, cbar.y + 0.5, cbar.w - 1, cbar.h - 1);

      const [vmin, vmax] = this._valueRange();
      const cbLabels = [
        { y: cbar.y, val: vmax },
        { y: cbar.y + cbar.h * 0.5, val: (vmin + vmax) * 0.5 },
        { y: cbar.y + cbar.h, val: vmin },
      ];
      ctx.fillStyle = tc.muted;
      ctx.font = '10.5px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';
      for (const cl of cbLabels) {
        ctx.fillText(formatSci(cl.val, 3), cbar.x + cbar.w + 6, cl.y);
      }

      ctx.restore();
    }

    exportPNG({ filename = 'astrosat_image.png', theme = null, scale = 2 } = {}) {
      if (!this.img) return;
      const w = Math.max(760, this.width);
      const h = Math.max(620, this.height);
      const off = document.createElement('canvas');
      off.width = w * scale;
      off.height = h * scale;
      const octx = off.getContext('2d');
      octx.setTransform(scale, 0, 0, scale, 0, 0);
      this.render(octx, w, h, theme || this.opts.theme);
      off.toBlob((blob) => {
        if (blob) triggerDownload(blob, filename);
      }, 'image/png');
    }
  }

  // =========================================================================
  // 3. Crown Jewel: Multi-Band Phase-Folded Pulse Alignment Plotter (2 Cycles)
  // =========================================================================
  class PulseAlignmentPlot {
    constructor(canvas, tooltipEl, opts = {}) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.tooltipEl = tooltipEl || null;
      this.opts = Object.assign({
        theme: 'dark',
        themeColors: null,
        mode: 'stacked', // 'stacked' | 'overlay'
        showErrors: true,
        showPeakGuide: true,
        title: 'Multi-Band Phase-Folded Pulse Profile',
        subtitle: '',
      }, opts);

      this.bands = []; // [{ label, color, phase[], values[], errors[], peakPhase, z2, chi2Red, pf, lagCycles, lagMs, counts }]
      this.refPeakPhase = null;
      this.hoverPhase = null;
      this.hoverPos = null;
      this.width = 700;
      this.height = 500;

      this._bindEvents();
      this._ro = new ResizeObserver(() => this.resize());
      if (canvas.parentElement) this._ro.observe(canvas.parentElement);
      this.resize();
    }

    setTheme(theme) {
      this.opts.theme = theme;
      this.render();
    }

    setOptions(patch) {
      Object.assign(this.opts, patch);
      this._adjustHeight();
      this.render();
    }

    setData({ bands = [], refPeakPhase = null, title, subtitle }) {
      this.bands = bands;
      this.refPeakPhase = refPeakPhase;
      if (title !== undefined) this.opts.title = title;
      if (subtitle !== undefined) this.opts.subtitle = subtitle;
      this._adjustHeight();
      this.render();
    }

    clear() {
      this.bands = [];
      this.refPeakPhase = null;
      this.hoverPhase = null;
      if (this.tooltipEl) this.tooltipEl.classList.add('hidden');
      this.render();
    }

    _adjustHeight() {
      const parent = this.canvas.parentElement;
      if (!parent) return;
      if (this.opts.mode === 'stacked' && this.bands.length > 0) {
        const targetH = Math.max(420, 88 + this.bands.length * 128);
        parent.style.height = `${targetH}px`;
      } else {
        parent.style.height = '470px';
      }
      this.resize();
    }

    resize() {
      const parent = this.canvas.parentElement;
      if (!parent) return;
      const rect = parent.getBoundingClientRect();
      if (rect.width < 10 || rect.height < 10) return;
      const dpr = Math.min(window.devicePixelRatio || 1, 2.5);
      this.width = Math.round(rect.width);
      this.height = Math.round(rect.height);
      this.canvas.width = Math.round(this.width * dpr);
      this.canvas.height = Math.round(this.height * dpr);
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.render();
    }

    _margins() {
      return { top: 44, right: 20, bottom: 44, left: 72 };
    }

    _bindEvents() {
      const c = this.canvas;
      c.addEventListener('mousemove', (ev) => {
        if (!this.bands.length) return;
        const rect = c.getBoundingClientRect();
        const mx = ev.clientX - rect.left;
        const my = ev.clientY - rect.top;
        const m = this._margins();
        const pw = this.width - m.left - m.right;
        const ph = this.height - m.top - m.bottom;
        if (mx < m.left || mx > m.left + pw || my < m.top || my > m.top + ph) {
          this.hoverPhase = null;
          if (this.tooltipEl) this.tooltipEl.classList.add('hidden');
          this.render();
          return;
        }
        this.hoverPhase = ((mx - m.left) / pw) * 2.0; // [0, 2)
        this.hoverPos = { mx, my };
        this._updateTooltip();
        this.render();
      });

      c.addEventListener('mouseleave', () => {
        this.hoverPhase = null;
        if (this.tooltipEl) this.tooltipEl.classList.add('hidden');
        this.render();
      });
    }

    _updateTooltip() {
      const el = this.tooltipEl;
      if (!el || this.hoverPhase == null || !this.bands.length) return;
      const phi = this.hoverPhase;
      const phiMod = ((phi % 1) + 1) % 1;
      let html = `<div class="tt-head">Phase φ = ${phi.toFixed(3)} (cycle ${phiMod.toFixed(3)})</div>`;
      for (const b of this.bands) {
        const n = b.values.length;
        if (!n) continue;
        const idx = Math.min(n - 1, Math.floor(phiMod * n));
        const val = b.values[idx];
        const err = b.errors ? b.errors[idx] : null;
        const errStr = (err != null && Number.isFinite(err)) ? ` ± ${formatSci(err, 3)}` : '';
        html += `<div class="tt-row"><span><span class="tt-swatch" style="background:${b.color}"></span>${escapeXml(b.label)}</span><span>${formatSci(val, 4)}${errStr}</span></div>`;
      }
      el.innerHTML = html;
      el.classList.remove('hidden');
      const { mx, my } = this.hoverPos;
      const tw = el.offsetWidth || 230;
      const th = el.offsetHeight || 90;
      let left = mx + 14;
      let top = my + 12;
      if (left + tw > this.width - 8) left = Math.max(8, mx - tw - 14);
      if (top + th > this.height - 8) top = Math.max(8, my - th - 12);
      el.style.left = `${left}px`;
      el.style.top = `${top}px`;
    }

    // Build 2-cycle stepped histogram coordinates (phi in [0, 2))
    _twoCycleBins(values, errors) {
      const n = values.length;
      const out = [];
      for (let cyc = 0; cyc < 2; cyc++) {
        for (let i = 0; i < n; i++) {
          out.push({
            x0: cyc + i / n,
            x1: cyc + (i + 1) / n,
            xc: cyc + (i + 0.5) / n,
            y: values[i],
            e: errors ? errors[i] : 0,
          });
        }
      }
      return out;
    }

    _normBand(b) {
      let lo = Infinity, hi = -Infinity;
      for (const v of b.values) {
        if (!Number.isFinite(v)) continue;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
      const span = hi - lo || 1;
      return {
        values: b.values.map((v) => (v - lo) / span),
        errors: (b.errors || []).map((e) => e / span),
      };
    }

    render(targetCtx = null, width = this.width, height = this.height, overrideTheme = null) {
      const ctx = targetCtx || this.ctx;
      const tc = resolveThemeColors(overrideTheme || this.opts.theme, this.opts.theme, this.opts.themeColors);
      ctx.save();
      ctx.fillStyle = tc.bg;
      ctx.fillRect(0, 0, width, height);

      if (!this.bands.length) {
        ctx.restore();
        return;
      }

      const m = this._margins();
      const pw = Math.max(40, width - m.left - m.right);
      const ph = Math.max(40, height - m.top - m.bottom);
      const xToPx = (phi) => m.left + (phi / 2.0) * pw;

      // Header title & ephemeris subtitle
      ctx.fillStyle = tc.text;
      ctx.font = '600 13px system-ui, -apple-system, sans-serif';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'top';
      ctx.fillText(this.opts.title, m.left, 8);
      if (this.opts.subtitle) {
        ctx.fillStyle = tc.muted;
        ctx.font = '11px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
        ctx.fillText(this.opts.subtitle, m.left, 25);
      }

      const xTicks = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0];

      if (this.opts.mode === 'stacked') {
        const nb = this.bands.length;
        const gap = 8;
        const panelH = Math.max(48, (ph - gap * (nb - 1)) / nb);

        for (let bi = 0; bi < nb; bi++) {
          const b = this.bands[bi];
          const pyTop = m.top + bi * (panelH + gap);
          const bins = this._twoCycleBins(b.values, b.errors);
          let yMin = Infinity, yMax = -Infinity;
          for (const bn of bins) {
            if (!Number.isFinite(bn.y)) continue;
            const e = this.opts.showErrors && Number.isFinite(bn.e) ? bn.e : 0;
            if (bn.y - e < yMin) yMin = bn.y - e;
            if (bn.y + e > yMax) yMax = bn.y + e;
          }
          if (!Number.isFinite(yMin)) { yMin = 0; yMax = 1; }
          const dy = (yMax - yMin) || Math.max(1, Math.abs(yMin) * 0.1);
          yMin -= dy * 0.08;
          yMax += dy * 0.16;
          const yToPx = (y) => pyTop + panelH - ((y - yMin) / (yMax - yMin || 1)) * panelH;

          ctx.fillStyle = tc.plotBg;
          ctx.fillRect(m.left, pyTop, pw, panelH);

          // Grid lines
          ctx.strokeStyle = tc.grid;
          ctx.lineWidth = 1;
          for (const tx of xTicks) {
            const px = Math.round(xToPx(tx)) + 0.5;
            ctx.beginPath(); ctx.moveTo(px, pyTop); ctx.lineTo(px, pyTop + panelH); ctx.stroke();
          }
          const yt = niceTicks(yMin, yMax, 3);
          for (const ty of yt.ticks) {
            const py = Math.round(yToPx(ty)) + 0.5;
            if (py < pyTop + 4 || py > pyTop + panelH - 4) continue;
            ctx.beginPath(); ctx.moveTo(m.left, py); ctx.lineTo(m.left + pw, py); ctx.stroke();
          }

          // Subtle cycle divider at phi = 1.0
          ctx.save();
          ctx.strokeStyle = tc.axis;
          ctx.globalAlpha = 0.4;
          ctx.setLineDash([2, 3]);
          const px1 = xToPx(1.0);
          ctx.beginPath(); ctx.moveTo(px1, pyTop); ctx.lineTo(px1, pyTop + panelH); ctx.stroke();
          ctx.restore();

          // Error bars + stepped pulse profile
          ctx.save();
          ctx.beginPath();
          ctx.rect(m.left, pyTop, pw, panelH);
          ctx.clip();

          if (this.opts.showErrors) {
            ctx.strokeStyle = b.color;
            ctx.globalAlpha = 0.5;
            ctx.lineWidth = 1;
            ctx.beginPath();
            for (const bn of bins) {
              if (!Number.isFinite(bn.e) || bn.e <= 0) continue;
              const px = xToPx(bn.xc);
              ctx.moveTo(px, yToPx(bn.y - bn.e));
              ctx.lineTo(px, yToPx(bn.y + bn.e));
            }
            ctx.stroke();
            ctx.globalAlpha = 1;
          }

          ctx.strokeStyle = b.color;
          ctx.lineWidth = 1.75;
          ctx.beginPath();
          for (let i = 0; i < bins.length; i++) {
            const bn = bins[i];
            const px0 = xToPx(bn.x0), px1b = xToPx(bn.x1), py = yToPx(bn.y);
            if (i === 0) ctx.moveTo(px0, py);
            else ctx.lineTo(px0, py);
            ctx.lineTo(px1b, py);
          }
          ctx.stroke();
          ctx.restore();

          // Panel border
          ctx.strokeStyle = tc.axis;
          ctx.lineWidth = 1;
          ctx.strokeRect(m.left + 0.5, pyTop + 0.5, pw - 1, panelH - 1);

          // Y tick labels
          ctx.fillStyle = tc.muted;
          ctx.font = '10.5px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
          ctx.textAlign = 'right';
          ctx.textBaseline = 'middle';
          for (const ty of yt.ticks) {
            const py = yToPx(ty);
            if (py < pyTop + 8 || py > pyTop + panelH - 8) continue;
            ctx.fillText(formatTick(ty, yt.step), m.left - 7, py);
          }

          // Clean band annotation inside top-left of panel
          ctx.font = '600 11.5px system-ui, -apple-system, sans-serif';
          ctx.textAlign = 'left';
          ctx.textBaseline = 'top';
          ctx.fillStyle = b.color;
          ctx.fillText(b.label, m.left + 9, pyTop + 6);

          const meta = [];
          if (Number.isFinite(b.z2)) meta.push(`Z²=${formatSci(b.z2, 4)}`);
          if (Number.isFinite(b.pf)) meta.push(`PF=${(b.pf * 100).toFixed(1)}%`);
          if (Number.isFinite(b.lagCycles)) {
            const sign = b.lagCycles >= 0 ? '+' : '';
            meta.push(`Δφ=${sign}${b.lagCycles.toFixed(4)}`);
          }
          if (meta.length) {
            ctx.font = '10.5px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
            ctx.textAlign = 'right';
            ctx.fillStyle = tc.muted;
            ctx.fillText(meta.join('  ·  '), m.left + pw - 9, pyTop + 6);
          }
        }

        // Full-height vertical dashed peak alignment lines across ALL stacked panels!
        if (this.opts.showPeakGuide && Number.isFinite(this.refPeakPhase)) {
          const p0 = ((this.refPeakPhase % 1) + 1) % 1;
          ctx.save();
          ctx.strokeStyle = tc.refLine;
          ctx.lineWidth = 1.25;
          ctx.setLineDash([5, 4]);
          for (const phi of [p0, p0 + 1]) {
            const px = xToPx(phi);
            ctx.beginPath();
            ctx.moveTo(px, m.top);
            ctx.lineTo(px, m.top + ph);
            ctx.stroke();
          }
          ctx.restore();
        }
      } else {
        // Normalized Overlay Mode
        const pyTop = m.top;
        const yMin = -0.1, yMax = 1.15;
        const yToPx = (y) => pyTop + ph - ((y - yMin) / (yMax - yMin)) * ph;

        ctx.fillStyle = tc.plotBg;
        ctx.fillRect(m.left, pyTop, pw, ph);

        ctx.strokeStyle = tc.grid;
        ctx.lineWidth = 1;
        for (const tx of xTicks) {
          const px = Math.round(xToPx(tx)) + 0.5;
          ctx.beginPath(); ctx.moveTo(px, pyTop); ctx.lineTo(px, pyTop + ph); ctx.stroke();
        }
        for (const ty of [0, 0.25, 0.5, 0.75, 1.0]) {
          const py = Math.round(yToPx(ty)) + 0.5;
          ctx.beginPath(); ctx.moveTo(m.left, py); ctx.lineTo(m.left + pw, py); ctx.stroke();
          ctx.fillStyle = tc.muted;
          ctx.font = '11px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
          ctx.textAlign = 'right';
          ctx.textBaseline = 'middle';
          ctx.fillText(ty.toFixed(2), m.left - 7, py);
        }

        if (this.opts.showPeakGuide && Number.isFinite(this.refPeakPhase)) {
          const p0 = ((this.refPeakPhase % 1) + 1) % 1;
          ctx.save();
          ctx.strokeStyle = tc.refLine;
          ctx.lineWidth = 1.25;
          ctx.setLineDash([5, 4]);
          for (const phi of [p0, p0 + 1]) {
            const px = xToPx(phi);
            ctx.beginPath(); ctx.moveTo(px, pyTop); ctx.lineTo(px, pyTop + ph); ctx.stroke();
          }
          ctx.restore();
        }

        ctx.save();
        ctx.beginPath(); ctx.rect(m.left, pyTop, pw, ph); ctx.clip();
        for (const b of this.bands) {
          const norm = this._normBand(b);
          const bins = this._twoCycleBins(norm.values, norm.errors);
          if (this.opts.showErrors) {
            ctx.strokeStyle = b.color;
            ctx.globalAlpha = 0.35;
            ctx.lineWidth = 1;
            ctx.beginPath();
            for (const bn of bins) {
              if (!Number.isFinite(bn.e) || bn.e <= 0) continue;
              const px = xToPx(bn.xc);
              ctx.moveTo(px, yToPx(bn.y - bn.e));
              ctx.lineTo(px, yToPx(bn.y + bn.e));
            }
            ctx.stroke();
            ctx.globalAlpha = 1;
          }
          ctx.strokeStyle = b.color;
          ctx.lineWidth = 1.8;
          ctx.beginPath();
          for (let i = 0; i < bins.length; i++) {
            const bn = bins[i];
            const px0 = xToPx(bn.x0), px1b = xToPx(bn.x1), py = yToPx(bn.y);
            if (i === 0) ctx.moveTo(px0, py);
            else ctx.lineTo(px0, py);
            ctx.lineTo(px1b, py);
          }
          ctx.stroke();
        }
        ctx.restore();

        ctx.strokeStyle = tc.axis;
        ctx.strokeRect(m.left + 0.5, pyTop + 0.5, pw - 1, ph - 1);

        // Overlay Legend
        let lx = m.left + pw - 8;
        ctx.font = '11px system-ui, -apple-system, sans-serif';
        ctx.textAlign = 'left';
        ctx.textBaseline = 'middle';
        for (let i = this.bands.length - 1; i >= 0; i--) {
          const b = this.bands[i];
          const tw = ctx.measureText(b.label).width;
          lx -= tw;
          if (lx < m.left + 200) break;
          ctx.fillStyle = tc.muted;
          ctx.fillText(b.label, lx, 16);
          lx -= 14;
          ctx.fillStyle = b.color;
          ctx.fillRect(lx, 13, 9, 6);
          lx -= 12;
        }
      }

      // Synchronized hover phase cursor across ALL panels
      if (!targetCtx && this.hoverPhase != null) {
        const px = xToPx(this.hoverPhase);
        ctx.save();
        ctx.strokeStyle = tc.crosshair;
        ctx.lineWidth = 1.1;
        ctx.setLineDash([4, 3]);
        ctx.beginPath();
        ctx.moveTo(px, m.top);
        ctx.lineTo(px, m.top + ph);
        ctx.stroke();
        ctx.restore();
      }

      // Shared Bottom X-Axis Ticks & Label
      ctx.fillStyle = tc.muted;
      ctx.font = '11px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      for (const tx of xTicks) {
        ctx.fillText(tx.toFixed(2), xToPx(tx), m.top + ph + 7);
      }

      ctx.fillStyle = tc.text;
      ctx.font = '500 12px system-ui, -apple-system, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'bottom';
      ctx.fillText('Pulse Phase φ (2 cycles shown)', m.left + pw / 2, height - 4);

      ctx.save();
      ctx.translate(15, m.top + ph / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      ctx.fillText(this.opts.mode === 'overlay' ? 'Normalized Intensity (0–1)' : 'Counts / Bin', 0, 0);
      ctx.restore();

      ctx.restore();
    }

    exportSVG({ filename = 'pulse_alignment.svg', theme = null, width = 980 } = {}) {
      if (!this.bands.length) return;
      const th = theme || this.opts.theme;
      const tc = resolveThemeColors(th, this.opts.theme, this.opts.themeColors);
      const m = this._margins();
      const height = this.opts.mode === 'stacked' ? Math.max(460, 96 + this.bands.length * 134) : 520;
      const pw = width - m.left - m.right;
      const ph = height - m.top - m.bottom;
      const xToPx = (phi) => m.left + (phi / 2.0) * pw;
      const xTicks = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0];

      const parts = [];
      parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`);
      parts.push(`<rect width="${width}" height="${height}" fill="${tc.bg}"/>`);
      parts.push(`<text x="${m.left}" y="20" fill="${tc.text}" font-family="sans-serif" font-size="13.5" font-weight="600">${escapeXml(this.opts.title)}</text>`);
      if (this.opts.subtitle) {
        parts.push(`<text x="${m.left}" y="35" fill="${tc.muted}" font-family="monospace" font-size="11">${escapeXml(this.opts.subtitle)}</text>`);
      }

      if (this.opts.mode === 'stacked') {
        const nb = this.bands.length;
        const gap = 8;
        const panelH = (ph - gap * (nb - 1)) / nb;
        for (let bi = 0; bi < nb; bi++) {
          const b = this.bands[bi];
          const pyTop = m.top + bi * (panelH + gap);
          const bins = this._twoCycleBins(b.values, b.errors);
          let yMin = Infinity, yMax = -Infinity;
          for (const bn of bins) {
            const e = this.opts.showErrors ? (bn.e || 0) : 0;
            if (bn.y - e < yMin) yMin = bn.y - e;
            if (bn.y + e > yMax) yMax = bn.y + e;
          }
          const dy = (yMax - yMin) || 1;
          yMin -= dy * 0.08; yMax += dy * 0.16;
          const yToPx = (y) => pyTop + panelH - ((y - yMin) / (yMax - yMin || 1)) * panelH;

          parts.push(`<rect x="${m.left}" y="${pyTop.toFixed(1)}" width="${pw}" height="${panelH.toFixed(1)}" fill="${tc.plotBg}" stroke="${tc.axis}"/>`);
          for (const tx of xTicks) {
            const px = xToPx(tx).toFixed(1);
            parts.push(`<line x1="${px}" y1="${pyTop.toFixed(1)}" x2="${px}" y2="${(pyTop + panelH).toFixed(1)}" stroke="${tc.grid}"/>`);
          }
          const yt = niceTicks(yMin, yMax, 3);
          for (const ty of yt.ticks) {
            const py = yToPx(ty);
            if (py < pyTop + 8 || py > pyTop + panelH - 8) continue;
            parts.push(`<line x1="${m.left}" y1="${py.toFixed(1)}" x2="${m.left + pw}" y2="${py.toFixed(1)}" stroke="${tc.grid}"/>`);
            parts.push(`<text x="${m.left - 8}" y="${(py + 3.5).toFixed(1)}" fill="${tc.muted}" font-family="monospace" font-size="10.5" text-anchor="end">${escapeXml(formatTick(ty, yt.step))}</text>`);
          }

          if (this.opts.showErrors) {
            let errD = '';
            for (const bn of bins) {
              if (!Number.isFinite(bn.e) || bn.e <= 0) continue;
              const px = xToPx(bn.xc).toFixed(1);
              errD += `M${px},${yToPx(bn.y - bn.e).toFixed(1)}V${yToPx(bn.y + bn.e).toFixed(1)}`;
            }
            if (errD) parts.push(`<path d="${errD}" stroke="${b.color}" stroke-opacity="0.5" stroke-width="1" fill="none"/>`);
          }

          let d = '';
          for (let i = 0; i < bins.length; i++) {
            const bn = bins[i];
            const px0 = xToPx(bn.x0).toFixed(1), px1 = xToPx(bn.x1).toFixed(1), py = yToPx(bn.y).toFixed(1);
            d += (i === 0 ? `M${px0},${py}` : `L${px0},${py}`) + `H${px1}`;
          }
          parts.push(`<path d="${d}" stroke="${b.color}" stroke-width="1.7" fill="none"/>`);
          parts.push(`<text x="${m.left + 9}" y="${(pyTop + 15).toFixed(1)}" fill="${b.color}" font-family="sans-serif" font-size="11.5" font-weight="600">${escapeXml(b.label)}</text>`);
          if (Number.isFinite(b.z2)) {
            const sign = (b.lagCycles || 0) >= 0 ? '+' : '';
            const meta = `Z²=${formatSci(b.z2, 4)}  PF=${((b.pf || 0) * 100).toFixed(1)}%  Δφ=${sign}${(b.lagCycles || 0).toFixed(4)}`;
            parts.push(`<text x="${m.left + pw - 9}" y="${(pyTop + 15).toFixed(1)}" fill="${tc.muted}" font-family="monospace" font-size="10.5" text-anchor="end">${escapeXml(meta)}</text>`);
          }
        }
      } else {
        const yMin = -0.1, yMax = 1.15;
        const yToPx = (y) => m.top + ph - ((y - yMin) / (yMax - yMin)) * ph;
        parts.push(`<rect x="${m.left}" y="${m.top}" width="${pw}" height="${ph}" fill="${tc.plotBg}" stroke="${tc.axis}"/>`);
        for (const ty of [0, 0.25, 0.5, 0.75, 1.0]) {
          const py = yToPx(ty);
          parts.push(`<line x1="${m.left}" y1="${py.toFixed(1)}" x2="${m.left + pw}" y2="${py.toFixed(1)}" stroke="${tc.grid}"/>`);
          parts.push(`<text x="${m.left - 8}" y="${(py + 3.5).toFixed(1)}" fill="${tc.muted}" font-family="monospace" font-size="11" text-anchor="end">${ty.toFixed(2)}</text>`);
        }
        for (const b of this.bands) {
          const norm = this._normBand(b);
          const bins = this._twoCycleBins(norm.values, norm.errors);
          let d = '';
          for (let i = 0; i < bins.length; i++) {
            const bn = bins[i];
            const px0 = xToPx(bn.x0).toFixed(1), px1 = xToPx(bn.x1).toFixed(1), py = yToPx(bn.y).toFixed(1);
            d += (i === 0 ? `M${px0},${py}` : `L${px0},${py}`) + `H${px1}`;
          }
          parts.push(`<path d="${d}" stroke="${b.color}" stroke-width="1.8" fill="none"/>`);
        }
      }

      if (this.opts.showPeakGuide && Number.isFinite(this.refPeakPhase)) {
        const p0 = ((this.refPeakPhase % 1) + 1) % 1;
        for (const phi of [p0, p0 + 1]) {
          const px = xToPx(phi).toFixed(1);
          parts.push(`<line x1="${px}" y1="${m.top}" x2="${px}" y2="${m.top + ph}" stroke="${tc.refLine}" stroke-dasharray="5,4" stroke-width="1.25"/>`);
        }
      }

      for (const tx of xTicks) {
        parts.push(`<text x="${xToPx(tx).toFixed(1)}" y="${m.top + ph + 17}" fill="${tc.muted}" font-family="monospace" font-size="11" text-anchor="middle">${tx.toFixed(2)}</text>`);
      }
      parts.push(`<text x="${(m.left + pw / 2).toFixed(1)}" y="${height - 8}" fill="${tc.text}" font-family="sans-serif" font-size="12" font-weight="500" text-anchor="middle">Pulse Phase φ (2 cycles shown)</text>`);
      parts.push(`<text x="16" y="${(m.top + ph / 2).toFixed(1)}" fill="${tc.text}" font-family="sans-serif" font-size="12" font-weight="500" text-anchor="middle" transform="rotate(-90 16 ${(m.top + ph / 2).toFixed(1)})">${this.opts.mode === 'overlay' ? 'Normalized Intensity (0–1)' : 'Counts / Bin'}</text>`);
      parts.push(`</svg>`);

      triggerDownload(new Blob([parts.join('\n')], { type: 'image/svg+xml;charset=utf-8' }), filename);
    }

    exportPNG({ filename = 'pulse_alignment.png', theme = null, scale = 2 } = {}) {
      if (!this.bands.length) return;
      const w = Math.max(900, this.width);
      const h = Math.max(480, this.height);
      const off = document.createElement('canvas');
      off.width = w * scale;
      off.height = h * scale;
      const octx = off.getContext('2d');
      octx.setTransform(scale, 0, 0, scale, 0, 0);
      this.render(octx, w, h, theme || this.opts.theme);
      off.toBlob((blob) => {
        if (blob) triggerDownload(blob, filename);
      }, 'image/png');
    }
  }

  window.AstroCharts = {
    getPalette,
    COLORMAPS,
    stretchValue,
    formatSci,
    niceTicks,
    logTicks,
    canvasToFitsPixel,
    fitsPixelToCanvas,
    fitsPixelToWorld,
    SeriesPlot,
    ImagePlot2D,
    PulseAlignmentPlot,
  };
})();
