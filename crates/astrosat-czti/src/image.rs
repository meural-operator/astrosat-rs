//! Stage 3b: `cztimage` — Coded-mask sky image reconstruction via Balanced Cross-Correlation
//! (`balcc`) and 2D FFT Cross-Correlation (`fft`), with CZTI ICRS & local camera WCS parameters.

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use serde::{Deserialize, Serialize};

use crate::types::{QuadAspect, MASK_HEIGHT_MM, PIXELS_PER_QUAD, QUAD_DIM, XLDET_MM};

/// Sky image reconstruction algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ImageAlgorithm {
    /// Balanced Cross-Correlation (`balcc`, CZTI v3.0 default).
    BalCC,
    /// 2D FFT cross-correlation (`fft`) with oversampled mask pattern.
    Fft,
}

impl ImageAlgorithm {
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "balcc" => Some(ImageAlgorithm::BalCC),
            "fft" => Some(ImageAlgorithm::Fft),
            _ => None,
        }
    }
}

/// Configuration for `cztimage`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImageConfig {
    pub algorithm: ImageAlgorithm,
    /// Mask/DPI oversampling factor for FFT mode (power of 2, default: `16`, or `4` for fast preview).
    pub oversampling_factor: usize,
    /// Quadrants to include (default: `[0, 1, 2, 3]`).
    pub quads: Vec<usize>,
}

impl Default for ImageConfig {
    fn default() -> Self {
        Self {
            algorithm: ImageAlgorithm::BalCC,
            oversampling_factor: 16,
            quads: vec![0, 1, 2, 3],
        }
    }
}

/// CZTI WCS coordinates (both ICRS `RA---TAN`/`DEC--TAN` and local camera `TANX`/`TANY`),
/// matching `Image::evaluate_wcs_coordinates` in `level2validation.cpp`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CztiWcs {
    pub ctype1: String,
    pub ctype2: String,
    pub crpix1: f64,
    pub crpix2: f64,
    pub crval1: f64,
    pub crval2: f64,
    pub cdelt1: f64,
    pub cdelt2: f64,
    pub crota1: f64,
    pub crota2: f64,
    pub cd1_1: f64,
    pub cd1_2: f64,
    pub cd2_1: f64,
    pub cd2_2: f64,
    // Alternate local camera WCS ("B")
    pub crpix1b: f64,
    pub crpix2b: f64,
    pub crval1b: f64,
    pub crval2b: f64,
    pub cdelt1b: f64,
    pub cdelt2b: f64,
    pub cd1_1b: f64,
    pub cd1_2b: f64,
    pub cd2_1b: f64,
    pub cd2_2b: f64,
}

impl CztiWcs {
    /// Compute CZTI WCS parameters from quadrant aspect (`RA`, `DEC`, `TWIST` in degrees),
    /// image grid dimensions `(nx, ny)`, and explicit angular pixel scale in degrees/pixel.
    pub fn from_aspect_with_scale(
        aspect: &QuadAspect,
        nx: usize,
        ny: usize,
        scale_deg: f64,
    ) -> Self {
        let scale = scale_deg;
        let crota = -aspect.twist_deg;
        let crota_rad = crota.to_radians();
        let cd1_1 = scale * crota_rad.cos();
        let cd1_2 = scale * crota_rad.sin();
        let cd2_1 = -cd1_2;
        let cd2_2 = cd1_1;

        let crpix1 = ((nx / 2) + 1) as f64;
        let crpix2 = ((ny / 2) + 1) as f64;

        Self {
            ctype1: "RA---TAN".to_string(),
            ctype2: "DEC--TAN".to_string(),
            crpix1,
            crpix2,
            crval1: aspect.ra_deg,
            crval2: aspect.dec_deg,
            cdelt1: scale,
            cdelt2: scale,
            crota1: crota,
            crota2: crota,
            cd1_1,
            cd1_2,
            cd2_1,
            cd2_2,
            crpix1b: crpix1,
            crpix2b: crpix2,
            crval1b: 0.0,
            crval2b: 0.0,
            cdelt1b: scale,
            cdelt2b: scale,
            cd1_1b: scale,
            cd1_2b: 0.0,
            cd2_1b: 0.0,
            cd2_2b: scale,
        }
    }

    /// Compute CZTI WCS parameters from quadrant aspect (`RA`, `DEC`, `TWIST` in degrees),
    /// image grid dimensions `(nx, ny)`, and effective oversampling factor.
    pub fn from_aspect(aspect: &QuadAspect, nx: usize, ny: usize, oversampling: usize) -> Self {
        let ov = oversampling.max(1) as f64;
        let scale = (XLDET_MM / (MASK_HEIGHT_MM * 16.0 * ov)).atan().to_degrees();
        Self::from_aspect_with_scale(aspect, nx, ny, scale)
    }
}

/// Reconstructed sky image plus WCS metadata.
#[derive(Debug, Clone)]
pub struct SkyImage {
    pub width: usize,
    pub height: usize,
    /// Row-major `data[y * width + x]`.
    pub data: Vec<f32>,
    pub wcs: CztiWcs,
}

impl SkyImage {
    /// Locate the peak pixel `(x, y, value)` in the reconstructed sky image.
    pub fn peak(&self) -> Option<(usize, usize, f32)> {
        let mut best: Option<(usize, usize, f32)> = None;
        for y in 0..self.height {
            for x in 0..self.width {
                let v = self.data[y * self.width + x];
                if v.is_finite() && best.is_none_or(|b| v > b.2) {
                    best = Some((x, y, v));
                }
            }
        }
        best
    }
}

/// Compute Balanced Cross-Correlation (`balcc`) sky image from a 64x64 quadrant DPI and a shadow
/// library (`shadows`: `grid_dim * grid_dim` shadow templates, each of length 4096),
/// parallelized over sky pixels with `rayon`.
///
/// Exact formula from `level2validation.cpp::generate_balanced_cross_correlation_image`:
/// `I[j, k] = (sum_{dpi!=0} s_p * dpi_p / sum_{dpi!=0} s_p) - (sum_{dpi!=0} (1-s_p) * dpi_p / sum_{dpi!=0} (1-s_p))`
pub fn balanced_cross_correlation(
    dpi_64x64: &[Vec<f32>],
    shadows: &[&[f32]],
    grid_dim: usize,
    aspect: &QuadAspect,
) -> SkyImage {
    let mut dpi1d = [0.0f32; PIXELS_PER_QUAD];
    for (y, row) in dpi_64x64.iter().take(QUAD_DIM).enumerate() {
        for (x, &v) in row.iter().take(QUAD_DIM).enumerate() {
            dpi1d[y * QUAD_DIM + x] = v;
        }
    }

    // CALDB `SHADOWLIB` rows are ordered with `ThetaX` in `[-3.0, +3.0]` deg as outer index
    // (`tx_i = row / grid_dim`) and `ThetaY` in `[-3.0, +3.0]` deg as inner index
    // (`ty_i = row % grid_dim`). Store row-major `data[ty_i * grid_dim + tx_i]` so FITS X=ThetaX
    // and FITS Y=ThetaY.
    let n_sky = grid_dim * grid_dim;
    let data: Vec<f32> = (0..n_sky)
        .into_par_iter()
        .map(|out_idx| {
            let ty_i = out_idx / grid_dim;
            let tx_i = out_idx % grid_dim;
            let shadow_idx = tx_i * grid_dim + ty_i;
            let Some(&shadow) = shadows.get(shadow_idx) else {
                return 0.0;
            };
            let mut cc_src = 0.0f64;
            let mut cc_back = 0.0f64;
            let mut norm = 0.0f64;
            let mut back_norm = 0.0f64;
            let n = shadow.len().min(PIXELS_PER_QUAD);
            for p in 0..n {
                let d = dpi1d[p] as f64;
                if d != 0.0 {
                    let s = shadow[p] as f64;
                    let b = 1.0 - s;
                    norm += s;
                    back_norm += b;
                    cc_src += s * d;
                    cc_back += b * d;
                }
            }
            if norm > 0.0 && back_norm > 0.0 {
                (cc_src / norm - cc_back / back_norm) as f32
            } else {
                0.0
            }
        })
        .collect();

    // Standard CZTI SHADOWLIB spans [-3.0, +3.0] deg across `grid_dim` samples (0.10 deg/pix for 61x61)
    let scale_deg = if grid_dim > 1 {
        6.0 / ((grid_dim - 1) as f64)
    } else {
        0.1
    };
    let wcs = CztiWcs::from_aspect_with_scale(aspect, grid_dim, grid_dim, scale_deg);
    SkyImage {
        width: grid_dim,
        height: grid_dim,
        data,
        wcs,
    }
}

/// Compute balanced cross-correlation directly from a 64x64 binary mask pattern by shifting the
/// mask across `(-half..=half)` pixel lags (useful when `SHADOWLIB` is unavailable or in tests).
pub fn balanced_cross_correlation_from_mask(
    dpi_64x64: &[Vec<f32>],
    mask_64x64: &[Vec<f32>],
    half_span: usize,
    aspect: &QuadAspect,
) -> SkyImage {
    let dim = 2 * half_span + 1;
    let data: Vec<f32> = (0..dim * dim)
        .into_par_iter()
        .map(|idx| {
            let gy = (idx / dim) as isize - (half_span as isize);
            let gx = (idx % dim) as isize - (half_span as isize);
            let mut cc_src = 0.0f64;
            let mut cc_back = 0.0f64;
            let mut norm = 0.0f64;
            let mut back_norm = 0.0f64;

            for (y, dpi_row) in dpi_64x64.iter().enumerate().take(QUAD_DIM) {
                let my = y as isize + gy;
                if !(0..QUAD_DIM as isize).contains(&my) {
                    continue;
                }
                for (x, &dv) in dpi_row.iter().enumerate().take(QUAD_DIM) {
                    let mx = x as isize + gx;
                    if !(0..QUAD_DIM as isize).contains(&mx) {
                        continue;
                    }
                    let d = dv as f64;
                    if d != 0.0 {
                        let Some(&ms) = mask_64x64
                            .get(my as usize)
                            .and_then(|r| r.get(mx as usize))
                        else {
                            continue;
                        };
                        let s = ms as f64;
                        let b = 1.0 - s;
                        norm += s;
                        back_norm += b;
                        cc_src += s * d;
                        cc_back += b * d;
                    }
                }
            }
            if norm > 0.0 && back_norm > 0.0 {
                (cc_src / norm - cc_back / back_norm) as f32
            } else {
                0.0
            }
        })
        .collect();

    let wcs = CztiWcs::from_aspect(aspect, dim, dim, 1);
    SkyImage {
        width: dim,
        height: dim,
        data,
        wcs,
    }
}

/// 2D FFT forward or inverse transform on a row-major `n x n` complex buffer using `rustfft`.
fn fft_2d(buf: &mut [Complex<f64>], n: usize, inverse: bool) {
    let mut planner = FftPlanner::<f64>::new();
    let fft = if inverse {
        planner.plan_fft_inverse(n)
    } else {
        planner.plan_fft_forward(n)
    };

    // Row-wise FFTs (parallelized with rayon)
    buf.par_chunks_mut(n).for_each(|row| {
        fft.process(row);
    });

    // Transpose
    let mut transposed = vec![Complex::new(0.0, 0.0); n * n];
    for y in 0..n {
        for x in 0..n {
            transposed[x * n + y] = buf[y * n + x];
        }
    }

    // Column-wise FFTs (now contiguous rows in `transposed`)
    transposed.par_chunks_mut(n).for_each(|col| {
        fft.process(col);
    });

    // Transpose back (and normalize if inverse)
    let scale = if inverse { 1.0 / ((n * n) as f64) } else { 1.0 };
    for y in 0..n {
        for x in 0..n {
            buf[y * n + x] = transposed[x * n + y] * scale;
        }
    }
}

/// Compute 2D FFT cross-correlation sky image (`cztimage` `algorithm=fft`) using `rustfft`,
/// applying explicit Fourier phase shifts `(xshift, yshift)` and taking complex magnitude
/// (`common/fft.h::generate_correlation_matrix2D` & `cztimage_v2.cpp:482-504`).
pub fn fft_cross_correlation_with_shift(
    dpi_64x64: &[Vec<f32>],
    mask_64x64: &[Vec<f32>],
    oversampling: usize,
    xshift: f64,
    yshift: f64,
    aspect: &QuadAspect,
) -> SkyImage {
    let ov = oversampling.clamp(1, 64);
    let n = (QUAD_DIM * ov).next_power_of_two();
    let mut a = vec![Complex::new(0.0, 0.0); n * n];
    let mut b = vec![Complex::new(0.0, 0.0); n * n];

    // Oversample DPI and balanced mask (2*M - 1) onto the n x n grid
    for y in 0..QUAD_DIM {
        for x in 0..QUAD_DIM {
            let dv = dpi_64x64
                .get(y)
                .and_then(|r| r.get(x))
                .copied()
                .unwrap_or(0.0) as f64;
            let mv = (2.0
                * (mask_64x64
                    .get(y)
                    .and_then(|r| r.get(x))
                    .copied()
                    .unwrap_or(0.0) as f64))
                - 1.0;
            for oy in 0..ov {
                let yy = y * ov + oy;
                if yy >= n {
                    continue;
                }
                for ox in 0..ov {
                    let xx = x * ov + ox;
                    if xx >= n {
                        continue;
                    }
                    a[yy * n + xx] = Complex::new(dv, 0.0);
                    b[yy * n + xx] = Complex::new(mv, 0.0);
                }
            }
        }
    }

    fft_2d(&mut a, n, false);
    fft_2d(&mut b, n, false);

    // Multiply A * conj(B) * phase_ramp(xshift, yshift) matching `fft.h:120-200`
    let shift_factor = std::f64::consts::TAU * (ov as f64);
    let xs = xshift * shift_factor;
    let ys = yshift * shift_factor;
    let half_n = n / 2;
    let n_f = n as f64;

    a.par_chunks_mut(n)
        .zip(b.par_chunks(n))
        .enumerate()
        .for_each(|(i, (row_a, row_b))| {
            let ii = if i < half_n {
                i as isize
            } else if i > half_n {
                (i as isize) - (n as isize)
            } else {
                0
            };
            // Row index `i` is Y frequency; column index `j` is X frequency.
            // `rustfft` forward uses exp(-i 2pi k n / N) (conjugate of NR `four1(..., +1)`),
            // so the shift phase ramp carries a minus sign `exp(-i * (fx + fy))`.
            let fy = (ii as f64 / n_f) * ys;
            for j in 0..n {
                let jj = if j < half_n {
                    j as isize
                } else if j > half_n {
                    (j as isize) - (n as isize)
                } else {
                    0
                };
                let fx = (jj as f64 / n_f) * xs;
                let phi = -(fx + fy);
                let phase = Complex::new(phi.cos(), phi.sin());
                row_a[j] = (row_a[j] * row_b[j].conj()) * phase;
            }
        });

    fft_2d(&mut a, n, true);

    // Center zero-lag (fftshift), take magnitude (`fft.h:213-214`), and extract central sky window
    // matching `fft.h:288-306` (`swap_correlation_matrix`: `sy = (jmid - cy) mod N`, `sx = (imid - cx) mod N`).
    let crop = (30 * ov + 1).min(n);
    let half_crop = crop / 2;
    let mut data = vec![0.0f32; crop * crop];

    for cy in 0..crop {
        let dy = (half_crop as isize) - (cy as isize);
        let sy = dy.rem_euclid(n as isize) as usize;
        for cx in 0..crop {
            let dx = (half_crop as isize) - (cx as isize);
            let sx = dx.rem_euclid(n as isize) as usize;
            data[cy * crop + cx] = a[sy * n + sx].norm() as f32;
        }
    }

    let wcs = CztiWcs::from_aspect(aspect, crop, crop, ov);
    SkyImage {
        width: crop,
        height: crop,
        data,
        wcs,
    }
}

/// Compute 2D FFT cross-correlation sky image (`cztimage` `algorithm=fft`) using `rustfft`.
///
/// Automatically applies the legacy `yshift = 1.0` Fourier phase ramp for `Q2`/`Q3`
/// (`cztimage_v2.cpp:482-504`).
pub fn fft_cross_correlation(
    dpi_64x64: &[Vec<f32>],
    mask_64x64: &[Vec<f32>],
    oversampling: usize,
    aspect: &QuadAspect,
) -> SkyImage {
    let yshift = if aspect.quad_id == 2 || aspect.quad_id == 3 {
        1.0
    } else {
        0.0
    };
    fft_cross_correlation_with_shift(dpi_64x64, mask_64x64, oversampling, 0.0, yshift, aspect)
}

/// Combine multiple quadrant sky images into a single summed all-quadrant sky image
/// (`cztimage_v2.cpp:532-550` `allQuadImage` for `GtiType::Common`).
pub fn combine_sky_images(images: &[SkyImage], aspect: &QuadAspect) -> SkyImage {
    let Some(first) = images.first() else {
        return SkyImage {
            width: 0,
            height: 0,
            data: Vec::new(),
            wcs: CztiWcs::from_aspect(aspect, 0, 0, 1),
        };
    };
    let mut data = vec![0.0f32; first.width * first.height];
    for img in images {
        if img.width == first.width && img.height == first.height {
            for (dst, &src) in data.iter_mut().zip(img.data.iter()) {
                *dst += src;
            }
        }
    }
    SkyImage {
        width: first.width,
        height: first.height,
        data,
        wcs: first.wcs.clone(),
    }
}
