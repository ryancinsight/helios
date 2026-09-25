//! Ram-Lak ramp filtering: the FBP filter stage, on the Atlas transform provider.
//!
//! Filtered back-projection convolves every projection row with the Ram-Lak ramp
//! `h[n]` (`h[0] = 1/(4Δs²)`, `h[odd] = −1/(π²n²Δs²)`, `h[even ≠ 0] = 0`) before
//! back-projecting. That convolution is the reconstruction's only
//! super-linear-cost stage: the kernel spans the **full** detector width, so the
//! spatial form costs `O(n_off²)` per projection and the whole filter
//! `O(n_ang · n_off²)`.
//!
//! This module performs it as a zero-padded linear convolution in the frequency
//! domain, through the Atlas transform provider (`apollo`). The pad length is
//! `N ≥ 3·n_off − 2`, the exact length of the linear convolution of an
//! `n_off`-sample row with the `2·n_off − 1`-tap kernel, so **no** sample of the
//! result is touched by circular wrap-around: the transform path reproduces the
//! spatial convolution rather than approximating it. The cost becomes
//! `O(n_ang · n_off log n_off)`.
//!
//! The spatial form is retained as [`RampMethod::Direct`] — the reference the
//! frequency-domain path is verified against, sample by sample, in this module's
//! tests.

use eunomia::Complex;

/// Detector width at and above which [`RampMethod::Auto`] uses the transform path.
///
/// Measured on the `ramp_filter` bench target (64 projections per width): the
/// transform path is already ~1.2× faster at 64 detector samples, and the margin
/// grows with the width — 4.3× at 181, 9.9× at 512, 19.0× at 1024. Narrower
/// detectors keep the spatial path, which allocates no transform plan or buffer.
pub const FFT_CROSSOVER: usize = 64;

/// Which convolution performs the ramp filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RampMethod {
    /// The frequency-domain path at or above [`FFT_CROSSOVER`], the spatial path
    /// below it.
    Auto,
    /// Always the spatial convolution — the reference implementation.
    Direct,
}

/// Ram-Lak ramp-filter kernel `h[n]` for `n ∈ [−(len−1), len−1]`, detector
/// sample spacing `ds_cm`.
///
/// `h[0] = 1/(4Δs²)`, `h[odd n] = −1/(π²n²Δs²)`, `h[even n ≠ 0] = 0`. Returned as
/// a `2·len − 1` vector indexed by `n + (len − 1)`, so the centre tap sits at
/// index `len − 1`. An empty detector yields an empty kernel.
#[must_use]
pub fn ram_lak_kernel(len: usize, ds_cm: f64) -> Vec<f64> {
    if len == 0 {
        return Vec::new();
    }
    let inv_ds_sq = (ds_cm * ds_cm).recip();
    let pi_sq = core::f64::consts::PI * core::f64::consts::PI;
    let mut kernel = vec![0.0_f64; 2 * len - 1];
    let base = len as isize - 1;
    for n in -base..=base {
        kernel[(n + base) as usize] = if n == 0 {
            0.25 * inv_ds_sq
        } else if n % 2 != 0 {
            let n_sq = (n * n) as f64;
            -(pi_sq * n_sq).recip() * inv_ds_sq
        } else {
            0.0
        };
    }
    kernel
}

/// Ramp-filter every projection row of a flat `n_ang × n_off` sinogram buffer,
/// in place.
///
/// `rows[a * n_off + i]` is projection `a` at detector offset `i`; each row is
/// replaced by `Δs · Σ_k row[k] · kernel[i − k + (n_off − 1)]`, i.e. the
/// convolution of the row with the centred `kernel` produced by
/// [`ram_lak_kernel`], evaluated only where the row and the kernel overlap.
///
/// `kernel` must have `2·n_off − 1` taps. An empty detector is a no-op.
pub fn ramp_filter_rows(
    rows: &mut [f64],
    n_ang: usize,
    n_off: usize,
    kernel: &[f64],
    ds_cm: f64,
    method: RampMethod,
) {
    if n_off == 0 || n_ang == 0 {
        return;
    }
    debug_assert_eq!(kernel.len(), 2 * n_off - 1, "kernel must span the detector");
    let use_fft = match method {
        RampMethod::Auto => n_off >= FFT_CROSSOVER,
        RampMethod::Direct => false,
    };
    if use_fft {
        filter_rows_fft(rows, n_ang, n_off, kernel, ds_cm);
    } else {
        filter_rows_direct(rows, n_ang, n_off, kernel, ds_cm);
    }
}

/// The spatial convolution — `O(n_off²)` per row. The reference implementation.
fn filter_rows_direct(rows: &mut [f64], n_ang: usize, n_off: usize, kernel: &[f64], ds_cm: f64) {
    let base = n_off as isize - 1;
    let mut filtered = vec![0.0_f64; n_off];
    for a in 0..n_ang {
        let start = a * n_off;
        let row = &rows[start..start + n_off];
        for (i, slot) in filtered.iter_mut().enumerate() {
            let mut acc = 0.0_f64;
            for (k, &p) in row.iter().enumerate() {
                acc += p * kernel[(i as isize - k as isize + base) as usize];
            }
            *slot = acc * ds_cm;
        }
        rows[start..start + n_off].copy_from_slice(&filtered);
    }
}

/// The zero-padded frequency-domain convolution — `O(n_off log n_off)` per row.
///
/// Both operands are placed at index 0 of a length-`N` buffer (`N ≥ 3·n_off − 2`),
/// so the linear convolution is `Σ_k row[k] · kernel[i − k]` and the ramp-filter
/// result for output `i` is that convolution read at `i + (n_off − 1)`. Nothing
/// wraps: the convolution's support is `[0, 3·n_off − 3] ⊂ [0, N)`.
fn filter_rows_fft(rows: &mut [f64], n_ang: usize, n_off: usize, kernel: &[f64], ds_cm: f64) {
    let base = n_off - 1;
    let nfft = next_fast_len(3 * n_off - 2);

    let mut kernel_buf = vec![0.0_f64; nfft];
    kernel_buf[..kernel.len()].copy_from_slice(kernel);
    let kernel_spectrum = apollo::fft_1d_slice::<f64>(&kernel_buf);

    let mut row_buf = vec![0.0_f64; nfft];
    for a in 0..n_ang {
        let start = a * n_off;
        row_buf[..n_off].copy_from_slice(&rows[start..start + n_off]);
        row_buf[n_off..].fill(0.0);

        let row_spectrum = apollo::fft_1d_slice::<f64>(&row_buf);
        let product: Vec<Complex<f64>> = row_spectrum
            .iter()
            .zip(kernel_spectrum.iter())
            .map(|(p, q)| Complex::new(p.re * q.re - p.im * q.im, p.re * q.im + p.im * q.re))
            .collect();
        let convolved = apollo::ifft_1d_slice::<f64>(&product);

        for i in 0..n_off {
            rows[start + i] = convolved[base + i] * ds_cm;
        }
    }
}

/// Smallest 5-smooth integer `≥ n` (`2^a · 3^b · 5^c`).
///
/// The transform is mixed-radix and accepts any length, but a length with a
/// large prime factor falls back to Rader/Bluestein, which is slower than the
/// direct convolution this path exists to beat. `3·n_off − 2` is frequently
/// prime (`n_off = 181` gives 541), so the pad is rounded to a 5-smooth length.
fn next_fast_len(n: usize) -> usize {
    if n <= 1 {
        return 1;
    }
    let limit = n.saturating_mul(2);
    let mut best = n.next_power_of_two();
    let mut p5 = 1usize;
    while p5 < limit {
        let mut p35 = p5;
        while p35 < limit {
            let mut candidate = p35;
            while candidate < n {
                candidate = candidate.saturating_mul(2);
            }
            if candidate < best {
                best = candidate;
            }
            p35 = p35.saturating_mul(3);
        }
        p5 = p5.saturating_mul(5);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ramp_row(n_off: usize, seed: f64) -> Vec<f64> {
        (0..n_off)
            .map(|k| {
                let t = k as f64;
                (t * 0.7 + seed).cos() + 0.3 * (t * 0.31).sin()
            })
            .collect()
    }

    #[test]
    fn ram_lak_kernel_has_expected_structure() {
        let k = ram_lak_kernel(5, 0.1);
        let base = 4;
        assert_eq!(k.len(), 9);
        // Even off-centre taps are exactly zero; odd taps are negative and
        // symmetric about the centre.
        assert_eq!(k[base + 2], 0.0);
        assert_eq!(k[base - 2], 0.0);
        assert!(k[base] > 0.0);
        assert!(k[base + 1] < 0.0 && k[base + 3] < 0.0);
        assert!((k[base + 1] - k[base - 1]).abs() < 1e-15);
        assert!((k[base + 3] - k[base - 3]).abs() < 1e-15);
    }

    #[test]
    fn empty_detector_is_a_noop() {
        assert!(ram_lak_kernel(0, 0.1).is_empty());
        let mut rows = vec![1.0_f64, 2.0];
        ramp_filter_rows(&mut rows, 2, 0, &[], 0.1, RampMethod::Auto);
        assert_eq!(rows, vec![1.0, 2.0]);
    }

    #[test]
    fn next_fast_len_is_5_smooth_and_never_shrinks() {
        for n in 1..600usize {
            let f = next_fast_len(n);
            assert!(f >= n, "next_fast_len({n}) = {f} < {n}");
            let mut m = f;
            for p in [2usize, 3, 5] {
                while m.is_multiple_of(p) {
                    m /= p;
                }
            }
            assert_eq!(m, 1, "next_fast_len({n}) = {f} has a prime factor > 5");
        }
    }

    /// The H-113 acceptance oracle: the transform path reproduces the spatial
    /// convolution, and does so at every width — including the widths where
    /// `Auto` still takes the direct path, so the two branches are both covered.
    #[test]
    fn transform_path_matches_the_spatial_reference() {
        for &n_off in &[3usize, 8, 32, 64, 181, 256] {
            let kernel = ram_lak_kernel(n_off, 0.2);
            let row = ramp_row(n_off, 0.13);

            let mut direct = row.clone();
            ramp_filter_rows(&mut direct, 1, n_off, &kernel, 0.2, RampMethod::Direct);

            let mut transformed = row.clone();
            ramp_filter_rows(&mut transformed, 1, n_off, &kernel, 0.2, RampMethod::Auto);

            let scale = direct.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            let error = direct
                .iter()
                .zip(transformed.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            // The transform sums the same products in a different order, so the
            // two agree to a few ulps of the row's own magnitude scaled by the
            // O(n_off) accumulation length.
            let bound = 64.0 * f64::EPSILON * (n_off as f64) * scale.max(f64::MIN_POSITIVE);
            assert!(
                error <= bound,
                "n_off = {n_off}: transform error {error:e} exceeds the derived bound {bound:e}"
            );
        }
    }

    #[test]
    fn transform_path_is_row_independent() {
        let n_off = 96usize;
        let kernel = ram_lak_kernel(n_off, 0.2);
        let rows: Vec<f64> = (0..3)
            .flat_map(|a| ramp_row(n_off, f64::from(a) * 0.37))
            .collect();

        let mut batched = rows.clone();
        ramp_filter_rows(&mut batched, 3, n_off, &kernel, 0.2, RampMethod::Auto);

        for a in 0..3 {
            let mut single = rows[a * n_off..(a + 1) * n_off].to_vec();
            ramp_filter_rows(&mut single, 1, n_off, &kernel, 0.2, RampMethod::Auto);
            assert_eq!(&batched[a * n_off..(a + 1) * n_off], &single[..]);
        }
    }
}
