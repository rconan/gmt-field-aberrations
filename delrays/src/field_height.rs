//! # Field-height vectors of misaligned segments
//!
//! Estimation of the field-height (field-shift) vectors `ς_j` of the two mirrors (j = M1, M2)
//! of a misaligned GMT segment from its pupil astigmatism (Noll 5, 6) and pupil coma
//! (Noll 7, 8) measured at a few field points, typically 3
//! (section "Field-height vectors of misaligned segments" of the double-Zernike paper).
//!
//! The wavefront of segment `k` at field point `ζ` is
//! `W = Σ_j W_j(ζ - ς_j, ρ + s_k)` with `W_j` the polynomial wave aberration of mirror `j`.
//! The pupil coefficients are linear in `ς_j` to first order; [estimate] solves the
//! linearized system by Gauss-Newton, with, for the outer segments, an extra unknown
//! `δω_040`, the field-independent spherical aberration of a segment despace.
//! [closed_form] is the third-order solution (Eq. "closed form" of the paper), useful
//! as a check, but biased by 10 to 80% for the outer segments of the GMT.
//!
//! Conventions:
//!  * field points and `ς_j` are in units of the field radius (10 arcmin for the GMT data),
//!    `ζ = (ζ_x, ζ_y)` with the azimuth of `ζ` the azimuth of the source;
//!  * pupil coefficients are Noll-normalized `[b5, b6, b7, b8]` over the segment
//!    (radius 4.1825 m), in the units of the model coefficients (nm for [PointModel::gmt]);
//!  * segments are numbered 1 to 7, 7 being the centre segment.
//!
//! ```ignore
//! use delrays::{field_height::{estimate, PointModel}, inversion::Layout};
//! let model = PointModel::gmt();
//! let fit = estimate(&model, &Layout::entrance(), 1, &field, &measured, &nominal);
//! println!("{fit}");
//! ```

use std::{f64::consts::PI, fmt, str::FromStr};

use num_complex::Complex64 as C64;

use crate::inversion::Layout;

/// Field radius of the GMT data \[arcmin\]
pub const FIELD_RADIUS_ARCMIN: f64 = 10.;

/// Number of radial and azimuthal nodes of the pupil quadrature
const N_RADIAL: usize = 16;
const N_AZIMUTHAL: usize = 48;

/// Gauss-Legendre nodes and weights on [-1, 1]
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0f64; n];
    let mut w = vec![0f64; n];
    for i in 0..n {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        loop {
            let (mut p0, mut p1) = (1f64, 0f64);
            for k in 0..n {
                let p2 = p1;
                p1 = p0;
                p0 = ((2 * k + 1) as f64 * z * p1 - k as f64 * p2) / (k + 1) as f64;
            }
            let dp = n as f64 * (z * p0 - p1) / (z * z - 1.);
            let dz = p0 / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                x[i] = z;
                w[i] = 2. / ((1. - z * z) * dp * dp);
                break;
            }
        }
    }
    (x, w)
}

/// Quadrature of the unit disc with the Noll 5 to 8 projectors folded into the weights
#[derive(Debug, Clone)]
struct Quadrature {
    x: Vec<f64>,
    y: Vec<f64>,
    /// `b_J = Σ_i proj[J-5][i] W(x_i, y_i)`
    proj: [Vec<f64>; 4],
}

impl Quadrature {
    fn new() -> Self {
        let (xr, wr) = gauss_legendre(N_RADIAL);
        let (mut x, mut y, mut qw) = (vec![], vec![], vec![]);
        for (xr, wr) in xr.iter().zip(&wr) {
            let r = 0.5 * (xr + 1.);
            for k in 0..N_AZIMUTHAL {
                let t = 2. * PI * k as f64 / N_AZIMUTHAL as f64;
                x.push(r * t.cos());
                y.push(r * t.sin());
                qw.push(0.5 * wr * r * 2. / N_AZIMUTHAL as f64);
            }
        }
        let proj = [
            |x: f64, y: f64| 6f64.sqrt() * 2. * x * y,
            |x: f64, y: f64| 6f64.sqrt() * (x * x - y * y),
            |x: f64, y: f64| 8f64.sqrt() * (3. * (x * x + y * y) - 2.) * y,
            |x: f64, y: f64| 8f64.sqrt() * (3. * (x * x + y * y) - 2.) * x,
        ]
        .map(|z| {
            x.iter()
                .zip(&y)
                .zip(&qw)
                .map(|((&x, &y), &w)| z(x, y) * w)
                .collect()
        });
        Self { x, y, proj }
    }
    fn project(&self, w: &[f64]) -> [f64; 4] {
        self.proj
            .each_ref()
            .map(|p| p.iter().zip(w).map(|(p, w)| p * w).sum())
    }
}

/// Wave-aberration model of the two mirrors:
/// `W_j = Σ_t ω_jt (ζ_j·ζ_j)^p (ρ_s·ρ_s)^n (ζ_j·ρ_s)^m`, j = M1, M2
#[derive(Debug, Clone)]
pub struct PointModel {
    terms: Vec<(u32, u32, u32)>,
    w: [Vec<f64>; 2],
    quad: Quadrature,
}

impl PointModel {
    /// New model from the terms `(p, n, m)` and the coefficients of M1 and M2
    pub fn new(terms: Vec<(u32, u32, u32)>, w_m1: Vec<f64>, w_m2: Vec<f64>) -> Self {
        assert_eq!(terms.len(), w_m1.len());
        assert_eq!(terms.len(), w_m2.len());
        Self {
            terms,
            w: [w_m1, w_m2],
            quad: Quadrature::new(),
        }
    }
    /// GMT model: 136 terms (`k ≤ 9`, `l ≤ 12`, `m ≤ 5`) per mirror, fitted to the ray-traced
    /// wavefront of M1 alone (prime focus) and of the telescope, entrance-pupil coordinates, nm
    pub fn gmt() -> Self {
        include_str!("../data/gmt_per_surface.txt")
            .parse()
            .expect("invalid built-in GMT model")
    }
    /// Keeps only the given terms
    pub fn truncated(&self, keep: &[(u32, u32, u32)]) -> Self {
        let idx: Vec<usize> = keep
            .iter()
            .filter_map(|t| self.terms.iter().position(|s| s == t))
            .collect();
        Self::new(
            idx.iter().map(|&i| self.terms[i]).collect(),
            idx.iter().map(|&i| self.w[0][i]).collect(),
            idx.iter().map(|&i| self.w[1][i]).collect(),
        )
    }
    /// Coefficient `ω_jklm` of mirror `j` (0: M1, 1: M2) for the term `(p, n, m)`
    pub fn coefficient(&self, mirror: usize, term: (u32, u32, u32)) -> Option<f64> {
        self.terms
            .iter()
            .position(|&t| t == term)
            .map(|i| self.w[mirror][i])
    }
    /// Noll 5 to 8 of the segment at `offset` (units of the segment radius), at `field`,
    /// with the field-height vectors `alpha` of M1 and M2
    pub fn noll(&self, offset: [f64; 2], field: [f64; 2], alpha: [[f64; 2]; 2]) -> [f64; 4] {
        let q = &self.quad;
        let mut w = vec![0f64; q.x.len()];
        for (j, a) in alpha.iter().enumerate() {
            let (ux, uy) = (field[0] - a[0], field[1] - a[1]);
            let uu = ux * ux + uy * uy;
            for (i, (x, y)) in q.x.iter().zip(&q.y).enumerate() {
                let (px, py) = (x + offset[0], y + offset[1]);
                let pp = px * px + py * py;
                let up = ux * px + uy * py;
                w[i] += self
                    .terms
                    .iter()
                    .zip(&self.w[j])
                    .map(|(&(p, n, m), c)| {
                        c * uu.powi(p as i32) * pp.powi(n as i32) * up.powi(m as i32)
                    })
                    .sum::<f64>();
            }
        }
        q.project(&w)
    }
    /// Noll 5 to 8 of a field-independent spherical aberration `(ρ_s·ρ_s)^2` of unit amplitude
    fn spherical(&self, offset: [f64; 2]) -> [f64; 4] {
        let q = &self.quad;
        let w: Vec<f64> =
            q.x.iter()
                .zip(&q.y)
                .map(|(x, y)| ((x + offset[0]).powi(2) + (y + offset[1]).powi(2)).powi(2))
                .collect();
        q.project(&w)
    }
}

/// Error parsing a [PointModel] file
#[derive(Debug, Clone, PartialEq)]
pub struct ModelParseError(String);
impl fmt::Display for ModelParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "point model parse error: {}", self.0)
    }
}
impl std::error::Error for ModelParseError {}

/// Parses lines `p n m ω_M1 ω_M2`; `#` starts a comment
impl FromStr for PointModel {
    type Err = ModelParseError;
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let (mut terms, mut w1, mut w2) = (vec![], vec![], vec![]);
        for (k, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap().trim();
            if line.is_empty() {
                continue;
            }
            let err = |e: String| ModelParseError(format!("line {}: {e}", k + 1));
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() != 5 {
                return Err(err(format!("expected 5 columns, found {}", f.len())));
            }
            let int = |s: &str| s.parse::<u32>().map_err(|e| err(e.to_string()));
            let real = |s: &str| s.parse::<f64>().map_err(|e| err(e.to_string()));
            terms.push((int(f[0])?, int(f[1])?, int(f[2])?));
            w1.push(real(f[3])?);
            w2.push(real(f[4])?);
        }
        if terms.is_empty() {
            return Err(ModelParseError("no terms".into()));
        }
        Ok(Self::new(terms, w1, w2))
    }
}

/// Segment offset `s` in units of the segment radius (segments 1 to 7)
pub fn segment_offset(layout: &Layout, segment: usize) -> [f64; 2] {
    assert!((1..=7).contains(&segment), "segment must be 1 to 7");
    if segment == 7 {
        [0., 0.]
    } else {
        let phi = layout.azimuth[segment - 1];
        [layout.s * phi.cos(), layout.s * phi.sin()]
    }
}

/// Field-height vectors of a segment
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Estimate {
    /// `[ς_M1, ς_M2]`, each `(x, y)` in units of the field radius
    pub alpha: [[f64; 2]; 2],
    /// Change of the field-independent spherical aberration (outer segments, despace), nm
    pub dw040: f64,
    /// Rms of the fit residuals, nm
    pub rms: f64,
    /// Number of Gauss-Newton iterations
    pub iterations: usize,
}

impl Estimate {
    /// `ς_M1 - ς_M2`, the best determined combination (driven by the field-constant coma)
    pub fn differential(&self) -> [f64; 2] {
        [
            self.alpha[0][0] - self.alpha[1][0],
            self.alpha[0][1] - self.alpha[1][1],
        ]
    }
}

impl fmt::Display for Estimate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mas = FIELD_RADIUS_ARCMIN * 60e3;
        for (name, a) in ["M1", "M2"].iter().zip(&self.alpha) {
            writeln!(
                f,
                "ς_{name} = ({:+10.3}, {:+10.3}) mas",
                a[0] * mas,
                a[1] * mas
            )?;
        }
        writeln!(f, "δω040 = {:+.4} nm", self.dw040)?;
        write!(
            f,
            "residual rms = {:.3e} nm ({} iterations)",
            self.rms, self.iterations
        )
    }
}

#[allow(clippy::needless_range_loop)]
/// Least-squares solution of `a x = r` by Householder QR (`a` row-major, `m x n`, `m ≥ n`)
fn lstsq(mut a: Vec<Vec<f64>>, mut r: Vec<f64>) -> Vec<f64> {
    let (m, n) = (a.len(), a[0].len());
    // column scaling
    let scale: Vec<f64> = (0..n)
        .map(|j| {
            let s = (0..m).map(|i| a[i][j] * a[i][j]).sum::<f64>().sqrt();
            if s > 0. { s } else { 1. }
        })
        .collect();
    a.iter_mut()
        .for_each(|row| row.iter_mut().zip(&scale).for_each(|(v, s)| *v /= s));
    for k in 0..n {
        let norm = (k..m).map(|i| a[i][k] * a[i][k]).sum::<f64>().sqrt();
        if norm == 0. {
            continue;
        }
        let alpha = if a[k][k] > 0. { -norm } else { norm };
        let mut v: Vec<f64> = (k..m).map(|i| a[i][k]).collect();
        v[0] -= alpha;
        let vv: f64 = v.iter().map(|x| x * x).sum();
        if vv == 0. {
            continue;
        }
        for j in k..n {
            let d = 2. * (k..m).map(|i| v[i - k] * a[i][j]).sum::<f64>() / vv;
            (k..m).for_each(|i| a[i][j] -= d * v[i - k]);
        }
        let d = 2. * (k..m).map(|i| v[i - k] * r[i]).sum::<f64>() / vv;
        (k..m).for_each(|i| r[i] -= d * v[i - k]);
    }
    let mut x = vec![0f64; n];
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| a[k][j] * x[j]).sum();
        x[k] = if a[k][k] != 0. {
            (r[k] - s) / a[k][k]
        } else {
            0.
        };
    }
    x.iter().zip(&scale).map(|(x, s)| x / s).collect()
}

/// Estimates the field-height vectors of `segment` (1 to 7) from the pupil astigmatism and
/// coma `measured` at the `field` points, given the same measurements on the aligned
/// telescope, `nominal`
///
/// The part of `nominal` that the model does not reproduce is carried over unchanged,
/// so that only the misalignment signal is fitted.
pub fn estimate(
    model: &PointModel,
    layout: &Layout,
    segment: usize,
    field: &[[f64; 2]],
    measured: &[[f64; 4]],
    nominal: &[[f64; 4]],
) -> Estimate {
    assert_eq!(field.len(), measured.len());
    assert_eq!(field.len(), nominal.len());
    let s = segment_offset(layout, segment);
    let outer = segment != 7;
    let sph = model.spherical(s);
    let n = if outer { 5 } else { 4 };
    let zero = [[0f64; 2]; 2];
    let offset: Vec<[f64; 4]> = field
        .iter()
        .zip(nominal)
        .map(|(&h, b)| {
            let m = model.noll(s, h, zero);
            std::array::from_fn(|k| b[k] - m[k])
        })
        .collect();
    let forward = |th: &[f64]| -> Vec<f64> {
        let alpha = [[th[0], th[1]], [th[2], th[3]]];
        let dw = if outer { th[4] } else { 0. };
        field
            .iter()
            .zip(&offset)
            .flat_map(|(&h, o)| {
                let m = model.noll(s, h, alpha);
                (0..4).map(move |k| m[k] + o[k] + dw * sph[k])
            })
            .collect()
    };
    let y: Vec<f64> = measured.iter().flatten().copied().collect();
    let mut th = vec![0f64; n];
    let mut iterations = 0;
    for _ in 0..10 {
        iterations += 1;
        let f0 = forward(&th);
        let r: Vec<f64> = y.iter().zip(&f0).map(|(y, f)| y - f).collect();
        let e = 1e-6;
        let cols: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let mut p = th.clone();
                let mut q = th.clone();
                p[i] += e;
                q[i] -= e;
                forward(&p)
                    .iter()
                    .zip(forward(&q))
                    .map(|(a, b)| (a - b) / (2. * e))
                    .collect()
            })
            .collect();
        let jac: Vec<Vec<f64>> = (0..y.len())
            .map(|row| cols.iter().map(|c| c[row]).collect())
            .collect();
        let dth = lstsq(jac, r);
        th.iter_mut().zip(&dth).for_each(|(t, d)| *t += d);
        // converged below the finite-difference noise: 1e-10 field radius, 1e-8 nm
        if dth
            .iter()
            .enumerate()
            .all(|(i, d)| d.abs() <= if i < 4 { 1e-10 } else { 1e-8 })
        {
            break;
        }
    }
    let f0 = forward(&th);
    let rms =
        (y.iter().zip(&f0).map(|(y, f)| (y - f).powi(2)).sum::<f64>() / y.len() as f64).sqrt();
    Estimate {
        alpha: [[th[0], th[1]], [th[2], th[3]]],
        dw040: if outer { th[4] } else { 0. },
        rms,
        iterations,
    }
}

/// Third-order coefficients `ω_j131` and `ω_j222` of M1 and M2
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThirdOrder {
    pub w131: [f64; 2],
    pub w222: [f64; 2],
}

impl From<&PointModel> for ThirdOrder {
    fn from(model: &PointModel) -> Self {
        let get = |j, t| model.coefficient(j, t).unwrap_or(0.);
        Self {
            w131: [get(0, (0, 1, 1)), get(1, (0, 1, 1))],
            w222: [get(0, (0, 0, 2)), get(1, (0, 0, 2))],
        }
    }
}

/// Third-order closed-form estimate from the changes `delta` = measured - nominal of
/// Noll 5 to 8 at 2 or more `field` points
///
/// With `η = ζ_x - i ζ_y`, `σ = s_x - i s_y` and the complex coefficients
/// `A = (b6 - i b5)/√6`, `C = (b8 - i b7)/√8`: `ΔA = a1 η + a0` (least squares),
/// `C̄ = ⟨ΔC⟩`, `δ = -3 Re[(a0 - 4σC̄)/σ²]`, `β_222 = -6 a1`, `β_131 = -24 C̄ + 4σδ`
/// and `[ω_M1,222 ω_M2,222; ω_M1,131 ω_M2,131] [α_M1; α_M2] = [β_222; β_131]`.
pub fn closed_form(
    w: &ThirdOrder,
    layout: &Layout,
    segment: usize,
    field: &[[f64; 2]],
    delta: &[[f64; 4]],
) -> Estimate {
    let s = segment_offset(layout, segment);
    let sigma = C64::new(s[0], -s[1]);
    let eta: Vec<C64> = field.iter().map(|h| C64::new(h[0], -h[1])).collect();
    let da: Vec<C64> = delta
        .iter()
        .map(|b| C64::new(b[1], -b[0]) / 6f64.sqrt())
        .collect();
    let dc: Vec<C64> = delta
        .iter()
        .map(|b| C64::new(b[3], -b[2]) / 8f64.sqrt())
        .collect();
    let n = eta.len() as f64;
    // least squares ΔA = a1 η + a0
    let (se, see) = (
        eta.iter().sum::<C64>(),
        eta.iter().map(|e| e.norm_sqr()).sum::<f64>(),
    );
    let (sa, sea) = (
        da.iter().sum::<C64>(),
        eta.iter().zip(&da).map(|(e, a)| e.conj() * a).sum::<C64>(),
    );
    let det = see * n - se.norm_sqr();
    let a1 = (sea * n - se.conj() * sa) / det;
    let a0 = (sa - a1 * se) / n;
    let cm = dc.iter().sum::<C64>() / n;
    let dw040 = if sigma.norm() > 0. {
        (-3. * (a0 - 4. * sigma * cm) / (sigma * sigma)).re
    } else {
        0.
    };
    let b222 = -6. * a1;
    let b131 = -24. * cm + 4. * sigma * dw040;
    let det = w.w222[0] * w.w131[1] - w.w222[1] * w.w131[0];
    let a_m1 = (b222 * w.w131[1] - b131 * w.w222[1]) / det;
    let a_m2 = (b131 * w.w222[0] - b222 * w.w131[0]) / det;
    Estimate {
        alpha: [[a_m1.re, -a_m1.im], [a_m2.re, -a_m2.im]],
        dw040,
        rms: f64::NAN,
        iterations: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[allow(dead_code)]
    mod fixtures {
        include!("../data/field_height_fixtures.rs");
    }
    use fixtures::*;

    fn h3() -> [[f64; 2]; 3] {
        [90f64, 210., 330.].map(|a| {
            let a = a.to_radians();
            [0.8 * a.cos(), 0.8 * a.sin()]
        })
    }

    #[test]
    fn quadrature_projects_noll() {
        let q = Quadrature::new();
        // pure Noll 6 and Noll 7 wavefronts
        let w: Vec<f64> = q
            .x
            .iter()
            .zip(&q.y)
            .map(|(x, y)| {
                2. * 6f64.sqrt() * (x * x - y * y) - 8f64.sqrt() * (3. * (x * x + y * y) - 2.) * y
            })
            .collect();
        let b = q.project(&w);
        for (b, e) in b.iter().zip([0., 2., -1., 0.]) {
            assert!((b - e).abs() < 1e-13, "{b} vs {e}");
        }
    }

    #[test]
    fn aligned_model_matches_ray_trace() {
        let model = PointModel::gmt();
        let layout = Layout::entrance();
        for (seg, nominal) in [(1, B_NOM_1), (7, B_NOM_7)] {
            for (h, b) in h3().iter().zip(nominal) {
                let m = model.noll(segment_offset(&layout, seg), *h, [[0.; 2]; 2]);
                for (m, b) in m.iter().zip(b) {
                    assert!((m - b).abs() < 1e-5, "segment {seg}: {m} vs {b}");
                }
            }
        }
    }

    #[test]
    fn synthetic_recovery() {
        let model = PointModel::gmt();
        let layout = Layout::entrance();
        let field = h3();
        let alpha = [[1.2e-3, -0.7e-3], [-0.4e-3, 0.9e-3]];
        for seg in [2, 7] {
            let s = segment_offset(&layout, seg);
            let nominal = field.map(|h| model.noll(s, h, [[0.; 2]; 2]));
            let measured = field.map(|h| model.noll(s, h, alpha));
            let fit = estimate(&model, &layout, seg, &field, &measured, &nominal);
            println!("segment {seg}:\n{fit}");
            for (a, e) in fit.alpha.iter().flatten().zip(alpha.iter().flatten()) {
                assert!((a - e).abs() < 1e-10, "{a} vs {e}");
            }
            assert!(fit.dw040.abs() < 1e-8);
        }
    }

    #[test]
    fn ray_traced_misalignments() {
        // reference: Python implementation on the same ray-traced data
        let model = PointModel::gmt();
        let layout = Layout::entrance();
        let field = h3();
        let cases = [
            (1, B_NOM_1, B_M2_DX_1, A_M2_DX_1, D_M2_DX_1),
            (1, B_NOM_1, B_M2_DZ_1, A_M2_DZ_1, D_M2_DZ_1),
            (1, B_NOM_1, B_M1_RX_1, A_M1_RX_1, D_M1_RX_1),
            (7, B_NOM_7, B_M2_DX_7, A_M2_DX_7, D_M2_DX_7),
            (7, B_NOM_7, B_M1_RX_7, A_M1_RX_7, D_M1_RX_7),
        ];
        for (seg, nominal, measured, alpha, dw) in cases {
            let fit = estimate(&model, &layout, seg, &field, &measured, &nominal);
            println!("segment {seg}:\n{fit}");
            for (a, e) in fit.alpha.iter().flatten().zip(alpha.iter().flatten()) {
                assert!((a - e).abs() < 1e-9 * (1. + e.abs() * 1e3), "{a} vs {e}");
            }
            assert!((fit.dw040 - dw).abs() < 1e-6, "{} vs {dw}", fit.dw040);
        }
    }

    #[test]
    fn closed_form_matches_third_order_model() {
        // on data generated by the third-order model, the closed form is exact to first order
        let full = PointModel::gmt();
        let model = full.truncated(&[(0, 2, 0), (0, 1, 1), (0, 0, 2)]);
        let layout = Layout::entrance();
        let field = h3();
        let alpha = [[1.0e-6, -2.0e-6], [3.0e-6, 0.5e-6]];
        for seg in [1, 4, 7] {
            let s = segment_offset(&layout, seg);
            let delta = field.map(|h| {
                let (a, b) = (model.noll(s, h, alpha), model.noll(s, h, [[0.; 2]; 2]));
                std::array::from_fn::<f64, 4, _>(|k| a[k] - b[k])
            });
            let fit = closed_form(&ThirdOrder::from(&model), &layout, seg, &field, &delta);
            for (a, e) in fit.alpha.iter().flatten().zip(alpha.iter().flatten()) {
                assert!((a - e).abs() < 1e-11, "segment {seg}: {a} vs {e}");
            }
        }
    }

    #[test]
    fn parse_model() {
        assert!("0 1 1 1.0".parse::<PointModel>().is_err());
        assert!("# empty".parse::<PointModel>().is_err());
        let m: PointModel = "0 1 1 1.0 2.0 # coma".parse().unwrap();
        assert_eq!(m.coefficient(1, (0, 1, 1)), Some(2.));
    }
}
