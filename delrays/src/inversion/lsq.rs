//! # Least-squares inversion of double-Zernike coefficients
//!
//! Generalizes [super::DoubleZernikes::invert] to any set of pupil modes (e.g. Noll 1 to 10),
//! any set of field modes (e.g. Noll 1 to 21, radial order 5) and any polynomial model
//! `W = Σ_t ω_t (ζ·ζ)^p (ρ_s·ρ_s)^n (ζ·ρ_s)^m` of an aligned telescope
//! (paragraph "Adding pupil modes 1-10" of the double-Zernike paper).
//!
//! The design matrix is the exact projection of each term onto the pupil Noll modes over
//! each segment and onto the field Noll modes over the unit field disc (Gauss-Legendre
//! quadrature, exact for the polynomial degrees involved). Terms that produce none of the
//! fitted modes are reported as not observable; the others are solved by a singular value
//! decomposition of the column-normalized matrix, and each coefficient is flagged as
//! identifiable when it lies in the row space of the matrix. Coefficients that are not
//! identifiable are given their minimum-norm values and must not be used.
//!
//! ```ignore
//! use delrays::inversion::{Layout, lsq::{terms, Coefficients, Design}};
//! let data = Coefficients::from(&segments).scaled(1e9);              // pupil 1-10, field 1-21
//! let design = Design::new(terms(5, 8, 3), data.pupil(), data.field(), &Layout::entrance());
//! let solution = design.solve(&data)?;
//! println!("{solution}");
//! ```

use std::{f64::consts::PI, fmt, str::FromStr};

use faer::Mat;

use super::Layout;

/// A term `(p, n, m)` of the wave-aberration polynomial, i.e. `ω_klm` with `k = 2p + m`, `l = 2n + m`
pub type Term = (u32, u32, u32);

/// Number of segments (6 outer + 1 centre)
pub const N_SEGMENT: usize = 7;

/// The terms with `k ≤ kmax`, `l ≤ lmax` and `m ≤ mmax`, ordered by `m`, `p`, `n`
pub fn terms(kmax: u32, lmax: u32, mmax: u32) -> Vec<Term> {
    let mut out = vec![];
    for m in 0..=mmax {
        for p in 0..=kmax / 2 {
            for n in 0..=lmax / 2 {
                if 2 * p + m <= kmax && 2 * n + m <= lmax {
                    out.push((p, n, m));
                }
            }
        }
    }
    out
}

/// Label `w{k}{l}{m}` of a term
pub fn label(t: Term) -> String {
    let (p, n, m) = t;
    format!("w{}{}{}", 2 * p + m, 2 * n + m, m)
}

/// Radial order and azimuthal frequency `(n, m)` of the Noll mode `j`
pub fn noll_nm(j: usize) -> (u32, u32) {
    assert!(j >= 1, "Noll indices start at 1");
    let mut n = 0usize;
    while (n + 1) * (n + 2) / 2 < j {
        n += 1;
    }
    let k = j - n * (n + 1) / 2 - 1;
    let ms: Vec<usize> = (n % 2..=n)
        .step_by(2)
        .flat_map(|m| if m == 0 { vec![0] } else { vec![m, m] })
        .collect();
    (n as u32, ms[k] as u32)
}

/// Noll-normalized Zernike polynomial `j` at polar coordinates `(r, θ)`
fn noll(j: usize, r: f64, theta: f64) -> f64 {
    let (n, m) = noll_nm(j);
    zernike::zernike(j as u32, n, m, r, theta)
}

/// Quadrature of the unit disc, exact for polynomials of degree `degree`;
/// returns `(x, y, r, θ, w)` with `Σ w = 1` (weights normalized by π)
fn disc_quadrature(degree: usize) -> Vec<[f64; 5]> {
    let n_r = degree / 2 + 2;
    let n_t = degree + 2;
    let (xr, wr) = gauss_legendre(n_r);
    let mut nodes = vec![];
    for (x, w) in xr.iter().zip(&wr) {
        let r = 0.5 * (x + 1.);
        for k in 0..n_t {
            let t = 2. * PI * k as f64 / n_t as f64;
            nodes.push([
                r * t.cos(),
                r * t.sin(),
                r,
                t,
                0.5 * w * r * 2. / n_t as f64,
            ]);
        }
    }
    nodes
}

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

/// Segment offsets in units of the segment radius (segments 1 to 7)
fn offsets(layout: &Layout) -> [[f64; 2]; N_SEGMENT] {
    std::array::from_fn(|k| {
        if k == 6 {
            [0., 0.]
        } else {
            let phi = layout.azimuth[k];
            [layout.s * phi.cos(), layout.s * phi.sin()]
        }
    })
}

/// Real Noll double-Zernike coefficients of the 7 segments, for any pupil and field modes
#[derive(Debug, Clone, PartialEq)]
pub struct Coefficients {
    pupil: Vec<usize>,
    field: Vec<usize>,
    /// `data[ip * field.len() + jf] = [segment 1, ..., segment 7]`
    data: Vec<[f64; N_SEGMENT]>,
}

impl Coefficients {
    /// New coefficients; `data` is ordered by pupil mode, then field mode
    pub fn new(pupil: Vec<usize>, field: Vec<usize>, data: Vec<[f64; N_SEGMENT]>) -> Self {
        assert_eq!(data.len(), pupil.len() * field.len());
        Self { pupil, field, data }
    }
    /// Pupil Noll modes
    pub fn pupil(&self) -> &[usize] {
        &self.pupil
    }
    /// Field Noll modes
    pub fn field(&self) -> &[usize] {
        &self.field
    }
    /// Coefficients of pupil mode `jp` and field mode `jf` for the 7 segments
    pub fn get(&self, jp: usize, jf: usize) -> Option<&[f64; N_SEGMENT]> {
        let ip = self.pupil.iter().position(|&j| j == jp)?;
        let i_f = self.field.iter().position(|&j| j == jf)?;
        self.data.get(ip * self.field.len() + i_f)
    }
    /// Multiplies all the coefficients by `factor` (e.g. `1e9` for m to nm)
    pub fn scaled(mut self, factor: f64) -> Self {
        self.data.iter_mut().flatten().for_each(|v| *v *= factor);
        self
    }
    /// Keeps the given pupil and field modes, if present
    pub fn select(&self, pupil: &[usize], field: &[usize]) -> Option<Self> {
        let mut data = vec![];
        for &jp in pupil {
            for &jf in field {
                data.push(*self.get(jp, jf)?);
            }
        }
        Some(Self::new(pupil.to_vec(), field.to_vec(), data))
    }
}

impl From<&crate::zernikes::SegmentsDoubleZernikes> for Coefficients {
    /// All the pupil and field modes of the first segment (unscaled, i.e. in m)
    fn from(segments: &crate::zernikes::SegmentsDoubleZernikes) -> Self {
        let first = &segments[0];
        let pupil: Vec<usize> = first.iter().map(|mode| mode.jnm.0).collect();
        let field: Vec<usize> = first[0].coef.iter().map(|mode| mode.jnm.0).collect();
        let mut data = vec![[0f64; N_SEGMENT]; pupil.len() * field.len()];
        for (k, seg) in segments.iter().take(N_SEGMENT).enumerate() {
            for (ip, pupil_mode) in seg.iter().enumerate() {
                for (i_f, field_mode) in pupil_mode.coef.iter().enumerate() {
                    data[ip * field.len() + i_f][k] = field_mode.coef;
                }
            }
        }
        Self::new(pupil, field, data)
    }
}

/// Error of the least-squares inversion
#[derive(Debug, Clone, PartialEq)]
pub enum LsqError {
    /// The data lack a pupil or field mode of the design
    MissingMode { pupil: usize, field: usize },
    /// Parse error of a double-Zernike table
    Parse(String),
}
impl fmt::Display for LsqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingMode { pupil, field } => {
                write!(f, "missing pupil mode {pupil} x field mode {field}")
            }
            Self::Parse(e) => write!(f, "double-Zernike table parse error: {e}"),
        }
    }
}
impl std::error::Error for LsqError {}

/// A row of a double-Zernike table: field mode and segment values
type Row = (usize, [f64; N_SEGMENT]);

/// Parses a table printed by `double_zernike`, with any pupil and field modes
///
/// A pupil mode header `(j, n, m)` is followed by rows `( j,  n,  m): b_1 ... b_7`;
/// every pupil mode must list the same field modes. Other lines are ignored.
impl FromStr for Coefficients {
    type Err = LsqError;
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let first = |s: &str| -> Option<usize> {
            let s = s.trim().strip_prefix('(')?.strip_suffix(')')?;
            s.split(',').next()?.trim().parse().ok()
        };
        let mut blocks: Vec<(usize, Vec<Row>)> = vec![];
        for (n, line) in text.lines().enumerate() {
            let line = line.trim();
            if let Some((head, values)) = line.split_once(':') {
                let (Some(block), Some(jf)) = (blocks.last_mut(), first(head)) else {
                    continue;
                };
                let v = values
                    .split_whitespace()
                    .map(|v| v.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| LsqError::Parse(format!("line {}: {e}", n + 1)))?;
                let v: [f64; N_SEGMENT] = v.try_into().map_err(|v: Vec<f64>| {
                    LsqError::Parse(format!(
                        "line {}: expected {N_SEGMENT} values, found {}",
                        n + 1,
                        v.len()
                    ))
                })?;
                block.1.push((jf, v));
            } else if let Some(jp) = first(line) {
                blocks.push((jp, vec![]));
            }
        }
        blocks.retain(|(_, rows)| !rows.is_empty());
        let Some((_, rows)) = blocks.first() else {
            return Err(LsqError::Parse("no data".into()));
        };
        let field: Vec<usize> = rows.iter().map(|(jf, _)| *jf).collect();
        if blocks
            .iter()
            .any(|(_, rows)| rows.iter().map(|(jf, _)| *jf).ne(field.iter().copied()))
        {
            return Err(LsqError::Parse(
                "pupil modes with different field modes".into(),
            ));
        }
        let pupil = blocks.iter().map(|(jp, _)| *jp).collect();
        let data = blocks
            .into_iter()
            .flat_map(|(_, rows)| rows.into_iter().map(|(_, v)| v))
            .collect();
        Ok(Self::new(pupil, field, data))
    }
}

/// Design matrix of a polynomial model for given pupil and field modes
#[derive(Debug, Clone)]
pub struct Design {
    terms: Vec<Term>,
    pupil: Vec<usize>,
    field: Vec<usize>,
    /// rows: pupil mode, field mode, segment; columns: terms
    a: Mat<f64>,
}

impl Design {
    /// Exact projections of each term onto the pupil Noll modes `pupil` over each segment
    /// and onto the field Noll modes `field` over the unit field disc
    pub fn new(terms: Vec<Term>, pupil: &[usize], field: &[usize], layout: &Layout) -> Self {
        let order = |modes: &[usize]| modes.iter().map(|&j| noll_nm(j).0).max().unwrap_or(0);
        let kmax = terms.iter().map(|&(p, _, m)| 2 * p + m).max().unwrap_or(0);
        let lmax = terms.iter().map(|&(_, n, m)| 2 * n + m).max().unwrap_or(0);
        let pq = disc_quadrature((lmax + order(pupil)) as usize);
        let fq = disc_quadrature((kmax + order(field)) as usize);
        let zp: Vec<Vec<f64>> = pupil
            .iter()
            .map(|&j| pq.iter().map(|q| noll(j, q[2], q[3]) * q[4]).collect())
            .collect();
        let zf: Vec<Vec<f64>> = field
            .iter()
            .map(|&j| fq.iter().map(|q| noll(j, q[2], q[3]) * q[4]).collect())
            .collect();
        let (np, nf) = (pupil.len(), field.len());
        let mut a = Mat::<f64>::zeros(np * nf * N_SEGMENT, terms.len());
        for (k, s) in offsets(layout).iter().enumerate() {
            for (c, &(p, n, m)) in terms.iter().enumerate() {
                // g[ip][f] = pupil projection at field node f
                let mut g = vec![vec![0f64; fq.len()]; np];
                for (f, h) in fq.iter().enumerate() {
                    let hh = h[0] * h[0] + h[1] * h[1];
                    let hp = hh.powi(p as i32);
                    for (i, q) in pq.iter().enumerate() {
                        let (rx, ry) = (q[0] + s[0], q[1] + s[1]);
                        let w = hp
                            * (rx * rx + ry * ry).powi(n as i32)
                            * (h[0] * rx + h[1] * ry).powi(m as i32);
                        for ip in 0..np {
                            g[ip][f] += zp[ip][i] * w;
                        }
                    }
                }
                for ip in 0..np {
                    for i_f in 0..nf {
                        let v: f64 = zf[i_f].iter().zip(&g[ip]).map(|(z, g)| z * g).sum();
                        a[((ip * nf + i_f) * N_SEGMENT + k, c)] = v;
                    }
                }
            }
        }
        Self {
            terms,
            pupil: pupil.to_vec(),
            field: field.to_vec(),
            a,
        }
    }
    /// Model terms
    pub fn terms(&self) -> &[Term] {
        &self.terms
    }
    /// Design-matrix entry for a term, pupil mode, field mode and segment (1 to 7)
    pub fn entry(&self, term: Term, jp: usize, jf: usize, segment: usize) -> Option<f64> {
        let c = self.terms.iter().position(|&t| t == term)?;
        let ip = self.pupil.iter().position(|&j| j == jp)?;
        let i_f = self.field.iter().position(|&j| j == jf)?;
        Some(self.a[((ip * self.field.len() + i_f) * N_SEGMENT + segment - 1, c)])
    }
    /// Forward model: coefficients of the given aberration coefficients (one per term)
    pub fn forward(&self, omega: &[f64]) -> Coefficients {
        assert_eq!(omega.len(), self.terms.len());
        let data = (0..self.pupil.len() * self.field.len())
            .map(|i| {
                std::array::from_fn(|k| {
                    omega
                        .iter()
                        .enumerate()
                        .map(|(c, w)| self.a[(i * N_SEGMENT + k, c)] * w)
                        .sum()
                })
            })
            .collect();
        Coefficients::new(self.pupil.clone(), self.field.clone(), data)
    }
    /// Least-squares solution for `data`, which must contain the pupil and field modes of the design
    pub fn solve(&self, data: &Coefficients) -> Result<Solution, LsqError> {
        let mut b = vec![];
        for &jp in &self.pupil {
            for &jf in &self.field {
                let v = data.get(jp, jf).ok_or(LsqError::MissingMode {
                    pupil: jp,
                    field: jf,
                })?;
                b.extend_from_slice(v);
            }
        }
        let (rows, cols) = (self.a.nrows(), self.a.ncols());
        let norms: Vec<f64> = (0..cols)
            .map(|c| {
                (0..rows)
                    .map(|r| self.a[(r, c)].powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .collect();
        let max_norm = norms.iter().cloned().fold(0., f64::max);
        let observable: Vec<bool> = norms.iter().map(|&c| c > 1e-10 * max_norm).collect();
        let live: Vec<usize> = (0..cols).filter(|&c| observable[c]).collect();
        let an = Mat::<f64>::from_fn(rows, live.len(), |r, c| {
            self.a[(r, live[c])] / norms[live[c]]
        });
        let svd = an.thin_svd().expect("SVD failed");
        let (u, s, v) = (svd.U(), svd.S().column_vector(), svd.V());
        let s0 = s[0];
        let rank = (0..live.len()).filter(|&i| s[i] > 1e-9 * s0).count();
        let mut coefficients = vec![0f64; cols];
        let mut identifiable = vec![false; cols];
        for (c, &col) in live.iter().enumerate() {
            let mut x = 0.;
            let mut leverage = 0.;
            for i in 0..rank {
                let ub: f64 = (0..rows).map(|r| u[(r, i)] * b[r]).sum();
                x += v[(c, i)] * ub / s[i];
                leverage += v[(c, i)].powi(2);
            }
            coefficients[col] = x / norms[col];
            identifiable[col] = leverage > 1. - 1e-6;
        }
        let residual_rms = ((0..rows)
            .map(|r| {
                let m: f64 = (0..cols).map(|c| self.a[(r, c)] * coefficients[c]).sum();
                (b[r] - m).powi(2)
            })
            .sum::<f64>()
            / rows as f64)
            .sqrt();
        Ok(Solution {
            terms: self.terms.clone(),
            coefficients,
            observable,
            identifiable,
            rank,
            condition: s0 / s[rank - 1],
            residual_rms,
        })
    }
}

/// Least-squares solution
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    /// Model terms
    pub terms: Vec<Term>,
    /// Aberration coefficients `Ω_klm`, in the units of the data (minimum-norm values
    /// for the coefficients that are not identifiable)
    pub coefficients: Vec<f64>,
    /// Whether the term produces any of the fitted modes
    pub observable: Vec<bool>,
    /// Whether the coefficient is determined by the data
    pub identifiable: Vec<bool>,
    /// Rank of the design matrix restricted to the observable terms
    pub rank: usize,
    /// Condition number of the column-normalized design matrix on its range
    pub condition: f64,
    /// Rms of the fit residuals
    pub residual_rms: f64,
}

impl Solution {
    /// Coefficient of a term, if it is identifiable
    pub fn get(&self, term: Term) -> Option<f64> {
        let i = self.terms.iter().position(|&t| t == term)?;
        self.identifiable[i].then_some(self.coefficients[i])
    }
    /// Coefficient by label (e.g. `"w131"`), if it is identifiable
    pub fn get_label(&self, name: &str) -> Option<f64> {
        let i = self.terms.iter().position(|&t| label(t) == name)?;
        self.identifiable[i].then_some(self.coefficients[i])
    }
    /// Number of observable terms
    pub fn n_observable(&self) -> usize {
        self.observable.iter().filter(|&&o| o).count()
    }
}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let p = f.precision().unwrap_or(4);
        writeln!(
            f,
            "{} terms, {} observable, rank {}, condition {:.1e}, residual rms {:.2e}",
            self.terms.len(),
            self.n_observable(),
            self.rank,
            self.condition,
            self.residual_rms
        )?;
        for (i, &t) in self.terms.iter().enumerate() {
            let flag = if !self.observable[i] {
                "  (not observable)"
            } else if !self.identifiable[i] {
                "  (not identifiable)"
            } else {
                ""
            };
            if self.observable[i] && self.identifiable[i] {
                writeln!(f, "{:>6} {:>16.p$}", label(t), self.coefficients[i])?;
            } else {
                writeln!(f, "{:>6} {:>16}{flag}", label(t), "-")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LSQ: &str = include_str!("../../data/double_zernike_lsq.txt");

    #[test]
    fn noll_indices() {
        for order in 1..7u32 {
            let (j, n, m) = zernike::jnm(order + 1);
            for ((j, n), m) in j.into_iter().zip(n).zip(m) {
                assert_eq!(noll_nm(j as usize), (n, m), "j={j}");
            }
        }
    }

    #[test]
    fn term_order() {
        let t = terms(5, 6, 3);
        assert_eq!(t.len(), 31);
        assert_eq!(terms(5, 8, 3).len(), 41);
        assert_eq!(t[0], (0, 0, 0));
        assert_eq!(label((1, 2, 1)), "w351");
    }

    #[test]
    fn design_entries() {
        // reference: Python quadrature (inverse.py / general.py)
        let pupil: Vec<usize> = (1..=10).collect();
        let field: Vec<usize> = (1..=21).collect();
        let d = Design::new(
            vec![(0, 1, 1), (1, 2, 0)],
            &pupil,
            &field,
            &Layout::entrance(),
        );
        let cases = [
            ((0, 1, 1), 8, 2, 1, 5.892556509887809e-02),
            ((1, 2, 0), 4, 4, 7, 8.333333333333312e-02),
            ((0, 1, 1), 7, 3, 3, 5.892556509887775e-02),
        ];
        for (t, jp, jf, seg, e) in cases {
            let v = d.entry(t, jp, jf, seg).unwrap();
            assert!((v - e).abs() < 1e-13, "{t:?} {jp} {jf} {seg}: {v} vs {e}");
        }
    }

    #[test]
    fn synthetic_recovery() {
        let pupil: Vec<usize> = (1..=10).collect();
        let field: Vec<usize> = (1..=21).collect();
        let design = Design::new(terms(5, 8, 3), &pupil, &field, &Layout::entrance());
        let omega: Vec<f64> = (0..design.terms().len())
            .map(|i| ((i * 7 + 3) as f64).sin())
            .collect();
        let data = design.forward(&omega);
        let sol = design.solve(&data).unwrap();
        assert_eq!(sol.rank, 41);
        for (x, w) in sol.coefficients.iter().zip(&omega) {
            assert!((x - w).abs() < 1e-8, "{x} vs {w}");
        }
    }

    #[test]
    fn exact_projections_of_ray_traced_wavefront() {
        // 136-term GMT wavefront projected on pupil 1-10 x field 1-21;
        // reference: Python least squares (pupil10 analysis), model k<=5, l<=8, m<=3
        type Case = (
            &'static [usize],
            usize,
            usize,
            [(&'static str, f64, bool); 11],
        );
        let cases: [Case; 3] = [
            (
                &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                41,
                41,
                [
                    ("w040", -3.896996523238e-01, true),
                    ("w131", -2.574823947710e+01, true),
                    ("w222", -2.972012319641e+03, true),
                    ("w240", 7.568543551435e+00, true),
                    ("w242", -4.908084254144e+00, true),
                    ("w331", 7.445795769662e+01, true),
                    ("w060", 3.133893266729e-02, true),
                    ("w151", -4.364575865439e-01, true),
                    ("w260", 2.043224409600e-04, true),
                    ("w262", -7.005965265024e-04, true),
                    ("w351", 1.316491274544e-02, true),
                ],
            ),
            (
                &[4, 5, 6, 7, 8],
                35,
                35,
                [
                    ("w040", -3.899791095000e-01, true),
                    ("w131", -2.574823355332e+01, true),
                    ("w222", -2.972012318908e+03, true),
                    ("w240", 7.568544269092e+00, true),
                    ("w242", -4.908084508545e+00, true),
                    ("w331", 7.445795664410e+01, true),
                    ("w060", 3.137861546349e-02, true),
                    ("w151", -4.364588233414e-01, true),
                    ("w260", 2.042218613501e-04, true),
                    ("w262", -7.005573471467e-04, true),
                    ("w351", 1.316516656724e-02, true),
                ],
            ),
            (
                &[5, 6, 7, 8],
                32,
                29,
                [
                    ("w040", f64::NAN, false),
                    ("w131", -2.574823355609e+01, true),
                    ("w222", -2.972012319112e+03, true),
                    ("w240", f64::NAN, false),
                    ("w242", -4.908083926228e+00, true),
                    ("w331", 7.445795613405e+01, true),
                    ("w060", f64::NAN, false),
                    ("w151", -4.364588285070e-01, true),
                    ("w260", f64::NAN, false),
                    ("w262", -7.006749592379e-04, true),
                    ("w351", 1.316524315199e-02, true),
                ],
            ),
        ];
        let all: Coefficients = LSQ.parse().unwrap();
        assert_eq!(all.pupil().len(), 10);
        assert_eq!(all.field().len(), 21);
        for (pupil, n_obs, rank, expected) in cases {
            let data = all.select(pupil, all.field()).unwrap();
            let design = Design::new(terms(5, 8, 3), pupil, data.field(), &Layout::entrance());
            let sol = design.solve(&data).unwrap();
            println!("pupil {pupil:?}: {sol:.6}");
            assert_eq!(sol.n_observable(), n_obs, "pupil {pupil:?}");
            assert_eq!(sol.rank, rank, "pupil {pupil:?}");
            for (name, value, ok) in expected {
                let got = sol.get_label(name);
                assert_eq!(got.is_some(), ok, "{name} identifiability, pupil {pupil:?}");
                if let Some(x) = got {
                    assert!(
                        (x - value).abs() < 1e-7 * (1. + value.abs()),
                        "{name}, pupil {pupil:?}: {x} vs {value}"
                    );
                }
            }
        }
    }

    #[test]
    fn parse_errors() {
        assert!("".parse::<Coefficients>().is_err());
        assert!(
            "(5, 2, 2)\n ( 1,  0,  0): 1 2 3"
                .parse::<Coefficients>()
                .is_err()
        );
        let two = "(5, 2, 2)\n (1, 0, 0): 1 2 3 4 5 6 7\n(6, 2, 2)\n (2, 1, 1): 1 2 3 4 5 6 7";
        assert!(two.parse::<Coefficients>().is_err());
    }
}
