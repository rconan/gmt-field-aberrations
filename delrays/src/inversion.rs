//! # Aberration coefficients from segment double-Zernike coefficients
//!
//! Closed-form inversion of the double-Zernike coefficients of pupil astigmatism
//! (Noll 5, 6) and pupil coma (Noll 7, 8), fitted with field Zernikes Noll 1 to 8,
//! into the wave-aberration coefficients
//! `Ω_klm = Σ_j ω_jklm` of an aligned two-mirror telescope
//! (section "Aberration coefficients from double-Zernike coefficients" of the
//! double-Zernike paper, Eqs. 31 to 36).
//!
//! Conventions:
//!  * segments are ordered as in `double_zernike`: 1 to 6 are the outer segments,
//!    7 is the centre segment;
//!  * `ρ` is normalized to the segment radius and `ζ` to the radius of the field mesh;
//!  * the outputs are in the units of the inputs (the binary prints nm);
//!  * the estimates are complex: the imaginary parts are a consistency check and
//!    should be small compared to the real parts.
//!
//! ```ignore
//! use delrays::inversion::{DoubleZernikes, Layout};
//! let dz = DoubleZernikes::from(&segments).scaled(1e9); // m -> nm
//! let omega = dz.invert(&Layout::entrance());
//! println!("{omega}");
//! ```

use std::{f64::consts::PI, fmt, ops::Index, str::FromStr};

use num_complex::Complex64 as C64;

pub mod lsq;

/// Number of segments (6 outer + 1 centre)
pub const N_SEGMENT: usize = 7;
/// Noll indices of the pupil modes: astigmatism (5, 6) and coma (7, 8)
pub const PUPIL_MODES: [usize; 4] = [5, 6, 7, 8];
/// Number of field modes (Noll 1 to 8)
pub const N_FIELD_MODE: usize = 8;

/// Real Noll double-Zernike coefficients
///
/// `b[segment][pupil Noll j - 5][field Noll j - 1]`, segment 7 (index 6) is the centre segment
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct DoubleZernikes(pub [[[f64; N_FIELD_MODE]; 4]; N_SEGMENT]);

impl DoubleZernikes {
    /// Multiplies all the coefficients by `factor` (e.g. `1e9` for m to nm)
    pub fn scaled(mut self, factor: f64) -> Self {
        self.0
            .iter_mut()
            .flatten()
            .flatten()
            .for_each(|b| *b *= factor);
        self
    }
    /// Complex coefficients of segment `id` (1 to 7)
    pub fn segment(&self, id: usize) -> SegmentCoefs {
        SegmentCoefs::from_noll(&self.0[id - 1])
    }
    /// Closed-form inversion into the aberration coefficients `Ω_klm` (Eqs. 34 to 36)
    pub fn invert(&self, layout: &Layout) -> Omega {
        cascade(self, layout)
    }
}

/// Pupil aberration
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pupil {
    /// Astigmatism `(n_ρ, m_ρ) = (2, 2)`
    Astigmatism,
    /// Coma `(n_ρ, m_ρ) = (3, 1)`
    Coma,
}

/// Field modes `(n_ζ, m_ζ)` of the complex coefficients
const FIELD_NM: [(u32, i32); 8] = [
    (0, 0),
    (1, 1),
    (1, -1),
    (2, 0),
    (2, 2),
    (2, -2),
    (3, 1),
    (3, -1),
];

/// Complex coefficients `a_{(n_ρ,m_ρ)(n_ζ,m_ζ)}` of one segment (Tables 2 and 3)
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SegmentCoefs {
    astigmatism: [C64; 8],
    coma: [C64; 8],
}

impl SegmentCoefs {
    /// Converts the real Noll coefficients of a segment, `b[pupil j - 5][field j - 1]`,
    /// into the complex coefficients (Eqs. 31 and 32)
    pub fn from_noll(b: &[[f64; N_FIELD_MODE]; 4]) -> Self {
        // Eq. (31): complex pupil coefficient from the (cos, sin) Noll pair, per field mode
        let pupil = |nr: f64, jc: usize, js: usize| -> [C64; 8] {
            let norm = (2. * (nr + 1.)).sqrt();
            std::array::from_fn(|k| C64::new(b[jc - 5][k], -b[js - 5][k]) / norm)
        };
        // Eq. (32): complex field coefficients from the (cos, sin) field Noll pairs
        let field = |c: [C64; 8]| -> [C64; 8] {
            let j = |j: usize| c[j - 1];
            let rot = |nh: f64, jc: usize, js: usize, sign: f64| {
                (j(jc) - sign * C64::i() * j(js)) / (2. * (nh + 1.)).sqrt()
            };
            [
                j(1),
                rot(1., 2, 3, 1.),
                rot(1., 2, 3, -1.),
                j(4) / 3f64.sqrt(),
                rot(2., 6, 5, 1.),
                rot(2., 6, 5, -1.),
                rot(3., 8, 7, 1.),
                rot(3., 8, 7, -1.),
            ]
        };
        Self {
            astigmatism: field(pupil(2., 6, 5)),
            coma: field(pupil(3., 8, 7)),
        }
    }
    /// Returns `a_{(n_ρ,m_ρ)(n_ζ,m_ζ)}`
    ///
    /// Panics if `(n_ζ, m_ζ)` is not one of the field modes up to radial order 3
    /// reached by Noll 1 to 8
    pub fn get(&self, pupil: Pupil, nh: u32, mh: i32) -> C64 {
        let k = FIELD_NM
            .iter()
            .position(|&nm| nm == (nh, mh))
            .unwrap_or_else(|| panic!("field mode ({nh},{mh}) is not in Noll 1 to 8"));
        match pupil {
            Pupil::Astigmatism => self.astigmatism[k],
            Pupil::Coma => self.coma[k],
        }
    }
}

/// GMT segment layout
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Layout {
    /// Offset of the outer segments in units of the segment radius
    pub s: f64,
    /// Azimuth of outer segments 1 to 6 \[rd\]
    pub azimuth: [f64; 6],
}

impl Layout {
    /// Segment offset: 8.71 m / 4.1825 m
    pub const S: f64 = 8.71 / 4.1825;
    /// Entrance-pupil layout: segment k at azimuth (3 - 2(k - 1))π/6
    pub fn entrance() -> Self {
        Self {
            s: Self::S,
            azimuth: std::array::from_fn(|k| (3. - 2. * k as f64) * PI / 6.),
        }
    }
    /// Exit-pupil layout: the entrance layout rotated by 180°
    pub fn exit() -> Self {
        let mut layout = Self::entrance();
        layout.azimuth.iter_mut().for_each(|a| *a += PI);
        layout
    }
    /// `σ_k = s exp(-i φ_k)` for the 6 outer segments
    fn sigma(&self) -> [C64; 6] {
        self.azimuth.map(|phi| C64::from_polar(self.s, -phi))
    }
}

impl Default for Layout {
    fn default() -> Self {
        Self::entrance()
    }
}

/// Aberration coefficients `Ω_klm = Σ_j ω_jklm`
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Omega {
    pub w040: C64,
    pub w131: C64,
    pub w222: C64,
    pub w240: C64,
    pub w242: C64,
    pub w331: C64,
    pub w060: C64,
    pub w151: C64,
    pub w260: C64,
    pub w262: C64,
    pub w351: C64,
}

impl Omega {
    /// Coefficient labels, in the order of [Omega::to_array]
    pub const LABELS: [&'static str; 11] = [
        "w040", "w131", "w222", "w240", "w242", "w331", "w060", "w151", "w260", "w262", "w351",
    ];
    /// The 11 coefficients in the order of [Omega::LABELS]
    pub fn to_array(&self) -> [C64; 11] {
        [
            self.w040, self.w131, self.w222, self.w240, self.w242, self.w331, self.w060, self.w151,
            self.w260, self.w262, self.w351,
        ]
    }
    /// Real parts in the order of [Omega::LABELS]
    pub fn real(&self) -> [f64; 11] {
        self.to_array().map(|w| w.re)
    }
    /// Iterator over `(label, coefficient)`
    pub fn iter(&self) -> impl Iterator<Item = (&'static str, C64)> {
        Self::LABELS.into_iter().zip(self.to_array())
    }
}

impl Index<&str> for Omega {
    type Output = C64;
    fn index(&self, label: &str) -> &C64 {
        match label {
            "w040" => &self.w040,
            "w131" => &self.w131,
            "w222" => &self.w222,
            "w240" => &self.w240,
            "w242" => &self.w242,
            "w331" => &self.w331,
            "w060" => &self.w060,
            "w151" => &self.w151,
            "w260" => &self.w260,
            "w262" => &self.w262,
            "w351" => &self.w351,
            _ => panic!("unknown aberration coefficient {label}"),
        }
    }
}

impl fmt::Display for Omega {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let p = f.precision().unwrap_or(4);
        writeln!(f, "{:>6} {:>14} {:>11}", "", "Re", "Im")?;
        for (label, w) in self.iter() {
            writeln!(f, "{label:>6} {:>14.p$} {:>+11.2e}", w.re, w.im)?;
        }
        Ok(())
    }
}

/// Solves the complex 2x2 system `m x = y` (Cramer's rule)
fn solve2(m: [[f64; 2]; 2], y: [C64; 2]) -> [C64; 2] {
    let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    [
        (y[0] * m[1][1] - y[1] * m[0][1]) / det,
        (y[1] * m[0][0] - y[0] * m[1][0]) / det,
    ]
}

/// Closed-form inversion (Eqs. 34 to 36)
fn cascade(dz: &DoubleZernikes, layout: &Layout) -> Omega {
    use Pupil::{Astigmatism as A, Coma as C};

    let s = layout.s;
    let s2 = s * s;
    let sig = layout.sigma();
    let segs: [SegmentCoefs; N_SEGMENT] = std::array::from_fn(|i| dz.segment(i + 1));
    let centre = |p: Pupil, nh: u32, mh: i32| segs[6].get(p, nh, mh);
    // ⟨a σ^-n⟩ (n > 0) or ⟨a conj(σ)^-1⟩ (conj) averaged over the 6 outer segments
    let mean = |p: Pupil, nh: u32, mh: i32, f: &dyn Fn(C64) -> C64| -> C64 {
        segs[..6]
            .iter()
            .zip(&sig)
            .map(|(seg, &sg)| seg.get(p, nh, mh) * f(sg))
            .sum::<C64>()
            / 6.
    };
    let inv = |n: i32| move |sg: C64| sg.powi(-n);

    // Eq. (34)
    let w351 = 96. * mean(C, 3, 1, &inv(2));
    let w262 = 36. * mean(C, 2, 2, &inv(3));
    let w151 = 16. * mean(C, 1, 1, &inv(2)) - 2. / 3. * w351;
    let w331 = 288. * centre(C, 3, -1) - 6. / 5. * w351;
    let w131 = 48. * centre(C, 1, -1) - 2. / 3. * w331 - 6. / 5. * w151 - 4. / 5. * w351;
    let w242 = 48. * mean(C, 2, -2, &|sg: C64| 1. / sg.conj()) - (8. / 5. + 4. * s2) * w262;
    let w222 = 36. * centre(A, 2, -2) - 3. / 4. * w242 - 3. / 5. * w262;

    // Eq. (35): field defocus
    let [w240, w260] = solve2(
        [
            [1. / 18., 1. / 8. + s2 / 6.],
            [1. / 36., 1. / 20. + s2 / 8.],
        ],
        [
            mean(A, 2, 0, &inv(2)) - w242 / 36. - (1. / 16. + s2 / 12.) * w262,
            mean(C, 2, 0, &inv(1)) - w242 / 72. - (1. / 40. + s2 / 16.) * w262,
        ],
    );

    // Eq. (36): field piston
    let [w040, w060] = solve2(
        [[1. / 3., 3. / 4. + s2], [1. / 6., 3. / 10. + 3. * s2 / 4.]],
        [
            mean(A, 0, 0, &inv(2))
                - w240 / 6.
                - w242 / 12.
                - (3. / 8. + s2 / 2.) * w260
                - (3. / 16. + s2 / 4.) * w262,
            mean(C, 0, 0, &inv(1))
                - w240 / 12.
                - w242 / 24.
                - (3. / 20. + 3. * s2 / 8.) * w260
                - (3. / 40. + 3. * s2 / 16.) * w262,
        ],
    );

    Omega {
        w040,
        w131,
        w222,
        w240,
        w242,
        w331,
        w060,
        w151,
        w260,
        w262,
        w351,
    }
}

/// Error parsing the text output of `double_zernike`
#[derive(Debug, Clone, PartialEq)]
pub struct ParseError(String);
impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "double-Zernike table parse error: {}", self.0)
    }
}
impl std::error::Error for ParseError {}

/// Parses the table printed by the `double_zernike` binary
///
/// A pupil mode header `(j, n, m)` is followed by rows `( j,  n,  m): b_1 ... b_7`
/// of field modes; other lines are ignored, as are pupil modes other than Noll 5 to 8
/// and field modes beyond Noll 8. All 4 x 8 rows must be present.
impl FromStr for DoubleZernikes {
    type Err = ParseError;
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let jnm = |s: &str| -> Option<usize> {
            let s = s.trim().strip_prefix('(')?.strip_suffix(')')?;
            s.split(',').next()?.trim().parse().ok()
        };
        let mut b = Self::default();
        let mut found = [[false; N_FIELD_MODE]; 4];
        let mut pupil: Option<usize> = None;
        for (n, line) in text.lines().enumerate() {
            let line = line.trim();
            if let Some((head, values)) = line.split_once(':') {
                let (Some(jp), Some(jf)) = (pupil, jnm(head)) else {
                    continue;
                };
                if !(1..=N_FIELD_MODE).contains(&jf) {
                    continue;
                }
                let values = values
                    .split_whitespace()
                    .map(|v| v.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| ParseError(format!("line {}: {e}", n + 1)))?;
                if values.len() != N_SEGMENT {
                    return Err(ParseError(format!(
                        "line {}: expected {N_SEGMENT} segment values, found {}",
                        n + 1,
                        values.len()
                    )));
                }
                for (seg, v) in b.0.iter_mut().zip(values) {
                    seg[jp - 5][jf - 1] = v;
                }
                found[jp - 5][jf - 1] = true;
            } else if let Some(j) = jnm(line) {
                pupil = PUPIL_MODES.contains(&j).then_some(j);
            }
        }
        for (i, row) in found.iter().enumerate() {
            if let Some(k) = row.iter().position(|f| !f) {
                return Err(ParseError(format!(
                    "missing field mode {} of pupil mode {}",
                    k + 1,
                    i + 5
                )));
            }
        }
        Ok(b)
    }
}

impl From<&crate::zernikes::SegmentsDoubleZernikes> for DoubleZernikes {
    /// Extracts pupil Noll 5 to 8 and field Noll 1 to 8 (unscaled, i.e. in m)
    fn from(segments: &crate::zernikes::SegmentsDoubleZernikes) -> Self {
        let mut b = Self::default();
        for (seg, field_coefs) in b.0.iter_mut().zip(segments.iter()) {
            for pupil_mode in field_coefs.iter() {
                let jp = pupil_mode.jnm.0;
                if !PUPIL_MODES.contains(&jp) {
                    continue;
                }
                for field_mode in pupil_mode.coef.iter() {
                    let jf = field_mode.jnm.0;
                    if (1..=N_FIELD_MODE).contains(&jf) {
                        seg[jp - 5][jf - 1] = field_mode.coef;
                    }
                }
            }
        }
        b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ENTRANCE: &str = include_str!("../data/double_zernike_entrance.txt");
    const SYNTHETIC: &str = include_str!("../data/double_zernike_synthetic.txt");

    #[test]
    fn noll_to_complex_round_trip() {
        // a pure Noll 6 (cos 2θ) pupil x Noll 2 (field x) term gives
        // a_{(2,2)(1,±1)} = 1/(√6 · 2)
        let mut b = [[0f64; 8]; 4];
        b[1][1] = 1.;
        let a = SegmentCoefs::from_noll(&b);
        let expected = 1. / (6f64.sqrt() * 2.);
        assert!((a.get(Pupil::Astigmatism, 1, 1).re - expected).abs() < 1e-15);
        assert!((a.get(Pupil::Astigmatism, 1, -1).re - expected).abs() < 1e-15);
        assert_eq!(a.get(Pupil::Coma, 1, 1), C64::new(0., 0.));
    }

    #[test]
    fn ceo_entrance() {
        // reference: Python implementation on the same (3-decimal) table, nm
        let reference = [
            ("w351", 0.019500597, -0.000399389),
            ("w262", -0.000863031, -0.000095892),
            ("w151", -0.436313331, 0.000360396),
            ("w331", 76.872599283, 0.000479267),
            ("w131", -26.540401524, -0.000432476),
            ("w242", -5.006241945, 0.001539632),
            ("w222", -2972.436800723, -0.001097188),
            ("w240", 7.446557199, -0.000979904),
            ("w260", 0.000003552, 0.000057802),
            ("w040", -0.327163741, -0.000060647),
            ("w060", 0.026705641, 0.000002845),
        ];
        let dz: DoubleZernikes = ENTRANCE.parse().unwrap();
        let omega = dz.invert(&Layout::entrance());
        println!("{omega:.6}");
        for (label, re, im) in reference {
            let w = omega[label];
            assert!((w.re - re).abs() < 1e-8, "{label}: {} vs {re}", w.re);
            assert!((w.im - im).abs() < 1e-8, "{label}: {} vs {im}", w.im);
        }
    }

    #[test]
    fn synthetic_exact() {
        // coefficients computed by quadrature from the 11-term wave aberration alone
        // (independent of Tables 2 and 3): the inversion must be exact
        let truth = [
            ("w040", -0.3911),
            ("w131", -25.7604),
            ("w222", -2972.0177),
            ("w240", 7.5680),
            ("w242", -4.9098),
            ("w331", 74.5477),
            ("w060", 0.0319),
            ("w151", -0.4364),
            ("w260", 0.0002),
            ("w262", -0.0006),
            ("w351", 0.0131),
        ];
        let dz: DoubleZernikes = SYNTHETIC.parse().unwrap();
        let omega = dz.invert(&Layout::entrance());
        println!("{omega:.6}");
        for (label, t) in truth {
            let w = omega[label];
            assert!((w.re - t).abs() < 1e-9, "{label}: {} vs {t}", w.re);
            assert!(w.im.abs() < 1e-9, "{label}: Im {}", w.im);
        }
    }

    #[test]
    fn parse_errors() {
        assert!(
            "(5, 2, 2)\n ( 1,  0,  0): 1 2 3"
                .parse::<DoubleZernikes>()
                .is_err()
        );
        assert!("".parse::<DoubleZernikes>().is_err());
    }
}
