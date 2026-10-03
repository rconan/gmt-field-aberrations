//! # Zernike basis

use faer::Mat;
use serde::{Deserialize, Serialize};

use crate::{Mesh, Set};

mod r#as;
pub mod fmt;
pub use r#as::AsZernikes;

/// Decomposition of Zernike coeffients field map in to Zernike modes
pub type FieldZernikeCoefficients = Set<Mode<Set<Mode>>>;

impl FieldZernikeCoefficients {
    /// Returns the mean of the field maps of Zernike coefficient
    pub fn mean(&self, field: &Mesh) -> Vec<f64> {
        let weights = field.lump_mass_matrix_weights();
        let area = field.area();
        self.iter()
            .map(|mode| {
                mode.mode
                    .iter()
                    .zip(&weights)
                    .map(|(c, w)| c * w)
                    .sum::<f64>()
                    / area
            })
            .collect()
    }
}

/// GMT segment double zernike expansion
pub type SegmentsDoubleZernikes = Set<FieldZernikeCoefficients>;

/// A Zernike mode or a Zernike coefficient field map
///
/// For a Zernike mode, the coefficient is a scalar (default) or
/// if the mode is a field map of a particular Zernike coefficient
/// then the coefficient is another [Set] of [Mode] with scalar coefficients
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Mode<C = f64> {
    pub jnm: (usize, usize, usize),
    pub(crate) mode: Vec<f64>,
    pub coef: C,
}

/// Zernike basis traits
pub trait Zernike {
    /// Zernike coefficients type
    type Coefs;
    /// Creates a new Zernike basis
    fn new<'a>(xy: impl Iterator<Item = &'a [f64]>, n_radial_order: usize) -> Set<Mode> {
        let (mut r, o): (Vec<f64>, Vec<f64>) = xy
            .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])))
            .unzip();
        let max_r = r
            .iter()
            .cloned()
            .max_by(|x, y| x.partial_cmp(y).unwrap())
            .unwrap();
        r.iter_mut().for_each(|r| *r /= max_r);

        let (j, n, m) = zernike::jnm(n_radial_order as u32);
        j.into_iter()
            .zip(n.into_iter())
            .zip(m.into_iter())
            .map(|((j, n), m)| {
                let mode = r
                    .iter()
                    .zip(&o)
                    .map(|(&r, &o)| zernike::zernike(j, n, m, r, o))
                    .collect::<Vec<_>>();
                Mode {
                    jnm: (j as usize, n as usize, m as usize),
                    mode,
                    coef: 0f64,
                }
            })
            .collect()
    }
    /// Computes the pseudo-inverse of the Zernike basis
    fn pseudo_inverse(&self) -> Mat<f64>;
    fn weighted_pseudo_inverse(&self, weights: &[f64]) -> Mat<f64>;
    /// Returns the indices of the Zernike modes
    fn jnm(&self) -> Vec<(usize, usize, usize)>;
    /// Returns the coefficients of the Zernike modes
    fn coefficients(&self) -> Vec<Self::Coefs>;
    /// Reduces the Zernike basis to the given modes
    fn reduce_into(self, j: impl IntoIterator<Item = usize>) -> Self;
    /// Return the modes self-projection matrix
    fn speye(&self) -> Mat<f64>;
}

impl<C: Clone> Zernike for Set<Mode<C>> {
    type Coefs = C;
    fn pseudo_inverse(&self) -> Mat<f64> {
        let ncols = self.len();
        let nrows = self[0].mode.len();
        let zerns: Vec<_> = self.iter().flat_map(|mode| mode.mode.to_vec()).collect();
        let mat = faer::MatRef::from_column_major_slice(&zerns, nrows, ncols);
        let svd = mat.svd().unwrap();
        svd.pseudoinverse()
    }
    fn weighted_pseudo_inverse(&self, weights: &[f64]) -> Mat<f64> {
        let ncols = self.len();
        let nrows = self[0].mode.len();
        let zerns: Vec<_> = self
            .iter()
            .flat_map(|mode| {
                mode.mode
                    .iter()
                    .zip(weights)
                    .map(|(m, w)| m * w)
                    .collect::<Vec<_>>()
            })
            .collect();
        let mat = faer::MatRef::from_column_major_slice(&zerns, nrows, ncols);
        let svd = mat.svd().unwrap();
        svd.pseudoinverse()
    }

    fn jnm(&self) -> Vec<(usize, usize, usize)> {
        self.iter().map(|mode| mode.jnm.clone()).collect()
    }

    fn coefficients(&self) -> Vec<Self::Coefs> {
        self.iter().map(|mode| mode.coef.clone()).collect()
    }

    fn reduce_into(self, j: impl IntoIterator<Item = usize>) -> Self {
        j.into_iter()
            .flat_map(|j| self.iter().filter(move |mode| mode.jnm.0 == j))
            .cloned()
            .collect()
    }

    fn speye(&self) -> Mat<f64> {
        let ncols = self.len();
        let nrows = self[0].mode.len();
        let zerns: Vec<_> = self.iter().flat_map(|mode| mode.mode.to_vec()).collect();
        let mat = faer::MatRef::from_column_major_slice(&zerns, nrows, ncols);
        mat.transpose() * mat
    }
}

/// A Zernike basis and a field location
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct FieldZernike {
    za: (f64, f64),
    modes: Set<Mode>,
}

impl FieldZernike {
    /// Associates a Zernike basis to a field point `(zenith,azimuth)`
    pub fn new(za: (f64, f64), modes: Set<Mode>) -> Self {
        Self { za, modes }
    }
}
