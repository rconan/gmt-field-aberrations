use faer::MatRef;

use super::{FieldZernike, FieldZernikeCoefficients, Mode, Zernike};
use crate::{Mesh, Set, opd::Opds};

/// Decomposition into Zernike modes
pub trait AsZernikes {
    /// Type of the modal expansion
    type Into;
    /// Projects object on Zernike modes
    ///
    /// The projection is a weighted least squares using the [Mesh] lump mass matrix weights
    fn as_zernikes(&self, n_radial_order: usize, mesh: &Mesh) -> Self::Into;
}

impl AsZernikes for Opds {
    type Into = Set<Mode>;
    /// Projects OPDs on Zernike modes
    fn as_zernikes(&self, n_radial_order: usize, mesh: &Mesh) -> Self::Into {
        let weights: Vec<_> = mesh
            .lump_mass_matrix_weights()
            .into_iter()
            .map(f64::sqrt)
            .collect();
        let mut zerns = Set::<Mode>::new(self.0.iter().map(|opd| &opd.xyz[..2]), n_radial_order);
        let delta: Vec<_> = self
            .0
            .iter()
            .zip(&weights)
            .map(|(opd, w)| opd.delta * w)
            .collect();
        let mat = MatRef::<f64>::from_column_major_slice(&delta, delta.len(), 1);
        let a = zerns.weighted_pseudo_inverse(&weights) * mat;
        a.col(0).iter().zip(zerns.iter_mut()).for_each(|(c, z)| {
            z.coef = *c;
        });
        zerns
    }
}

impl AsZernikes for Set<FieldZernike> {
    type Into = FieldZernikeCoefficients;
    /// Projects field maps of Zernike coefficients on Zernike modes
    fn as_zernikes(&self, n_radial_order: usize, mesh: &Mesh) -> Self::Into {
        let weights: Vec<_> = mesh
            .lump_mass_matrix_weights()
            .into_iter()
            .map(f64::sqrt)
            .collect();
        let xy: Vec<_> = self
            .iter()
            .map(|field| {
                let (z, a) = field.za;
                let (s, c) = a.sin_cos();
                [z * c, z * s]
            })
            .collect();
        let mut zerns = Set::<Mode>::new(xy.iter().map(|xy| xy.as_slice()), n_radial_order);
        let n_mode = self.0[0].modes.len();
        let pinv = zerns.weighted_pseudo_inverse(&weights);
        let mut modes = vec![];
        for i in 0..n_mode {
            let c: Vec<_> = self
                .iter()
                .map(|field| field.modes[i].coef)
                .zip(&weights)
                .map(|(c, w)| c * w)
                .collect();
            let mat = MatRef::<f64>::from_column_major_slice(&c, c.len(), 1);
            let a = &pinv * mat;

            a.col(0).iter().zip(zerns.iter_mut()).for_each(|(c, z)| {
                z.coef = *c;
            });
            let Mode { jnm, .. } = self[0].modes[i].clone();
            modes.push(Mode {
                jnm,
                mode: c,
                coef: zerns.clone(),
            });
        }
        modes.into_iter().collect()
    }
}
