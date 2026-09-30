//! Rays optical path differences

use std::fmt::Display;

use faer::MatRef;
use serde::{Deserialize, Serialize};

use crate::{
    Set,
    zernikes::{Mode, Zernike},
};

/// Single ray optical path difference
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Opd {
    pub(crate) xyz: [f64; 3],
    pub(crate) delta: f64,
}

impl Opd {
    /// Checks if two OPDs share the same coordinates
    pub fn match_coordinates(&self, other: &Opd) -> bool {
        self.xyz
            .iter()
            .zip(other.xyz.iter())
            .all(|(lhs, rhs)| lhs == rhs)
    }
}

/// Rays optical path differences
pub type Opds = Set<Opd>;

impl Opds {
    /// Project OPDs on Zernike modes
    pub fn as_zernikes(&self, n_radial_order: usize) -> Set<Mode> {
        let mut zerns = Set::<Mode>::new(self.0.iter().map(|opd| &opd.xyz[..2]), n_radial_order);
        let delta: Vec<_> = self.0.iter().map(|opd| opd.delta).collect();
        let mat = MatRef::<f64>::from_column_major_slice(&delta, delta.len(), 1);
        let a = zerns.pseudo_inverse() * mat;
        a.col(0).iter().zip(zerns.iter_mut()).for_each(|(c, z)| {
            z.coef = *c;
        });
        zerns
    }
    /// Substract two OPDs
    ///
    /// Returns None if coordinates do not match
    pub fn sub(&self, rhs: &Opds) -> Option<Opds> {
        self.0
            .iter()
            .zip(rhs.0.iter())
            .map(|(lhs, rhs)| {
                if lhs.match_coordinates(rhs) {
                    Some(Opd {
                        xyz: rhs.xyz.clone(),
                        delta: rhs.delta - lhs.delta,
                    })
                } else {
                    None
                }
            })
            .collect()
    }
}

/// OPDs statistics
#[derive(Default, Debug, Clone)]
pub struct Stats {
    n_sample: usize,
    min_max: (f64, f64),
    mean: f64,
    var: f64,
}
impl Display for Stats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "sample #: {};min/max: [{:.3e},{:.3e}]; mean: {:.3e}; std: {:.3e}",
            self.n_sample,
            self.min_max.0,
            self.min_max.1,
            self.mean,
            self.var.sqrt()
        )
    }
}
impl From<&Opds> for Stats {
    fn from(opds: &Opds) -> Self {
        let n_sample = opds.0.len();
        let min = opds
            .0
            .iter()
            .min_by(|a, b| a.delta.partial_cmp(&b.delta).unwrap())
            .unwrap()
            .delta;
        let max = opds
            .0
            .iter()
            .max_by(|a, b| a.delta.partial_cmp(&b.delta).unwrap())
            .unwrap()
            .delta;
        let mean = opds.0.iter().map(|opd| opd.delta).sum::<f64>() / n_sample as f64;
        let var = opds
            .0
            .iter()
            .map(|opd| opd.delta - mean)
            .map(|x| x * x)
            .sum::<f64>()
            / n_sample as f64;
        Self {
            n_sample,
            min_max: (min, max),
            mean,
            var,
        }
    }
}

impl From<&Set<Opds>> for Stats {
    fn from(set: &Set<Opds>) -> Self {
        let (n_sample, mins, maxs, means, vars) = set.iter().fold(
            (0usize, vec![], vec![], 0f64, 0f64),
            |(mut n, mut mins, mut maxs, mut means, mut vars), opds| {
                let Stats {
                    n_sample,
                    min_max: (min, max),
                    mean,
                    var,
                } = opds.into();
                n += n_sample;
                mins.push(min);
                maxs.push(max);
                means += mean * n_sample as f64;
                vars += var * n_sample as f64;
                (n, mins, maxs, means, vars)
            },
        );
        let min = mins
            .into_iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let max = maxs
            .into_iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        Self {
            n_sample,
            min_max: (min, max),
            mean: means / n_sample as f64,
            var: vars / n_sample as f64,
        }
    }
}
