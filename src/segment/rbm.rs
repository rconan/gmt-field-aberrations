use std::ops::Add;

use serde::{Deserialize, Serialize};
use skyangle::Conversion;

use super::{
    dof::Dof,
    units::{Rxyz, Txyz},
};

/// Rigid body motions
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
pub struct Rbm {
    /// Translations
    pub t_xyz: [Txyz; 3],
    /// Rotations
    pub r_xyz: [Rxyz; 3],
}
impl Rbm {
    pub fn t_x(v: Txyz) -> Self {
        let v0 = Txyz::default();
        Self {
            t_xyz: [v, v0, v0],
            r_xyz: [Default::default(); 3],
        }
    }
    pub fn t_y(v: Txyz) -> Self {
        let v0 = Txyz::default();
        Self {
            t_xyz: [v0, v, v0],
            r_xyz: [Default::default(); 3],
        }
    }
    pub fn t_z(v: Txyz) -> Self {
        let v0 = Txyz::default();
        Self {
            t_xyz: [v0, v0, v],
            r_xyz: [Default::default(); 3],
        }
    }
    pub fn r_x(v: Rxyz) -> Self {
        let v0 = Rxyz::default();
        Self {
            t_xyz: [Default::default(); 3],
            r_xyz: [v, v0, v0],
        }
    }
    pub fn r_y(v: Rxyz) -> Self {
        let v0 = Rxyz::default();
        Self {
            t_xyz: [Default::default(); 3],
            r_xyz: [v0, v, v0],
        }
    }
}
impl From<Dof> for Rbm {
    fn from(value: Dof) -> Self {
        match value {
            Dof::Tx(txyz) => Rbm::t_x(txyz),
            Dof::Ty(txyz) => Rbm::t_y(txyz),
            Dof::Rx(sky_angle) => Rbm::r_x(sky_angle),
            Dof::Ry(sky_angle) => Rbm::r_y(sky_angle),
        }
    }
}
impl Add for Rbm {
    type Output = Rbm;

    fn add(self, rhs: Self) -> Self::Output {
        let Rbm { t_xyz: t, r_xyz: r } = self;
        let Rbm {
            t_xyz: tt,
            r_xyz: rr,
        } = rhs;
        let ta: Vec<_> = t
            .into_iter()
            .zip(tt.into_iter())
            .map(|(x, y)| x.as_f64() + y.as_f64())
            .map(|z| Txyz::Nm(1e9 * z))
            .collect();
        let ra: Vec<_> = r
            .into_iter()
            .zip(rr.into_iter())
            .map(|(x, y)| x.to_radians() + y.to_radians())
            .map(|z| Rxyz::MilliArcsec(z.to_mas()))
            .collect();
        Rbm {
            t_xyz: ta.try_into().unwrap(),
            r_xyz: ra.try_into().unwrap(),
        }
    }
}
