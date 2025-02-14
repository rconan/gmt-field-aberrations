use std::{ops::Neg, vec};

use super::units::{Rxyz, Txyz};

#[derive(Debug, Clone, Copy)]
pub enum Dof {
    Tx(Txyz),
    Ty(Txyz),
    Tz(Txyz),
    Rx(Rxyz),
    Ry(Rxyz),
}
impl Dof {
    pub fn r_xy() -> vec::IntoIter<Dof> {
        vec![Dof::Rx(Rxyz::Arcsecond(1.)), Dof::Ry(Rxyz::Arcsecond(1.))].into_iter()
    }
    pub fn into_iter() -> vec::IntoIter<Dof> {
        vec![
            Dof::Tx(Txyz::Mu(1.)),
            Dof::Ty(Txyz::Mu(1.)),
            Dof::Tz(Txyz::Mu(1.)),
            Dof::Rx(Rxyz::Arcsecond(1.)),
            Dof::Ry(Rxyz::Arcsecond(1.)),
        ]
        .into_iter()
    }
    pub fn as_f64(&self) -> f64 {
        match self {
            Dof::Tx(txyz) | Dof::Ty(txyz) | Dof::Tz(txyz) => txyz.as_f64(),
            Dof::Rx(sky_angle) | Dof::Ry(sky_angle) => sky_angle.to_radians(),
        }
    }
    pub fn len() -> usize {
        4
    }
}
impl Neg for Dof {
    type Output = Self;

    fn neg(self) -> Self::Output {
        match self {
            Dof::Tx(txyz) => Dof::Tx(-txyz),
            Dof::Ty(txyz) => Dof::Ty(-txyz),
            Dof::Tz(txyz) => Dof::Tz(-txyz),
            Dof::Rx(sky_angle) => Dof::Rx(-sky_angle),
            Dof::Ry(sky_angle) => Dof::Ry(-sky_angle),
        }
    }
}
