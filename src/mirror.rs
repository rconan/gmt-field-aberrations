use std::ops::{Add, Sub};

use crseo::Gmt;
use serde::{Deserialize, Serialize};
use skyangle::SkyAngle;

use crate::segment::{rbm::Rbm, units::Txyz};

/// GMT mirror selection
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum Mirror {
    M1(Rbm),
    M2(Rbm),
    Both(Rbm, Rbm),
}
impl Mirror {
    pub fn apply_rbms(&self, sid: i32, gmt: &mut Gmt) {
        let to_meters =
            |t_xyz: &[Txyz; 3]| t_xyz.into_iter().map(|x| x.as_f64()).collect::<Vec<f64>>();

        let to_radians = |r_xyz: &[SkyAngle<f64>; 3]| {
            r_xyz
                .into_iter()
                .map(|x| x.to_radians())
                .collect::<Vec<f64>>()
        };
        match self {
            Mirror::M1(Rbm { t_xyz, r_xyz }) => {
                gmt.m1_segment_state(sid, &to_meters(t_xyz), &to_radians(r_xyz))
            }
            Mirror::M2(Rbm { t_xyz, r_xyz }) => {
                gmt.m2_segment_state(sid, &to_meters(t_xyz), &to_radians(r_xyz))
            }
            Mirror::Both(
                Rbm {
                    t_xyz: m1_t_xyz,
                    r_xyz: m1_r_xyz,
                },
                Rbm {
                    t_xyz: m2_t_xyz,
                    r_xyz: m2_r_xyz,
                },
            ) => {
                gmt.m1_segment_state(sid, &to_meters(m1_t_xyz), &to_radians(m1_r_xyz));
                gmt.m2_segment_state(sid, &to_meters(m2_t_xyz), &to_radians(m2_r_xyz))
            }
        }
    }
}
impl Add for Mirror {
    type Output = Mirror;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Mirror::M1(x), Mirror::M1(y)) => Mirror::M1(x + y),
            (Mirror::M1(x), Mirror::M2(y)) => Mirror::Both(x, y),
            (Mirror::M1(x), Mirror::Both(y, z)) => Mirror::Both(x + y, z),
            (Mirror::M2(x), Mirror::M1(y)) => Mirror::Both(y, x),
            (Mirror::M2(x), Mirror::M2(y)) => Mirror::M2(x + y),
            (Mirror::M2(x), Mirror::Both(y, z)) => Mirror::Both(y, x + z),
            (Mirror::Both(x, y), Mirror::M1(z)) => Mirror::Both(x + z, y),
            (Mirror::Both(x, y), Mirror::M2(z)) => Mirror::Both(x, y + z),
            (Mirror::Both(x, y), Mirror::Both(z, w)) => Mirror::Both(x + z, y + w),
        }
    }
}
impl Sub for Mirror {
    type Output = Mirror;

    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Mirror::M1(x), Mirror::M1(y)) => Mirror::M1(x - y),
            (Mirror::M1(x), Mirror::M2(y)) => Mirror::Both(x, -y),
            (Mirror::M1(x), Mirror::Both(y, z)) => Mirror::Both(x - y, -z),
            (Mirror::M2(x), Mirror::M1(y)) => Mirror::Both(-y, x),
            (Mirror::M2(x), Mirror::M2(y)) => Mirror::M2(x - y),
            (Mirror::M2(x), Mirror::Both(y, z)) => Mirror::Both(-y, x - z),
            (Mirror::Both(x, y), Mirror::M1(z)) => Mirror::Both(x - z, y),
            (Mirror::Both(x, y), Mirror::M2(z)) => Mirror::Both(x, y - z),
            (Mirror::Both(x, y), Mirror::Both(z, w)) => Mirror::Both(x - z, y - w),
        }
    }
}
