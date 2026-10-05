//! # GMT optical prescription

use crseo::{
    Builder, CrseoError, FromBuilder,
    raytracing::{Conic, Rays},
};

/// M1 optical parameters
pub mod m1 {
    /// Conic constant
    pub const CONIC: f64 = 0.9982857;
    /// Mirror curvature
    pub const CURVATURE: f64 = 36.;
}
/// M2 optical parameters
pub mod m2 {
    /// Conic constant
    pub const CONIC: f64 = 0.71692784;
    /// Mirror curvature
    pub const CURVATURE: f64 = -4.1639009;
    /// Mirror height
    pub const HEIGHT: f64 = 20.26247614;
}
/// Nominal focal plane height (rounded)
///
/// The paraxial focus of the prescription is 5.3 µm below, at [paraxial_focus_z];
/// a reference sphere centred at -5.83 m adds 1.08 nm of defocus and, because the
/// chief ray is inclined at the focus, 0.006 nm of field-linear coma (Ω131)
pub const FOCAL_PLANE_Z: f64 = -5.83;

/// Paraxial image, by M2, of an on-axis point at height `z_object` \[m\]
///
/// Mirror equation for M2 with distances from its vertex along z:
/// `1/(z_i - z_2) + 1/(z_o - z_2) = 2/r_2`
pub fn m2_paraxial_image(z_object: f64) -> f64 {
    m2::HEIGHT + 1. / (2. / m2::CURVATURE - 1. / (z_object - m2::HEIGHT))
}
/// Height of the paraxial focus of the telescope \[m\]: image by M2 of the M1 focus at `r_1/2`
pub fn paraxial_focus_z() -> f64 {
    m2_paraxial_image(0.5 * m1::CURVATURE)
}
/// Height of the exit pupil \[m\]: image by M2 of the M1 vertex (the stop)
pub fn exit_pupil_z() -> f64 {
    m2_paraxial_image(0.)
}
/// Nominal focal surface radius (CEO value)
///
/// Within 2.1 mm of the medial surface of the prescription, [focal_surface_radius]
/// with [FocalSurface::Medial]; the radius of the reference surface changes the
/// field-cubic coma ω331 by 0.0056 nm and ω240 by 1.06e-4 nm per nm of field defocus
pub const FOCAL_PLANE_RADIUS: f64 = 2.197173;

/// Third-order focal surfaces
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FocalSurface {
    Petzval,
    Sagittal,
    Medial,
    Tangential,
}

/// Seidel astigmatism sum `S_III` of M1 and M2 for a unit marginal ray height at M1
/// and a unit chief ray angle, stop at M1 (Welford's conventions: light along +z,
/// `n' = -n` at a mirror, curvature `c = -1/r` with `r` the radius used by CEO)
fn seidel_astigmatism() -> f64 {
    // (curvature, conic constant, distance to the next mirror)
    let mirrors = [
        (-1. / m1::CURVATURE, -m1::CONIC, -m2::HEIGHT),
        (-1. / m2::CURVATURE, -m2::CONIC, 0.),
    ];
    let (mut h, mut u, mut hb, mut ub, mut n) = (1f64, 0f64, 0f64, 1f64, 1f64);
    let mut s3 = 0.;
    for (c, k, d) in mirrors {
        let np = -n;
        let up = (n * u - h * c * (np - n)) / np;
        let ubp = (n * ub - hb * c * (np - n)) / np;
        let ab = n * (ub + hb * c);
        let du = up / np - u / n;
        s3 += -ab * ab * h * du + k * c.powi(3) * (np - n) * h * h * hb * hb;
        (u, ub, n) = (up, ubp, np);
        h += d * u;
        hb += d * ub;
    }
    s3
}

/// Vertex radius of curvature \[m\] of a third-order focal surface of the telescope:
/// `1/R = 2/r_1 + 2/r_2 - j S_III`, `j = 0, 1, 2, 3` for Petzval, sagittal, medial and
/// tangential, with `S_III` from [seidel_astigmatism] (equal to `2|Ω222|/(R²θ²)`)
pub fn focal_surface_radius(surface: FocalSurface) -> f64 {
    let petzval = 2. / m1::CURVATURE - 2. / m2::CURVATURE;
    let j = match surface {
        FocalSurface::Petzval => 0.,
        FocalSurface::Sagittal => 1.,
        FocalSurface::Medial => 2.,
        FocalSurface::Tangential => 3.,
    };
    1. / (petzval - j * seidel_astigmatism())
}

use crate::trace::RayTracing;

/// GMT M1 and M2 optical model
pub struct Gmt {
    m1: Conic,
    m2: Conic,
}
impl Gmt {
    /// Creates a new [Gmt] instance
    pub fn new() -> Result<Self, CrseoError> {
        let m1 = Conic::builder()
            .conic_cst(1. - m1::CONIC)
            .curvature_radius(m1::CURVATURE)
            .build()?;
        let m2 = Conic::builder()
            .conic_cst(1. - m2::CONIC)
            .curvature_radius(m2::CURVATURE)
            .origin([0f64, 0f64, m2::HEIGHT])
            .build()?;
        Ok(Self { m1, m2 })
    }
}
impl RayTracing for Gmt {
    fn ray_tracing(&mut self, rays: &mut Rays) {
        self.m1.trace(rays);
        self.m2.trace(rays);
        // reference sphere centred on the medial focal surface through the paraxial focus,
        // both computed from the prescription (not the rounded FOCAL_PLANE_Z/RADIUS)
        rays.to_sphere(
            paraxial_focus_z(),
            focal_surface_radius(FocalSurface::Medial),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paraxial_focus_and_exit_pupil() {
        // high-precision values from the same prescription
        assert!((paraxial_focus_z() - (-5.83000531118439)).abs() < 1e-11);
        assert!((exit_pupil_z() - 17.9421102695571).abs() < 1e-11);
        // CEO's on-axis reference-sphere radius is z_E - FOCAL_PLANE_Z
        assert!((exit_pupil_z() - FOCAL_PLANE_Z - 23.772110269559725).abs() < 1e-11);
    }

    #[test]
    fn focal_surfaces() {
        // radii at which the ray-traced field defocus of each surface vanishes
        let expected = [
            (FocalSurface::Petzval, 1.86611),
            (FocalSurface::Sagittal, 2.01728),
            (FocalSurface::Medial, 2.19510),
            (FocalSurface::Tangential, 2.40729),
        ];
        for (surface, r) in expected {
            let radius = focal_surface_radius(surface);
            assert!((radius - r).abs() < 1e-5, "{surface:?}: {radius} vs {r}");
        }
        // the CEO value is within 2.1 mm of the medial surface
        assert!((focal_surface_radius(FocalSurface::Medial) - FOCAL_PLANE_RADIUS).abs() < 2.2e-3);
    }
}
