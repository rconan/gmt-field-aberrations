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
/// Focal surface radius
pub const FOCAL_PLANE_RADIUS: f64 = 2.197173;

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
        // reference sphere centred at the paraxial focus (not the rounded FOCAL_PLANE_Z)
        rays.to_sphere(paraxial_focus_z(), FOCAL_PLANE_RADIUS);
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
}
