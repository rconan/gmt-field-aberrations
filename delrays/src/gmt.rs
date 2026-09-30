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
/// Focal plane height
pub const FOCAL_PLANE_Z: f64 = -5.83;
/// Focal surface radius
pub const FOCAL_PLANE_RADIUS: f64 = 2.197173;

use crate::RayTracing;

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
        rays.to_sphere(FOCAL_PLANE_Z, FOCAL_PLANE_RADIUS);
    }
}
