use crseo::{
    Builder, CrseoError, FromBuilder,
    raytracing::{Conic, Rays},
};

pub mod m1 {
    pub const CONIC: f64 = 0.9982857;
    pub const CURVATURE: f64 = 36.;
}
pub mod m2 {
    pub const CONIC: f64 = 0.71692784;
    pub const CURVATURE: f64 = -4.1639009;
    pub const HEIGHT: f64 = 20.26247614;
}
pub const FOCAL_PLANE_Z: f64 = -5.83;
pub const FOCAL_PLANE_RADIUS: f64 = 2.197173;

use crate::RayTracing;

pub struct Gmt {
    m1: Conic,
    m2: Conic,
}
impl Gmt {
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
