use crseo::{
    Builder, CrseoError, FromBuilder,
    raytracing::{Conic, Rays},
};

use crate::RayTracing;

pub struct Gmt {
    m1: Conic,
    m2: Conic,
}
impl Gmt {
    pub fn new() -> Result<Self, CrseoError> {
        let m1 = Conic::builder()
            .conic_cst(1. - 0.9982857)
            .curvature_radius(36.)
            .build()?;
        let m2 = Conic::builder()
            .conic_cst(1. - 0.71692784)
            .curvature_radius(-4.1639009)
            .origin([0f64, 0f64, 20.26247614])
            .build()?;
        Ok(Self { m1, m2 })
    }
}
impl RayTracing for &mut Gmt {
    fn ray_tracing(self, rays: &mut Rays) {
        self.m1.trace(rays);
        self.m2.trace(rays);
        rays.to_sphere(-5.83, 2.197173);
    }
}
