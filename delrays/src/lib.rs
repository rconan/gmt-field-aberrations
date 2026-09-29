use crseo::{
    FromBuilder,
    raytracing::{Rays, RaysBuilder},
};
use triangle_rs::Delaunay;

mod delaunay;
mod gmt;
mod opd;
mod set;
pub use delaunay::Mesh;
pub use gmt::Gmt;
pub use opd::{Opd, Opds};
pub use set::Set;

#[derive(Debug, thiserror::Error)]
pub enum DelraysError {
    #[error("crseo failure")]
    Crseo(#[from] crseo::CrseoError),
}

pub trait RayTracing {
    fn ray_tracing(self, rays: &mut Rays);
}
pub trait Trace {
    fn from_mesh(mesh: &Delaunay) -> RaysBuilder;
    fn trace(&mut self, object: impl RayTracing);
    fn opds(&mut self) -> Opds;
}
impl Trace for Rays {
    fn from_mesh(mesh: &Delaunay) -> RaysBuilder {
        Rays::builder().xy(mesh.vertex_iter().flatten().cloned().collect())
    }
    fn trace(&mut self, object: impl RayTracing) {
        object.ray_tracing(self);
    }
    fn opds(&mut self) -> Opds {
        let opds = self.optical_path_difference();
        let xyz = self.coordinates();
        xyz.chunks(3)
            .zip(opds.into_iter())
            .map(|(xyz, delta)| Opd {
                xyz: xyz.try_into().unwrap(),
                delta,
            })
            .collect()
    }
}
