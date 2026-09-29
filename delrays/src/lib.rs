use crseo::{
    Builder, CrseoError, raytracing::{Rays, RaysBuilder},
};
use triangle_rs::Delaunay;

mod delaunay;
mod gmt;
mod opd;
mod set;
pub use delaunay::Mesh;
pub use gmt::Gmt;
pub use opd::{Opd, Opds, Stats};
pub use set::Set;

#[derive(Debug, thiserror::Error)]
pub enum DelraysError {
    #[error("crseo failure")]
    Crseo(#[from] crseo::CrseoError),
}

pub trait RayTracing {
    fn ray_tracing(&mut self, rays: &mut Rays);
}
pub trait Trace {
    type From;
    type OpdData;
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError>
    where
        Self: Sized;
    fn trace<T: RayTracing>(&mut self, object: &mut T);
    fn opds(&mut self) -> Self::OpdData;
}
impl Trace for Rays {
    type From = Delaunay;
    type OpdData = Opds;
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError> {
        builder
            .xy(mesh.vertex_iter().flatten().cloned().collect())
            .build()
    }
    fn trace<T: RayTracing>(&mut self, object: &mut T) {
        object.ray_tracing(self);
    }
    fn opds(&mut self) -> Self::OpdData {
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
impl Trace for Set<Rays> {
    type From = Set<Delaunay>;
    type OpdData = Set<Opds>;
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError>
    where
        Self: Sized,
    {
        mesh.iter()
            .map(|del| Rays::from_mesh(del, builder.clone()))
            .collect()
    }

    fn trace<T: RayTracing>(&mut self, object: &mut T) {
        self.iter_mut().for_each(|rays| rays.trace(object));
    }

    fn opds(&mut self) -> Self::OpdData {
        self.iter_mut().map(|rays| rays.opds()).collect()
    }
}
