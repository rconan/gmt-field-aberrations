use crseo::{
    Builder, CrseoError,
    raytracing::{Rays, RaysBuilder},
};
use triangle_rs::Delaunay;

use crate::{
    Mesh, Set, opd::{Opd, Opds},
};

/// Ray tracing through an optical system
pub trait RayTracing {
    /// Draws rays through the system
    fn ray_tracing(&mut self, rays: &mut Rays);
}
/// [Rays] extensions
pub trait Trace {
    /// Type of the object used to initialize the [Rays] coordinates
    type From;
    /// Type of the object that contains the [Rays] optical path differences
    type OpdData;
    /// Creates a new object from a mesh and a [Rays] builder
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError>
    where
        Self: Sized;
    /// Traces the [Rays]
    fn trace<T: RayTracing>(&mut self, object: &mut T) -> &mut Self;
    /// Retrieves the [Rays] OPD
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
    fn trace<T: RayTracing>(&mut self, object: &mut T) -> &mut Self {
        object.ray_tracing(self);
        self
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
    type From = Set<Mesh>;
    type OpdData = Set<Opds>;
    fn from_mesh(mesh: &Self::From, builder: RaysBuilder) -> Result<Self, CrseoError>
    where
        Self: Sized,
    {
        mesh.iter()
            .map(|del| Rays::from_mesh(del, builder.clone()))
            .collect()
    }

    fn trace<T: RayTracing>(&mut self, object: &mut T) -> &mut Self {
        self.iter_mut().for_each(|rays| {
            rays.trace(object);
        });
        self
    }

    fn opds(&mut self) -> Self::OpdData {
        self.iter_mut().map(|rays| rays.opds()).collect()
    }
}
