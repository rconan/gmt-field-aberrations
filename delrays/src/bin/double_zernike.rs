use std::io::{Write, stdout};

use crseo::{CrseoError, FromBuilder, raytracing::Rays};
use delrays::{
    Gmt, Mesh, Set, Trace,
    zernikes::{AsZernikes, FieldZernike, Zernike},
};
use skyangle::Conversion;
use triangle_rs::Delaunay;

const PUPIL_MODES: [usize; 2] = [5, 6];

fn main() -> anyhow::Result<()> {
    let field_mesh = Delaunay::disc(20f64, 2., None);
    field_mesh.plot();
    // println!("{field_mesh}");

    let mut lock = stdout().lock();

    let mut results = vec![];
    for id in 1..=7 {
        write!(lock, "{id}").unwrap();
        lock.flush().unwrap();

        let delaunay = Delaunay::gmt_segment_with_id(id)?;

        let mut gmt = Gmt::new()?;

        let iter = field_mesh
            .vertex_iter()
            .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])));

        let field_zernikes = iter
            .map(|(zen, azi)| {
                let rays_builder = Rays::builder().zenith(zen.from_arcmin()).azimuth(azi);
                let mut rays: Rays = Rays::from_mesh(&delaunay, rays_builder)?;
                let modes = rays
                    .trace(&mut gmt)
                    .opds()
                    .as_zernikes(4)
                    .reduce_into(PUPIL_MODES);
                Ok(FieldZernike::new((zen, azi), modes))
            })
            .collect::<Result<Set<FieldZernike>, CrseoError>>()?;
        let field_zern_coefs = field_zernikes.as_zernikes(4);
        results.push(field_zern_coefs);
    }
    println!();

    'pupil: for j in 0.. {
        'field: for k in 0.. {
            for (i, field_zern_coefs) in results.iter().enumerate() {
                let Some(fzc) = field_zern_coefs.get(j) else {
                    break 'pupil;
                };
                let coef = fzc.coef.clone().reduce_into(1..=8);
                let jnm = coef.jnm();
                let coef = coef.coefficients();
                let Some((jnm, c)) = jnm.get(k).zip(coef.get(k)) else {
                    break 'field;
                };
                if i == 0 {
                    if k == 0 {
                        println!("{:?}", fzc.jnm);
                    }
                    print!(" {:2?}: ", jnm);
                }
                print!("{:+6.0}", c * 1e9);
            }
            println!()
        }
    }
    Ok(())
}
