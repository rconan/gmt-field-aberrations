use std::{
    io::{Write, stdout},
    thread,
    time::Instant,
};

use crseo::{CrseoError, FromBuilder, raytracing::Rays};
use delrays::{
    Gmt, Mesh, Set, Trace,
    zernikes::{AsZernikes, FieldZernike, Zernike},
};
use skyangle::Conversion;

const PUPIL_MODES: [usize; 2] = [5, 6];

fn main() -> anyhow::Result<()> {
    // field_mesh.plot();
    // println!("{field_mesh}");

    let mut lock = stdout().lock();

    let now = Instant::now();
    let mut results = vec![];
    for id in 1..=7 {
        write!(lock, "{id}").unwrap();
        lock.flush().unwrap();
        results.push(thread::spawn(move || {
            let field_mesh = Mesh::disc(20f64, 2., None);
            let delaunay = Mesh::gmt_segment_with_id(id).unwrap();

            let mut gmt = Gmt::new().unwrap();

            let iter = field_mesh
                .vertex_iter()
                .map(|xy| (xy[0].hypot(xy[1]), xy[1].atan2(xy[0])));

            let field_zernikes = iter
                .map(|(zen, azi)| {
                    let rays_builder = Rays::builder().zenith(zen.from_arcmin()).azimuth(azi);
                    let mut rays: Rays = Rays::from_mesh(&delaunay, rays_builder)?;
                    let modes = rays
                        .trace(&mut gmt)
                        .opds_centered(&delaunay)
                        .unwrap()
                        .as_zernikes(4, &delaunay)
                        .reduce_into(PUPIL_MODES);
                    Ok(FieldZernike::new((zen, azi), modes))
                })
                .collect::<Result<Set<FieldZernike>, CrseoError>>()
                .unwrap();
            field_zernikes.as_zernikes(4, &field_mesh)
        }));
    }
    println!();

    let results: Vec<_> = results.into_iter().map(|res| res.join().unwrap()).collect();
    println!("Elapsed time: {:.3?}", now.elapsed());

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
