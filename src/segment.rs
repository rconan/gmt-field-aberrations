use crseo::{gmt, utilities::MaskFilter, Builder, FromBuilder, Source};
use geotrans::{Segment, SegmentTrait, Transform, M1};
use skyangle::SkyAngle;

use crate::{field::FieldError, Mirror, Projection, Rbm, Txyz, ZernikeBasis};

type Result<T> = std::result::Result<T, FieldError>;

pub fn field_zernike(
    sid: i32,
    z: f32,
    a: f32,
    n_radial_order: u32,
    mirror: &Mirror,
    zeroed: bool,
) -> Result<Projection> {
    // segment ID
    // let sid = 1;
    // segment origin
    let segment = Segment::<M1>::new(sid)?;
    let [x0, y0, _] = Transform::to([0f64; 3], segment);
    log::info!("segment origin: {:.3?}", (x0, y0));

    let mut gmt = gmt!();

    // segment chief rays (@ origin)
    let xyz_chief = {
        let mut src = Source::builder()
            .zenith_azimuth(vec![z], vec![a])
            .rays_coordinates(vec![x0], vec![y0])
            .build()?;
        let rays = &mut src.as_raw_mut_ptr().rays;
        unsafe {
            gmt.m1.trace(rays);
            gmt.m2.trace(rays);
            rays.to_sphere1(-5.830, 2.197173);
        };

        src.rays().coordinates()
    };
    log::info!("segment chief ray at exit pupil: {:.3?}", &xyz_chief);

    // selecting the GMT segment
    gmt.keep(&[sid]);

    let to_meters = |t_xyz: &[Txyz; 3]| t_xyz.into_iter().map(|x| x.as_f64()).collect::<Vec<f64>>();

    let to_radians = |r_xyz: &[SkyAngle<f64>; 3]| {
        r_xyz
            .into_iter()
            .map(|x| x.to_radians())
            .collect::<Vec<f64>>()
    };

    let zernp = if zeroed {
        // reference source (unperturbed GMT)
        let mut src_ref = Source::builder().zenith_azimuth(vec![z], vec![a]).build()?;
        src_ref.through(&mut gmt).xpupil();

        // probing source (perturbed GMT)
        match mirror {
            Mirror::M1(Rbm { t_xyz, r_xyz }) => {
                gmt.m1_segment_state(sid, &to_meters(t_xyz), &to_radians(r_xyz))
            }
            Mirror::M2(Rbm { t_xyz, r_xyz }) => {
                gmt.m2_segment_state(sid, &to_meters(t_xyz), &to_radians(r_xyz))
            }
        }
        let mut src = Source::builder().zenith_azimuth(vec![z], vec![a]).build()?;
        src.through(&mut gmt).xpupil();

        // building the zernike basis
        let xyz = src.rays().coordinates();
        let xy: Vec<_> = (&src_ref.rays().mask(), &src.rays().mask())
            .filter(xyz.chunks(3))
            .map(|c| [c[0] - xyz_chief[0], c[1] - xyz_chief[1]])
            .collect();
        log::info!("segment rays #: {}", xy.len());
        let zern = ZernikeBasis::new(n_radial_order, &xy);

        // opd wrt to reference
        let opd: Vec<_> = src
            .phase()
            .iter()
            .zip(src_ref.phase().iter())
            .map(|(p, p0)| p - p0)
            .collect();
        let opd: Vec<_> = (&src_ref.rays().mask(), &src.rays().mask())
            .filter(opd.iter())
            .cloned()
            .collect();
        // projection onto Zernike
        let mut zernp = Projection::new(zern);

        zernp.project(&*opd);
        zernp
    } else {
        // probing source (perturbed GMT)
        match mirror {
            Mirror::M1(Rbm { t_xyz, r_xyz }) => {
                gmt.m1_segment_state(sid, &to_meters(t_xyz), &to_radians(r_xyz))
            }
            Mirror::M2(Rbm { t_xyz, r_xyz }) => {
                gmt.m2_segment_state(sid, &to_meters(t_xyz), &to_radians(r_xyz))
            }
        }
        let mut src = Source::builder().zenith_azimuth(vec![z], vec![a]).build()?;
        src.through(&mut gmt).xpupil();

        // building the zernike basis
        let xyz = src.rays().coordinates();
        let xy: Vec<_> = src
            .rays()
            .mask()
            .filter(xyz.chunks(3))
            .map(|c| [c[0] - xyz_chief[0], c[1] - xyz_chief[1]])
            .collect();
        log::info!("segment rays #: {}", xy.len());
        let zern = ZernikeBasis::new(n_radial_order, &xy);

        // opd wrt to reference
        let opd = src.phase().iter();
        let opd: Vec<_> = src.rays().mask().filter(opd).cloned().collect();
        // projection onto Zernike
        let mut zernp = Projection::new(zern);
        zernp.project(&*opd);
        zernp
    };
    Ok(zernp)
}
