use std::error::Error;

use gmt_field_aberrations::Probes;
use image::{ImageBuffer, Rgb};
use triangle_rs::Delaunay;

pub fn tripcolor(tri: Delaunay, probes: &Probes, j: usize) -> Result<(), Box<dyn Error>> {
    let mut img = ImageBuffer::new(600, 600);
    for v in tri.triangle_iter() {
        let (p0, v0) = pv(probes, v[0], j);
        let (p1, v1) = pv(probes, v[1], j);
        let (p2, v2) = pv(probes, v[2], j);

        rasterize_triangle(&mut img, p0, p1, p2, v0, v1, v2);
    }
    img.save("test.png")?;
    Ok(())
}
pub fn pv(probes: &Probes, index: usize, j: usize) -> ((f64, f64), f64) {
    let v_max = probes
        .projections
        .iter()
        .map(|p| p.coefficients()[j])
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let v_min = probes
        .projections
        .iter()
        .map(|p| p.coefficients()[j])
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    dbg!((v_min, v_max));

    let (field, projection) = probes.get(index).expect(&format!(
        "Probe #{index} not found (Probe #:{})",
        probes.len()
    ));
    let p = {
        let (sa, ca) = (field.at().azimuth.to_radians() as f64).sin_cos();
        (ca, sa)
    };
    let v = projection.coefficients()[j];
    (p, (v - v_min) / (v_max - v_min))
}

pub fn rasterize_triangle(
    img: &mut ImageBuffer<Rgb<u8>, Vec<u8>>,
    p0: (f64, f64),
    p1: (f64, f64),
    p2: (f64, f64),
    v0: f64,
    v1: f64,
    v2: f64,
) {
    let width = img.width() as f64;
    let height = img.height() as f64;

    // Convert normalized coordinates to pixel coordinates
    let px0 = (0.5 * (1. + p0.0) * width) as i32;
    let py0 = (0.5 * (1. + p0.1) * height) as i32;
    let px1 = (0.5 * (1. + p1.0) * width) as i32;
    let py1 = (0.5 * (1. + p1.1) * height) as i32;
    let px2 = (0.5 * (1. + p2.0) * width) as i32;
    let py2 = (0.5 * (1. + p2.1) * height) as i32;

    // Get bounding box
    let min_x = px0.min(px1).min(px2).max(0);
    let max_x = px0.max(px1).max(px2).min(img.width() as i32 - 1);
    let min_y = py0.min(py1).min(py2).max(0);
    let max_y = py0.max(py1).max(py2).min(img.height() as i32 - 1);

    // Rasterize: check each pixel in bounding box
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            // Compute barycentric coordinates
            if let Some((w0, w1, w2)) = barycentric(
                x as f64, y as f64, px0 as f64, py0 as f64, px1 as f64, py1 as f64, px2 as f64,
                py2 as f64,
            ) {
                // Point is inside triangle, interpolate value
                let value = w0 * v0 + w1 * v1 + w2 * v2;
                let color = value_to_color(value);
                img.put_pixel(x as u32, y as u32, color);
            }
        }
    }
}

fn barycentric(
    px: f64,
    py: f64,
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
) -> Option<(f64, f64, f64)> {
    let v0x = bx - ax;
    let v0y = by - ay;
    let v1x = cx - ax;
    let v1y = cy - ay;
    let v2x = px - ax;
    let v2y = py - ay;

    let den = v0x * v1y - v1x * v0y;
    if den.abs() < 1e-10 {
        return None;
    }

    let v = (v2x * v1y - v1x * v2y) / den;
    let w = (v0x * v2y - v2x * v0y) / den;
    let u = 1.0 - v - w;

    // Check if point is inside triangle
    if u >= 0.0 && v >= 0.0 && w >= 0.0 {
        Some((u, v, w))
    } else {
        None
    }
}

fn value_to_color(value: f64) -> Rgb<u8> {
    // Clamp value to [0, 1]
    // let value = value.max(0.0).min(1.0);
    if value > 1f64 || value < 0f64 {
        panic!("value={value}");
    }

    // Simple blue-to-red colormap
    // let r = (value * 255.0) as u8;
    // let b = ((1.0 - value) * 255.0) as u8;
    let rgb = colorous::CIVIDIS.eval_continuous(value).into_array();
    Rgb(rgb)
}
