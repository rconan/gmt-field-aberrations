use triangle_rs::{Builder, Delaunay};

pub trait Mesh {
    fn disc(diameter: f64, perimeter_pitch: f64) -> Delaunay {
        let mut builder = Builder::new();
        let triangle_area = perimeter_pitch.powi(2) * 3f64.sqrt() / 4f64;
        let n_rim = (std::f64::consts::PI * diameter / perimeter_pitch).round();
        let outer_rim: Vec<_> = (0..n_rim as usize)
            .flat_map(|i| {
                let o = 2. * std::f64::consts::PI * i as f64 / n_rim;
                let (s, c) = o.sin_cos();
                let radius = 0.5 * diameter;
                vec![radius * c, radius * s]
            })
            .collect();
        builder.add_polygon(&outer_rim);
        builder
            .set_switches(&format!("Qpqa{}", triangle_area))
            .build()
    }
}

impl Mesh for Delaunay {}
