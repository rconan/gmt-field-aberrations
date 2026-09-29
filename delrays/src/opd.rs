use std::fmt::Display;

use serde::{Deserialize, Serialize};

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Opd {
    pub(crate) xyz: [f64; 3],
    pub(crate) delta: f64,
}
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Opds(Vec<Opd>);

impl FromIterator<Opd> for Opds {
    fn from_iter<T: IntoIterator<Item = Opd>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

#[derive(Default, Debug, Clone)]
pub struct Stats {
    n_sample: usize,
    min_max: (f64, f64),
    mean: f64,
    var: f64,
}
impl Display for Stats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "sample #: {};min/max: [{:.3e},{:.3e}]; mean: {:.3e}; std: {:.3e}",
            self.n_sample,
            self.min_max.0,
            self.min_max.1,
            self.mean,
            self.var.sqrt()
        )
    }
}
impl From<&Opds> for Stats {
    fn from(opds: &Opds) -> Self {
        let n_sample = opds.0.len();
        let min = opds
            .0
            .iter()
            .min_by(|a, b| a.delta.partial_cmp(&b.delta).unwrap())
            .unwrap()
            .delta;
        let max = opds
            .0
            .iter()
            .max_by(|a, b| a.delta.partial_cmp(&b.delta).unwrap())
            .unwrap()
            .delta;
        let mean = opds.0.iter().map(|opd| opd.delta).sum::<f64>() / n_sample as f64;
        let var = opds
            .0
            .iter()
            .map(|opd| opd.delta - mean)
            .map(|x| x * x)
            .sum::<f64>()
            / n_sample as f64;
        Self {
            n_sample,
            min_max: (min, max),
            mean,
            var,
        }
    }
}
