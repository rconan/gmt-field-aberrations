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
