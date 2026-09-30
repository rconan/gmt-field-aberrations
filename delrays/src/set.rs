//! Generic collection

use std::ops::{Deref, DerefMut};

use serde::{Deserialize, Serialize};

/// Wrapper around a [Vec]
#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct Set<T>(pub(crate) Vec<T>);

impl<T> Deref for Set<T> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> DerefMut for Set<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl<T> FromIterator<T> for Set<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<T> IntoIterator for Set<T> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}
