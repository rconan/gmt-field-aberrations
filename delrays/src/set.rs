use std::ops::Deref;

pub struct Set<T>(pub(crate) Vec<T>);

impl<T> Deref for Set<T> {
    type Target=[T];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
