use ndarray::{Array1,Array2};
pub mod poisson;
pub mod triangle;
mod calc;

pub trait LinearFemElement<KE=Array2<f32>,FE=Array1<f32>,K=Array2<f32>,F=Array1<f32>> {
    fn create_Ke(&self) -> KE;
    fn create_fe(&self) -> FE;
    fn patch_to_K(&self, Ke: &KE, K: &mut K);
    fn patch_to_fe(&self, fe: &FE, f: &mut F);
}

