pub mod poisson;
pub mod triangle;

pub trait LinearFem<Elememt> {
    fn create_Ke(elem:&Elememt) -> ndarray::Array2<f32>;
    fn create_fe(elem:&Elememt, f_i: &[f32]) -> ndarray::Array1<f32>;
    fn patch_to_K(elem:&Elememt, Ke: &ndarray::Array2<f32>, K: &mut ndarray::Array2<f32>);
    fn patch_to_fe(elem:&Elememt, fe: &ndarray::Array1<f32>, f: &mut ndarray::Array1<f32>);
}

pub struct ElasticEquation;
