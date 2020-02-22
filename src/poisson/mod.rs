use crate::LinearFem;
use crate::triangle::{Triangle,TriangleQuad};
pub struct Poisson;

impl LinearFem<Triangle> for Poisson{
    fn create_Ke(elem:&Triangle) -> ndarray::Array2<f32>{
        unimplemented!();
    }
    fn create_fe(elem:&Triangle, f_i: &[f32]) -> ndarray::Array1<f32>{
        unimplemented!();
    }
    fn patch_to_K(elem:&Triangle, Ke: &ndarray::Array2<f32>, K: &mut ndarray::Array2<f32>){
        unimplemented!();
    }
    fn patch_to_fe(elem:&Triangle, fe: &ndarray::Array1<f32>, f: &mut ndarray::Array1<f32>){
        unimplemented!();
    }
}

impl LinearFem<TriangleQuad> for Poisson{
    fn create_Ke(elem:&TriangleQuad) -> ndarray::Array2<f32>{
        unimplemented!();
    }
    fn create_fe(elem:&TriangleQuad, f_i: &[f32]) -> ndarray::Array1<f32>{
        unimplemented!();
    }
    fn patch_to_K(elem:&TriangleQuad, Ke: &ndarray::Array2<f32>, K: &mut ndarray::Array2<f32>){
        unimplemented!();
    }
    fn patch_to_fe(elem:&TriangleQuad, fe: &ndarray::Array1<f32>, f: &mut ndarray::Array1<f32>){
        unimplemented!();
    }
}
