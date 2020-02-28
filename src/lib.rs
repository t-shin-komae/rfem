use ndarray::{Array1,Array2};
pub mod poisson;
pub mod triangle;
pub mod line;
mod calc;
pub mod io;

pub trait LinearFemElement<KE=Array2<f64>,FE=Array1<f64>,K=Array2<f64>,F=Array1<f64>> {
    fn create_Ke(&self) -> KE;
    fn create_fe(&self) -> FE;
    fn patch_to_K(&self, Ke: &KE, K: &mut K);
    fn patch_to_fe(&self, fe: &FE, f: &mut F);
}

pub fn dirichlet(K:&mut ndarray::Array2<f64>,f:&mut ndarray::Array1<f64>,nth:usize,val:f64){
    *f -= &(val*(&K.column(nth)));
    f[nth] = val;
    K.column_mut(nth).fill(0.) ;
    K.row_mut(nth).fill(0.);
    K[[nth,nth]] = 1.0;
}