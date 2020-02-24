use ndarray::{Array1, Array2};
use rfem::poisson::PoissonTriangleElement;
use rfem::triangle::*;
use rfem::LinearFemElement;
fn main() {
    let triangles = make_triangles();
    let mut k_mat: Array2<f32> = Array2::zeros((100, 100));
    let mut f_vec: Array1<f32> = Array1::zeros((100,));
    let f = [1.;100];
    let poisson_elements_iter = triangles.into_iter().map(|tri| PoissonTriangleElement::from_all_f(tri,&f));
    for po in poisson_elements_iter {
        let ke:Array2<f32> = po.create_Ke();
        let fe:Array1<f32> = po.create_fe();
        po.patch_to_K(&ke, &mut k_mat);
        po.patch_to_fe(&fe, &mut f_vec);
    }
    //dirichlet(&mut k_mat, &mut f_vec, , );
}

fn make_triangles() -> Vec<Triangle> {
    unimplemented!();
}



fn dirichlet(K:&mut ndarray::Array2<f32>,f:&mut ndarray::Array1<f32>,nth:usize,val:f32){
    *f -= &(val*(&K.column(nth)));
    f[nth] = val;
    K.column_mut(nth).fill(0.) ;
    K.row_mut(nth).fill(0.);
    K[[nth,nth]] = 1.0;
}
