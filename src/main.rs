use ndarray::{Array1, Array2};
use rfem::poisson::Poisson;
use rfem::triangle::*;
use rfem::LinearFem;
fn main() {
    let triangles = make_triangles();
    let mut k_mat: Array2<f32> = Array2::zeros((100, 100));
    let mut f_vec: Array1<f32> = Array1::zeros((100,));
    for tri in triangles.iter() {
        let ke = Poisson::create_Ke(tri);
        let fe = Poisson::create_fe(tri, &[1., 1., 1.]);
        Poisson::patch_to_K(tri, &ke, &mut k_mat);
        Poisson::patch_to_fe(tri, &fe, &mut f_vec);
    }
}

fn make_triangles() -> Vec<Triangle> {
    unimplemented!();
}
