#[cfg(test)]
extern crate openblas_src;
use crate::triangle::{Triangle, TriangleQuad};
use crate::LinearFemElement;
/// Triangle element that has several convenient method for creating Matrix
pub struct PoissonTriangleElement<'a> {
    /// refrence to a triangle element
    pub triangle: &'a Triangle,
    /// the value at each nodes
    pub f_i: [f64; 3],
}

impl<'a> PoissonTriangleElement<'a> {
    pub fn new(tri: &'a Triangle, f_i: &[f64; 3]) -> Self {
        Self {
            triangle: tri,
            f_i: *f_i,
        }
    }
    pub fn from_all_f(tri: &'a Triangle, f: &[f64]) -> Self {
        let ids = tri.get_all_ids();
        let f_i = [f[ids[0]], f[ids[1]], f[ids[2]]];
        Self {
            triangle: tri,
            f_i: f_i,
        }
    }
}
use crate::calc::Vector2d;
use ndarray::{Array1, Array2};
impl<'a> LinearFemElement for PoissonTriangleElement<'a> {
    fn create_Ke(&self) -> Array2<f64> {
        let xi_xi = ndarray::arr2(&[[0.5, -0.5, 0.], [-0.5, 0.5, 0.], [0., 0., 0.]]);
        let xi_eta_eta_xi = ndarray::arr2(&[[1., -0.5, -0.5], [-0.5, 0., 0.5], [-0.5, 0.5, 0.]]);
        let eta_eta = ndarray::arr2(&[[0.5, 0., -0.5], [0., 0., 0.], [-0.5, 0., 0.5]]);
        // ((a,b),(c,d))J^-1 (J^-1)^T |J|
        let edge_vec1 = Vector2d::start_end(&self.triangle.get_node(0), &self.triangle.get_node(1));
        let edge_vec2 = Vector2d::start_end(&self.triangle.get_node(0), &self.triangle.get_node(2));
        let J = [&edge_vec1, &edge_vec2];
        let det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        let a = edge_vec2.square_sum() / det;
        let b = - edge_vec1.dot(&edge_vec2) / det;
        let d = edge_vec1.square_sum() / det;
        a * xi_xi + b * xi_eta_eta_xi + d * eta_eta
    }
    fn create_fe(&self) -> Array1<f64> {
        let matrix = ndarray::arr2(&[[2., 1., 1.], [1., 2., 1.], [1., 1., 2.]]) / 24.;
        let edge_vec1 = Vector2d::start_end(&self.triangle.get_node(0), &self.triangle.get_node(1));
        let edge_vec2 = Vector2d::start_end(&self.triangle.get_node(0), &self.triangle.get_node(2));
        let J = [&edge_vec1, &edge_vec2];
        let det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        let f_vec = ndarray::arr1(&self.f_i) * det;
        f_vec.dot(&matrix)
    }
    fn patch_to_K(&self, Ke: &Array2<f64>, K: &mut Array2<f64>) {
        for ((i, j), k_ij) in Ke.indexed_iter() {
            let id_i = self.triangle.get_id(i);
            let id_j = self.triangle.get_id(j);
            K[[id_i, id_j]] += k_ij;
        }
    }
    fn patch_to_fe(&self, fe: &Array1<f64>, f: &mut Array1<f64>) {
        for (i, f_i) in fe.indexed_iter() {
            let id_i = self.triangle.get_id(i);
            f[id_i] += f_i;
        }
    }
}
