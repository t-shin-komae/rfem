//! <script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>
//! <script type="text/x-mathjax-config">
//!  MathJax.Hub.Config({
//!  tex2jax: {
//!  inlineMath: [["\\(","\\)"] ],
//!  displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
//!  }
//!  });
//! </script>
#[cfg(test)]
extern crate openblas_src;
use crate::triangle::{Triangle, TriangleQuad};
use crate::LinearFemElement;
/// Triangle element that has several convenient method for creating Matrix
pub struct PoissonTriangleElement<'a> {
    /// refrence to a triangle element
    pub triangle: &'a Triangle,
    /// the `f` value at each nodes
    pub f_i: [f64; 3],
}

impl<'a> PoissonTriangleElement<'a> {
    /// Construct a new PoissonTriangleElement.
    /// `f_i` is a value of f at the triangle's nodes.
    pub fn new(tri: &'a Triangle, f_i: &[f64; 3]) -> Self {
        Self {
            triangle: tri,
            f_i: *f_i,
        }
    }
    /// Construct a new PoissonTriangleElement from all `f` values at all nodes.
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
/// Let \\(\Omega_e\\) be a triangle, \\( \phi_i ~~~(i=1,2,3)\\) be a interpolation function, and \\(J_e\\) be a jacobi matrix for the transformation from \\(x,y\\) coordinate to \\(\xi,\eta\\) coordinate(local coordinate).
impl<'a> LinearFemElement for PoissonTriangleElement<'a> {
//! <script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>
//! <script type="text/x-mathjax-config">
//!  MathJax.Hub.Config({
//!  tex2jax: {
//!  inlineMath: [["\\(","\\)"] ],
//!  displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
//!  }
//!  });
//! </script>
    /// This method create the following matrix.
    /// $$
    /// \begin{align}
    /// \mathbf{K^e} &= \int_{\Omega_e}  \nabla \phi_i^e \nabla \phi_j^e \mathrm{d}\Omega\\\\
    ///              &= \int_0^1 \int_0^{1-\eta}
    /// \begin{pmatrix}
    /// \frac{\partial \phi_i^e}{\partial \xi} &
    /// \frac{\partial \phi_i^e}{\partial \eta} 
    /// \end{pmatrix}
    /// J_e^{-1} (J_e^{-1})^T
    /// \begin{pmatrix}
    /// \frac{\partial \phi_j^e}{\partial \xi}\\\\
    /// \frac{\partial \phi_j^e}{\partial \eta} 
    /// \end{pmatrix}
    /// |J_e|
    ///  \mathrm{d}\xi \mathrm{d}\eta
    /// \end{align}
    /// $$
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
    /// This method create the following matrix.
    /// $$
    /// \begin{align}
    /// \mathbf{f^e} &= \mathbf{f_i} \int_{\Omega_e}  \phi_i^e \phi_j^e \mathrm{d}\Omega\\\\
    ///              &= \mathbf{f_i} \int_0^1 \int_0^{1-\eta}
    /// \\phi_i^e \phi_j^e
    /// |J_e|
    ///  \mathrm{d}\xi \mathrm{d}\eta
    /// \end{align}
    /// $$
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
