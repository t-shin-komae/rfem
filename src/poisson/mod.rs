extern crate openblas_src;
use crate::triangle::{Triangle, TriangleQuad};
use crate::LinearFemElement;
use crate::LinearFemBoundaryElement;
pub struct PoissonTriangleElement<'a> {
    pub triangle: &'a Triangle,
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

use crate::line::Line;
pub struct PoissonLineBoundaryElement<'a>{
    pub line: &'a Line,
    pub q_i: [f64; 2],
}
impl<'a> PoissonLineBoundaryElement<'a> {
    pub fn new(line: &'a Line, q_i: &[f64; 2]) -> Self {
        Self {
            line:line,q_i:*q_i
        }
    }
}

impl<'a> LinearFemBoundaryElement for PoissonLineBoundaryElement<'a>{
    fn create_qe(&self)-> Array1<f64>{
        let psi_i_psi_j = ndarray::arr2(&[[1./3.,1./6.],[1./6.,1./3.]]);
        ndarray::arr1(&self.q_i).dot(&psi_i_psi_j)*self.line.jacobian()
    }
    fn patch_to_qe(&self,qe:&Array1<f64>, q:&mut Array1<f64>){
        for (i, q_i) in qe.indexed_iter() {
            let id_i = self.line.get_id(i);
            q[id_i] += q_i;
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use crate::dirichlet;
    use crate::io::gmsh::process_2d_trianglation_file;
    use ndarray::prelude::*;
    use ndarray::{Array1, Array2};
    use ndarray_linalg::Solve;
    #[test]
    fn test_laplace() {
        fn exact_solution(x: f64, y: f64) -> f64 {
            x*x - y*y
        }
        let (physicalnames, points, lines, triangles) =
            process_2d_trianglation_file("triangulation.msh").unwrap();
        let f_at_points = vec![0.; points.len()];
        let mut f_vec = Array1::<f64>::zeros(points.len().f());
        let mut K = Array2::<f64>::zeros((points.len(), points.len()).f());
        let poisson_elements_iter = triangles
            .iter()
            .map(|tri| PoissonTriangleElement::from_all_f(&tri.0, &f_at_points));
        for po in poisson_elements_iter {
            let ke: Array2<f64> = po.create_Ke();
            let fe: Array1<f64> = po.create_fe();
            po.patch_to_K(&ke, &mut K);
            po.patch_to_fe(&fe, &mut f_vec);
        }
        for line_with_tags in lines.iter() {
            let (line, tags) = line_with_tags;
            let [x,y] = line.nodes[0];
            if physicalnames[tags[0] - 1].name == "ZERO" {
                dirichlet(&mut K, &mut f_vec, line.ids[0], x*x-y*y);
            } else if physicalnames[tags[0] - 1].name == "HIGH" {
                dirichlet(
                    &mut K,
                    &mut f_vec,
                    line.ids[0],
                    x*x-y*y
                );
            } else if physicalnames[tags[0] - 1].name == "LOW" {
                dirichlet(
                    &mut K,
                    &mut f_vec,
                    line.ids[0],
                    x*x-y*y
                );
            }
        }
        let u = K.solve_into(f_vec).unwrap();
        for (u_at_p,point) in u.iter().zip(points.iter()){
            assert!((u_at_p - exact_solution(point[0], point[1])).abs() < 1e-4 );
        }
    }


    #[test]
    fn test_poisson() {
        fn exact_solution(x: f64, y: f64) -> f64 {
            x*x - 0.5*y*y
        }
        let (_, points, lines, triangles) =
            process_2d_trianglation_file("triangulation.msh").unwrap();
        let f_at_points = vec![-1.; points.len()];
        let mut f_vec = Array1::<f64>::zeros(points.len().f());
        let mut K = Array2::<f64>::zeros((points.len(), points.len()).f());
        let poisson_elements_iter = triangles
            .iter()
            .map(|tri| PoissonTriangleElement::from_all_f(&tri.0, &f_at_points));
        for po in poisson_elements_iter {
            let ke: Array2<f64> = po.create_Ke();
            let fe: Array1<f64> = po.create_fe();
            po.patch_to_K(&ke, &mut K);
            po.patch_to_fe(&fe, &mut f_vec);
        }
        for line_with_tags in lines.iter() {
            let (line, _) = line_with_tags;
            let [x,y] = line.nodes[0];
            dirichlet(&mut K, &mut f_vec, line.ids[0], exact_solution(x, y));
        }
        let u = K.solve_into(f_vec).unwrap();
        for (u_at_p,point) in u.iter().zip(points.iter()){
            assert!((u_at_p - exact_solution(point[0], point[1])).abs() < 1e-4 );
        }
    }
    #[test]
    fn test_poisson_with_boundary_condition() {
        fn exact_solution(x: f64, y: f64) -> f64 {
            x*x - 0.5*y*y
        }
        let (physicalnames, points, lines, triangles) =
            process_2d_trianglation_file("triangulation.msh").unwrap();
        let f_at_points = vec![-1.; points.len()];
        let mut f_vec = Array1::<f64>::zeros(points.len().f());
        let mut K = Array2::<f64>::zeros((points.len(), points.len()).f());
        let poisson_elements_iter = triangles
            .iter()
            .map(|tri| PoissonTriangleElement::from_all_f(&tri.0, &f_at_points));
        for po in poisson_elements_iter {
            let ke: Array2<f64> = po.create_Ke();
            let fe: Array1<f64> = po.create_fe();
            po.patch_to_K(&ke, &mut K);
            po.patch_to_fe(&fe, &mut f_vec);
        }
        for line_with_tags in lines.iter() {
            let (line, tags) = line_with_tags;
            let [x_1,y_1] = line.nodes[0];
            let [x_2,y_2] = line.nodes[1];
            if physicalnames[tags[0] - 1].name == "LOW" {
                let p_line = PoissonLineBoundaryElement::new(&line,&[-1.,-1.]);
                let qe = p_line.create_qe();
                p_line.patch_to_qe(&qe, &mut f_vec);
            }
        }
        for line_with_tags in lines.iter() {
            let (line, tags) = line_with_tags;
            let [x_1,y_1] = line.nodes[0];
            let [x_2,y_2] = line.nodes[1];
            if physicalnames[tags[0] - 1].name != "LOW" {
                dirichlet(&mut K, &mut f_vec, line.ids[0], exact_solution(x_1, y_1));
                dirichlet(&mut K, &mut f_vec, line.ids[1], exact_solution(x_2, y_2));
            }
        }
        let u = K.solve_into(f_vec).unwrap();
        let mut file_node = std::fs::File::create("output_node.dat").unwrap();
        use std::io::Write;
        for (x,point) in u.iter().zip(points.iter()){
            writeln!(file_node,"{}\t{}\t{}\t{}\t{}",point[0],point[1],x,exact_solution(point[0], point[1]),x-exact_solution(point[0], point[1])).unwrap();
        }
        for (u_at_p,point) in u.iter().zip(points.iter()){
            assert!((u_at_p - exact_solution(point[0], point[1])).abs() < 1e-4 );
        }
    }
}
