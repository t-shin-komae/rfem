//! <script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>
//! <script type="text/x-mathjax-config">
//!  MathJax.Hub.Config({
//!  tex2jax: {
//!  inlineMath: [["\\(","\\)"] ],
//!  displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
//!  }
//!  });
//! </script>
//! This module has elements for poisson equation as below.
//! $$\Delta u = -f(x,y)$$
pub mod triangle;
pub mod line;
pub use triangle::PoissonTriangleElement;
pub use line::PoissonLineBoundaryElement;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dirichlet;
    use crate::io::gmsh::process_2d_trianglation_file;
    use ndarray::prelude::*;
    use ndarray::{Array1, Array2};
    use ndarray_linalg::Solve;
    use crate::LinearFemElement;
    use crate::LinearFemBoundaryElement;
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
            let [x_1,y_1] = line.get_node(0);
            let [x_2,y_2] = line.get_node(1);
            if physicalnames[tags[0] - 1].name == "LOW" {
                let p_line = PoissonLineBoundaryElement::new(&line,&[-1.,-1.]);
                let qe = p_line.create_qe();
                p_line.patch_to_qe(&qe, &mut f_vec);
            }
        }
        for line_with_tags in lines.iter() {
            let (line, tags) = line_with_tags;
            let [x_1,y_1] = line.get_node(0);
            let [x_2,y_2] = line.get_node(1);
            if physicalnames[tags[0] - 1].name != "LOW" {
                dirichlet(&mut K, &mut f_vec, line.get_id(0), exact_solution(x_1, y_1));
                dirichlet(&mut K, &mut f_vec, line.get_id(1), exact_solution(x_2, y_2));
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