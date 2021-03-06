extern crate openblas_src;
use ndarray::prelude::*;
use rfem::io::gmsh::process_2d_trianglation_file;
use rfem::poisson::PoissonTriangleElement;
use rfem::LinearFemElement;
use rfem::dirichlet;
use ndarray_linalg::Solve;
use std::fs::File;
use std::io::Write;
fn main() {
    // Δu = -f
    let (physicalnames, points, lines, triangles) =
        process_2d_trianglation_file("test.msh").unwrap();
    let f_at_points = vec![1.; points.len()];
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
    for line_with_tags in lines.iter(){
        let (line,tags) = line_with_tags;
        let [x,y] = line.get_node(0);
        dirichlet(&mut K, &mut f_vec, line.get_id(0), 0.);
    }
    println!("start solve");
    let u = K.solve_into(f_vec).unwrap();
    println!("end solve");
    let mut file = File::create("output.dat").unwrap();
    let mut file_node = File::create("output_node.dat").unwrap();
    for x in (0..101).map(|i|i as f64/100.){
        for y in (0..101).map(|i| i as f64 / 100.){
            let mut u_at_point = 0.;
            for (tri,_) in triangles.iter(){
                if tri.is_include(&[x,y]){
                    u_at_point = tri.physical_from_ids(x,y,u.as_slice().unwrap());
                    break;
                }
            }
            let exact = exact_solution(x, y);
            writeln!(file,"{}\t{}\t{}\t{}\t{}",x,y,u_at_point,exact,u_at_point-exact).unwrap();
        }
        writeln!(file,"").unwrap();
    }
    for (x,point) in u.iter().zip(points.iter()){
        writeln!(file_node,"{}\t{}\t{}\t{}\t{}",point[0],point[1],x,exact_solution(point[0], point[1]),x-exact_solution(point[0], point[1])).unwrap();
    }
}

use std::f64::consts::PI;
fn exact_solution(x: f64, y: f64) -> f64 {
    // http://www.kz.tsukuba.ac.jp/~abe/ohp-cfd/cfd-2.pdf

    x*x-0.5*y*y
}