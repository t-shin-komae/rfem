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
    for line_with_tags in lines.iter(){
        let (line,tags) = line_with_tags;
        let [x,y] = line.nodes[0];
        if physicalnames[tags[0]-1].name == "ZERO"{
            dirichlet(&mut K, &mut f_vec, line.ids[0], 0.);
        }else if physicalnames[tags[0]-1].name == "HIGH"{
            dirichlet(&mut K, &mut f_vec, line.ids[0], (1.-x)*x);
        }else if physicalnames[tags[0]-1].name == "LOW"{
            dirichlet(&mut K, &mut f_vec, line.ids[0], -(1.-x)*x);
        }
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
    let mut sol = 0.;
    for m in (1..100).step_by(2) {
        let m = m as f64;
        sol += -2.0 * f64::sinh(m * PI * (y - 0.5)) * f64::sin(m * PI * x)
            / (m.powf(3.) * f64::sinh(m * PI / 2.));
    }
    4. * sol / (PI.powf(3.))
}
