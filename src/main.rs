extern crate openblas_src;
use ndarray::prelude::*;
use rfem::io::gmsh::process_2d_trianglation_file;
use rfem::poisson::PoissonTriangleElement;
use rfem::triangle::Triangle;
use rfem::LinearFemElement;
use rfem::dirichlet;
use ndarray_linalg::Solve;
use std::fs::File;
use std::io::{self,Read,Write};
fn main() {
    let (physicalnames, points, lines, triangles) =
        process_2d_trianglation_file("triangulation.msh").unwrap();
    let f_at_points = vec![0.; points.len()];
    let mut f_vec = Array1::<f32>::zeros(points.len().f());
    let mut K = Array2::<f32>::zeros((points.len(), points.len()).f());
    let poisson_elements_iter = triangles
        .into_iter()
        .map(|tri| PoissonTriangleElement::from_all_f(tri.0, &f_at_points));
    for po in poisson_elements_iter {
        let ke: Array2<f32> = po.create_Ke();
        let fe: Array1<f32> = po.create_fe();
        po.patch_to_K(&ke, &mut K);
        po.patch_to_fe(&fe, &mut f_vec);
    }
    for line_with_tags in lines.iter(){
        let (line,tags) = line_with_tags;
        if physicalnames[tags[0]-1].name == "ZERO"{
            dirichlet(&mut K, &mut f_vec, line.ids[0], 0.);
        }else if physicalnames[tags[0]-1].name == "HIGH"{
            dirichlet(&mut K, &mut f_vec, line.ids[0], (1.-line.nodes[0][0])*line.nodes[0][0]);
        }else if physicalnames[tags[0]-1].name == "LOW"{
            dirichlet(&mut K, &mut f_vec, line.ids[0], -(1.-line.nodes[0][0])*line.nodes[0][0]);
        }
    }
    println!("start solve");
    let x = K.solve_into(f_vec).unwrap();
    println!("end solve");
    let mut file = File::create("output.dat").unwrap();
    for (x,point) in x.iter().zip(points.iter()){
        writeln!(file,"{}\t{}\t{}",point[0],point[1],x);
    }
}
