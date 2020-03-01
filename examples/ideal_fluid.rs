extern crate openblas_src;
use ndarray::prelude::*;
use rfem::io::gmsh::process_2d_trianglation_file;
use rfem::poisson::PoissonTriangleElement;
use rfem::poisson::PoissonLineBoundaryElement;
use rfem::LinearFemElement;
use rfem::LinearFemBoundaryElement;
use rfem::dirichlet;
use ndarray_linalg::Solve;
use std::fs::File;
use std::io::Write;
fn main() {
    let (physicalnames, points, lines, triangles) =
        process_2d_trianglation_file("triangulation_for_fluid.msh").unwrap();
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
        let [x_1,y_1] = line.get_node(0);
        let [x_2,y_2] = line.get_node(1);
        if physicalnames[tags[0] - 1].name == "WALL" {
            let p_line = PoissonLineBoundaryElement::new(&line,&[0.,0.]);
            let qe = p_line.create_qe();
            p_line.patch_to_qe(&qe, &mut f_vec);
        }
        if physicalnames[tags[0] - 1].name == "IN" {
            let p_line = PoissonLineBoundaryElement::new(&line,&[-1.,-1.]);
            let qe = p_line.create_qe();
            p_line.patch_to_qe(&qe, &mut f_vec);
        }
    }
    for line_with_tags in lines.iter(){
        let (line,tags) = line_with_tags;
        if physicalnames[tags[0]-1].name == "OUT"{
            dirichlet(&mut K, &mut f_vec, line.get_id(0), 0.);
            dirichlet(&mut K, &mut f_vec, line.get_id(1), 0.);
        }
    }
    println!("start solve");
    let u = K.solve_into(f_vec).unwrap();
    println!("end solve");
    let mut file = File::create("output.dat").unwrap();
    for x in (0..21).map(|i|i as f64/20.){
        for y in (0..21).map(|i| i as f64 / 20.){
            if (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5) < 0.04001{
                continue;
            }
            let mut diff_x = 0.;
            let mut diff_y = 0.;
            for (tri,_) in triangles.iter(){
                if tri.is_include(&[x,y]){
                    diff_x = tri.physical_diff_x_from_ids(x,y,u.as_slice().unwrap());
                    diff_y = tri.physical_diff_y_from_ids(x,y,u.as_slice().unwrap());
                    break;
                }
            }
            writeln!(file,"{}\t{}\t{}\t{}",x,y,diff_x,diff_y).unwrap();
        }
    }
    // for (x,point) in u.iter().zip(points.iter()){
    //     writeln!(file_node,"{}\t{}\t{}",point[0],point[1],x).unwrap();
    // }
}
