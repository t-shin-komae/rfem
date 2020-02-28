use gmsh_v2_parser::*;
use std::fs::File;
use std::path::Path;
use std::error::Error;
use std::io::prelude::*;
use std::io::BufReader;
use crate::triangle::Triangle;
use crate::line::Line;

pub fn process_2d_trianglation_file<P:AsRef<Path>>(filepath:P) -> Result<(Vec<PhysicalName>,Vec<Point2D>,Vec<(Line,Tags)>,Vec<(Triangle,Tags)>),Box<dyn Error>>{
    let f = File::open(filepath)?;
    let mut buf_reader = BufReader::new(f);
    let mut contents = String::new();
    buf_reader.read_to_string(&mut contents)?;
    let lines:Vec<&str> = contents.lines().collect();
    let (_,physical_name,nodes,elements) = parse(&lines)?;
    let points:Vec<Point2D> = nodes.into_iter().map(|node|convert_node_to_point(node)).collect();
    let lines_with_tags = extract_line2d(&points, &elements);
    let triangles_with_tags = extract_triangle2d(&points, &elements);
    Ok((physical_name,points,lines_with_tags,triangles_with_tags))
}

pub type Tags = [usize;3];
type Point2D = [f64;2];

fn convert_node_to_point(node:Node) -> Point2D {
    [node.coord[0],node.coord[1]]
}

fn extract_triangle2d(points:&[Point2D],elements:&[Element]) -> Vec<(Triangle,Tags)>{
    let mut tris = Vec::new();
    for element in elements.iter(){
        match element.element{
            ElementType::Triangle(point_ids_start_1) => {
                let mut point_ids = [0;3];
                point_ids.iter_mut().zip(point_ids_start_1.iter()).for_each(|(start_zero,start_one)|*start_zero=*start_one-1);
                tris.push((Triangle::from_ids(point_ids,points),element.tags))
            },
            _ =>{}
        }
    }
    tris
}
fn extract_line2d(points:&[Point2D],elements:&[Element]) -> Vec<(Line,Tags)>{
    let mut lines = Vec::new();
    for element in elements.iter(){
        match element.element{
            ElementType::Line(point_ids_start_1) => {
                let mut point_ids = [0;2];
                point_ids.iter_mut().zip(point_ids_start_1.iter()).for_each(|(start_zero,start_one)|*start_zero=*start_one-1);
                lines.push((Line::from_ids(point_ids,points),element.tags))
            },
            _ =>{}
        }
    }
    lines
}