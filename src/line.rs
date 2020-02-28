#[derive(Debug)]
pub struct Line {
    pub nodes: [[f64; 2]; 2],
    pub ids: [usize; 2],
}

impl Line{
    pub fn new(nodes:[[f64;2];2],ids:[usize;2]) -> Self{
        Self{nodes,ids}
    }
    pub fn from_ids(ids:[usize;2],points:&[[f64;2]])-> Self{
        Self{
            nodes:[points[ids[0]],points[ids[1]]],
            ids:ids
        }
    }
}