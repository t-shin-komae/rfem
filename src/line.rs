use crate::calc::Vector2d;
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
    pub fn get_all_ids(&self) -> [usize;2] {
        self.ids
    }
    pub fn get_node(&self,index:usize) -> [f64;2]{
        self.nodes[index]
    }
    pub fn get_id(&self,index:usize) -> usize{
        self.ids[index]
    }

    pub fn jacobian(&self) -> f64{
        Vector2d::start_end(&self.nodes[0],&self.nodes[1]).len()
    }
}