pub struct Triangle {
    nodes: [[f32; 2]; 3],
    ids: [usize; 3],
}

impl Triangle {
    pub fn new(nodes: [[f32; 2]; 3], ids: [usize; 3]) -> Self {
        unimplemented!();
    }
    pub fn from_ids(ids: [usize; 3], points: &[[f32; 2]]) -> Self {
        unimplemented!();
    }
    pub fn physical_quantity(&self, x: f32, y: f32, physics: [f32; 3]) -> f32 {
        unimplemented!();
    }
    pub fn physical_from_ids(&self, x: f32, y: f32, physics: &[f32]) -> f32 {
        unimplemented!();
    }
    pub fn local_coord(&self, x: f32, y: f32) -> [f32; 2] {
        unimplemented!();
    }
    pub fn is_include(&self, point: [f32; 2]) -> bool {
        unimplemented!();
    }
    fn phi(xi: f32, eta: f32, physics: [f32; 3]) -> f32 {
        // physics[0]*phi1(xi,eta)+
        // physics[1]*phi2(xi,eta)+
        // physics[2]*phi3(xi,eta)
        unimplemented!();
    }
    #[inline]
    fn phi1(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    #[inline]
    fn phi2(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    #[inline]
    fn phi3(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
}

pub struct TriangleQuad {
    nodes: [[f32; 2]; 6],
    ids: [usize; 6],
}
impl TriangleQuad {
    pub fn new(nodes: [[f32; 2]; 6], ids: [usize; 6]) -> Self {
        unimplemented!();
    }
    pub fn from_ids(ids: [usize; 6], points: &[[f32; 2]]) -> Self {
        unimplemented!();
    }
    pub fn local_coord(&self, x: f32, y: f32) -> [f32; 2] {
        unimplemented!();
    }
    pub fn physical_quantity(&self, x: f32, y: f32, physics: [f32; 6]) -> f32 {
        unimplemented!();
    }
    pub fn physical_from_ids(&self, x: f32, y: f32, physics: &[f32]) -> f32 {
        unimplemented!();
    }
    pub fn is_include(&self, point: [f32; 2]) -> bool {
        unimplemented!();
    }
    fn phi(xi: f32, eta: f32, physics: [f32; 6]) -> f32 {
        unimplemented!();
    }
    fn phi1(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    fn phi2(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    fn phi3(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    fn phi4(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    fn phi5(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
    fn phi6(xi: f32, eta: f32) -> f32 {
        unimplemented!();
    }
}
