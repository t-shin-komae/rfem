//! <script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS_CHTML"></script>
//! <script type="text/x-mathjax-config">
//!  MathJax.Hub.Config({
//!  tex2jax: {
//!  inlineMath: [["\\(","\\)"] ],
//!  displayMath: [ ['$$','$$'], ["\\[","\\]"] ]
//!  }
//!  });
//! </script>
/// Triangle element type in 2D plane.
/// 
/// This structure stores three *nodes*(vertexes) and their *ids*.
/// The *ids* are *zero based numbering*.
pub struct Triangle {
    nodes: [[f64; 2]; 3],
    ids: [usize; 3],
}

impl Triangle {
    /// Construct a new `Triangle` element from 3 nodes and their correspoinding ids.
    pub fn new(nodes: [[f64; 2]; 3], ids: [usize; 3]) -> Self {
        Self{nodes,ids}
    }
    /// Construct a new `Triangle` element from 3 ids and nodes slice.
    /// This methods take 3 nodes from the nodes slice.
    /// # Example
    /// ```
    /// let nodes = [
    ///     [0.,0.],
    ///     [1.,0.],
    ///     [0.,1.],
    ///     [1.,1.],
    ///     [2.,1.],
    ///     [2.,2.],
    /// ];
    /// // This triangle has these nodes, [1.,0.],[0.,1.],[1.,0.]
    /// let triangle = Triangle::from_ids([2,3,1],&nodes);
    /// ```
    /// # Panic
    /// if nodes.len() < max(ids) then panic
    /// ```
    /// let nodes = [
    ///     [0.,0.],
    ///     [1.,0.],
    ///     [0.,1.],
    ///     [1.,1.],
    ///     [2.,1.],
    ///     [2.,2.],
    /// ];
    /// // panic because the nodes has 6 element but designate 6th index, so internally, array out of bounds occurs
    /// let triangle = Triangle::from_ids([2,6,1],&nodes);
    /// ```
    pub fn from_ids(ids: [usize; 3], nodes: &[[f64; 2]]) -> Self {
        Self{
            nodes:[nodes[ids[0]],nodes[ids[1]],nodes[ids[2]]],
            ids:ids
        }
    }
    /// Get node in index th.
    #[inline]
    pub fn get_node(&self,index:usize) -> [f64;2]{
        self.nodes[index]
    }
    /// Get id in index th.
    #[inline]
    pub fn get_id(&self,index:usize) -> usize{
        self.ids[index]
    }
    /// Get all ids in triangle.
    #[inline]
    pub fn get_all_ids(&self) -> [usize;3]{
        self.ids
    }
    /// Calculate physical quantity in triangle.
    /// `x` and `y` should be the point in triangle.
    /// `physics` is a physical quantity at each nodes.
    pub fn physical_quantity(&self, x: f64, y: f64, physics: &[f64; 3]) -> f64 {
        let [xi,eta] = self.local_coord(x,y);
        Self::phi(xi,eta,physics)
    }
    /// Calculate the partial derivative of physical quantity with respect to x in triangle.
    /// `x` and `y` should be the point in triangle.
    /// `physics` is a physical quantity at each nodes.
    pub fn physical_quantity_diff_x(&self, x: f64, y: f64, physics: &[f64; 3]) -> f64 {
        // The differencial is constant in a one order triangle element
        let [xi,eta] = self.local_coord(x,y);
        let [xi_delta,eta_delta] = self.local_coord(x+0.1,y);
        (Self::phi(xi_delta,eta_delta,physics)-Self::phi(xi,eta,physics))/0.1
    }
    /// Calculate the partial derivative of physical quantity with respect to y in triangle.
    /// `x` and `y` should be the point in triangle.
    /// `physics` is a physical quantity at each nodes.
    pub fn physical_quantity_diff_y(&self, x: f64, y: f64, physics: &[f64; 3]) -> f64 {
        let [xi,eta] = self.local_coord(x,y);
        let [xi_delta,eta_delta] = self.local_coord(x,y+0.1);
        (Self::phi(xi_delta,eta_delta,physics)-Self::phi(xi,eta,physics))/0.1
    }
    /// Convenient wrapper for `physical_quantity_diff_x`. see also `from_ids`
    pub fn physical_diff_x_from_ids(&self, x: f64, y: f64, physics: &[f64]) -> f64 {
        self.physical_quantity_diff_x(x,y,&[physics[self.ids[0]],physics[self.ids[1]],physics[self.ids[2]]])
    }
    /// Convenient wrapper for `physical_quantity_diff_y`. see also `from_ids`
    pub fn physical_diff_y_from_ids(&self, x: f64, y: f64, physics: &[f64]) -> f64 {
        self.physical_quantity_diff_y(x,y,&[physics[self.ids[0]],physics[self.ids[1]],physics[self.ids[2]]])
    }
    /// Convenient wrapper for `physical_quantity`. see also `from_ids`
    pub fn physical_from_ids(&self, x: f64, y: f64, physics: &[f64]) -> f64 {
        let [xi,eta] = self.local_coord(x,y);
        Self::phi(xi,eta,&[physics[self.ids[0]],physics[self.ids[1]],physics[self.ids[2]]])
    }
    /// Calculate the local coordinate of the designated point `(x,y)`
    pub fn local_coord(&self, x: f64, y: f64) -> [f64; 2] {
        // x_e = xi*(p2-p1) + eta*(p3-p1)
        // (xi,eta) = J^-1 x_e
        use crate::calc::Vector2d;
        let nodes = self.nodes;
        let v1 = Vector2d::start_end(&nodes[0],&nodes[1]); // vector (p2-p1)
        let v2 = Vector2d::start_end(&nodes[0],&nodes[2]); // vector (p3-p1)
        let x_e = Vector2d::start_end(&nodes[0],&[x,y]);   // vector (x-p1)
        let det = v1[0]*v2[1]-v1[1]*v2[0];
        let xi = 1./det * Vector2d([v2[1],-v2[0]]).dot(&x_e);
        let eta = 1./det * Vector2d([-v1[1],v1[0]]).dot(&x_e);
        [xi,eta]
    }
    /// Check if or if not the `point` is included in this triangle.
    pub fn is_include(&self, point: &[f64; 2]) -> bool {
        let [xi,eta] = self.local_coord(point[0],point[1]);
        0.<=xi && 0.<=eta && xi+eta <= 1.00001
    }
    fn phi(xi: f64, eta: f64, physics: &[f64; 3]) -> f64 {
        physics[0]*Self::phi1(xi,eta)+
        physics[1]*Self::phi2(xi,eta)+
        physics[2]*Self::phi3(xi,eta)
    }
    #[inline]
    fn phi1(xi: f64, eta: f64) -> f64 {
        1.-xi-eta
    }
    #[inline]
    fn phi2(xi: f64, _eta: f64) -> f64 {
        xi
    }
    #[inline]
    fn phi3(_xi: f64, eta: f64) -> f64 {
        eta
    }
}

pub struct TriangleQuad {
    nodes: [[f64; 2]; 6],
    ids: [usize; 6],
}
impl TriangleQuad {
    pub fn new(nodes: [[f64; 2]; 6], ids: [usize; 6]) -> Self {
        unimplemented!();
    }
    pub fn from_ids(ids: [usize; 6], points: &[[f64; 2]]) -> Self {
        unimplemented!();
    }
    pub fn local_coord(&self, x: f64, y: f64) -> [f64; 2] {
        unimplemented!();
    }
    pub fn physical_quantity(&self, x: f64, y: f64, physics: [f64; 6]) -> f64 {
        unimplemented!();
    }
    pub fn physical_from_ids(&self, x: f64, y: f64, physics: &[f64]) -> f64 {
        unimplemented!();
    }
    pub fn is_include(&self, point: [f64; 2]) -> bool {
        unimplemented!();
    }
    fn phi(xi: f64, eta: f64, physics: [f64; 6]) -> f64 {
        unimplemented!();
    }
    fn phi1(xi: f64, eta: f64) -> f64 {
        unimplemented!();
    }
    fn phi2(xi: f64, eta: f64) -> f64 {
        unimplemented!();
    }
    fn phi3(xi: f64, eta: f64) -> f64 {
        unimplemented!();
    }
    fn phi4(xi: f64, eta: f64) -> f64 {
        unimplemented!();
    }
    fn phi5(xi: f64, eta: f64) -> f64 {
        unimplemented!();
    }
    fn phi6(xi: f64, eta: f64) -> f64 {
        unimplemented!();
    }
}


#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn local_coor_and_include_test(){
        let tri = Triangle::new([[0.,0.],[0.2,1.],[1.,0.2]],[0,1,2]);
        let [xi,eta] = tri.local_coord(0.6,0.6);
        println!("{:?}",[xi,eta]);
        assert!((xi-0.5).abs()<core::f64::EPSILON);
        assert!((eta-0.5).abs()<core::f64::EPSILON);
        println!("{:?}",tri.local_coord(0.36,0.36));
        println!("{:?}",tri.local_coord(0.0,-0.3));
        assert!(tri.is_include(&[0.3,0.6]));
        assert!(!tri.is_include(&[0.0,-0.3]));
    }
    #[test]
    fn physical_quantity(){
        let tri = Triangle::new([[0.,0.],[0.2,1.],[1.,0.2]],[0,1,2]);
        let physical_quantity = [1.,2.,3.];
        let physical_quantity_at_point = tri.physical_quantity(0.1,0.5,&physical_quantity);
        assert!((physical_quantity_at_point-1.5)<core::f64::EPSILON);
        println!("{:?}",tri.physical_quantity(0.6,0.6,&physical_quantity));
        assert!((tri.physical_quantity(0.6,0.6,&physical_quantity)-2.5)<2.5*core::f64::EPSILON);
    }
}