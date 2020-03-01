use crate::calc::Vector2d;
#[derive(Debug)]
/// Line element type in 2D plane.
/// 
/// This structure stores two *nodes* and their *ids*.
/// The *ids* are *zero based numbering*.
pub struct Line {
    nodes: [[f64; 2]; 2],
    ids: [usize; 2],
}

impl Line{
    /// Construct a new `Line` element from 2 nodes and their correspoinding ids.
    pub fn new(nodes:[[f64;2];2],ids:[usize;2]) -> Self{
        Self{nodes,ids}
    }
    /// Construct a new `Line` element from 2 ids and nodes slice.
    /// This methods take 2 nodes from the nodes slice.
    /// # Example
    /// ```
    /// let nodes = [
    ///     [0.,0.],
    ///     [1.,0.],
    ///     [0.,1.],
    ///     [1.,1.],
    /// ];
    /// // This line has these nodes, [1.,0.],[0.,1.]
    /// let line = Line::from_ids([2,3],&nodes);
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
    /// let line = Line::from_ids([2,6],&nodes);
    /// ```
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