use crate::LinearFemBoundaryElement;
use crate::line::Line;
use crate::calc::Vector2d;
use ndarray::{Array1, Array2};
pub struct PoissonLineBoundaryElement<'a>{
    pub line: &'a Line,
    pub q_i: [f64; 2],
}
impl<'a> PoissonLineBoundaryElement<'a> {
    pub fn new(line: &'a Line, q_i: &[f64; 2]) -> Self {
        Self {
            line:line,q_i:*q_i
        }
    }
}

impl<'a> LinearFemBoundaryElement for PoissonLineBoundaryElement<'a>{
    fn create_qe(&self)-> Array1<f64>{
        let psi_i_psi_j = ndarray::arr2(&[[1./3.,1./6.],[1./6.,1./3.]]);
        ndarray::arr1(&self.q_i).dot(&psi_i_psi_j)*self.line.jacobian()
    }
    fn patch_to_qe(&self,qe:&Array1<f64>, q:&mut Array1<f64>){
        for (i, q_i) in qe.indexed_iter() {
            let id_i = self.line.get_id(i);
            q[id_i] += q_i;
        }
    }
}