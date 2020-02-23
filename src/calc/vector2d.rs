pub(crate) struct Vector2d(pub(crate)[f32;2]);

impl Vector2d{
    pub(crate) fn start_end(start:&[f32;2],end:&[f32;2])->Self{
        Self([end[0]-start[0],end[1]-start[1]])
    }

    pub(crate) fn dot(&self,other:&Self)->f32{
        self[0]*other[0]+self[1]*other[1]
    }
    pub(crate) fn square_sum(&self)->f32{
        self[0]*self[0]+self[1]*self[1]
    }
}

use core::ops::*;
impl Index<usize> for Vector2d{
    type Output = f32;
    fn index(&self,index:usize) -> &Self::Output{
        &self.0[index]
    }
}