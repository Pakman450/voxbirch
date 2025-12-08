use crate::voxel::VoxelGrid;

// Itani similarity index for binary voxel grids
// NOTE: assumes all grids have the same dimensions and 
// number of bits 
pub fn itani_bin(l_grids: Vec<VoxelGrid>) -> f32 {

    // num of bits
    let m = l_grids[0].data.len();

    // num of grids
    let n = l_grids.len();

    // vector to hold column wise sums
    let mut kq: Vec<f32> = vec![0.0; m]; 

    // compute columwise sums
    for grid in l_grids.iter() {
        for i in 0..m {
            kq[i] += grid.data[i] as f32;
        }
    }

    // compute itani index
    let mut numer = 0.0;
    let mut denom = 0.0;

    for i in 0..m {
        numer += kq[i] * (kq[i] - 1.0) / 2.0;
    }

    for i in 0..m {
        denom += (kq[i] * (kq[i] - 1.0) / 2.0)
            + (kq[i] * (n as f32 - kq[i]));
    }

    let itani = numer / denom;

    itani
}

// Itani similarity index for real numbers
// NOTE: assumes all grids have the same dimensions and 
// number of bits 

pub fn itani_real(l_grids: Vec<VoxelGrid>) -> f32 {


    // num of bits
    let m = l_grids[0].data.len();

    // num of grids
    let n = l_grids.len();

    // vector to hold column wise sums
    let mut xq: Vec<f32> = vec![0.0; m]; 

    // vector to hold column wise sums
    let mut xq_sq: Vec<f32> = vec![0.0; m]; 

    for grid in l_grids.iter() {
        for i in 0..m {
            // convert u8 to f32
            xq[i] += grid.data[i] as f32;
            xq_sq[i] += (grid.data[i] as f32) * (grid.data[i] as f32);
        }
    }

    // compute itani index
    let mut numer = 0.0;
    let mut denom1 = 0.0;
    let mut denom2 = 0.0;
    for i in 0..m {
        numer += (xq[i] * xq[i]) - xq_sq[i];
    }   

    numer = numer / 2.0;

    for i in 0..m {
        denom1 += xq_sq[i];
    }   

    denom1 = (n-1) as f32 * denom1;

    for i in 0..m {
        denom2 += (xq[i] * xq[i]) - xq_sq[i];
    }   

    denom2 = denom2 / 2.0;

    let itani = numer / (denom1 - denom2);

    itani
}