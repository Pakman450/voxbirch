use crate::voxel::VoxelGrid;

// Itani similarity index for binary voxel grids
// NOTE: assumes all grids have the same dimensions and 
// number of bits 
pub fn itani_bin(l_grids: Vec<VoxelGrid>) -> f32 {

    // num of bits
    let m = l_grids[0].data.len();

    // num of grids
    let n = l_grids.len() as f32;

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
            + (kq[i] * (n - kq[i]));
    }

    let itani = numer / denom;

    itani
}