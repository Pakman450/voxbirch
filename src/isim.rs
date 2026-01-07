use crate::voxel::VoxelGrid;
use log::{error, warn, info, debug, trace};

// Itani similarity index for binary voxel grids
// NOTE: assumes all grids have the same dimensions and 
// number of bits 
pub fn itani_bin(l_grids: &Vec<VoxelGrid>) -> f32 {

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

pub fn itani_real(l_grids: &Vec<VoxelGrid>) -> f32 {


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



// BUG: This is not correct. Need to implement diameter calculation properly.
// Here it will always return zero. The way binary and real calulcations 
// are different. Double check
pub fn itani_real_no_list(
    c_total: &Vec<f32>,
    n: usize // num_of_grids
) -> f32 {

    // num of bits
    let m = c_total.len();

    // vector to hold column wise sums
    let mut xq: Vec<f32> = vec![0.0; m];

    // vector to hold column wise sums
    let mut xq_sq: Vec<f32> = vec![0.0; m]; 


    for i in 0..m {
        xq[i] = c_total[i];
        xq_sq[i] = (c_total[i] as f32) * (c_total[i] as f32);
    }

    // print out trye if xq and xq_sq have nonzero values
    for i in 0..m {
        if xq[i] != 0.0 || xq_sq[i] !=
            0.0 {
                debug!("xq[{}] = {}, xq_sq[{}] = {}", i, xq[i], i, xq_sq[i]);
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

    debug!("jt_isim = {} : n_bits = {}, n_objects = {}, numer = {}, denom1 = {}, denom2 = {}", 
        itani,
        m,
        n,
        numer,
        denom1,
        denom2
        );

    itani
}


pub fn diameter_real(itani_real: f32) -> f32{
    1.0-itani_real
}

// radius function and its helper functions
// NOTE: not tested
// fn centroid(
//     ls: &Vec<f32>, 
//     nj: u32,
//     dims: [usize; 3]
// ) -> VoxelGrid {

//     let mut cj: Vec<f32> = vec![0.0; ls.len()];

//     let mut cj_grid = VoxelGrid::new(dims);

//     for i in 0..ls.len(){
//         cj[i] = ls[i]/nj as f32 + 0.5;
//         cj[i] = cj[i].floor();
//     }

//     cj_grid.data = cj.iter().map(|&x| x as u8).collect();

//     cj_grid
// }

// NOTE: not tested
pub fn radius_real(
    itani_r: f32, 
    xj: &Vec<VoxelGrid>, 
    cj: &VoxelGrid, 
    nj: u32 
) -> f32 {

    // Make a new vector with capacity for xj + cj
    let mut xj_cj = Vec::with_capacity(xj.len() + 1);

    // Copy values from xj
    xj_cj.extend_from_slice(xj);

    // Push the extra element
    xj_cj.push(cj.clone()); 

    let rj = 1.0 - ( 
        (((nj + 1) as f32 * itani_real(&xj_cj)) -
        ((nj-1) as f32 * itani_r)) / 2.0 
                    );

    rj
}