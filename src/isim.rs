use crate::voxel::VoxelGrid;
use log::{debug};

// Itani similarity index for binary voxel grids
// NOTE: assumes all grids have the same dimensions and 
// number of bits 
pub fn itani_bin_grids(l_grids: &Vec<VoxelGrid>) -> f32 {

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

pub fn itani_real_grids(l_grids: &Vec<VoxelGrid>) -> f32 {

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

// Generates the jt isim score by using 
// ls and ss vectors
pub fn jt_isim_real(
    ls_total: &Vec<f32>,
    ss_total: &Vec<f32>,
    n: usize // num_of_grids
) -> f32 {

    // num of bits
    let m = ls_total.len();

    // vector to hold column wise sums
    let mut xq: Vec<f32> = vec![0.0; m];

    // vector to hold column wise sums
    let mut xq_sq: Vec<f32> = vec![0.0; m]; 

    // compute columwise sums and squares
    for i in 0..m {
        xq[i] = ls_total[i];
        xq_sq[i] = ss_total[i];
    }

    // print out trye if xq and xq_sq have nonzero values
    // for i in 0..m {
    //     if xq[i] != 0.0 || xq_sq[i] !=
    //         0.0 {
    //             debug!("xq[{}] = {}, xq_sq[{}] = {}", i, xq[i], i, xq_sq[i]);
    //         }
    // }

    // compute itani index
    let mut numer = vec![0.0; m];

    for i in 0..m {
        numer[i] = ((xq[i] * xq[i]) - xq_sq[i]) / 2.0;
    }   

    let numer_sum: f32 = numer.iter().sum();
    let ss_sum: f32 = ss_total.iter().sum();
    let inners=(n as f32 - 1.0) * ss_sum;

    let result_score = numer_sum / (inners - numer_sum);

    debug!("jt_isim_real = {} : n_bits = {}, n_objects = {}, numer_sum = {}", 
        result_score,
        m,
        n,
        numer_sum
        );

    result_score

}

// Define the jt_isim function
pub fn jt_isim_binary(c_total: &Vec<f32>, n_objects: usize) -> f32 {
    // Sum of the elements in c_total (column-wise sum)
    let sum_kq: f32 = c_total.iter().copied().sum();

    // Sum of squares (dot product of c_total with itself)
    let sum_kqsq: f32 = c_total.iter().copied().map(|x| x * x).sum();

    // Compute the variable a
    let a = (sum_kqsq - sum_kq) / 2.0;

    let jt_isim_val = a / (a + (n_objects as f32) * sum_kq - sum_kqsq);

    println!("ji_isim score: {}", jt_isim_val);

    debug!("jt_isim = {} : a = {}, sum_kq = {}, sum_kqsq = {}, n_objects = {}", 
        jt_isim_val,
        a, 
        sum_kq,
        sum_kqsq, 
        n_objects
        );

    // Return the iSIM Jaccard-Tanimoto value
    jt_isim_val
}
