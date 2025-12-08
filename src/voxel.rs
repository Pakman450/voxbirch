use crate::file_io::VoxMol;


pub struct VoxelGrid {
    pub dims: [usize; 3],
    pub data: Vec<u8>,    
}

impl VoxelGrid {
    pub fn new(dims: [usize; 3]) -> Self {
        let size = dims[0] * dims[1] * dims[2];
        Self {
            dims,
            data: vec![0; size],
        }
    }

    // This function computes the 1D index in the voxel grid from 3D coordinates
    pub fn voxel_index(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.dims[0] * (y + self.dims[1] * z)
    }

}

// Example "voxelizer": fills grid with a single value
// FIX: currently only for positive numbers
// I have to Shift coordinates to positive first
// I should find the minimium coordinates first and shift accordingly
pub fn voxelize(
    l_mols: Vec<VoxMol>,
    dims: [usize; 3],
    rs: f32,
    x0: f32,
    y0: f32,
    z0: f32,
    method: String,
) -> Vec::<VoxelGrid> {

    // create a list of voxel grids based on number of molecules

    let mut grids = Vec::<VoxelGrid>::new();


    for mol in l_mols.iter() {

        let mut grid = VoxelGrid::new(dims);

        for atom_idx in 0..mol.num_atoms() {
            let x = mol.x[atom_idx];
            let y = mol.y[atom_idx];
            let z = mol.z[atom_idx];

            if x == -1000000000.0 {
                println!("Skipping invalid atom coordinates");
            }

            let ix = ((x - x0) / rs).floor() as usize;
            let iy = ((y - y0) / rs).floor() as usize;
            let iz = ((z - z0) / rs).floor() as usize;
        
            if x == -1000000000.0 {
                println!("ix: {}", ix);
                println!("iy: {}", iy);
                println!("iz: {}", iz);
            }
            let index = grid.voxel_index(ix, iy, iz);
            if method == "binary" {
                grid.data[index] = 1; // Increment voxel value
            } else {
                grid.data[index] += 1; // Mark voxel as occupied

            }

        }

        grids.push(grid);
    }

    let grid = &grids[0];  // for now just return the first grid
    // print out the number of occupied voxels
    let occupied_voxels = grid.data.iter().filter(|&&v| v > 0).count();
    println!("Number of occupied voxels: {}", occupied_voxels);

    // print out the sum of all voxel values
    let voxel_sum: usize = grid.data.iter().map(|&v| v as usize).sum();
    println!("Sum of all voxel values: {}", voxel_sum);

    grids
}




