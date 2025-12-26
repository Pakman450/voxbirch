use voxelizer::voxelize;
use voxelizer::read_mol2_file;
use voxelizer::itani_bin;
use voxelizer::itani_real;
use voxelizer::diameter_real;
use voxelizer::birch::VoxBirch;

use std::path::{Path};
use nalgebra::DMatrix;
use clap::Parser;
use rand::Rng;  

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the MOL2 file
    #[arg(short, long)]
    path: String,

    /// Dimensions of the voxel grid (x, y, z), comma separated
    #[arg(short, long, value_delimiter = ',', default_value = "100,100,100")]
    dims: Vec<usize>,

    /// resolution of the voxel grid
    #[arg(short, long, default_value_t = 0.125)]
    resolution: f32,

    /// target origin x0 y0 z0, comma separated
    #[arg(short, long, value_delimiter = ',', default_value = "0.0,0.0,0.0")]
    origin: Vec<f32>

}



fn main() {

    let args = Args::parse();

    // Argument unpacking
    let file_path: String = args.path;
    let dimx = args.dims[0];
    let dimy = args.dims[1];
    let dimz = args.dims[2];
    let x0 = args.origin[0];
    let y0 = args.origin[1];       
    let z0 = args.origin[2];
    let resolution = args.resolution;
    // let method = args.method;

    // Read MOL2 file
    let path = Path::new(&file_path);

    // Get molecule list
    let l_mols = read_mol2_file(path).expect("Failed to read MOL2 file");

    // Give user the recommended origin values for placing voxels.
    let mut l_vals = Vec::<f32>::new();

    for mol in &l_mols {
        l_vals.extend(mol.x.clone())
    }

    // Give user origin recommendation 
    if let Some(min_value) = &l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()) {
        println!("The minimum x value is: {}", min_value);
    } else {
        panic!("No valid minimum found (possibly due to NaN values).");
    }

    l_vals.clear();

    for mol in &l_mols {
        l_vals.extend(mol.y.clone())
    }

    // Give user origin recommendation 
    if let Some(min_value) = &l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()) {
        println!("The minimum y value is: {}", min_value);
    } else {
        panic!("No valid minimum found (possibly due to NaN values).");
    }

    l_vals.clear();


    for mol in &l_mols {
        l_vals.extend(mol.z.clone())
    }

    // Give user origin recommendation 
    if let Some(min_value) = &l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()) {
        println!("The minimum z value is: {}", min_value);
    } else {
        panic!("No valid minimum found (possibly due to NaN values).");
    }

    // Voxelization of molecules's xyz's
    let grids = voxelize(l_mols, [dimx, dimy, dimz], resolution, x0, y0, z0); 

    // Print some voxel grid info
    println!("Voxel Grid Dimensions: {:?}", grids[0].dims);

    let mut vb = VoxBirch::new(
        0.5, // threshold
        50 // branches
    );

    // Get the number of rows (which is the number of VoxelGrids)
    let num_rows = grids.len();
    
    // Get the number of columns (which is the length of the data in each VoxelGrid)
    let num_cols = grids[0].data.len();  // Assuming all VoxelGrids have the same length of data

    // Create the DMatrix with the correct size
    let mut input_matrix: DMatrix<f32> = DMatrix::zeros(
        num_rows, 
        num_cols
    );

    // Fill the matrix with the data from each VoxelGrid
    for (i, grids) in grids.iter().enumerate() {
        for (j, &value) in grids.data.iter().enumerate() {
            input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
        }
    }

    // start clustering
    vb.fit(&input_matrix, true);
        
}
