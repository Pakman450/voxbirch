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
    origin: Vec<f32>,

    /// target origin x0 y0 z0, comma separated
    #[arg(short, long, default_value = "binary")]
    method: String,

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
    let method = args.method;

    // Read MOL2 file
    let path = Path::new(&file_path);

    // Get molecule list
    let l_mols = read_mol2_file(path).expect("Failed to read MOL2 file");

    // Voxelization
    let grids = voxelize(l_mols, [dimx, dimy, dimz], resolution, x0, y0, z0, method.clone()); 

    // Print some voxel grid info
    println!("Voxel Grid Dimensions: {:?}", grids[0].dims);

    if method == "binary" {
        println!("Itani Similarity Score: {}", itani_bin(&grids));
    } else if method == "real" {
        let itani_r = itani_real(&grids);
        println!("Itani Similarity Score: {}", itani_r);
        println!("Diameter: {}", diameter_real(itani_r));
    }


    let mut vb = VoxBirch::new(0.5, 5);
    let mut randmat: DMatrix<f32> = DMatrix::zeros(10, 10);
        // Fill each element with a random float between 0.0 and 1.0
    let mut rng = rand::thread_rng();
    for i in 0..randmat.nrows() {
        for j in 0..randmat.ncols() {
            randmat[(i, j)] = rng.gen::<f32>(); // gen() generates a float in [0,1)
        }
    }
    vb.fit(&randmat, true);
        
}
