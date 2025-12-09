use voxelizer::voxelize;
use voxelizer::read_mol2_file;
use voxelizer::itani_bin;
use voxelizer::itani_real;
use voxelizer::diameter_real;

use std::path::{Path};

use clap::Parser;

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
        return;
    } else if method == "real" {
        let itani_r = itani_real(&grids);
        println!("Itani Similarity Score: {}", itani_r);
        println!("Diameter: {}", diameter_real(itani_r));
        return;
    }
    
}
