use voxelizer::voxelize;
use voxelizer::read_mol2_file;
use voxelizer::write_cluster_mol_ids;

use voxelizer::birch::VoxBirch;

use std::path::{Path};
use nalgebra::DMatrix;
use clap::Parser;
use std::time::{Instant};
use log::{error, warn, info, debug, trace};
use std::env;

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
    #[arg(short, long, default_value_t = 0.5)]
    resolution: f32,

    /// target origin x0 y0 z0, comma separated
    #[arg(short, long, value_delimiter = ',', default_value = "0.0,0.0,0.0")]
    origin: Vec<f32>,

    // threshold
    #[arg(short, long, default_value_t = 0.65)]
    threshold: f32,

    // max_branches
    #[arg(short, long, default_value_t = 50)]
    max_branches: usize,

    // verbosity level
    #[arg(short, long, default_value_t = false)]
    verbosity: bool,
}



fn main() {



    // Get the current time
    let start_time = Instant::now();



    let ascii_art = r#"
                                     ,---,.                            ,---,     
       ,---.                       ,'  .'  \  ,--,                   ,--.' |     
      /__./|   ,---.             ,---.' .' |,--.'|    __  ,-.        |  |  :     
 ,---.;  ; |  '   ,'\ ,--,  ,--, |   |  |: ||  |,   ,' ,'/ /|        :  :  :     
/___/ \  | | /   /   ||'. \/ .`| :   :  :  /`--'_   '  | |' | ,---.  :  |  |,--. 
\   ;  \ ' |.   ; ,. :'  \/  / ; :   |    ; ,' ,'|  |  |   ,'/     \ |  :  '   | 
 \   \  \: |'   | |: : \  \.' /  |   :     \'  | |  '  :  / /    / ' |  |   /' : 
  ;   \  ' .'   | .; :  \  ;  ;  |   |   . ||  | :  |  | ' .    ' /  '  :  | | | 
   \   \   '|   :    | / \  \  \ '   :  '; |'  : |__;  : | '   ; :__ |  |  ' | : 
    \   `  ; \   \  /./__;   ;  \|   |  | ; |  | '.'|  , ; '   | '.'||  :  :_:,' 
     :   \ |  `----' |   :/\  \ ;|   :   /  ;  :    ;---'  |   :    :|  | ,'     
      '---"          `---'  `--` |   | ,'   |  ,   /        \   \  / `--''       
                                 `----'      ---`-'          `----'              
    "#;

    // Print the ASCII art
    println!("{}", ascii_art);
    println!("Code Written by: Steven Pak\n");

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
    let threshold = args.threshold;
    let max_branches = args.max_branches;
    let verbose = args.verbosity;
    // let method = args.method;

    // Initialize the logger with appropriate level
    if verbose {
        // Set RUST_LOG to debug level if verbose
        env::set_var("RUST_LOG", "debug");
        env::set_var("RUST_LOG", "warn");
        env::set_var("RUST_LOG", "info");
        env::set_var("RUST_LOG", "error");
        env::set_var("RUST_LOG", "trace");
    } 
    env_logger::init();

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
        println!("The recommended minimum x value is: {}", min_value);
    } else {
        panic!("No valid minimum found (possibly due to NaN values).");
    }

    l_vals.clear();

    for mol in &l_mols {
        l_vals.extend(mol.y.clone())
    }

    // Give user origin recommendation 
    if let Some(min_value) = &l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()) {
        println!("The recommended minimum y value is: {}", min_value);
    } else {
        panic!("No valid minimum found (possibly due to NaN values).");
    }

    l_vals.clear();


    for mol in &l_mols {
        l_vals.extend(mol.z.clone())
    }

    // Give user origin recommendation 
    if let Some(min_value) = &l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()) {
        println!("The recommended minimum z value is: {}", min_value);
    } else {
        panic!("No valid minimum found (possibly due to NaN values).");
    }

    // Voxelization of molecules's xyz's
    let grids = voxelize(l_mols, [dimx, dimy, dimz], resolution, x0, y0, z0); 

    // Print some voxel grid info
    println!("Voxel Grid Dimensions: {:?}", grids[0].dims);
    println!("Voxel Grid Resolution: {}", resolution);
    println!("Voxel Grid Number of grids of the first grid: {}", grids[0].data.len());

    let mut vb = VoxBirch::new(
        threshold, // threshold
        max_branches // branches
    );

    println!("Thresold: {}", threshold);
    println!("Max brances: {}", max_branches);

    // Get the number of rows (which is the number of VoxelGrids)
    let num_rows = grids.len();
    
    println!("Number of molecules: {}", num_rows);

    // Get the number of columns (which is the length of the data in each VoxelGrid)
    let num_cols = grids[0].data.len();  // Assuming all VoxelGrids have the same length of data

    println!("Number of bits: {}", num_cols);

    // Create the DMatrix with the correct size
    let mut input_matrix: DMatrix<f32> = DMatrix::zeros(
        num_rows, 
        num_cols
    );

    let mut titles: Vec<String> = Vec::new();

    // Fill the matrix with the data from each VoxelGrid
    for (i, grids) in grids.iter().enumerate() {
        for (j, &value) in grids.data.iter().enumerate() {
            input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
        }
    }

    // Collect titles
    for grids in grids.iter() {
        titles.push(grids.title.clone());
    }

    // start clustering
    vb.fit(&input_matrix, titles, true);


    // Get results after clustering. 
    let cluster_mol_ids = vb.get_cluster_mol_ids();
    let path_cluster_ids: String = String::from("./clusters_mol_ids.txt");
    let write_to_path = Path::new(&path_cluster_ids);
    let _ = write_cluster_mol_ids(&write_to_path, &cluster_mol_ids);

    // Get the elapsed time
    let duration = start_time.elapsed();

    // Get total seconds and calculate hours, minutes, seconds
    let total_secs = duration.as_secs();
    let hours = total_secs / 3600;
    let minutes = (total_secs % 3600) / 60;
    let seconds = total_secs % 60;

    // Get milliseconds from remaining nanoseconds
    let milliseconds = duration.subsec_millis();

    // Format the time as a string
    println!(
        "\nFinished\nElapsed time: {}h {}m {}s {}ms",
        hours, minutes, seconds, milliseconds
    );
}
