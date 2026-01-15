use voxbirch::voxelize;
use voxbirch::read_mol2_file;
use voxbirch::write_cluster_mol_ids;
use voxbirch::birch::VoxBirch;
use voxbirch::get_recommended_info;

use std::path::{Path};
use nalgebra::DMatrix;
use clap::Parser;
use std::time::{Instant};
use std::env;
use env_logger::{Builder};
use std::io::Write;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the MOL2 file (required)
    #[arg(short, long, required = true)]
    path: String,

    /// Dimensions of the voxel grid (x, y, z), comma separated
    #[arg(short, long, value_delimiter = ',', default_value = "20,20,20")]
    dims: Vec<usize>,

    /// Resolution of the voxel grid in Angstroms
    #[arg(short, long, default_value_t = 2.0)]
    resolution: f32,

    /// Target origin x0 y0 z0 via comma separated string
    #[arg(short, long, value_delimiter = ',', default_value = "0.0,0.0,0.0")]
    origin: Vec<f32>,

    /// Threshold of similarity
    #[arg(short, long, default_value_t = 0.65)]
    threshold: f32,

    /// Number of max branches
    #[arg(short, long, default_value_t = 50)]
    max_branches: usize,

    /// Clustered mol ids output name
    #[arg(long, default_value = "./clustered_mol_ids.txt")]
    output_file_path: String,

    /// Do not condense voxel grids. Leaving this out condenses grids.
    #[arg(long)]
    no_condense: bool,

    /// Verbosity level
    #[arg(short, long, action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,

}



fn main() {

    // Get the current time
    let start_time = Instant::now();

    // Print the ASCII art
    voxbirch::ascii::print_ascii_art();
        
    // Argument unpacking
    let args = Args::parse();
    let file_path= args.path;
    let dimx = args.dims[0];
    let dimy = args.dims[1];
    let dimz = args.dims[2];
    let x0 = args.origin[0];
    let y0 = args.origin[1];       
    let z0 = args.origin[2];
    let resolution = args.resolution;
    let threshold = args.threshold;
    let max_branches = args.max_branches;
    let clustered_mol_id_string = args.output_file_path;
    let no_condense = args.no_condense;
    let verbosity = args.verbosity;

    // Initialize the logger with appropriate level
    if verbosity == 2 {
        // Set RUST_LOG to debug level if verbose
        env::set_var("RUST_LOG", "debug");
        env::set_var("RUST_LOG", "warn");
        env::set_var("RUST_LOG", "info");
        env::set_var("RUST_LOG", "error");
        env::set_var("RUST_LOG", "trace");
    } 

    // Initialize the logger with custom format
    Builder::from_default_env()
        .format(|buf, record| {
            // IMPORTANT: keep the style alive
            let level_style = buf.default_level_style(record.level());
            let level = level_style.value(record.level());

            let file = record.file().unwrap_or("unknown");
            let line = record.line().unwrap_or(0);

            writeln!(
                buf,
                "[{} {}:{} {}] {}",
                level,
                file,
                line,
                record.target(),
                record.args()
            )
        }).init();

    // Read MOL2 file
    let path = Path::new(&file_path);
    let l_mols = read_mol2_file(path).expect("Failed to read MOL2 file");

    // Print some input info
    println!("################################################");
    println!("MOL2 file path: {}", file_path);
    println!("Number of molecules read: {}", l_mols.len());
    println!("Voxel Grid Dimensions: {} x {} x {}", dimx, dimy, dimz);
    println!("Voxel Grid Resolution: {}", resolution);
    println!("Voxel Grid Origin: ({}, {}, {})", x0, y0, z0);
    println!("Thresold: {}", threshold);
    println!("Max branches: {}", max_branches);
    println!("Condense Voxel Grids: {}", !no_condense);
    println!("################################################");

    println!("\nGrabbing recommended voxelization parameters...");
    // Give user the recommended origin values for placing voxels.
    let (
        min_x, 
        min_y, 
        min_z, 
        need_x, 
        need_y, 
        need_z,
        need_x_user,
        need_y_user,
        need_z_user
    ) = get_recommended_info(&l_mols, resolution, x0, y0, z0);
    println!("The recommended origin: {},{},{}", min_x.floor(), min_y.floor(), min_z.floor());

    println!(
        "Minimal voxel grid dimensions to cover all molecules\n\t(from absolute origin {}): {},{},{}",
        format!("{:.3},{:.3},{:.3}", min_x, min_y, min_z), need_x, need_y, need_z
    );
    println!(
        "Required dims from provided origin ({:.3},{:.3},{:.3}): {},{},{}", 
        x0, y0, z0, 
        need_x_user, need_y_user, need_z_user
    );

    println!("\nVoxelizing...");
    // Voxelization of molecules's xyz's
    let grids = voxelize(&l_mols, [dimx, dimy, dimz], resolution, x0, y0, z0); 

    // Get the number of rows (which is the number of VoxelGrids)
    let num_rows = grids.len();
    let num_cols;

    if !no_condense {
        num_cols = grids[0].condensed_data.len(); 
    } else{
        num_cols = grids[0].data.len(); 
    }

    println!("Shape of data: ({} molecules, {} voxels)", num_rows, num_cols);

    // Initialize VoxBirch
    let mut vb = VoxBirch::new(
        threshold, 
        max_branches
    );

    // Create the DMatrix with the correct size
    let mut input_matrix: DMatrix<f32> = DMatrix::zeros(
        num_rows, 
        num_cols
    );

    let mut titles: Vec<String> = Vec::new();

    if !no_condense {
        // Fill the matrix with the data from each VoxelGrid
        for (i, grids) in grids.iter().enumerate() {
            for (j, &value) in grids.condensed_data.iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }
    } else{
        // Fill the matrix with the data from each VoxelGrid
        for (i, grids) in grids.iter().enumerate() {
            for (j, &value) in grids.data.iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }
    }

    if !input_matrix.iter().any(|&x| x != 0.0) {
        panic!("All of your voxels for all rows have 0.0s ");
    }

    // Collect titles
    for grids in grids.iter() {
        titles.push(grids.title.clone());
    }

    // start clustering
    vb.fit(&input_matrix, titles, true);


    // Get results after clustering. 
    let cluster_mol_ids = vb.get_cluster_mol_ids();
    let path_cluster_ids: String = String::from(clustered_mol_id_string);
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
