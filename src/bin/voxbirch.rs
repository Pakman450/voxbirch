use voxbirch::voxelize;
use voxbirch::read_mol2_file;
use voxbirch::{
    write_cluster_mol_ids,
    write_mol2s_via_cluster_ind
};
use voxbirch::birch::VoxBirch;
use voxbirch::get_recommended_info;
use voxbirch::{ calc_time_breakdown, init_logging, mem_logging };
use voxbirch::ArgsV;

use std::path::Path;
use nalgebra::{DMatrix, RowDVector};
use clap::Parser;
use std::time::Instant;
use std::fs::File;
use std::io::{Write, BufWriter};

fn main() {

    // Get the current time
    let start_time = Instant::now();
        
    // Argument unpacking
    let args = ArgsV::parse();
    let file_path= args.path.clone();
    let dimx = args.dims[0];
    let dimy = args.dims[1];
    let dimz = args.dims[2];
    let x0 = args.origin[0];
    let y0 = args.origin[1];       
    let z0 = args.origin[2];
    let resolution = args.resolution;
    let threshold = args.threshold;
    let max_branches = args.max_branches;
    let mut clustered_mol_id_string = args.clustered_ids_path.clone();
    let cluster_write_limit = args.cluster_write_limit;
    if args.clustered_ids_path.is_none() {
        clustered_mol_id_string = Some(String::from("./clustered_mol_ids.txt"));
    }
    let no_condense = args.no_condense;
    let atom_typing = args.atom_typing;
    
    init_logging(args.verbosity);
    let quiet = args.quiet;
    let rss = args.rss;

    if rss {
        mem_logging();
    }
    // Print some input info
    let mut stdout: Box<dyn Write> = if let Some(path) = &args.output_path
    {
        Box::new(BufWriter::new(File::create(path).unwrap()))
    } else if quiet {
        Box::new(std::io::sink())
    }else {
        Box::new(std::io::stdout().lock())
    };

    // Print the ASCII art
    voxbirch::ascii::print_ascii_art(& mut stdout, false);

    // Read MOL2 file
    let path = Path::new(&file_path);
    let (
        l_mols,
        all_atom_types
     ) = read_mol2_file(path, atom_typing).expect("Failed to read MOL2 file");

    writeln!(stdout,"################################################").unwrap();
    writeln!(stdout,"MOL2 file path: {}", file_path).unwrap();
    writeln!(stdout,"Number of molecules read: {}", l_mols.len()).unwrap();
    writeln!(stdout,"Voxel Grid Dimensions: {} x {} x {}", dimx, dimy, dimz).unwrap();
    writeln!(stdout,"Voxel Grid Resolution: {}", resolution).unwrap();
    writeln!(stdout,"Voxel Grid Origin: ({}, {}, {})", x0, y0, z0).unwrap();
    writeln!(stdout,"Thresold: {}", threshold).unwrap();
    writeln!(stdout,"Max Branches: {}", max_branches).unwrap();
    writeln!(stdout,"Enforce Atom Typing: {}", atom_typing).unwrap();
    writeln!(stdout,"Condense Voxel Grids: {}", !no_condense).unwrap();
    writeln!(stdout,"Quiet mode: {}", quiet).unwrap();
    writeln!(stdout,"################################################").unwrap();

    writeln!(stdout,"\nGrabbing recommended voxelization parameters...").unwrap();
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
    writeln!(stdout,"The recommended origin: {},{},{}", min_x.floor(), min_y.floor(), min_z.floor()).unwrap();

    writeln!(stdout,
        "Minimal voxel grid dimensions to cover all molecules\n\t(from absolute origin {}): {},{},{}",
        format!("{:.3},{:.3},{:.3}", min_x, min_y, min_z), need_x, need_y, need_z
    ).unwrap();
    writeln!(stdout,
        "Required dims from provided origin ({:.3},{:.3},{:.3}): {},{},{}", 
        x0, y0, z0, 
        need_x_user, need_y_user, need_z_user
    ).unwrap();

    let (
        _, 
        _, 
        _, 
        _, 
        _, 
        _,
        rec_need_x_user,
        rec_need_y_user,
        rec_need_z_user
    ) = get_recommended_info(&l_mols, resolution, min_x.floor(), min_y.floor(), min_z.floor());

    writeln!(stdout,
        "Required dims from recommended origin ({},{},{}): {},{},{}", 
        min_x.floor(), min_y.floor(), min_z.floor(),
        rec_need_x_user, rec_need_y_user, rec_need_z_user
    ).unwrap();

    writeln!(stdout,"\nVoxelizing...").unwrap();
    // Voxelization of molecules's xyz's
    let grids = voxelize(
        &l_mols, 
        [dimx, dimy, dimz], 
        resolution, 
        x0, y0, z0, 
        atom_typing, &all_atom_types, 
        &mut stdout
    ); 

    let voxelize_duration: std::time::Duration = start_time.elapsed();

    // Get the number of rows (which is the number of VoxelGrids)
    let num_rows = grids.as_ref().unwrap().len();

    // Get the number of cols (which is the number of voxels)
    let num_cols;

    if !no_condense {
        num_cols = grids.as_ref().unwrap()[0].condensed_data.len(); 
    } else{
        num_cols = grids.as_ref().unwrap()[0].data.len(); 
    }
    
    writeln!(stdout,"Shape of data: ({} molecules, {} voxels)", num_rows, num_cols).unwrap();

    // Create the DMatrix with the correct size
    let mut input_matrix: DMatrix<f32> = DMatrix::zeros(
        num_rows, 
        num_cols
    );

    if !no_condense {
        for (i, grid) in grids.as_ref().unwrap().iter().enumerate(){

            for (j, &value) in grid.condensed_data
                .iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }
    } else{
        for (i, grid) in grids.as_ref().unwrap().iter().enumerate(){

            for (j, &value) in grid.data
                .iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }
    }

    let mut vb = VoxBirch::new(
        threshold, 
        max_branches
    );

    if !input_matrix.iter().any(|&x| x != 0.0) {
        panic!("All of your voxels for all rows have 0.0s ");
    }

    // start clustering by inserting the entire data as a whole
    vb.cluster(
        &input_matrix, 
        grids.as_ref().unwrap().iter().map(|g| g.title.clone()).collect(),
        & mut stdout
    );

    let clustering_duration: std::time::Duration = start_time.elapsed();

    // Get results after clustering. 
    let cluster_mol_ids: Vec<Vec<(String, usize)>> = vb.get_cluster_mol_ids();
    let num_clusters = cluster_mol_ids.len();
    let path_cluster_ids: String = clustered_mol_id_string.unwrap();
    let write_to_path = Path::new(&path_cluster_ids);

    let _ = write_mol2s_via_cluster_ind(
        &cluster_mol_ids,
        &path,
        cluster_write_limit
    );

    writeln!(stdout,"Writing cluster mol ids to: {:?}", write_to_path).unwrap();
    let _ = write_cluster_mol_ids(&write_to_path, &cluster_mol_ids);

    // TODO: I need to write out mol2s for the 
    // the top ranking clusters.  

    // Get the breakdown of elapsed time
    let (
        vox_secs_entire,
        vox_hours,
        vox_minutes,
        vox_seconds,
        vox_milliseconds
    ) = calc_time_breakdown(&voxelize_duration);

    let (
        clust_secs_entire,
        clust_hours,
        clust_minutes,
        clust_seconds,
        clust_milliseconds
    ) = calc_time_breakdown(&clustering_duration);

    let total_duration: std::time::Duration = start_time.elapsed();
    let (
        _,
        tot_hours,
        tot_minutes,
        tot_seconds,
        tot_milliseconds
    ) = calc_time_breakdown(&total_duration);

    let milli_or_sec: &str = 
        if clust_secs_entire == clust_milliseconds as u64
        {
            "ms"
        } else {
            "sec"
        };

    
    let process_rate: f32 = if clust_secs_entire == vox_secs_entire 
    {       
         l_mols.len() as f32 / 0.0000000000001
    } else {
         l_mols.len() as f32 / (clust_secs_entire - vox_secs_entire) as f32

    };


    writeln!(stdout,
"\nFinished
Summary statistics:
Total number of clusters: {}
Elapsed time for Voxelization: {}h {}m {}s {}ms
Elapsed time for Clustering: {}h {}m {}s {}ms
Number of molecules processed per {} during Clustering: {}
Number of {} per molecule processed during Clustering: {}
Total Elapsed time: {}h {}m {}s {}ms",
        num_clusters,
        
        vox_hours, 
        vox_minutes,
        vox_seconds, 
        vox_milliseconds,
        
        clust_hours - vox_hours, 
        clust_minutes - vox_minutes, 
        clust_seconds - vox_seconds, 
        clust_milliseconds - vox_milliseconds,

        milli_or_sec,
        process_rate,
        
        milli_or_sec,
        1.0 / process_rate,

        tot_hours, 
        tot_minutes, 
        tot_seconds, 
        tot_milliseconds
    ).unwrap();

}
