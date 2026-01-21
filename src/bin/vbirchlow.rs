use voxbirch::voxelize_stream;
use voxbirch::read_mol2_file_stream;
use voxbirch::{
    write_cluster_mol_ids,
    write_mol2s_via_cluster_ind
};
use voxbirch::birch::VoxBirch;
use voxbirch::{get_recommended_info_stream,condense_data_stream};
use voxbirch::{ calc_time_breakdown, init_logging };
use voxbirch::ArgsV;
use voxbirch::VoxelGrid;

use std::path::Path;
use clap::Parser;
use std::time::Instant;
use std::fs::File;
use std::io::{Write, BufWriter};
use wincode::{deserialize};
use std::io::{BufReader, Read};

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
    // TODO: I am not using atom_typing nor condensation
    let no_condense = args.no_condense;
    let atom_typing = args.atom_typing;
    
    init_logging(args.verbosity);
    let quiet = args.quiet;

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
    voxbirch::ascii::print_ascii_art(& mut stdout,true);

    // Read MOL2 file
    let path = Path::new(&file_path);
    let (
        all_atom_types,
        num_mols
     ) = read_mol2_file_stream(path, atom_typing).expect("Failed to read MOL2 file");

    writeln!(stdout,"################################################").unwrap();
    writeln!(stdout,"MOL2 file path: {}", file_path).unwrap();
    writeln!(stdout,"Number of molecules read: {}", num_mols).unwrap();
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
    ) = get_recommended_info_stream(resolution, x0, y0, z0);
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
    ) = get_recommended_info_stream(resolution, min_x.floor(), min_y.floor(), min_z.floor());

    writeln!(stdout,
        "Required dims from recommended origin ({},{},{}): {},{},{}", 
        min_x.floor(), min_y.floor(), min_z.floor(),
        rec_need_x_user, rec_need_y_user, rec_need_z_user
    ).unwrap();
    writeln!(stdout,"\nVoxelizing...").unwrap();
    // Voxelization of molecules's xyz's
    let (
        num_rows,
        num_cols,
        num_condensed_cols,
        condensed_data_idx
    ) = voxelize_stream(
        [dimx, dimy, dimz], 
        resolution, 
        x0, y0, z0, 
        atom_typing,
        all_atom_types,
        &mut stdout
    ).expect("Failed to voxelize"); 

    if !no_condense {
        condense_data_stream(condensed_data_idx).expect("Error during condensation");
    }

    let voxelize_duration: std::time::Duration = start_time.elapsed();
    
    if !no_condense {
        writeln!(stdout,"Shape of data: ({} molecules, {} voxels)", num_rows, num_condensed_cols).unwrap();
    } else {
        writeln!(stdout,"Shape of data: ({} molecules, {} voxels)", num_rows, num_cols).unwrap();
    }


    let mut vb = VoxBirch::new(
        threshold, 
        max_branches
    );

    // --- STREAM READ ---
    let read_binary_file;

    if !no_condense {
        read_binary_file = File::open("./tmp/grids_condensed_stream.binary.tmp");
    } else {
        read_binary_file = File::open("./tmp/grids_stream.binary.tmp");
    }
    let mut reader = BufReader::new(read_binary_file.unwrap());

    let mut iter: u64 = 0;
    
    writeln!(stdout,"\n#############################\nFitting grids with the VoxBirch Clustering\n#############################\n\n").unwrap();
    loop {

        // Read length prefix (4 bytes)
        let mut len_buf = [0u8; 4];

        match reader.read_exact(&mut len_buf) {
            Ok(_) => {},
            Err(e) => {
                if e.kind() == std::io::ErrorKind::UnexpectedEof {
                    break;
                } else {
                    break;
                }
            }
        }

        let len = u32::from_le_bytes(len_buf) as usize;

        // Read exactly `len` bytes
        let mut record_buf = vec![0u8; len];
        let _ = reader.read_exact(&mut record_buf);

        // Deserialize the struct
        let grid: VoxelGrid = deserialize(&record_buf).unwrap();

        let vec_grid: Vec<f32> = if !no_condense {
            grid.condensed_data.into_iter().map(|x| x as f32).collect()
        } else {
            grid.data.into_iter().map(|x| x as f32).collect()
        };

        writeln!(stdout, "\nInserting {}: {}/{}", 
            &grid.title, 
            iter+1, num_rows).unwrap();
        // start clustering by inserting the entire data as a whole
        vb.insert(
            Some(&vec_grid), 
            &grid.title,
            iter,
            & mut stdout
        );

        iter += 1;

        vb.first_call = false;
    }

    
    
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

    let process_rate = if clust_secs_entire == vox_secs_entire {
        9999999999.0
    } else {
        iter as f64 / (clust_secs_entire - vox_secs_entire)as f64
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
