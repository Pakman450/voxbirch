pub mod voxel;
pub mod file_io;
pub mod isim;
pub mod birch;
pub mod ascii;
pub mod utils;
pub mod args;
pub mod memmap;

pub use voxel::{
    voxelize_stream,
    VoxelGrid,
    get_recommended_info_stream,
    condense_data_stream
};
pub use file_io::{
    read_mol2_file_stream,
    write_cluster_mol_ids, 
    write_mol2s_via_cluster_ind
};
pub use isim::{jt_isim_real, jt_isim_binary};
pub use utils::{calc_time_breakdown, init_logging, mem_logging};
pub use birch::VoxBirch;
pub use args::ArgsV;
pub use memmap::read_in_file;
 

// NOTE: This is 
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_in_one_mol() {
        use std::path::Path;
        let file_path = String::from("./test_files/one.mol2");

        let (
            all_atom_types,
            num_mols
        ) = 
            read_mol2_file_stream(Path::new(&file_path),true)
            .expect("Failed to read MOL2 file");
        assert_eq!(
            all_atom_types.len(),
            6
        );

        assert_eq!(num_mols, 1);
        let _ = std::fs::remove_dir_all("./tmp");

    }

    #[test]
    fn voxelize_mol(){
        use std::path::Path;
        use std::io::Write;

        let mut stdout: Box<dyn Write> = Box::new(std::io::sink());

        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            all_atom_types,
            num_mols
        ) = 
            read_mol2_file_stream(Path::new(&file_path),false)
            .expect("Failed to read MOL2 file"); 


        let (
            num_rows,
            num_cols,
            num_condensed_cols,
            condensed_data_idx
        ) = 
            voxelize_stream(
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0,
                false,
                all_atom_types,
                & mut stdout
            ).expect("Failed to voxelize");

        println!("{} {} {:?}", num_cols, num_condensed_cols , condensed_data_idx);
        assert_eq!(num_rows, 2);
        assert_eq!(num_mols, 2);
        let _ = std::fs::remove_dir_all("./tmp");


    }


    #[test]
    fn read_in_two_mol() {
        use std::path::Path;
        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            all_atom_types,
            num_mols
        ) = 
            read_mol2_file_stream(Path::new(&file_path),true)
            .expect("Failed to read MOL2 file");
        assert_eq!(
            all_atom_types.len(),
            7
        );

        assert_eq!(
            num_mols,
            2
        );
        let _ = std::fs::remove_dir_all("./tmp");

    }

    #[test]
    fn voxelize_two_mols(){
        use std::path::Path;
        use std::io::Write;

        let mut stdout: Box<dyn Write> = Box::new(std::io::sink());

        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            all_atom_types,
            num_mols
        ) = 
            read_mol2_file_stream(Path::new(&file_path), false)
            .expect("Failed to read MOL2 file"); 


        let (
            num_rows,
            num_cols,
            num_condensed_cols,
            condensed_data_idx
        ) = 
            voxelize_stream(
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0,
                false,
                all_atom_types,
                & mut stdout
            ).expect("Failed to voxelize");
        
        println!("{} {} {:?}", num_cols, num_condensed_cols , condensed_data_idx);
        assert_eq!(num_rows, 2);
        assert_eq!(num_mols, 2);
        let _ = std::fs::remove_dir_all("./tmp");


    }

    #[test]
    fn get_rec_info(){
        use std::path::Path;

        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            all_atom_types,
            num_mols
        ) = 
            read_mol2_file_stream(Path::new(&file_path), false)
            .expect("Failed to read MOL2 file"); 


        let x0 = 0.0;
        let y0 = 0.0;
        let z0 = 0.0;
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
        ) = get_recommended_info_stream(
            0.5, 
            x0, 
            y0, 
            z0
        );

        assert_eq!(min_x.floor(), 20.0);
        assert_eq!(min_y.floor(), -13.0);
        assert_eq!(min_z.floor(), 4.0);
        assert_eq!(
            format!(
                "({}): {},{},{}", 
                format!("{:.3},{:.3},{:.3}", min_x, min_y, min_z), 
                need_x, need_y, need_z
            ),
            "(20.012,-12.414,4.123): 12,35,11"
        );
        assert_eq!(
            format!(
                "({:.3},{:.3},{:.3}): {},{},{}",
                x0, y0, z0, 
                need_x_user, need_y_user, need_z_user
            ),
            "(0.000,0.000,0.000): 52,10,20"
        );
        let _ = std::fs::remove_dir_all("./tmp");

    }

    #[test]
    fn cluster_and_writeout(){
        use std::path::Path;
        use std::io::Write;
        use std::io::{BufReader, Read};
        use std::fs::File;
        use wincode::{deserialize};


        let args = ArgsV {
            path: "test_files/two.mol2".into(),
            dims: vec![20,20,20],
            resolution: 2.0,
            origin: vec![0.0,0.0,0.0],
            threshold: 0.65,
            max_branches: 50,
            clustered_ids_path: Some(String::from("./clustered_mol_ids.txt")),
            output_path: Some(String::from("./voxbirch.out")),
            cluster_write_limit: 100,
            no_condense: false,
            atom_typing: false,
            verbosity: 0,
            quiet: true,
            rss: false
        };

        let mut stdout: Box<dyn Write> = Box::new(std::io::sink());


        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join(&args.path);

        let (
            all_atom_types,
            num_mols
        ) = 
            read_mol2_file_stream(Path::new(&file_path),false)
            .expect("Failed to read MOL2 file"); 

        let (
            num_rows,
            num_cols,
            num_condensed_cols,
            condensed_data_idx
        ) = voxelize_stream(
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0,
                false,
                all_atom_types,
                & mut stdout
            ).expect("Failed to voxelize"); 

        assert_eq!(num_rows, 2);
        assert_eq!(num_cols, 180);

        let mut vb = VoxBirch::new(
            0.50, 
            10,
            num_cols as usize
        );

        let read_binary_file;

        read_binary_file = File::open("./tmp/grids_stream.binary.tmp");

        let mut reader = BufReader::new(read_binary_file.unwrap());

        let mut iter: u64 = 0;

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

            let vec_grid: Vec<f32> = grid.data.into_iter().map(|x| x as f32).collect();

            // start clustering by inserting the entire data as a whole
            vb.insert(
                Some(&vec_grid), 
                &grid.title,
                iter,
                & mut stdout
            );

            iter += 1;
        }

        let cluster_mol_ids: Vec<Vec<(String, usize)>> = vb.get_cluster_mol_ids();

        assert_eq!(cluster_mol_ids.len(),2);
        assert_eq!(cluster_mol_ids[0][0].0,"ZINC000004771104");
        assert_eq!(cluster_mol_ids[1][0].0,"ZINC000108479470");
        let _ = std::fs::remove_dir_all("./tmp");

    }
}