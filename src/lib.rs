pub mod voxel;
pub mod file_io;
pub mod isim;
pub mod birch;
pub mod ascii;
pub mod utils;
pub mod args;

pub use voxel::{
    voxelize,
    voxelize_stream,
    VoxelGrid,
    get_recommended_info,
    get_recommended_info_stream
};
pub use file_io::{
    read_mol2_file, 
    read_mol2_file_stream,
    write_cluster_mol_ids, 
    write_mol2s_via_cluster_ind
};
pub use isim::{jt_isim_real, jt_isim_binary};
pub use utils::{calc_time_breakdown, init_logging};
pub use birch::VoxBirch;
pub use args::ArgsV;
 

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_in_one_mol() {
        use std::path::Path;
        let file_path = String::from("./test_files/one.mol2");

        let (
            l_mols,
            _
        ) = 
            read_mol2_file(Path::new(&file_path),false)
            .expect("Failed to read MOL2 file");
        assert_eq!(
            l_mols.len(),
            1
        );

    }

    #[test]
    fn voxelize_mol(){
        use std::path::Path;
        use std::io::Write;

        let mut stdout: Box<dyn Write> = Box::new(std::io::sink());

        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            l_mols,
            _
        ) = 
            read_mol2_file(Path::new(&file_path),false)
            .expect("Failed to read MOL2 file"); 
        
        let num_atoms = l_mols[0].num_atoms();

        let grids = 
            voxelize(
                &l_mols, 
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0,
                false,
                &Vec::<String>::new(),
                & mut stdout
            ); 

        let condense_sum: u32 = grids.as_ref().unwrap()[0].condensed_data.iter().sum();
        let sum: u32 = grids.as_ref().unwrap()[0].data.iter().sum();
        let total_length = grids.as_ref().unwrap()[0].data.len();

        assert_eq!(num_atoms,condense_sum as usize);
        assert_eq!(num_atoms, sum as usize);
        assert_eq!(12*3*5, total_length);
    }


    #[test]
    fn read_in_two_mol() {
        use std::path::Path;
        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            l_mols,
            _
        ) = 
            read_mol2_file(Path::new(&file_path),false)
            .expect("Failed to read MOL2 file");
        assert_eq!(
            l_mols.len(),
            2
        );

    }

    #[test]
    fn voxelize_two_mols(){
        use std::path::Path;
        use std::io::Write;

        let mut stdout: Box<dyn Write> = Box::new(std::io::sink());

        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            l_mols,
            _
        )= 
            read_mol2_file(Path::new(&file_path), false)
            .expect("Failed to read MOL2 file"); 

        let grids = 
            voxelize(
                &l_mols, 
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0,
                false,
                &Vec::<String>::new(),
                & mut stdout
            ); 
        
        let l_titles = vec!["ZINC000004771104", "ZINC000108479470"];

        for (i,mol) in l_mols.iter().enumerate() {
            let num_atoms = mol.num_atoms();
            let condense_sum: u32 = grids.as_ref().unwrap()[i].condensed_data.iter().sum();
            let sum: u32 = grids.as_ref().unwrap()[i].data.iter().sum();
            let total_length = grids.as_ref().unwrap()[i].data.len();

            assert_eq!(num_atoms,condense_sum as usize);
            assert_eq!(num_atoms, sum as usize);
            assert_eq!(12*3*5, total_length);
            assert_eq!(grids.as_ref().unwrap()[i].title, l_titles[i]);
        }

    }

    #[test]
    fn get_rec_info(){
        use std::path::Path;

        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/two.mol2");

        let (
            l_mols,
            _
        ) = 
            read_mol2_file(Path::new(&file_path), false)
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
        ) = get_recommended_info(
            &l_mols, 
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
    }

    #[test]
    fn cluster_and_writeout(){
        use std::path::Path;
        use nalgebra::DMatrix;
        use std::io::Write;

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
            quiet: true
        };

        let mut stdout: Box<dyn Write> = Box::new(std::io::sink());


        let file_path = 
            Path::new(env!("CARGO_MANIFEST_DIR")).join(&args.path);

        let (
            l_mols,
            _
        ) = 
            read_mol2_file(Path::new(&file_path),false)
            .expect("Failed to read MOL2 file"); 

        let grids = 
            voxelize(
                &l_mols, 
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0,
                false,
                &Vec::<String>::new(),
                & mut stdout
            ); 

        let mut titles: Vec<String> = Vec::new();

        for grid in grids.as_ref().unwrap().iter() {
            titles.push(grid.title.clone());
        }

        let num_rows = grids.as_ref().unwrap().len();
        let num_cond_cols = grids.as_ref().unwrap()[0].condensed_data.len(); 
        let num_cols = grids.as_ref().unwrap()[0].data.len(); 
        let num_cond_cols_sec = grids.as_ref().unwrap()[1].condensed_data.len(); 
        let num_cols_sec = grids.as_ref().unwrap()[1].data.len(); 

        assert_eq!(num_rows, 2);
        assert_eq!(num_cond_cols, 12);
        assert_eq!(num_cols, num_cols_sec);
        assert_eq!(num_cond_cols_sec, 12);

        let mut vb = VoxBirch::new(
            0.50, 
            10
        );

        let mut input_matrix: DMatrix<f32> = DMatrix::zeros(
            num_rows, 
            num_cols
        );

        for (i, grid) in grids.unwrap().iter().enumerate(){

            for (j, &value) in grid.data
                .iter().enumerate() {
                input_matrix[(i, j)] = value as f32; // Convert u8 to f32 and assign
            }
        }

        vb.cluster(&input_matrix, titles, & mut stdout);

        let cluster_mol_ids: Vec<Vec<(String, usize)>> = vb.get_cluster_mol_ids();

        assert_eq!(cluster_mol_ids.len(),2);
        assert_eq!(cluster_mol_ids[0][0].0,"ZINC000004771104");
        assert_eq!(cluster_mol_ids[1][0].0,"ZINC000108479470");

    }
}