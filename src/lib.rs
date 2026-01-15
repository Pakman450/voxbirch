pub mod voxel;
pub mod file_io;
pub mod isim;
pub mod birch;
pub mod ascii;

pub use voxel::{voxelize, VoxelGrid, get_recommended_info};
pub use file_io::{read_mol2_file, write_cluster_mol_ids};
pub use isim::{jt_isim_real, jt_isim_binary};
pub use birch::VoxBirch;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_in_one_mol() {
        use std::path::Path;
        let file_path = String::from("./test_files/one.mol2");

        let l_mols = 
            read_mol2_file(Path::new(&file_path))
            .expect("Failed to read MOL2 file");
        assert_eq!(
            l_mols.len(),
            1
        );

    }

    #[test]
    fn voxelize_mol(){
        use std::path::Path;
        let file_path = String::from("./test_files/one.mol2");

        let l_mols = 
            read_mol2_file(Path::new(&file_path))
            .expect("Failed to read MOL2 file"); 
        
        let num_atoms = l_mols[0].num_atoms();

        let grids = 
            voxelize(
                &l_mols, 
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0
            ); 

        let condense_sum: u32 = grids[0].condensed_data.iter().sum();
        let sum: u32 = grids[0].data.iter().sum();
        let total_length = grids[0].data.len();

        assert_eq!(num_atoms,condense_sum as usize);
        assert_eq!(num_atoms, sum as usize);
        assert_eq!(12*3*5, total_length);
    }


    #[test]
    fn read_in_two_mol() {
        use std::path::Path;
        let file_path = String::from("./test_files/two.mol2");

        let l_mols = 
            read_mol2_file(Path::new(&file_path))
            .expect("Failed to read MOL2 file");
        assert_eq!(
            l_mols.len(),
            2
        );

    }

    #[test]
    fn voxelize_two_mols(){
        use std::path::Path;
        let file_path = String::from("./test_files/two.mol2");

        let l_mols = 
            read_mol2_file(Path::new(&file_path))
            .expect("Failed to read MOL2 file"); 

        let grids = 
            voxelize(
                &l_mols, 
                [12, 3, 5], 
                2.0, 0.0, 0.0, 0.0
            ); 

        for (i,mol) in l_mols.iter().enumerate() {
            let num_atoms = mol.num_atoms();
            let condense_sum: u32 = grids[i].condensed_data.iter().sum();
            let sum: u32 = grids[i].data.iter().sum();
            let total_length = grids[i].data.len();

            assert_eq!(num_atoms,condense_sum as usize);
            assert_eq!(num_atoms, sum as usize);
            assert_eq!(12*3*5, total_length);
        }

    }

    #[test]
    fn get_rec_info(){
        use std::path::Path;

        let file_path = String::from("./test_files/two.mol2");

        let l_mols = 
            read_mol2_file(Path::new(&file_path))
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
}