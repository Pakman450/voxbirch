use crate::file_io::VoxMol;
use std::io::Write;
use std::fs::{OpenOptions,File};
use wincode_derive::{SchemaWrite, SchemaRead};
use wincode::{deserialize, serialize_into};
use std::io::{BufReader, Read, Result};

#[derive(Clone, SchemaWrite, SchemaRead)]
pub struct VoxelGrid {
    pub title : String,
    pub dims: [usize; 3],
    pub data: Vec<u8>
}

impl VoxelGrid {
    pub fn new(
        dims: [usize; 3], 
        atom_typing: bool,
        num_atom_types: usize
    ) -> Self {

        let size: usize = if atom_typing 
        {
            // Note: + 1 means "other" atom_types
            dims[0] * dims[1] * dims[2] * (num_atom_types + 1)
        } else {
            dims[0] * dims[1] * dims[2]
        };

        Self {
            title: String::new(),
            dims,
            data: vec![0; size]
        }
    }

    // This function computes the 1D index in the voxel grid from 3D coordinates
    pub fn voxel_index(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.dims[0] * (y + self.dims[1] * z)
    }

    fn push_condensed_data(&mut self, datum: u8) {
        self.data.push(datum);
    }
}

fn get_atom_typing_index(
    ind: usize, 
    all_atom_types: &Vec<String>,
    atom_type: &String,
) -> usize {

    if let Some(shifted_idx) = all_atom_types.iter().position(|x| x == atom_type)
    {
        return ind + shifted_idx
    }

    ind + all_atom_types.len()

}

pub fn voxelize_stream(
    dims: [usize; 3],
    rs: f32,
    x0: f32,
    y0: f32,
    z0: f32,
    atom_typing: bool,
    all_atom_types: Vec<String>,
    stdout: & mut Box<dyn Write>
) -> Result<(
    u64, 
    u64,
    usize,
    Vec<u32>
)> {

    // create a list of voxel grids based on number of molecules

    let mut num_rows: u64 = 0;
    let mut num_cols: u64 = 0;

    let mut binary_file = OpenOptions::new()
        .create(true)
        .write(true)
        .open("./tmp/grids_stream.binary.tmp")?;

    // --- STREAM READ ---
    let read_binary_file = File::open("./tmp/mols_stream.binary.tmp")?;
    let mut reader = BufReader::new(read_binary_file);

    let mut condense_data_idx = Vec::<u32>::new();

    let mut occupied_voxels = 0;
    let mut voxel_sum = 0;
    let mut num_heavy_atoms = 0;
    let num_atom_types = all_atom_types.len();

    loop {
        // Read length prefix (4 bytes)
        let mut len_buf = [0u8; 4];
        match reader.read_exact(&mut len_buf) {
            Ok(_) => {},
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => break, // EOF
            Err(e) => return Err(e.into()),
        }

        let len = u32::from_le_bytes(len_buf) as usize;

        // Read exactly `len` bytes
        let mut record_buf = vec![0u8; len];
        reader.read_exact(&mut record_buf)?;

        // Deserialize the struct
        let mol: VoxMol = deserialize(&record_buf).unwrap();

        let mut grid = VoxelGrid::new(dims, atom_typing, num_atom_types);

        grid.title = mol.title.clone();

        for atom_idx in 0..mol.num_atoms() {
            let x = mol.x[atom_idx];
            let y = mol.y[atom_idx];
            let z = mol.z[atom_idx];


            let ix = ((x - x0) / rs).floor() as usize;
            let iy = ((y - y0) / rs).floor() as usize;
            let iz = ((z - z0) / rs).floor() as usize;
        
            let mut index = grid.voxel_index(ix, iy, iz);

            if atom_typing {
                index = get_atom_typing_index(
                 index, 
                 &all_atom_types,
                 &mol.atom_types.as_ref().unwrap()[atom_idx]
                );
            }
            // TODO:now I want the indexes around the origin of ix,iy,iz 
            // to spread the voxelization to neighboring voxels
            // for simplicity, we will just do a 3x3x3 cube around the voxel
            // In a real application, you might want to use a more sophisticated method
            // such as a Gaussian spread or a sphere of influence
            // but for now, let's just increment the neighboring voxels

            num_heavy_atoms += 1;
            grid.data[index] += 1;
            voxel_sum += 1;
            if grid.data[index] == 1 {
                occupied_voxels += 1;
                if !condense_data_idx.contains(&(index as u32)) {
                    condense_data_idx.push(index as u32);
                }
                
            }

        }
        num_cols = grid.data.len() as u64;
        // --- STREAM WRITE ---
        let mut buffer = Vec::new();
        let _ = serialize_into(&mut buffer, &grid); // fills buffer
        let len = buffer.len() as u32; 
        // Write bytes to file
        binary_file.write_all(&len.to_le_bytes())?;
        binary_file.write_all(&buffer)?;
        buffer.clear(); // reuse buffer
        num_rows += 1;
        
    }
    let total_voxels = if atom_typing {
            dims[0] * dims[1] * dims[2] * ( num_atom_types + 1 ) * num_rows as usize
    } else {
            dims[0] * dims[1] * dims[2] * num_rows as usize
    };

    writeln!(
        stdout,
        "Total number of atom types per voxel: {}",
        all_atom_types.len() + 1
    ).unwrap(); 
    writeln!(
        stdout,
        "Total number of voxels =\n\tTotal number of occupied voxels + Total number of unoccupied voxels:\n\t{} = {} + {}",
        total_voxels, occupied_voxels, total_voxels - occupied_voxels
    ).unwrap();  
    writeln!(
        stdout,
        "Total number of occupied voxels for one grid without condensation: {}", 
        num_cols
    ).unwrap();
    writeln!(
        stdout,
        "Total number of occupied voxels across all grids: = {}", 
        condense_data_idx.len()
    ).unwrap();
    writeln!(
        stdout,
        "Total number of heavy atoms voxelized: {}", num_heavy_atoms
    ).unwrap();
    writeln!(
        stdout,
        "Sum of all voxel values: {}", voxel_sum
    ).unwrap();

    Ok((num_rows, num_cols as u64, condense_data_idx.len(), condense_data_idx))
}


pub fn condense_data_stream(condensed_data_idx: Vec<u32>)
-> Result<()>{

    let mut binary_file = OpenOptions::new()
        .create(true)
        .write(true)
        .open("./tmp/grids_condensed_stream.binary.tmp")?;

    let read_binary_file = File::open("./tmp/grids_stream.binary.tmp");
    let mut reader = BufReader::new(read_binary_file.unwrap());
    
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
        let mut grid: VoxelGrid = deserialize(&record_buf).unwrap();

        let vec_data_grid = grid.data.clone();

        grid.data = Vec::<u8>::new();
        
        for &ind in & condensed_data_idx {

            grid.push_condensed_data(vec_data_grid[ind as usize]);
        }

        drop(vec_data_grid);

        // --- STREAM WRITE ---
        let mut buffer = Vec::new();
        let _ = serialize_into(&mut buffer, &grid); // fills buffer
        let len = buffer.len() as u32; 
        // Write bytes to file
        binary_file.write_all(&len.to_le_bytes())?;
        binary_file.write_all(&buffer)?;
        buffer.clear(); // reuse buffer

    }


    Ok(())
}

// pub fn voxelize(
//     l_mols: &Vec<VoxMol>,
//     dims: [usize; 3],
//     rs: f32,
//     x0: f32,
//     y0: f32,
//     z0: f32,
//     atom_typing: bool,
//     all_atom_types: &Vec<String>,
//     stdout: & mut Box<dyn Write>
// ) -> Result<Vec::<VoxelGrid>> {

//     // create a list of voxel grids based on number of molecules

//     let mut grids = Vec::<VoxelGrid>::new();
//     let mut num_cols: u64 = 0;

//     let mut occupied_voxels = 0;
//     let mut voxel_sum = 0;
//     let mut num_heavy_atoms = 0;
//     let mut condense_data_idx = Vec::<u32>::new();
//     let num_atom_types = all_atom_types.len();

//     for mol in l_mols.iter() {

//         let mut grid = VoxelGrid::new(dims, atom_typing, num_atom_types);

//         grid.title = mol.title.clone();

//         for atom_idx in 0..mol.num_atoms() {
//             let x = mol.x[atom_idx];
//             let y = mol.y[atom_idx];
//             let z = mol.z[atom_idx];


//             let ix = ((x - x0) / rs).floor() as usize;
//             let iy = ((y - y0) / rs).floor() as usize;
//             let iz = ((z - z0) / rs).floor() as usize;
        
//             let mut index = grid.voxel_index(ix, iy, iz);

//             if atom_typing {
//                 index = get_atom_typing_index(
//                  index, 
//                  &all_atom_types,
//                  &mol.atom_types.as_ref().unwrap()[atom_idx]
//                 );
//             }
//             // TODO:now I want the indexes around the origin of ix,iy,iz 
//             // to spread the voxelization to neighboring voxels
//             // for simplicity, we will just do a 3x3x3 cube around the voxel
//             // In a real application, you might want to use a more sophisticated method
//             // such as a Gaussian spread or a sphere of influence
//             // but for now, let's just increment the neighboring voxels

//             num_heavy_atoms += 1;
//             grid.data[index] += 1;
//             voxel_sum += 1;

//             if grid.data[index] == 1 {
//                 occupied_voxels += 1;
//                 if !condense_data_idx.contains(&(index as u32)) {
//                     condense_data_idx.push(index as u32);
//                 }
                
//             }

//         }

//         grids.push(grid);
//     }

//     num_cols = grids[0].data.len() as u64;

//     for grid in & mut grids{
//         let vec_data_grid = grid.data.clone();
        
//         grid.data.clear();

//         for &ind in & condense_data_idx {

//             grid.push_condensed_data(vec_data_grid[ind as usize]);
//         }
        
//     }

//     let total_voxels = if atom_typing {
//             dims[0] * dims[1] * dims[2] * ( num_atom_types + 1 ) * grids.len()
//     } else {
//             dims[0] * dims[1] * dims[2] * grids.len()
//     };

//     writeln!(
//         stdout,
//         "Total number of atom types per voxel: {}",
//         all_atom_types.len() + 1
//     ).unwrap(); 
//     writeln!(
//         stdout,
//         "Total number of voxels =\n\tTotal number of occupied voxels + Total number of unoccupied voxels:\n\t{} = {} + {}",
//         total_voxels, occupied_voxels, total_voxels - occupied_voxels
//     ).unwrap();    
//     writeln!(
//         stdout,
//         "Total number of occupied voxels for one grid without condensation: {}", 
//         num_cols
//     ).unwrap();
//     writeln!(
//         stdout,
//         "Total number of occupied voxels across all grids: = {}", 
//         condense_data_idx.len()
//     ).unwrap();
//     writeln!(
//         stdout,
//         "Total number of heavy atoms voxelized: {}", num_heavy_atoms
//     ).unwrap();
//     writeln!(
//         stdout,
//         "Sum of all voxel values: {}", voxel_sum
//     ).unwrap();

//     Ok(grids)
// }


// pub fn get_recommended_info(l_mols: &Vec<VoxMol>, resolution: f32, x0: f32, y0: f32, z0: f32) -> (
//     f32, f32, f32, usize, usize, usize, usize, usize, usize
// ) {
//     let mut l_vals = Vec::<f32>::new();

//     for mol in l_mols {
//         l_vals.extend(mol.x.clone())
//     }

//     let min_x = l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()).expect("No valid minimum found (possibly due to NaN values).");

//     l_vals.clear();

//     for mol in l_mols {
//         l_vals.extend(mol.y.clone())
//     }

//     let min_y = l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()).expect("No valid minimum found (possibly due to NaN values).");

//     l_vals.clear();
//     for mol in l_mols {
//         l_vals.extend(mol.z.clone())
//     }
//     let min_z = l_vals.iter().cloned().filter(|&x| !x.is_nan()).min_by(|a, b| a.partial_cmp(b).unwrap()).expect("No valid minimum found (possibly due to NaN values).");
//         // Compute maximum coordinates across all molecules so we can compute
//     // the minimal dimensions required to cover them from the recommended origin
//     let mut xs: Vec<f32> = Vec::new();
//     let mut ys: Vec<f32> = Vec::new();
//     let mut zs: Vec<f32> = Vec::new();

//     for mol in l_mols {
//         xs.extend(mol.x.iter().cloned());
//         ys.extend(mol.y.iter().cloned());
//         zs.extend(mol.z.iter().cloned());
//     }

//     let max_x = xs.iter().cloned().filter(|v| !v.is_nan()).max_by(|a, b| 
//         a.partial_cmp(b).unwrap()).expect("No valid maximum found (possibly due to NaN values)."
//     );
//     let max_y = ys.iter().cloned().filter(|v| !v.is_nan()).max_by(|a, b| 
//         a.partial_cmp(b).unwrap()).expect("No valid maximum found (possibly due to NaN values)."
//     );
//     let max_z = zs.iter().cloned().filter(|v| !v.is_nan()).max_by(|a, b| 
//         a.partial_cmp(b).unwrap()).expect("No valid maximum found (possibly due to NaN values)."
//     );

//     // span and required voxels (inclusive of boundary). Use floor(span / resolution) + 1
//     let span_x = (max_x - min_x).max(0.0);
//     let span_y = (max_y - min_y).max(0.0);
//     let span_z = (max_z - min_z).max(0.0);

//     let need_x = (span_x / resolution).floor() as usize + 1;
//     let need_y = (span_y / resolution).floor() as usize + 1;
//     let need_z = (span_z / resolution).floor() as usize + 1;

//     // Also compute required dims if using the user-provided origin (x0,y0,z0)
//     let span_x_user = (max_x - x0).max(0.0);
//     let span_y_user = (max_y - y0).max(0.0);
//     let span_z_user = (max_z - z0).max(0.0);

//     let need_x_user = (span_x_user / resolution).floor() as usize + 1;
//     let need_y_user = (span_y_user / resolution).floor() as usize + 1;
//     let need_z_user = (span_z_user / resolution).floor() as usize + 1;

//     (
//         min_x, min_y, min_z, // minimum coordinates
//         need_x, need_y, need_z, // required dims from absolute origin
//         need_x_user, need_y_user, need_z_user // required dims from user-provided origin
//     )
// }




pub fn get_recommended_info_stream(resolution: f32, x0: f32, y0: f32, z0: f32) ->
    (f32, f32, f32, usize, usize, usize, usize, usize, usize) {
        
    // --- STREAM READ ---
    // let mut min_xyz = ( 0.0, 0.0, 0.0 );
    // let mut max_xyz = ( 0.0, 0.0, 0.0 );
    
    let mut min_xyz = (f32::INFINITY, f32::INFINITY, f32::INFINITY);
    let mut max_xyz = (f32::NEG_INFINITY, f32::NEG_INFINITY, f32::NEG_INFINITY);

    let read_binary_file = File::open("./tmp/mols_stream.binary.tmp");
    let mut reader = BufReader::new(read_binary_file.unwrap());
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
        let mol: VoxMol = deserialize(&record_buf).unwrap();

        for atom_idx in 0..mol.num_atoms() {

            if mol.x[atom_idx] < min_xyz.0 {
                min_xyz.0 = mol.x[atom_idx];
            }
            if mol.y[atom_idx] < min_xyz.1 {
                min_xyz.1 = mol.y[atom_idx];
            }
            if mol.z[atom_idx] < min_xyz.2 {
                min_xyz.2 = mol.z[atom_idx];
            }

            if mol.x[atom_idx] > max_xyz.0 {
                max_xyz.0 = mol.x[atom_idx];
            }
            if mol.y[atom_idx] > max_xyz.1 {
                max_xyz.1 = mol.y[atom_idx];
            }
            if mol.z[atom_idx] > max_xyz.2 {
                max_xyz.2 = mol.z[atom_idx];
            }
        }
    }

    let max_x = max_xyz.0;
    let max_y = max_xyz.1;
    let max_z = max_xyz.2;
    let min_x = min_xyz.0;
    let min_y = min_xyz.1;
    let min_z = min_xyz.2;

    
    // span and required voxels (inclusive of boundary). Use floor(span / resolution) + 1
    let span_x = (max_x - min_x).max(0.0);
    let span_y = (max_y - min_y).max(0.0);
    let span_z = (max_z - min_z).max(0.0);

    let need_x = (span_x / resolution).floor() as usize + 1;
    let need_y = (span_y / resolution).floor() as usize + 1;
    let need_z = (span_z / resolution).floor() as usize + 1;

    // Also compute required dims if using the user-provided origin (x0,y0,z0)
    let span_x_user = (max_x - x0).max(0.0);
    let span_y_user = (max_y - y0).max(0.0);
    let span_z_user = (max_z - z0).max(0.0);

    let need_x_user = (span_x_user / resolution).floor() as usize + 1;
    let need_y_user = (span_y_user / resolution).floor() as usize + 1;
    let need_z_user = (span_z_user / resolution).floor() as usize + 1;
    
    (
        min_x, min_y, min_z, // minimum coordinates
        need_x, need_y, need_z, // required dims from absolute origin
        need_x_user, need_y_user, need_z_user // required dims from user-provided origin
    )
}



