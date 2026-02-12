use core::panic;
use std::path::Path;
use std::fs::{OpenOptions,File};
use std::io::{self, BufRead, BufReader, BufWriter, Seek, Result, Write};
use std::fs;
use wincode_derive::{SchemaWrite, SchemaRead};
use wincode::{serialize_into};
use clap::ValueEnum;


#[derive(PartialEq, Clone, ValueEnum, Debug)]
pub enum AtomTyping {
    NoType,
    ExplicitType,
    ElementalType
}


#[derive(Clone, Debug, SchemaWrite, SchemaRead)]
pub struct VoxMol {
    pub title: String,
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub atom_types: Option<Vec<String>>
}

impl VoxMol {
    pub fn new(atom_typing: &AtomTyping) -> Self {

        match atom_typing {
            AtomTyping::NoType => {
                return Self {
                    title: String::new(),
                    x: Vec::new(),
                    y: Vec::new(),
                    z: Vec::new(),
                    atom_types: None
                }
            },
            _ => {
                return Self {
                    title: String::new(),
                    x: Vec::new(),
                    y: Vec::new(),
                    z: Vec::new(),
                    atom_types: Some(Vec::new())
                }
            }
        }
    }

    pub fn num_atoms(&self) -> usize {
        if self.x.len() == self.y.len() && self.y.len() == self.z.len() {
            self.x.len()
        } else {
            panic!("Inconsistent atom coordinate lengths");
        }
    }
}


fn create_folder_n_exists(path: &str) -> Result<()> {

    if fs::metadata(path).is_err() {
        fs::create_dir_all(path)?;
    } else {
        println!("Directory '{}' already exists. Writing files to this directory.", path);
        println!("    Ignore if you starting fresh.");
    }

    Ok(())
}

//TODO: make sure that the tmp folder is being written
// the any pre exisiting tmp. if you don't do this,
// it will use the previous file. 

/// Reads a MOL2 file and extracts the title, x, y, z coordinates of heavy atoms, and optionally atom types.
/// Writes each molecule as a binary file in a streaming manner to avoid high memory usage.
/// Compatible with DOCK6 output MOL2 files, which may have a "Name:" as the
/// first header information before the "@<TRIPOS>MOLECULE" section. If "Name:" is present, 
/// it will use the line following "Name:" as the title of the molecule; otherwise, it will use the line following "@<TRIPOS>MOLECULE" as the title.
/// Returns a tuple containing a vector of all unique atom types and the total number of molecules processed.
/// # Arguments
/// * `path` - A reference to the path of the MOL2 file to be read
/// * `atom_typing` - A AtomTyping enum type indicating how to extract atom type information from the MOL2 file. See AtomTyping enum for more details.
///   # Returns
///     * `Result<(Vec<String>, u64)>` - A result containing a tuple of a vector of unique atom types and the total number of molecules, or an error if the file cannot be read
/// # Errors
/// * Returns an error if the specified file cannot be opened or read
/// # Example
/// ```
/// use std::path::Path;
/// let path = Path::new("molecules.mol2");
/// let atom_typing = AtomTyping::ExplicitType; // or AtomTyping::ElementalType, or AtomTyping::NoType
/// let (all_atom_types, num_mols) = read_mol2_file_stream(&path, &atom_typing).expect("Failed to read MOL2 file");
/// println!("Unique atom types: {:?}", all_atom_types);
/// println!("Number of molecules: {}", num_mols);
/// ``` 
pub fn read_mol2_file_stream(
    path: &Path, 
    atom_typing: &AtomTyping
) -> Result<(Vec<String>, u64)> {

    let file = File::open(path)?;
    let reader = io::BufReader::new(file);

    // Create Voxmol instance and a list to hold multiple molecules
    let mut mol = VoxMol::new(atom_typing);
    let mut num_mols: u64 = 0;
    // Flags to track sections in MOL2 file
    let mut in_atom_section: bool = false;
    let mut in_molecule_section: bool = false;
    let mut all_atom_types: Vec<String> = Vec::new();

    create_folder_n_exists("./processed_mols")?;

    let binary_file: Option<File>;

    binary_file = Some(OpenOptions::new()
        .create(true)
        .write(true)
        .open("./processed_mols/mols_stream.binary")?);

    // Read file line by line
    for line in reader.lines() {

        let line = line?;

        if line.starts_with("@<TRIPOS>MOLECULE") {
            in_molecule_section = true;
            continue;
        }

        if in_molecule_section{
            mol.title = line.trim().to_string();
            in_molecule_section = false;
            continue;
        }

        if line.starts_with("@<TRIPOS>ATOM") {
            in_atom_section = true;
            continue;
        }

        if line.starts_with("@<TRIPOS>BOND") {
            in_atom_section = false;
            continue;
        }

        if in_atom_section {
            // Parse atom line to extract x, y, z coordinates
            // Example line format (not actual MOL2 format):
            // 1 C1 0.000 0.000 0.000 C.3
            // Ignores Hs and only takes coordinate of heavy atoms.  

            let parts: Vec<&str> = line.split_whitespace().collect();

            if parts.len() >= 5 && parts[5] != "H" {
                let x: f32 = parts[2].parse().unwrap_or(0.0);
                let y: f32 = parts[3].parse().unwrap_or(0.0);
                let z: f32 = parts[4].parse().unwrap_or(0.0);
                mol.x.push(x);
                mol.y.push(y);
                mol.z.push(z);
                
                // Depending on the atom typing method, 
                // extract the appropriate atom type information and store it in the VoxMol struct.
                match atom_typing {
                    AtomTyping::ExplicitType => 
                    {
                        let at_extype: String =  parts[5]
                            .parse()
                            .unwrap_or(
                                "other".into()
                            );

                        if !all_atom_types.contains(&at_extype) {
                            all_atom_types.push(
                                at_extype.clone()
                            );
                        }

                        mol.atom_types
                        .as_mut()
                        .unwrap()
                        .push(
                            at_extype.clone()
                        );
                    },
                    AtomTyping::ElementalType => 
                    {
                        let at_extype: String =  parts[5]
                            .parse()
                            .unwrap_or(
                                "other".into()
                            );
                        let chars: Vec<char> = at_extype.chars().collect(); 

                        let first_char = chars[0];

                        // NOTE: I don't know know if this is the best way to 
                        // check if the first character is an alphabetic character. 
                        if !first_char.is_ascii_alphabetic() {
                            panic!("Unexpected atom type format: {}", at_extype);
                        }

                        if !all_atom_types.contains(&first_char.to_string()) {
                            all_atom_types.push(
                                first_char.to_string()
                            );
                        }

                        mol.atom_types
                        .as_mut()
                        .unwrap()
                        .push(
                            first_char.to_string()
                        );
                    },
                    AtomTyping::NoType => {}
                }

            }
        }

        if !in_atom_section && !mol.x.is_empty() {

            let mut buffer = Vec::new();
            let _ = serialize_into(&mut buffer, &mol); // fills buffer
            let len = buffer.len() as u32; 
            // Write bytes to file
            binary_file.as_ref().unwrap().write_all(&len.to_le_bytes())?;
            binary_file.as_ref().unwrap().write_all(&buffer)?;
            buffer.clear(); // reuse buffer

            // Reset for next molecule
            mol = VoxMol::new(atom_typing);
            num_mols += 1;
        }

    }   

    Ok((all_atom_types, num_mols))
    
}


/// Writes the cluster molecule IDs to a specified file path in a structured format.
/// If the input vector is already ordered, it will be written in that order.
/// Each molecule ID is written with its corresponding index within the cluster.
/// Each cluster will be separated by a blank line for readability.
/// # Arguments
/// * `path` - A reference to the path where the cluster molecule IDs will be written
/// * `cluster_mol_ids` - A reference to a vector of vectors containing tuples of molecule IDs and their corresponding indices
/// # Returns
/// * `Result<()>` - A result indicating success or failure of the write operation
/// # Errors
/// * Returns an error if the specified file cannot be created or written to
/// # Example
/// ```
/// use std::path::Path;
/// let cluster_mol_ids = vec![vec![("mol1".to_string(), 0)], vec![("mol2".to_string(), 1)]];
/// let path =  Path::new("cluster_mol_ids.txt");
/// write_cluster_mol_ids(&path, &cluster_mol_ids).expect("Failed to write cluster molecule IDs");
/// ```
/// The output file will have the following format:
/// index: 0
/// mol 0: mol1
/// mol 1: mol4
/// index: 1
/// mol 1: mol2
/// Each cluster is separated by a blank line for readability.
pub fn write_cluster_mol_ids(path: &Path, cluster_mol_ids: &Vec<Vec<(String, usize)>>) -> Result<()>  {

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (index, row) in cluster_mol_ids.iter().enumerate() {
        for (i, val) in row.iter().enumerate() {
            if i == 0 {
                write!(writer, "index: {}\n", index)?;
            }
            write!(writer, "mol {}: {}\n", val.1, val.0)?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

fn index_molecules(path: &Path) -> Result<Vec<u64>> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    let mut indices = Vec::new();
    let mut offset = 0u64;

    let mut if_contains_name: bool = false;
    let mut offset_before: u64 = 0u64;

    loop {
        let mut buf = String::new();
        let bytes_read = reader.read_line(&mut buf)?;
        if bytes_read == 0 {
            break;
        }

        if buf.contains("Name:") {
            if_contains_name = true;
            offset_before = offset;
        }

        if buf.contains("@<TRIPOS>MOLECULE") {
            if !if_contains_name {
                indices.push( offset );
            } else {
                indices.push( offset_before );
            }
            if_contains_name = false;
        } 

        offset += bytes_read as u64;
    }

    Ok(indices)
}

fn write_molecule(
    reader: & mut BufReader<File>, 
    writer: & mut BufWriter<File>
) -> std::io::Result<()> {

    loop {
        let mut buf = String::new();
        let bytes = reader.read_line(&mut buf)?;
        if bytes == 0 { break; }

        if buf.ends_with("ROOT\n") {
            writeln!(writer,"{}", buf).expect("Error writing line in write_molecule");
            break;
        }
        write!(writer,"{}", buf).expect("Error writing line in write_molecule");
    }

    Ok(())
}

/// Writes the clustered molecule IDs to separate MOL2 files based on their cluster indices.
/// Each cluster will be written to a separate file named "cluster_{index}.mol2" in the "./molecules" directory.
/// # Arguments
/// * `cluster_mol_ids` - A reference to a vector of vectors containing tuples of molecule IDs and their corresponding indices for each cluster
/// * `path` - A reference to the path of the original MOL2 file containing all molecules
/// * `cluster_write_limit` - A limit on the number of clusters to write (if set to 0, all clusters will be written)
/// # Returns
/// * `Result<()>` - A result indicating success or failure of the write operation
/// # Errors
/// * Returns an error if the specified file cannot be opened, read, or written to
/// # Example
/// ```
/// use std::path::Path;
/// let cluster_mol_ids = vec![vec![("mol1".to_string(), 0)], vec![("mol2".to_string(), 1)]];
/// let path =  Path::new("molecules.mol2");
/// write_mol2s_via_cluster_ind(&cluster_mol_ids, &path, 0).expect("Failed to write clustered MOL2 files");
/// ```
/// This will create two files: "./molecules/cluster_00.mol2" containing "mol1" and "./molecules/cluster_01.mol2" containing "mol2
/// Note: The function assumes that the original MOL2 file is structured in a way that allows for sequential reading of molecules based on their byte offsets, and that each molecule ends with a line containing "ROOT". The `cluster_write_limit` parameter can be used to limit the number of clusters written to files, which can be useful for testing or when dealing with a large number of clusters.
pub fn write_mol2s_via_cluster_ind( 
    cluster_mol_ids :&Vec<Vec<(String, usize)>>,
    path: &Path,
    cluster_write_limit: usize
) -> io::Result<()>
{

    create_folder_n_exists("./molecules")?;

    let last_index = cluster_mol_ids.len() - 1;
    let num_digits = last_index.to_string().len();

    let molecular_file = File::open(path)?;
    let mut reader = BufReader::new(molecular_file);

    let l_mol_byte_indices = index_molecules(&path).expect("Can't open file");

    for (index, row) in cluster_mol_ids.iter().enumerate() {

        if cluster_write_limit == index {
            break;
        }

        let path_file: String = format!(
            "./molecules/cluster_{:0width$}.mol2",
            index,
            width = num_digits
        );

        let writer_file = File::create(&path_file)?;
        let mut writer: BufWriter<File> = BufWriter::new(writer_file); 

        for val in row.iter() {

            reader.seek(io::SeekFrom::Start(l_mol_byte_indices[val.1]))?;

            write_molecule(& mut reader, & mut writer).expect("Error during writing mols");
        }
    }

    Ok(())

}


// NOTE: save this as reference just in case
// pub fn write_mol2s_via_cluster_ind_v2( 
//     cluster_mol_ids :&Vec<Vec<(String, usize)>>,
//     path: &Path,
//     cluster_write_limit: usize
// ) -> io::Result<()>
// {
//     fs::create_dir_all("./molecules")?;
//     let lines: Vec<String> = {
//         let file = File::open(path)?;
//         let reader = io::BufReader::new(file);
//         reader.lines().collect::<io::Result<_>>()?
//     };

//     let last_index = cluster_mol_ids.len() - 1;
//     let num_digits = last_index.to_string().len();
//     let mut l_mol_lines_indices = Vec::<usize>::new();

//     for (index, row) in cluster_mol_ids.iter().enumerate() {
        
//         if cluster_write_limit == index {
//             break;
//         }

//         let path_file: String = format!(
//             "./molecules/cluster_{:0width$}.mol2",
//             index,
//             width = num_digits
//         );
//         let file = File::create(path_file)?;
//         let mut writer = BufWriter::new(file); 

//         // Vector to store matching lines or strings
//         let mut l_titles: Vec<String> = Vec::new();

//         for (line_num, line) in lines.iter().enumerate() {
//             if line.contains("Name:") {
//                 let parts: Vec<&str> = 
//                     line.split_whitespace().collect();
//                 l_titles.push(parts[2].to_string() );
//                 l_mol_lines_indices.push(line_num);
//             }
//         }

//         for val in row.iter() {

//             let start_line_num: usize = l_mol_lines_indices[val.1];

//             for line in &lines[start_line_num..] {
//                 if line.ends_with("ROOT") {
//                     write!(writer,"{}\n\n", line)?;
//                     break;
//                 }
//                 write!(writer,"{}\n", line)?;
//                 continue;
//             }
//         }
//     }

//     Ok(())

// }
