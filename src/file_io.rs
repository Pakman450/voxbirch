use core::panic;
use std::path::Path;
use std::fs::{OpenOptions,File};
use std::io::{self, BufRead, BufWriter, Result, Write};
use std::fs;
use wincode_derive::{SchemaWrite, SchemaRead};
use wincode::{serialize_into};

#[derive(Clone, Debug, SchemaWrite, SchemaRead)]
pub struct VoxMol {
    pub title: String,
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub atom_types: Option<Vec<String>>
}

impl VoxMol {
    pub fn new(atom_typing: bool) -> Self {



        if atom_typing {
            return Self {
                title: String::new(),
                x: Vec::new(),
                y: Vec::new(),
                z: Vec::new(),
                atom_types: Some(Vec::new())
            }
        } else {
            return Self {
                title: String::new(),
                x: Vec::new(),
                y: Vec::new(),
                z: Vec::new(),
                atom_types: None
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

//

//TODO: make sure that the tmp folder is being written
// the any pre exisiting tmp. if you don't do this,
// it will use the previous file. 
pub fn read_mol2_file_stream(
    path: &Path, 
    atom_typing: bool
) -> Result<(Vec<String>, u64)> {

    let file = File::open(path)?;
    let reader = io::BufReader::new(file);

    // Create Voxmol instance and a list to hold multiple molecules
    let mut mol = VoxMol::new(atom_typing);
    let mut num_mols: u64 = 0;
    // Flags to track sections in MOL2 file
    let mut in_atom_section: bool = false;

    let mut all_atom_types: Vec<String> = Vec::new();

    fs::create_dir_all("./tmp")?;

    let mut binary_file: Option<File> = None;

    binary_file = Some(OpenOptions::new()
        .create(true)
        .write(true)
        .open("./tmp/mols_stream.binary.tmp")?);


    // Read file line by line
    for line in reader.lines() {

        let line = line?;

        if line.starts_with("@<TRIPOS>ATOM") {
            in_atom_section = true;
            continue;
        }

        if line.starts_with("@<TRIPOS>BOND") {
            in_atom_section = false;
            continue;
        }

        if line.contains("Name:") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                mol.title = parts[2].to_string();
            }
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
                
                if atom_typing {
                    let at_type =  parts[5]
                        .parse()
                        .unwrap_or(
                            "other".into()
                        );

                    if !all_atom_types.contains(&at_type) {
                        all_atom_types.push(
                            at_type.clone()
                        );
                    }

                    mol.atom_types
                    .as_mut()
                    .unwrap()
                    .push(
                        at_type.clone()
                    );
                }
            }
        }

        if !in_atom_section && !mol.x.is_empty() {

            // --- STREAM WRITE ---
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



pub fn read_mol2_file(
    path: &Path, 
    atom_typing: bool
) -> Result<(Vec<VoxMol>,Vec<String>)> {

    let file = File::open(path)?;
    let reader = io::BufReader::new(file);

    // Create Voxmol instance and a list to hold multiple molecules
    let mut mol = VoxMol::new(atom_typing);
    let mut l_mols: Vec<VoxMol> = Vec::new();

    // Flags to track sections in MOL2 file
    let mut in_atom_section: bool = false;

    let mut all_atom_types: Vec<String> = Vec::new();


    // Read file line by line
    for line in reader.lines() {

        let line = line?;

        if line.starts_with("@<TRIPOS>ATOM") {
            in_atom_section = true;
            continue;
        }

        if line.starts_with("@<TRIPOS>BOND") {
            in_atom_section = false;
            continue;
        }

        if line.contains("Name:") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                mol.title = parts[2].to_string();
            }
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
                
                if atom_typing {
                    let at_type =  parts[5]
                        .parse()
                        .unwrap_or(
                            "other".into()
                        );

                    if !all_atom_types.contains(&at_type) {
                        all_atom_types.push(
                            at_type.clone()
                        );
                    }

                    mol.atom_types
                    .as_mut()
                    .unwrap()
                    .push(
                        at_type.clone()
                    );
                }
            }
        }

        if !in_atom_section && !mol.x.is_empty() {
            // Finished reading one molecule's atoms
            l_mols.push(mol.clone());
            // Reset for next molecule
            mol = VoxMol::new(atom_typing);
        }

    }   

    Ok((l_mols,all_atom_types))
    
}

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

pub fn write_mol2s_via_cluster_ind( 
    cluster_mol_ids :&Vec<Vec<(String, usize)>>,
    path: &Path,
    cluster_write_limit: usize
) -> io::Result<()>
{
    fs::create_dir_all("./molecules")?;
    let lines: Vec<String> = {
        let file = File::open(path)?;
        let reader = io::BufReader::new(file);
        reader.lines().collect::<io::Result<_>>()?
    };

    let last_index = cluster_mol_ids.len() - 1;
    let num_digits = last_index.to_string().len();

    for (index, row) in cluster_mol_ids.iter().enumerate() {
        
        if cluster_write_limit == index {
            break;
        }

        let path_file: String = format!(
            "./molecules/cluster_{:0width$}.mol2",
            index,
            width = num_digits
        );
        let file = File::create(path_file)?;
        let mut writer = BufWriter::new(file); 

        // Vector to store matching lines or strings
        let mut l_titles: Vec<String> = Vec::new();

        for line in &lines {
            if line.contains("Name:") {
                let parts: Vec<&str> = 
                    line.split_whitespace().collect();
                l_titles.push(parts[2].to_string() );
            }
        }

        for val in row.iter() {

            // Flags to track sections in MOL2 file
            let mut in_mol_section: bool = false;
            let mut index_counter: usize = 0;

            for line in &lines {

                if line.contains("Name:") {
                    
                    let parts: Vec<&str> = 
                        line.split_whitespace().collect();
                    if parts.len() >= 2 && 
                        val.0 == parts[2].to_string() 
                        && val.1 == index_counter
                    {
                        in_mol_section = true;
                        
                    }
                    index_counter += 1;
                }

                if in_mol_section {
                    if line.ends_with("ROOT") {
                        write!(writer,"{}\n\n", line)?;
                        break;
                    }
                    write!(writer,"{}\n", line)?;

                    continue;
                }
            }

        }
    }

    Ok(())

}