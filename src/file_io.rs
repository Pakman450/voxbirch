use std::path::Path;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write, Result};

#[derive(Clone)]
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
pub fn read_mol2_file(
    path: &Path, atom_typing: bool
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
            l_mols.push(mol);
            // Reset for next molecule
            mol = VoxMol::new(atom_typing);
        }

    }   

    Ok((l_mols,all_atom_types))
    
}

pub fn write_cluster_mol_ids(path: &Path, cluster_mol_ids: &Vec<Vec<String>>) -> Result<()>  {

    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    for (index, row) in cluster_mol_ids.iter().enumerate() {
        for (i, val) in row.iter().enumerate() {
            if i == 0 {
                write!(writer, "index: {}\n", index)?;
            }
            write!(writer, "mol: {}\n", val)?;
        }
        writeln!(writer)?;
    }

    Ok(())


}