

#[pyo3::pymodule]
mod vb_py {
    use pyo3::prelude::*;
    use voxbirch::read_mol2_file_stream;
    use std::path::Path;

    #[pyfunction]
    fn read_mol2(     
        path_str: String, 
        atom_typing: bool
    ) -> PyResult<(Vec<String>, u64)> {

        let path = Path::new(&path_str);
        Ok(read_mol2_file_stream(path, atom_typing)?)

    }
}
