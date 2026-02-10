

#[pyo3::pymodule]
mod vb_py {
    use pyo3::prelude::*;
    use pyo3::ffi;
    use pyo3::{Python, PyErr};
    use pyo3::exceptions::PyTypeError;

    use voxbirch::{
        read_mol2_file_stream,
        get_recommended_info_stream,
        voxelize_stream,
        condense_data_stream
    };
    use voxbirch::birch::VoxBirch;
    use voxbirch::VoxelGrid;

    use std::path::Path;
    use std::io::{Write,Read};
    use std::fs::{self};
    use std::panic;

    use wincode::{deserialize};

    #[pyclass(unsendable)]
    struct Voxbirch {
        inner: VoxBirch,  // just holds the Rust struct
    }

    #[pymethods]
    impl Voxbirch {
        #[new]
        fn new(
            threshold: f32,
            max_branches: usize,
            num_features: usize
        ) -> Self {
            Voxbirch {
                inner: VoxBirch::new(
                    threshold, 
                    max_branches, 
                    num_features
                )
            }
        }

        fn insert(& mut self, data: Vec<f32>, title: String, iter: u64) {
            let mut stdsink: Box<dyn Write> = Box::new(std::io::sink());

            let gil_state = unsafe { ffi::PyEval_SaveThread() };
            
            let result = self.inner.insert(
                Some(&data), 
                &title, 
                iter, 
                &mut stdsink
            );  

            unsafe { ffi::PyEval_RestoreThread(gil_state) };

            match result {
                Ok(v) => v,
                Err(e) => panic!("Error during insertion: {}", e)
            };
        }
        
        fn cluster(& mut self) -> PyResult<()> {
            let read_binary_file =
                std::fs::File::open("./tmp/grids_condensed_stream.binary.tmp")?;
            let mut reader = std::io::BufReader::new(read_binary_file);

            let mut iter: u64 = 0;

            loop {
                // --- read record (pure Rust) ---
                let mut len_buf = [0u8; 4];
                if reader.read_exact(&mut len_buf).is_err() {
                    break;
                }

                let len = u32::from_le_bytes(len_buf) as usize;
                let mut record_buf = vec![0u8; len];
                reader.read_exact(&mut record_buf)?;

                let grid: VoxelGrid = deserialize(&record_buf).unwrap();
                let vec_grid: Vec<f32> =
                    grid.data.into_iter().map(|x| x as f32).collect();

                // --- PRINT (GIL is held here) ---
                let _ = Python::try_attach(|_| -> PyResult<()> {
                    println!("Inserting {}: {}", grid.title, iter + 1);
                    Ok(())
                });

                // --- HEAVY WORK (NO GIL, but Restored after) ---
                self.insert(vec_grid, grid.title, iter);

                iter += 1;
            }

            Ok(())
        }

    }

    #[pyfunction]
    #[pyo3(signature = (path_str))]
    fn read_mol2(     
        path_str: String
    ) -> PyResult<(Vec<String>, u64)> {

        let path = Path::new(&path_str);
        // --- RELEASE GIL ---
        let gil_state = unsafe { ffi::PyEval_SaveThread() };

        let result = read_mol2_file_stream(path, true);

        // --- RE-ACQUIRE GIL ---
        unsafe { ffi::PyEval_RestoreThread(gil_state) };

        let (all_atom_types, num_mols) = match result {
            Ok(v) => v,
            Err(e) => panic!("Error during reading in file: {}", e)
        };

        Ok((all_atom_types, num_mols))
    }

    #[pyfunction]
    #[pyo3(signature = (resolution=1.4142135))]
    fn get_optim_ori_dims(   
        resolution: f32 
    ) -> PyResult<(
        [f32; 3],
        [usize; 3]
    )> {


        if !fs::metadata("./tmp").is_ok() {
            return Err(
                PyErr::new::<PyTypeError, _>(
                    "Temporary directory './tmp' not found.\nPlease ensure the previous step of reading the mol file was completed successfully\nand that the './tmp' directory is present."
                ));
        }

        let gil_state = unsafe { ffi::PyEval_SaveThread() };

        let result = panic::catch_unwind(|| {
            let (
                min_x, min_y, min_z,
                need_x, need_y, need_z,
                need_x_user, need_y_user, need_z_user
            ) = get_recommended_info_stream(resolution, 0.0, 0.0, 0.0);

            let (
                _, _, _, _, _, _,
                rec_need_x_user, rec_need_y_user, rec_need_z_user
            ) = get_recommended_info_stream(
                resolution,
                min_x.floor(),
                min_y.floor(),
                min_z.floor()
            );

            (
                min_x, min_y, min_z,
                need_x, need_y, need_z,
                need_x_user, need_y_user, need_z_user,
                rec_need_x_user, rec_need_y_user, rec_need_z_user
            )
        });

        unsafe { ffi::PyEval_RestoreThread(gil_state) };

        let (
            min_x, min_y, min_z,
            need_x, need_y, need_z,
            need_x_user, need_y_user, need_z_user,
            rec_need_x_user, rec_need_y_user, rec_need_z_user
        ) = match result {
            Ok(v) => v,
            Err(e) => {
                let msg = if let Some(s) = e.downcast_ref::<&str>() {
                    *s 
                } else if let Some(s) = e.downcast_ref::<String>() {
                    s.as_str()
                } else {
                    "Unknown panic"
                };
                panic!("Rust panic caught: {}", msg);
            }
        };

        
        Ok(
            (
                [min_x.floor(), min_y.floor(), min_z.floor()],
                [rec_need_x_user, rec_need_y_user, rec_need_z_user]
            )
        )
    }

    #[pyfunction]
    #[pyo3(signature = (dims=[100,100,100], resolution=1.4142135, origin=[0.0,0.0,0.0], all_atom_types=Vec::new()))]
    fn voxelize(
        dims: [usize; 3],
        resolution: f32,
        origin: [f32; 3],
        all_atom_types: Vec<String>
    ) -> PyResult<(
        u64,
        u64,
        usize,
        Vec<u32>
    )> {

        let atom_typing = !all_atom_types.is_empty();

        if !atom_typing {
            return Err(PyErr::new::<PyTypeError, _>(
                "Provide `all_atom_types` for voxelization."
            ));
        }

        let mut stdsink: Box<dyn Write> = Box::new(std::io::sink());
        
        // --- RELEASE GIL ---
        let gil_state = unsafe { ffi::PyEval_SaveThread() };

        let result = voxelize_stream(
                dims, 
                resolution, 
                origin[0], origin[1], origin[2],
                atom_typing,
                all_atom_types,
                &mut stdsink
            );

        // --- RE-ACQUIRE GIL ---
        unsafe { ffi::PyEval_RestoreThread(gil_state) };

        let (
            num_rows,
            num_cols,
            num_condensed_cols,
            condensed_data_idx
            ) = match result {
            Ok(v) => v,
            Err(e) => panic!("Error during reading in file: {}", e)
        };

        Ok((
            num_rows,
            num_cols,
            num_condensed_cols,
            condensed_data_idx
        ))
    }

    #[pyfunction]
    #[pyo3(signature = (condensed_data_idx=None))]
    fn condense_grids(
        condensed_data_idx: Option<Vec<u32>>
    ) -> PyResult<()> {

        // --- RELEASE GIL ---
        let gil_state = unsafe { ffi::PyEval_SaveThread() };

        match &condensed_data_idx {
            Some(condensed_data_idx) => {
                condense_data_stream(
                    condensed_data_idx.to_vec()
                ).expect("Error during condensation");
            },
            None => {
                println!("No condensed data index provided, skipping condensation.");
            }
        }

        // --- RE-ACQUIRE GIL ---
        unsafe { ffi::PyEval_RestoreThread(gil_state) };

        Ok(())
    }
}
