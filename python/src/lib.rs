

#[pyo3::pymodule]
mod vb_py {
    use pyo3::prelude::*;
    use pyo3::ffi;
        
    use voxbirch::{
        read_mol2_file_stream,
        get_recommended_info_stream,
        voxelize_stream,
        condense_data_stream
    };
    use voxbirch::birch::VoxBirch;
    use voxbirch::VoxelGrid;

    use std::path::Path;
    use std::io::{Write,BufReader,Read};
    use std::fs::{File, self};
    use wincode::{deserialize};



// NOTE: THIS WORKS FOR PRINTING TO JUPYTER NOTEBOOK,
// BUT NOT FOR STREAMING OUTPUT TO A FILE (E.G. FOR LOGGING PURPOSES)
#[pyfunction]
fn long_running(py: Python) -> PyResult<()> {
    let sys = py.import("sys")?;
    let stdout = sys.getattr("stdout")?;

    for i in 0..5 {
        stdout.call_method1("write", (format!("Step {}\n", i),))?;
        stdout.call_method0("flush")?;

        std::thread::sleep(std::time::Duration::from_secs(1));
    }

    Ok(())
}

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

        fn insert(&mut self, data: Vec<f32>, title: String, iter: u64) {
            let mut stdout: Box<dyn Write> = Box::new(std::io::stdout().lock());
            self.inner.insert(
                Some(&data), 
                &title, 
                iter, 
                &mut stdout
            );
        }
        
        fn cluster(&mut self) -> PyResult<()> {
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
                Python::try_attach(|py| -> PyResult<()> {
                    let sys = py.import("sys")?;
                    let stdout = sys.getattr("stdout")?;
                    stdout.call_method1(
                        "write",
                        (format!("\nInserting {}: {}", grid.title, iter + 1),),
                    )?;
                    stdout.call_method0("flush")?;
                    Ok(())
                }).unwrap();

                // --- RELEASE GIL ---
                let gil_state = unsafe { ffi::PyEval_SaveThread() };

                // --- HEAVY RUST WORK (NO GIL) ---
                self.insert(vec_grid, grid.title, iter);

                // --- RE-ACQUIRE GIL ---
                unsafe { ffi::PyEval_RestoreThread(gil_state) };

                iter += 1;
            }

            Ok(())
        }

    }

    #[pyfunction]
    #[pyo3(signature = (path_str, atom_typing=true))]
    fn read_mol2(     
        path_str: String, 
        atom_typing: bool
    ) -> PyResult<(Vec<String>, u64)> {

        let path = Path::new(&path_str);

        Ok(read_mol2_file_stream(path, atom_typing)?)
    }

    #[pyfunction]
    #[pyo3(signature = (resolution=1.4142135, x0=0.0, y0=0.0, z0=0.0))]
    fn get_optim_ori_dims(   
        resolution: f32,
        x0: f32,
        y0: f32,
        z0: f32,  
    ) -> PyResult<(
        [f32; 3],
        [usize; 3]
    )> {


        if !fs::metadata("./tmp").is_ok() {
            panic!("Temporary directory './tmp' does not exist.\n Please ensure that the previous step of reading the mol file was completed successfully\n and that the './tmp' directory is present.");
        }

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
                resolution, 
                x0, 
                y0, 
                z0
            );
        println!("The recommended origin: {},{},{}", min_x.floor(), min_y.floor(), min_z.floor());

        println!(
            "Minimal voxel grid dimensions to cover all molecules\n\t(from absolute origin {}): {},{},{}",
            format!("{:.3},{:.3},{:.3}", min_x, min_y, min_z), need_x, need_y, need_z
        );
        println!(
            "Required dims from provided origin ({:.3},{:.3},{:.3}): {},{},{}", 
            x0, y0, z0, 
            need_x_user, need_y_user, need_z_user
        );
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
        ) = get_recommended_info_stream(
            resolution, 
            min_x.floor(), 
            min_y.floor(), 
            min_z.floor()
        );


        println!(
            "Required dims from recommended origin ({},{},{}): {},{},{}", 
            min_x.floor(), min_y.floor(), min_z.floor(),
            rec_need_x_user, rec_need_y_user, rec_need_z_user
        );

        Ok(
            (
                [min_x.floor(), min_y.floor(), min_z.floor()],
                [rec_need_x_user, rec_need_y_user, rec_need_z_user]
            )
        )
    }

    #[pyfunction]
    #[pyo3(signature = (dims=[100,100,100], resolution=1.4142135, origin=[0.0,0.0,0.0], atom_typing=false, all_atom_types=Vec::new()))]
    fn voxelize(
        dims: [usize; 3],
        resolution: f32,
        origin: [f32; 3],
        atom_typing: bool,
        all_atom_types: Vec<String>
    ) -> PyResult<(
        u64,
        u64,
        usize,
        Vec<u32>
    )> {
        let mut stdout: Box<dyn Write> = Box::new(std::io::stdout().lock());

        let (
            num_rows,
            num_cols,
            num_condensed_cols,
            condensed_data_idx
        ) = voxelize_stream(
            dims, 
            resolution, 
            origin[0], origin[1], origin[2],
            atom_typing,
            all_atom_types,
            &mut stdout
        ).expect("Failed to voxelize"); 

        println!("Shape of data: ({} molecules, {} voxels)", num_rows, num_cols);
        println!("Shape of condensed data: ({} molecules, {} voxels)", num_rows, num_condensed_cols);

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

        Ok(())
    }
}
