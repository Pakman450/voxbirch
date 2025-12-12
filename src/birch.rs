use nalgebra::{DMatrix, RowVector, VecStorage, U1, Dyn};

use crate::voxel::VoxelGrid;
use std::any::TypeId;
use std::rc::Rc;
use std::cell::RefCell;

#[derive(Debug, Clone)]
struct BFNode {
    threshold: f32,
    max_branches: usize,
    is_leaf: bool,
    n_features: usize,
    d_type: TypeId,
    next_leaf: Option<Rc<RefCell<BFNode>>>,
    prev_leaf: Option<Rc<RefCell<BFNode>>>,
    subclusters: Option<Vec<BFSubcluster>>,
    init_centroids: Option<DMatrix<f32>>,
    centroids: Option<DMatrix<f32>>,
}

impl BFNode {
    pub fn new(
        threshold: f32,
        max_branches: usize,
        is_leaf: bool,
        n_features: usize,
        d_type: TypeId,

        ) -> Self {
        BFNode {
            threshold,
            max_branches,
            is_leaf,
            n_features,
            d_type,
            next_leaf: None,
            prev_leaf: None,
            subclusters: None,
            init_centroids: Some(DMatrix::<f32>::zeros(max_branches+1, n_features)),
            centroids: Some(DMatrix::<f32>::zeros(max_branches+1, n_features)),
        }
    }
    
    pub fn append_subcluster (&mut self, subcluster: BFSubcluster) {

        
        let n_samples: usize = match &self.subclusters {
            Some(subclusters) => subclusters.len(),
            None => 0,
        };


        // --- take the centroid BEFORE moving subcluster ---
        let centroid = subcluster.centroid.as_ref().unwrap().clone();
        // let centroid_row = centroid.transpose(); // 1×D

        match &mut self.subclusters {
            Some(subclusters) => subclusters.push(subcluster),
            None => self.subclusters = Some(vec![subcluster]),
        };

        // self.init_centroids.as_mut().unwrap().copy_from(&centroid);

        // NOTE. THIS IS DIFFERENT FROM ORIGINAL CODE
        self.init_centroids.as_mut().unwrap().row_mut(n_samples).copy_from(&centroid.row(0));

        // self.init_centroids.as_mut().unwrap().rows_mut(0, n_samples+1).copy_from(&centroid.slice(
        //     (0,0), 
        //     (n_samples+1, centroid.ncols())
        // ));


        // I am gonna keep this here for now.
        // This code copies all centroids at once.
        // I wonder why would you do this? 
        // maybe there is a reason to overwrite all centroids at once.
        self.centroids.as_mut().unwrap().rows_mut(0, n_samples+1).copy_from(&self.init_centroids.as_ref().unwrap().slice(
            (0,0), 
            (n_samples+1, self.init_centroids.as_ref().unwrap().ncols())
        ));

        // println!("init_centroids {} {}", self.init_centroids.clone().unwrap(), n_samples)


        // NOTE: The above code block is equivalent to the line below, but less efficient.
        // But there must be reason to do it this way.
        // self.centroids.as_mut().unwrap().row_mut(n_samples).copy_from(&self.init_centroids.as_ref().unwrap().row(n_samples));

        // println!("Appending subcluster. Total subclusters: {}",&self.init_centroids.as_ref().unwrap().slice(
        //     (0,0), 
        //     (n_samples+1, self.init_centroids.as_ref().unwrap().ncols())
        // ));
        
        // println!("Appending subcluster. Total subclusters: {:?}", self.centroids);
    }

    pub fn insert_bf_subcluster(
        &mut self,
        subcluster: BFSubcluster,
        set_bits: f32,
        parent: Option<Rc<RefCell<BFNode>>>,
        singly: bool
    ) -> bool {

        if self.subclusters.is_none() {
            self.append_subcluster(subcluster.clone());
            return false;      
        }

        let threshold = self.threshold;
        let max_branches = self.max_branches;

        // Find the closest subcluster among all subclusters
        // so we can insert our new subcluster there

        // println!("self.centroids{}", self.centroids.as_ref().unwrap());
        // println!("subcluster.centroid {}", &subcluster.centroid.as_ref().unwrap().transpose());


        let a = self.centroids.as_ref().unwrap() * &subcluster.centroid.as_ref().unwrap().transpose();

        let row_sums: Vec<f32> = (0..self.centroids.as_ref().unwrap().nrows())
            .map(|i| self.centroids.as_ref().unwrap().row(i).sum())
            .collect();

        let mut sim_matrix = a.clone();

        for i in 0..sim_matrix.nrows() {
            let denom = row_sums[i] + set_bits;
            for j in 0..sim_matrix.ncols() {
                sim_matrix[(i,j)] = sim_matrix[(i,j)] / (denom - sim_matrix[(i,j)]);
            }
        }

        let mut max_val = f32::MIN;
        let mut closest_index = (0,0);

        for i in 0..sim_matrix.nrows() {
            for j in 0..sim_matrix.ncols() {
                if sim_matrix[(i,j)] > max_val {
                    max_val = sim_matrix[(i,j)];
                    closest_index = (i,j);
                }
            }
        }

        let closest_subcluster = &self.subclusters.as_ref().unwrap()[closest_index.0];


        println!("closest_index = {:?}", closest_index);
        println!("closest_subcluster {:?} ", closest_subcluster);
        // println!("sim_matrx{}", sim_matrix);
        // Placeholder implementation
        true
    }
}

#[derive(Debug, Clone)]
struct BFSubcluster {
    nj: u32,
    ls: Option<Vec<f32>>,
    mols: Vec::<u32>,
    cj: Option<Vec<f32>>,
    child: Option<Rc<RefCell<BFNode>>>,
    parent: Option<Rc<RefCell<BFNode>>>,
    centroid: Option<DMatrix<f32>>,
}

impl BFSubcluster {
    pub fn new(linear_sum: Option<Vec<f32>>, mol_indices: Vec<u32>, max_branches: usize, n_features: usize) -> Self {

        // add linear sum into centroid__
        // centroid__ = 
        // Some(DMatrix::from_row_slice(1, n, &v))
        if linear_sum == None {
            BFSubcluster {
                nj: 0,
                ls: Some(linear_sum.clone().unwrap()),
                mols: Vec::<u32>::new(),
                cj: None,
                child: None,
                parent: None,
                centroid:  Some(DMatrix::<f32>::zeros(max_branches+1, n_features)),
            }
            
        } else {

            let mut centroid_zeros = DMatrix::<f32>::zeros(max_branches + 1, n_features);

            // Convert linear_sum Vec<f32> into a dynamic RowVector
            let row_vec = RowVector::<f32, Dyn, VecStorage<f32, U1, Dyn>>::from_vec(linear_sum.clone().unwrap());

            // Set the first row
            centroid_zeros.row_mut(0).copy_from(&row_vec);


            BFSubcluster {
                nj: 1,
                ls: Some(linear_sum.clone().unwrap()),
                mols: mol_indices,
                cj: Some(linear_sum.clone().unwrap()),
                child: None,
                parent: None,
                centroid: Some(centroid_zeros),
            }

        }
        
    }


}

#[derive(Debug)]
pub struct VoxBirch {
    threshold: f32,
    max_branches: usize,
    index_tracker: u32,
    first_call: bool,
    root: Option<Rc<RefCell<BFNode>>>,
    dummy_leaf: Option<Rc<RefCell<BFNode>>>,
}


impl VoxBirch {
    pub fn new(threshold: f32, max_branches: usize) -> Self {
        VoxBirch {
            threshold,
            max_branches,
            index_tracker: 0,
            first_call: true,
            root: None,
            dummy_leaf: None
        }
    }


    // Fit function. only takes in one grid. 
    pub fn fit(&mut self, grids : &DMatrix<f32>, singly: bool) -> &VoxBirch {


        println!("Fitting BIRCH model to voxel grids...");

        let n_features = grids.ncols();

        // get data types
        // let d_type = grids.nrows();


        if self.first_call {

            self.root = 
                Some(Rc::new(RefCell::new(BFNode::new(
                    self.threshold,
                    self.max_branches,
                    true,
                    n_features,
                    TypeId::of::<f32>(),
                ))));


            self.dummy_leaf =                        
                Some(Rc::new(RefCell::new(BFNode::new(
                    self.threshold,
                    self.max_branches,
                    true,
                    n_features,
                    TypeId::of::<f32>(),
                ))));


            // self.root = root;
            // self.dummy_leaf = dummy_leaf;


            self.dummy_leaf
                .as_ref()
                .unwrap()
                .borrow_mut()
                .next_leaf = self.root.clone();

            self.root
                .as_ref()
                .unwrap()
                .borrow_mut()
                .prev_leaf = self.dummy_leaf.clone();

        }

        // check if matrix is sparse
        // NOTE: original BitBirch checks for sparse matrices
        // This is not applicable here since we are convert 3D mols to voxel grids
        // if not spars.issparse() {
        //     panic!("Currently only sparse matrices are supported.");
        // }

        for iter in 0..grids.nrows() {
            
            let grid: Option<Vec<f32>> = Some(grids.row(iter).iter().copied().collect());
            let set_bits: f32 = grids.row(iter).sum();
            let mol_indices: Vec<u32> = vec![self.index_tracker];
            let subcluster = BFSubcluster::new(grid.clone(), mol_indices, self.max_branches, grid.unwrap().len());

            let split = self.root.as_ref().unwrap().borrow_mut().insert_bf_subcluster(
                subcluster.clone(),
                set_bits,
                subcluster.parent,
                singly
            );

            // println!("{:?}", subcluster);


        }

        



        self
    }


}



