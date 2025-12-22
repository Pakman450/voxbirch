use nalgebra::{DMatrix, RowVector, VecStorage, U1, Dyn};

use crate::voxel::VoxelGrid;
use std::any::TypeId;
use std::rc::Rc;
use std::cell::RefCell;

#[derive(Debug)]
enum Parent {
    Node(BFNode),
    Subcluster(BFSubcluster)
}

fn max_separation(centroids: &DMatrix<f32>, max_branches: usize) -> (usize, usize, Vec<f32>, Vec<f32>){

    

    // Get the centroid of the set
    let n_samples: u32 = centroids.nrows().try_into().unwrap();
    // NOTE: row_sum adds all rows in an elementwise fashion per column. very confusing... 
    let linear_sum: Vec<f32> = centroids.row_sum().as_slice().to_vec();
    let centroid = calc_centroid( &linear_sum, n_samples, max_branches, centroids.ncols() );

    // Get the similarity of each molecule to the centroid
    // NOTE: column_sum adds all cols in an elementwise fashion per row. very confusing... 
    let pop_counts: Vec<f32> = centroids.column_sum().as_slice().to_vec();

    //This needs to calculate dot products for each row of centroids
    // with centroid. centrpoid has a shape of n_features of cols
    let a_centroid: Vec<f32> = (centroids * centroid.transpose()).as_slice().to_vec(); 

    let sims_med: Vec<f32> = a_centroid.iter()
        .zip(pop_counts.iter())
        .map(|(&a, &p)| a / (p + centroid.row(0).sum()  - a))
        .collect();


    // # Get the least similar molecule to the centroid
    // mol1 = np.argmin(sims_med)
    let mol1 = sims_med
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap().0;
    // # Get the similarity of each molecule to mol1
    // a_mol1 = np.dot(X, X[mol1])

    let a_mol1 : Vec<f32> = (centroids * centroids.row(mol1).transpose()).as_slice().to_vec();

    // sims_mol1 = a_mol1 / (pop_counts + pop_counts[mol1] - a_mol1)

    let sims_mol1: Vec<f32> = a_mol1.iter()
        .zip(pop_counts.iter())
        .map(|(&a, &p)| a / (p + pop_counts[mol1]  - a))
        .collect();

    // # Get the least similar molecule to mol1
    // mol2 = np.argmin(sims_mol1)

    let mol2 = sims_mol1
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap().0;

    // # Get the similarity of each molecule to mol2
    // a_mol2 = np.dot(X, X[mol2])
    let a_mol2 : Vec<f32> = (centroids * centroids.row(mol2).transpose()).as_slice().to_vec();


    // sims_mol2 = a_mol2 / (pop_counts + pop_counts[mol2] - a_mol2)
    let sims_mol2: Vec<f32> = a_mol2.iter()
        .zip(pop_counts.iter())
        .map(|(&a, &p)| a / (p + pop_counts[mol2]  - a))
        .collect();

    return (mol1, mol2, sims_mol1, sims_mol2)
}

// NOTE: Original did something different here. 
// This function doesn't take in max_branches or n_features.
// But because the original Code language can mutate any variable
// with any other data type. We need max_branches and n_featires 
// to retern a generated Dmatrix for self.centroid. 
fn calc_centroid( ls: &Vec<f32>, nj: u32, max_branches: usize, n_features: usize) -> DMatrix<f32> {
    
    debug_assert_eq!(ls.len(), n_features);

    let threshold = nj as f32 * 0.5;
    let mut centroid = DMatrix::<f32>::zeros(max_branches + 1, n_features);

    for (j, &x) in ls.iter().enumerate() {
        centroid[(0, j)] = if x >= threshold { 1.0 } else { 0.0 };
    }

    centroid
}

fn split_node(node_child: &Option<Rc<RefCell<BFNode>>>, threshold: f32, max_branches: usize, singly: bool )-> (BFSubcluster, BFSubcluster){

    let mut node = node_child.as_ref().unwrap().borrow_mut();

    let mut new_subcluster1 = BFSubcluster::new(None, Vec::<u32>::new(), max_branches, node.n_features);
    let mut new_subcluster2 = BFSubcluster::new(None, Vec::<u32>::new(), max_branches, node.n_features);

    // Assuming `node.n_features` is a field in your `Node` struct
    let n_features = node.n_features; 

    let  new_node1 = Rc::new(RefCell::new(BFNode::new(
                    threshold,
                    max_branches,
                    node.is_leaf,
                    n_features,
                    //dtype=node.init_centroids_.dtype
                    node.d_type,
                )));

    let  new_node2 = Rc::new(RefCell::new(BFNode::new(
                    threshold,
                    max_branches,
                    node.is_leaf,
                    n_features,
                    //dtype=node.init_centroids_.dtype
                    node.d_type,
                )));
    
    new_subcluster1.child = Some(Rc::clone(&new_node1));
    new_subcluster2.child = Some(Rc::clone(&new_node2));
         
    if node.is_leaf {
        if !node.prev_leaf.is_none() {
            node.prev_leaf
            .as_ref()
            .unwrap()
            .borrow_mut().next_leaf = Some(Rc::clone( &new_node1 ));
        }

        new_node1.borrow_mut().prev_leaf = node.prev_leaf.clone();
        new_node1.borrow_mut().next_leaf = Some(Rc::clone(&new_node2));

        new_node2.borrow_mut().prev_leaf = Some(Rc::clone(&new_node1));
        new_node2.borrow_mut().next_leaf = node.next_leaf.clone();

        if !node.next_leaf.is_none() {
            node.next_leaf
            .as_ref()
            .unwrap()
            .borrow_mut().prev_leaf = Some(Rc::clone( &new_node2 ));
        }

    }

    let (
        farthest_idx1,
        farthest_idx2,
        node1_dist,
        node2_dist
    ) = max_separation(&node.centroids.clone().unwrap(), max_branches);


    let mut node1_closer: Vec<bool> = node1_dist.iter()
        .zip(node2_dist.iter())
        .map(|(&d1, &d2)| d1 > d2)
        .collect();

    // FROM BITBIRCH: "Make sure node1 is closest to itself even if all distances are equal.
    // This can only happen when all node.centroids_ are duplicates leading to all
    // distances between centroids being zero.""
    node1_closer[farthest_idx1] = true; 
    
    for (idx, subcluster) in node.subclusters.as_mut().unwrap().iter_mut().enumerate() {
        if node1_closer[idx] {
            new_node1.borrow_mut().append_subcluster(subcluster.clone());
            new_subcluster1.update(subcluster, max_branches, n_features);
            if !singly {
                subcluster.parent = Some(Rc::new(RefCell::new(
                    Parent::Subcluster(new_subcluster1.clone())
                )));
            }
        } else {
            new_node2.borrow_mut().append_subcluster(subcluster.clone());
            new_subcluster2.update(subcluster, max_branches, n_features);
            if !singly {
                // let parent_subcluster = Parent::Subcluster(new_subcluster2.clone());
                subcluster.parent = Some(Rc::new(RefCell::new(
                    Parent::Subcluster(new_subcluster2.clone())
                )));
            }
        }
    }

    (new_subcluster1, new_subcluster2)
}   




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

        // This also returns the index for the last subcluster via index. 
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

    fn update_split_subclusters(
        & mut self,
        subcluster: Rc<RefCell<BFSubcluster>> ,
        new_subcluster1: BFSubcluster ,
        new_subcluster2: BFSubcluster ,
        singly : bool
    ) {
        // if !singly {
        //     new_subcluster1.parent = self.subclusters[0].parent
        //     new_subcluster2.parent = self.subclusters[0].parent
        // }

        // ind = self.subclusters.index()
        // self.subclusters_[ind] = new_subcluster1
        // self.init_centroids_[ind] = new_subcluster1.centroid_
        // self.centroids_[ind] = new_subcluster1.centroid_
        // self.append_subcluster(new_subcluster2)

    }

    pub fn insert_bf_subcluster(
        &mut self,
        subcluster: BFSubcluster,
        set_bits: f32,
        mut parent: BFSubcluster,
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

        // perform dot product between two matrices. not inner dot product
        let a = self.centroids.as_ref().unwrap() * &subcluster.centroid.as_ref().unwrap().transpose();

        let row_sums: Vec<f32> = (0..self.centroids.as_ref().unwrap().nrows())
            .map(|i| self.centroids.as_ref().unwrap().row(i).sum())
            .collect();

        let mut sim_matrix = a.clone();

        // generate sim matrix
        for i in 0..sim_matrix.nrows() {
            let denom = row_sums[i] + set_bits;
            for j in 0..sim_matrix.ncols() {
                sim_matrix[(i,j)] = sim_matrix[(i,j)] / (denom - sim_matrix[(i,j)]);
            }
        }

        let mut max_val = f32::MIN;
        let mut closest_index = (0,0);

        // Find index to the maximum value. 
        for i in 0..sim_matrix.nrows() {
            for j in 0..sim_matrix.ncols() {
                if sim_matrix[(i,j)] > max_val {
                    max_val = sim_matrix[(i,j)];
                    closest_index = (i,j);
                }
            }
        }

        let (row_idx, _) = closest_index; 

        // get a reference towards self.subclusters.as_mut().unwrap()[row_idx]
        let closest_subcluster =
            &mut self.subclusters.as_mut().unwrap()[row_idx];

        if !closest_subcluster.child.is_none() {
            println!("ee2");
            parent = closest_subcluster.clone();

            let split_child = 
                closest_subcluster.child.as_ref().unwrap().borrow_mut().insert_bf_subcluster(
                    subcluster.clone(),
                    set_bits,
                    parent.clone(),
                    singly
                );

            if !split_child {

                closest_subcluster.update(&subcluster, self.max_branches, self.n_features );

                self.init_centroids.as_mut().unwrap()
                    .row_mut(row_idx)
                    .copy_from(&self.subclusters.as_ref().unwrap()[row_idx].centroid.clone().unwrap());
                
                self.centroids.as_mut().unwrap()
                    .row_mut(row_idx)
                    .copy_from(&self.subclusters.as_ref().unwrap()[row_idx].centroid.clone().unwrap());
                return false
            } else {

                let (new_subcluster1, new_subcluster2) = split_node(
                    &closest_subcluster.child, // by reference. meaning this must be mutated
                    threshold,
                    max_branches,
                    singly
                );

                self.update_split_subclusters(
                    closest_subcluster,
                    new_subcluster1, 
                    new_subcluster2,
                    singly
                );

                return false
            }
        }

        println!("closest_index = {:?}", closest_index);
        // println!("closest_subcluster {:?} ", closest_subcluster);
        println!("sim_matrx{}", sim_matrix);
        // Placeholder implementation
        true
    }
}

#[derive(Debug, Clone)]
struct BFSubcluster {
    nj: u32,
    ls: Option<Vec<f32>>,
    mols: Option<Vec<u32>>,
    cj: Option<Vec<f32>>,
    child: Option<Rc<RefCell<BFNode>>>,
    parent: Option<Rc<RefCell<Parent>>>,
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
                mols: Some(Vec::<u32>::new()),
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
                mols: Some(mol_indices),
                cj: Some(linear_sum.clone().unwrap()),
                child: None,
                parent: None,
                centroid: Some(centroid_zeros),
            }

        }
        
    }

    pub fn update(& mut self, subcluster: &BFSubcluster, max_branches: usize, n_features: usize) {


        self.nj += subcluster.nj;


        // NOTE: the original wanted to `self.linear_sum_ += subcluster.linear_sum_`
        // This IS the same as elementwise increments between two ndarrays 
        if let (Some(a), Some(b)) = (self.ls.as_mut(), subcluster.ls.as_ref()) {
            assert_eq!(a.len(), b.len());

            for (x, y) in a.iter_mut().zip(b.iter()) {
                *x += *y;
            }
        }

        // NOTE: the original wanted to `self.mol_indices += subcluster.mol_indices`
        // This is not the same as elementwise increments. This is an extension between two
        // arrays. 
        if let (Some(a), Some(b)) = (self.mols.as_mut(), subcluster.mols.as_ref()) {
            assert_eq!(a.len(), b.len());

            a.extend_from_slice(b);
        }

        // NOTE: Original did something different here. 
        // This function doesn't take in max_branches or n_features.
        // But because the original Code language can mutate any variable
        // with any other data type. We need max_branches and n_featires 
        // to retern a generated Dmatrix for self.centroid. 
        self.centroid = Some(calc_centroid(&self.ls.as_ref().unwrap(), self.nj, max_branches, n_features));
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
    pub fn fit(&mut self, grids : &DMatrix<f32>, singly: bool) {


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
                // Here, BitBirch feeds in subcluster.parent_
                // But since Rust is a strongly typed language
                // we must send in BFSubcluster rather than. BFNode.
                subcluster.clone(),
                singly
            );

            // println!("{:?}", subcluster);


        }

        
    }


}



