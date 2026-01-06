use nalgebra::{DMatrix, RowVector, RowDVector, VecStorage, U1, Dyn};

use std::any::TypeId;
use std::rc::Rc;
use std::cell::RefCell;
use log::{error, warn, info, debug, trace};

#[derive(Debug)]
enum Parent {
    // Node(BFNode),
    Subcluster(BFSubcluster)
}

#[derive(Clone, Copy, Debug)]
enum MergeCriterion {
    Radius,
    Diameter,
    ToleranceTough,
    Tolerance,
}

fn element_wise_add(vec1: &Vec<f32>, vec2: &Vec<f32>) -> Vec<f32> {
    // Ensure the vectors have the same length
    if vec1.len() != vec2.len() {
        panic!("Vectors must be of the same length");
    }

    // Create a new Vec<f32> to hold the result
    let result: Vec<f32> = vec1.iter()
        .zip(vec2.iter()) // Pair elements from both vectors
        .map(|(a, b)| a + b) // Add them element-wise
        .collect(); // Collect the results into a new vector

    result
}

fn set_merge(merge_criterion: MergeCriterion, tolerance: f32) -> Box<dyn Fn(f32, &Vec<f32>, &Vec<f32>, usize, &Vec<f32>, &Vec<f32>, usize, usize) -> bool> {
    match merge_criterion {
        MergeCriterion::Radius => {
            Box::new(move |threshold, new_ls, new_centroid, new_n, old_ls, nom_ls, old_n, nom_n| {
                let jt_sim = jt_isim(&[new_ls.clone(), new_centroid.clone()].concat(), new_n + 1)
                    * (new_n + 1) as f32
                    - jt_isim(&new_ls, new_n) * (new_n - 1) as f32;
                jt_sim >= threshold * 2.0
            })
        }
        MergeCriterion::Diameter => {
            Box::new(move |threshold, new_ls, new_centroid, new_n, old_ls, nom_ls, old_n, nom_n| {
                let jt_radius = jt_isim(&new_ls, new_n);
                jt_radius >= threshold
            })
        }
        MergeCriterion::ToleranceTough => {
            Box::new(move |threshold, new_ls, new_centroid, new_n, old_ls, nom_ls, old_n, nom_n| {
                let jt_radius = jt_isim(&new_ls, new_n);
                if jt_radius < threshold {
                    return false;
                } else {
                    if old_n == 1 && nom_n == 1 {
                        return true;
                    } else if nom_n == 1 {
                        (jt_isim(&(element_wise_add(&old_ls, &nom_ls)), old_n + 1) * (old_n + 1) as f32
                            - jt_isim(&old_ls, old_n) * (old_n - 1) as f32)
                            / 2.0
                            >= jt_isim(&old_ls, old_n) - tolerance
                            && jt_radius >= threshold
                    } else {
                        (jt_isim(&(element_wise_add(&old_ls, &nom_ls)), old_n + nom_n) * (old_n + nom_n) as f32 * (old_n + nom_n - 1) as f32
                            - jt_isim(&old_ls, old_n) * old_n as f32 * (old_n - 1) as f32
                            - jt_isim(&nom_ls, nom_n) * nom_n as f32 * (nom_n - 1) as f32)
                            / (2.0 * (old_n * nom_n) as f32)
                            >= jt_isim(&old_ls, old_n) - tolerance
                            && jt_radius >= threshold
                    }
                }
            })
        }
        MergeCriterion::Tolerance => {
            Box::new(move |threshold, new_ls, new_centroid, new_n, old_ls, nom_ls, old_n, nom_n| {
                let jt_radius = jt_isim(&new_ls, new_n);
                if jt_radius < threshold {
                    return false;
                } else {
                    if old_n == 1 && nom_n == 1 {
                        return true;
                    } else if nom_n == 1 {
                        (jt_isim(&(element_wise_add(&old_ls, &nom_ls)), old_n + 1) * (old_n + 1) as f32
                            - jt_isim(&old_ls, old_n) * (old_n - 1) as f32)
                            / 2.0
                            >= jt_isim(&old_ls, old_n) - tolerance
                            && jt_radius >= threshold
                    } else {
                        return true;
                    }
                }
            })
        }
    }
}

// Define the jt_isim function
fn jt_isim(c_total: &Vec<f32>, n_objects: usize) -> f32 {
    // Sum of the elements in c_total (column-wise sum in Python)
    let sum_kq: f32 = c_total.iter().copied().sum();

    // Sum of squares (dot product of c_total with itself)
    let sum_kqsq: f32 = c_total.iter().copied().map(|x| x * x).sum();

    // Compute the variable a
    let a = (sum_kqsq - sum_kq) / 2.0;

    println!("jt_isim = {} : a = {}, sum_kq = {}, sum_kqsq = {}, n_objects = {}", 
        a / (a + (n_objects as f32) * sum_kq - sum_kqsq),
        a, 
        sum_kq,
        sum_kqsq, 
        n_objects
        );

    // Return the iSIM Jaccard-Tanimoto value
    a / (a + (n_objects as f32) * sum_kq - sum_kqsq)
}
fn max_separation(centroids: &DMatrix<f32>, max_branches: usize) -> (usize, usize, Vec<f32>, Vec<f32>){

    

    // Get the centroid of the set
    let n_samples: u32 = centroids.nrows().try_into().unwrap();
    // NOTE: row_sum adds all rows in an elementwise fashion per column. very confusing... 
    let linear_sum: Vec<f32> = centroids.row_sum().as_slice().to_vec();
    let has_nonzero2 = linear_sum
        .iter()
        .any(|&x| x != 0.0);
    println!("&linear_sum in max_seperation has nonzero? {}", has_nonzero2);

    println!("linear_sum in max_seperation: {:?}", linear_sum);

    // BUG: Here we know that calc_centroids returns zero.
    let mut centroid = calc_centroid( &linear_sum, n_samples, max_branches, centroids.ncols() );
    
    // If centroid is all zeros, pick the most central molecule via medoid calculation
    // as centroid. We should make this robust to all zero centroids.
    // Also, we need to make sure centroid is not all zeros.
    // Also, we need to make a medoid calculation function, which is shown below.
    let centroid_is_zero = centroid.row(0).iter().all(|&x| x == 0.0);
    if centroid_is_zero {
        let pop_counts: Vec<f32> = centroids.column_sum().as_slice().to_vec();

        let mut best_idx = 0;
        let mut best_score = -1.0;

        for i in 0..centroids.nrows() {
            let a_i: Vec<f32> =
                (centroids * centroids.row(i).transpose()).as_slice().to_vec();

            let sims: Vec<f32> = a_i.iter()
                .zip(pop_counts.iter())
                .map(|(&a, &p)| a / (p + pop_counts[i] - a))
                .collect();

            let avg_sim: f32 = sims.iter().sum::<f32>() / sims.len() as f32;

            if avg_sim > best_score {
                best_score = avg_sim;
                best_idx = i;
            }
        }
        centroid.row_mut(0).copy_from(&centroids.row(best_idx));

    }
    
    let has_nonzero1 = centroid
        .iter()
        .any(|&x| x != 0.0);
    println!("&centroid in max_seperation has nonzero? {}", has_nonzero1);
  
    // Get the similarity of each molecule to the centroid
    // NOTE: column_sum adds all cols in an elementwise fashion per row. very confusing... 
    let pop_counts: Vec<f32> = centroids.column_sum().as_slice().to_vec();

    println!("pop_counts: {:?}", pop_counts);

    // println!("centroids: {:?}", centroids);
    // println!("centroid: {:?}", centroid);

    //This needs to calculate dot products for each row of centroids
    // with centroid. centrpoid has a shape of n_features of cols
    let a_centroid: Vec<f32> = (centroids * centroid.transpose()).as_slice().to_vec(); 

    let has_nonzeroz = a_centroid
        .iter()
        .any(|&x| x != 0.0);
    println!("&a_centroid in max_seperation has nonzero? {}", has_nonzeroz);
   
    println!("a_centroid: {:?}", a_centroid);

    let sims_med: Vec<f32> = a_centroid.iter()
        .zip(pop_counts.iter())
        .map(|(&a, &p)| a / (p + centroid.row(0).sum()  - a))
        .collect();
    
    let has_nonzero = sims_med.iter().any(|&x| x != 0.0);
    println!("sims_med has nonzero? {}", has_nonzero);

    println!("sims_med: {:?}", sims_med);


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
    
    let has_nonzero = a_mol1.iter().any(|&x| x != 0.0);
    println!("a_mol1 has nonzero? {}", has_nonzero);

    println!("a_mol1: {:?}", a_mol1);

    let has_nonzero = pop_counts.iter().any(|&x| x != 0.0);
    println!("pop_counts has nonzero? {}", has_nonzero);

    println!("pop_counts: {:?}", pop_counts);
    
    // sims_mol1 = a_mol1 / (pop_counts + pop_counts[mol1] - a_mol1)

    let sims_mol1: Vec<f32> = a_mol1.iter()
        .zip(pop_counts.iter())
        .map(|(&a, &p)| a / (p + pop_counts[mol1]  - a))
        .collect();

    // # Get the least similar molecule to mol1
    // mol2 = np.argmin(sims_mol1)
    let has_nonzero = sims_mol1.iter().any(|&x| x != 0.0);
    println!("sims_mol1 has nonzero? {}", has_nonzero);

    println!("sims_mol1: {:?}", sims_mol1);

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
    
    println!("nj, threshold: {}, {}", nj, threshold);



    for (j, &x) in ls.iter().enumerate() {

        if x > threshold {
            println!("in calc_centroid, GREATER. x > threshold {}, {}",x, threshold);
        }
        
        centroid[(0, j)] = if x >= threshold { 1.0 } else { 0.0 };

        // centroid[(0, j)] = if x > 0.0 { 1.0 } else { 0.0 };

    }



    centroid
}

fn split_node(
    node_child: &Option<Rc<RefCell<BFNode>>>,
    threshold: f32,
    max_branches: usize,
    singly: bool
)-> (BFSubcluster, BFSubcluster){


    let mut node = node_child.as_ref().unwrap().borrow_mut();

    let mut new_subcluster1 = BFSubcluster::new(
        None, 
        Vec::<String>::new(), 
        max_branches, 
        node.n_features);
    let mut new_subcluster2 = BFSubcluster::new(
        None, 
        Vec::<String>::new(), 
        max_branches, 
        node.n_features
    );

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

    let has_nonzero = node.centroids
        .iter()
        .any(|mat| mat.iter().any(|&v| v != 0.0));
    println!("&node.centroids has nonzero? {}", has_nonzero);
    println!("node.centroids[0] {}", &node.centroids.clone().unwrap().row(0));
    let has_nonzero2 = &node.centroids.clone().unwrap().row(0)
        .iter()
        .any(|&x| x != 0.0);

    println!("&node.centroids.clone().unwrap().row(0) has zero?: {:?}", has_nonzero2);

    let (
        farthest_idx1,
        _,
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

        // println!("INIT {:?}", self.init_centroids);
        // println!("CENTROID {:?}", centroid);


        // NOTE. THIS IS DIFFERENT FROM ORIGINAL CODE
        // self.init_centroids.as_mut().unwrap().row_mut(n_samples).copy_from(&centroid.row(0));
        let num_rows = self.init_centroids.as_ref().unwrap().nrows();
        if n_samples < num_rows {
            self.init_centroids.as_mut().unwrap().row_mut(n_samples).copy_from(&centroid.row(0));
        } else {
            // Handle the error, possibly initialize more rows or log a message
            eprintln!("n_samples out of bounds. Matrix has {} rows, but tried to access row {}", num_rows, n_samples);
        }
        // self.init_centroids.as_mut().unwrap().rows_mut(0, n_samples+1).copy_from(&centroid.slice(
        //     (0,0), 
        //     (n_samples+1, centroid.ncols())
        // ));


        // I am gonna keep this here for now.
        // This code copies all centroids at once.
        // I wonder why would you do this? 
        // maybe there is a reason to overwrite all centroids at once.
        self.centroids.as_mut().unwrap().rows_mut(0, n_samples+1).copy_from(&self.init_centroids.as_ref().unwrap().view(
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
        row_idx: usize , 
        mut new_subcluster1: BFSubcluster ,
        mut new_subcluster2: BFSubcluster ,
        singly : bool
    ) {
        if !singly {
            new_subcluster1.parent = self.subclusters.as_ref().unwrap()[0].parent.clone();
            new_subcluster2.parent = self.subclusters.as_ref().unwrap()[0].parent.clone();
        }

        // ind = self.subclusters.index()
        self.subclusters.as_mut().unwrap()[row_idx] = new_subcluster1.clone();
        self.init_centroids.as_mut().unwrap().row_mut(row_idx).copy_from(
            &new_subcluster1.centroid.clone().unwrap().row(0)
        );
        self.centroids.as_mut().unwrap().row_mut(row_idx).copy_from(
            &new_subcluster1.centroid.clone().unwrap().row(0)
        );
        // self.centroids.unwrap()[row_idx] = new_subcluster1.centroid;
        self.append_subcluster(new_subcluster2.clone());

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

        println!("closeset_index: {}", row_idx);

        // println!("self.centroids: {:?}", self.centroids);


        if !self.subclusters.as_mut().unwrap()[row_idx].child.is_none() {

            parent = self.subclusters.as_mut().unwrap()[row_idx].clone();

            let split_child = 
                self.subclusters.as_mut().unwrap()[row_idx]
                .child.as_ref().unwrap().borrow_mut().insert_bf_subcluster(
                    subcluster.clone(),
                    set_bits,
                    parent.clone(),
                    singly
                );

            if !split_child {

                self.subclusters.as_mut().unwrap()[row_idx].update(&subcluster, self.max_branches, self.n_features );

                // NOTE: this saves the first row of centroids
                self.init_centroids.as_mut().unwrap()
                    .row_mut(row_idx)
                    .copy_from(&self.subclusters.as_ref().unwrap()[row_idx].centroid.clone().unwrap().row(0));
                
                self.centroids.as_mut().unwrap()
                    .row_mut(row_idx)
                    .copy_from(&self.subclusters.as_ref().unwrap()[row_idx].centroid.clone().unwrap().row(0));

                let row = self.centroids.as_ref().unwrap().row(row_idx);
                let has_nonzero = row.iter().any(|&x| x != 0.0);
                println!("1Row {} has nonzero? {}", row_idx, has_nonzero);
                // println!("1self.centroids[row_idx]: {:?}", self.centroids.as_mut().unwrap().row(row_idx).iter().take(10).collect::<Vec<_>>());
                return false

            } else {

                let (new_subcluster1, new_subcluster2) = split_node(
                    &self.subclusters.as_mut().unwrap()[row_idx].child, // by reference. meaning this must be mutated
                    threshold,
                    max_branches,
                    singly
                );

                self.update_split_subclusters(
                    row_idx,
                    new_subcluster1, 
                    new_subcluster2,
                    singly
                );

                if self.subclusters.as_ref().unwrap().len() > self.max_branches {
                    return true
                }
                let row = self.centroids.as_ref().unwrap().row(row_idx);
                let has_nonzero = row.iter().any(|&x| x != 0.0);
                println!("2Row {} has nonzero? {}", row_idx, has_nonzero);
                return false
            }

        } else {

            let merged = self.subclusters.as_mut().unwrap()[row_idx].merge_subcluster(
                subcluster.clone(), max_branches, threshold
            );
            
            println!("MERGED? {}", merged);


            if merged {

                let closest_subcluster = self.subclusters.as_mut().unwrap(); // Unwrap Option to get a mutable reference
                let row = & closest_subcluster[row_idx].centroid;

                // NOTE: this saves the first row of centroids
                self.centroids.as_mut().unwrap().row_mut(row_idx).copy_from(&row.clone().unwrap().row(0));
                self.init_centroids.as_mut().unwrap().row_mut(row_idx).copy_from(&row.clone().unwrap().row(0));
                let row = self.centroids.as_ref().unwrap().row(row_idx);
                let has_nonzero = row.iter().any(|&x| x != 0.0);
                println!("3Row {} has nonzero? {}", row_idx, has_nonzero);
                if !singly{
                    // closest_subcluster.parent = ps; 
                    self.subclusters.as_mut().unwrap()[row_idx] = parent;
                }
                return false

            // # not close to any other subclusters, and we still
            // # have space, so add.
            } else if  self.subclusters.as_ref().unwrap().len() < self.max_branches {
                self.append_subcluster(subcluster.clone());
                if !singly{
                    // closest_subcluster.parent = ps; 
                    self.subclusters.as_mut().unwrap()[row_idx] = parent;
                }
                let row = self.centroids.as_ref().unwrap().row(row_idx);
                let has_nonzero = row.iter().any(|&x| x != 0.0);
                println!("4Row {} has nonzero? {}", row_idx, has_nonzero);

                // println!("4self.centroids[row_idx]: {:?}", self.centroids.as_mut().unwrap().row(row_idx).iter().take(10).collect::<Vec<_>>());

                return false

            // # We do not have enough space nor is it closer to an
            // # other subcluster. We need to split.
            } else {
                self.append_subcluster(subcluster.clone());

                let row = self.centroids.as_ref().unwrap().row(row_idx);
                let has_nonzero = row.iter().any(|&x| x != 0.0);
                println!("5Row {} has nonzero? {}", row_idx, has_nonzero);

                // println!("5self.centroids[row_idx]: {:?}", self.centroids.as_mut().unwrap().row(row_idx).iter().take(10).collect::<Vec<_>>());

                return true

            }


        }

        // true // Rust compiler states it is an unreachable expression
    }
}

#[derive(Debug, Clone)]
struct BFSubcluster {
    nj: u32,
    ls: Option<Vec<f32>>,
    mols: Option<Vec<String>>,
    // cj: Option<Vec<f32>>,
    child: Option<Rc<RefCell<BFNode>>>,
    parent: Option<Rc<RefCell<Parent>>>,
    centroid: Option<DMatrix<f32>>,
}

impl BFSubcluster {
    pub fn new(linear_sum: Option<Vec<f32>>, mol_titles: Vec<String>, max_branches: usize, n_features: usize) -> Self {

        if linear_sum == None {

            BFSubcluster {
                nj: 0,
                ls: Some(vec![0.0; n_features]),
                mols: Some(Vec::<String>::new()),
                // cj: None,
                child: None,
                parent: None,
                centroid:  Some(DMatrix::<f32>::zeros(max_branches+1, n_features)),
            }
            
        } else {

            let mut centroid_zeros = DMatrix::<f32>::zeros(max_branches + 1, n_features);

            // Convert linear_sum Vec<f32> into a dynamic RowVector
            let row_vec = RowVector::<f32, Dyn, VecStorage<f32, U1, Dyn>>::from_vec(linear_sum.clone().unwrap());
            let has_nonzero = row_vec.iter().any(|&x| x != 0.0);

            println!("row_vec in new BFsubcluster, has nonzero? {}", has_nonzero);
            // Set the first row
            centroid_zeros.row_mut(0).copy_from(&row_vec);


            BFSubcluster {
                nj: 1,
                ls: Some(linear_sum.clone().unwrap()),
                mols: Some(mol_titles),
                // cj: Some(linear_sum.clone().unwrap()),
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

            a.extend_from_slice(b);
        }

        // NOTE: Original did something different here. 
        // This function doesn't take in max_branches or n_features.
        // But because the original Code language can mutate any variable
        // with any other data type. We need max_branches and n_featires 
        // to retern a generated Dmatrix for self.centroid. 

        let new_cent = calc_centroid(
            &self.ls.as_ref().unwrap(), 
            self.nj, 
            max_branches, 
            n_features
        );

        // if centroids is all zeros, fallback to presence-based centroid.
        // THIS IS NOT THE SAME AS ORIGINAL CODE
        // NOR IS IT THE SAME LOGIC AS max_seperation FUNCTION
        let centroid_is_zero = new_cent.row(0).iter().all(|&x| x == 0.0);
        if centroid_is_zero {
            // no previous centroid: fallback to presence-based centroid
            let mut fb = new_cent.clone();
            for (j, &v) in self.ls.as_ref().unwrap().iter().enumerate() {
                fb[(0, j)] = if v > 0.0 { 1.0 } else { 0.0 };
            }
            self.centroid = Some(fb);

        } else {
            self.centroid = Some(new_cent);
        }
    }

    fn merge_subcluster (
        & mut self, 
        nominee_cluster: BFSubcluster,
        max_branches: usize,
        threshold: f32 
    ) -> bool {
        
        let new_ls: Vec<f32> = self.ls.clone().unwrap().iter()
            .zip(nominee_cluster.ls.clone().unwrap().iter())
            .map(|(x, y)| *x + *y)
            .collect();

        let new_n = self.nj + nominee_cluster.nj;

        // BUG: Here we know that calc_centroids returns zero.
        let new_centroid = calc_centroid(&new_ls, new_n, max_branches, new_ls.len() );
        let has_nonzero1 = new_centroid.row(0).iter().any(|&x| x != 0.0);
        println!("&centroid in merge_subcluster has nonzero? {}", has_nonzero1);
        // TODO: this needs to be changed where set_merge can called anywhere.
        let merge_accept = set_merge(MergeCriterion::Diameter, 0.05);

        // Collect the row into a Vec<f32> instead of trying to call as_slice()
        let row_as_vec: Vec<f32> = new_centroid.row(0).iter().cloned().collect();

        if merge_accept(
            threshold, 
            &new_ls, 
            &row_as_vec, 
            new_n as usize, 
            &self.ls.as_ref().unwrap(), 
            &nominee_cluster.ls.unwrap(), 
            self.nj as usize, 
            nominee_cluster.nj as usize
        ) {
            self.nj = new_n;
            self.ls = Some(new_ls);
            self.centroid = Some(new_centroid); 
            // self.mol_indices = self.mol_indices + nominee_cluster.mol_indices;
            if let (Some(a), Some(b)) = (self.mols.as_mut(), nominee_cluster.mols.as_ref()) {
                println!("Merging mols {:?} with {:?}", a, b);
                a.extend_from_slice(b);
            }
            return true
        }
        false
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
    subcluster_centers: Option<Vec<RowDVector<f32>>>,
    n_features_out: usize
}


impl VoxBirch {
    pub fn new(threshold: f32, max_branches: usize) -> Self {
        VoxBirch {
            threshold,
            max_branches,
            index_tracker: 0,
            first_call: true,
            root: None,
            dummy_leaf: None,
            subcluster_centers: None,
            n_features_out: 0
        }
    }


    // Fit function. only takes in one grid. 
    pub fn fit(
        &mut self, 
        grids : &DMatrix<f32>, 
        titles: Vec<String>, 
        singly: bool) -> &mut VoxBirch 
        {

        // error!("1This is an error.");
        // warn!("2This is a warning.");
        // info!("3This is normal info.");
        // debug!("4This is debug info (visible only with --verbose).");
        // trace!("5Trace info (only with RUST_LOG=trace).");

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

            println!("grid: {}", iter);
             
            let grid: Option<Vec<f32>> = Some(grids.row(iter).iter().copied().collect());
            let set_bits: f32 = grids.row(iter).sum();
            let mol_indices: Vec<String> = vec![titles.get(iter).unwrap().to_string()];
            let subcluster = BFSubcluster::new(
                grid.clone(),
                mol_indices, 
                self.max_branches, 
                grid.unwrap().len());

            let split = self.root.as_ref().unwrap().borrow_mut().insert_bf_subcluster(
                subcluster.clone(),
                set_bits,
                // Here, BitBirch feeds in subcluster.parent_
                // But since Rust is a strongly typed language
                // we must send in BFSubcluster rather than. BFNode.
                subcluster.clone(),
                singly
            );

            if split {

                let (new_subcluster1, new_subcluster2) = split_node(
                    &self.root,
                    self.threshold,
                    self.max_branches,
                    singly
                );



                self.root = None;

                self.root = Some(Rc::new(RefCell::new(BFNode::new(
                    self.threshold,
                    self.max_branches,
                    false,
                    n_features,
                    TypeId::of::<f32>(),
                ))));

                self.root.as_ref()
                    .unwrap()
                    .borrow_mut()
                    .append_subcluster(new_subcluster1.clone());

                self.root.as_ref()
                    .unwrap()
                    .borrow_mut()
                    .append_subcluster(new_subcluster2.clone());       

                if !singly {
                    // Iterate over subclusters in new_subcluster1
                    if let Some(child_node1) = &new_subcluster1.child {
                        if let Some(subclusters1) = &child_node1.borrow().subclusters {
                            for subcluster in subclusters1 {
                                if let Some(parent) = &subcluster.parent {
                                    // Mutate the parent reference to point to new_subcluster1
                                    *parent.borrow_mut() = Parent::Subcluster(new_subcluster1.clone());
                                }
                            }
                        }
                    }

                    // Iterate over subclusters in new_subcluster2
                    if let Some(child_node2) = &new_subcluster2.child {
                        if let Some(subclusters2) = &child_node2.borrow().subclusters {
                            for subcluster in subclusters2 {
                                if let Some(parent) = &subcluster.parent {
                                    // Mutate the parent reference to point to new_subcluster2
                                    *parent.borrow_mut() = Parent::Subcluster(new_subcluster2.clone());
                                }
                            }
                        }
                    }

                }
            }
            self.index_tracker += 1;
        }

        // Get leaves by calling get_leaves()
        let leaves = self.get_leaves();

        // self.subcluster_centers = Some(leaves);
        // Concatenate the centroids of each leaf node
        let mut rows: Vec<RowDVector<f32>> = Vec::new();

        for leaf in leaves {
            let leaf_ref = leaf.borrow();
            if let Some(c) = &leaf_ref.centroids {
                for i in 0..c.nrows() {
                    rows.push(c.row(i).into_owned());
                }
            }
        }

        // Assuming the centroids are a matrix and need to be reshaped
        self.subcluster_centers = Some(rows.clone());

        self.n_features_out = self.subcluster_centers.as_ref().unwrap().len();
        
        self.first_call = false;
        return self
    }

    fn get_leaves(&self) -> Vec<Rc<RefCell<BFNode>>> {
        let mut leaves = Vec::new();
        let mut leaf_ptr = self.dummy_leaf.as_ref().unwrap().borrow().next_leaf.clone();
        
        while let Some(leaf) = leaf_ptr {
            leaves.push(leaf.clone());
            leaf_ptr = leaf.borrow().next_leaf.clone();
        }
        
        leaves
    }
    pub fn get_cluster_mol_ids(& self) ->  Vec<Vec<String>> {

        if self.first_call {
            panic!("The model has not been fitted yet.");
        }        

        let mut clusters_mol_id = Vec::<Vec::<String>>::new();
        
        for leaf in self.get_leaves() {
            let leaf_ref = leaf.borrow();

            if let Some(subclusters) = leaf_ref.subclusters.as_ref() {
                for subcluster in subclusters.iter() {
                    if let Some(mols) = subcluster.mols.as_ref() {
                        clusters_mol_id.push(mols.clone());
                    }
                }
            }
        }

        // Sort the clusters by the number of samples in the cluster
        // This is descending order 
        clusters_mol_id.sort_by(|a, b| b.len().cmp(&a.len()));

        return clusters_mol_id
    }
}



