use nalgebra::{DMatrix, RowVector, RowDVector, VecStorage, U1, Dyn};
use core::panic;
use std::{any::TypeId};
use std::rc::Rc;
use std::cell::RefCell;
use log::{error, warn, debug};
use std::io::{Write};

use crate::isim::{jt_isim_real, jt_isim_binary};


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

fn set_merge(
    merge_criterion: MergeCriterion, 
    tolerance: f32
) -> Box<dyn FnMut(
    f32, 
    &Vec<f32>, 
    &Vec<f32>,
    &Vec<f32>, 
    usize, 
    &Vec<f32>, 
    &Vec<f32>,
    &Vec<f32>, 
    usize, 
    usize,
    &mut dyn Write,
) -> bool> {
    match merge_criterion {
        MergeCriterion::Radius => {
            Box::new(move |threshold, new_ls, new_ss, new_centroid, new_n, old_ls, old_ss, nom_ls, old_n, nom_n, write_out| {
                let jt_sim = jt_isim_real(&[new_ls.clone(), new_centroid.clone()].concat(), &new_ss, new_n + 1)
                    * (new_n + 1) as f32
                    - jt_isim_real(&new_ls, &new_ss, new_n) * (new_n - 1) as f32;
                jt_sim >= threshold * 2.0
            })
        }
        MergeCriterion::Diameter => {
            Box::new(move |threshold, new_ls, new_ss, new_centroid, new_n, old_ls, old_ss, nom_ls, old_n, nom_n, write_out| {
                let jt_radius = jt_isim_real(
                    &new_ls, 
                    &new_ss, 
                    new_n,
                );

                if jt_radius >= threshold {
                    writeln!(write_out, "Merging due to diameter criterion: jt_radius = {}", jt_radius).unwrap();
                }

                jt_radius >= threshold
            })
        }
        MergeCriterion::ToleranceTough => {
            Box::new(move |threshold, new_ls,new_ss, new_centroid, new_n, old_ls, old_ss, nom_ls, old_n, nom_n, write_out| {
                let jt_radius = jt_isim_real(&new_ls, &new_ss, new_n);
                if jt_radius < threshold {
                    return false;
                } else {
                    if old_n == 1 && nom_n == 1 {
                        return true;
                    } else if nom_n == 1 {
                        (jt_isim_binary(&(element_wise_add(&old_ls, &nom_ls)), old_n + 1) * (old_n + 1) as f32
                            - jt_isim_binary(&old_ls, old_n) * (old_n - 1) as f32)
                            / 2.0
                            >= jt_isim_binary(&old_ls, old_n) - tolerance
                            && jt_radius >= threshold
                    } else {
                        (jt_isim_binary(&(element_wise_add(&old_ls, &nom_ls)), old_n + nom_n) * (old_n + nom_n) as f32 * (old_n + nom_n - 1) as f32
                            - jt_isim_binary(&old_ls, old_n) * old_n as f32 * (old_n - 1) as f32
                            - jt_isim_binary(&nom_ls, nom_n) * nom_n as f32 * (nom_n - 1) as f32)
                            / (2.0 * (old_n * nom_n) as f32)
                            >= jt_isim_binary(&old_ls, old_n) - tolerance
                            && jt_radius >= threshold
                    }
                }
            })
        }
        MergeCriterion::Tolerance => {
            Box::new(move |threshold, new_ls, new_ss, new_centroid, new_n, old_ls, old_ss, nom_ls, old_n, nom_n, write_out| {
                let jt_radius = jt_isim_binary(&new_ls, new_n);
                if jt_radius < threshold {
                    return false;
                } else {
                    if old_n == 1 && nom_n == 1 {
                        return true;
                    } else if nom_n == 1 {
                        (jt_isim_binary(&(element_wise_add(&old_ls, &nom_ls)), old_n + 1) * (old_n + 1) as f32
                            - jt_isim_binary(&old_ls, old_n) * (old_n - 1) as f32)
                            / 2.0
                            >= jt_isim_binary(&old_ls, old_n) - tolerance
                            && jt_radius >= threshold
                    } else {
                        return true;
                    }
                }
            })
        }
    }
}


fn max_seperation(centroids: &DMatrix<f32>, max_branches: usize) -> (usize, usize, Vec<f32>, Vec<f32>){

    // Get the centroid of the set
    let n_samples: u32 = centroids.nrows().try_into().unwrap();
    // NOTE: row_sum adds all rows in an elementwise fashion per column. very confusing... 
    let linear_sum: Vec<f32> = centroids.row_sum().as_slice().to_vec();

    let mut centroid = calc_centroid( 
        &linear_sum, 
        n_samples, 
        max_branches, 
        centroids.ncols() 
    );

    // iterate through each row of centroids and compute if it is all zeros at a per row fashion
    for i in 0..centroids.nrows() {
        let row = centroids.row(i);
        if row.iter().all(|&x| x == 0.0) {
            error!("Row {i} is all zeros");
        }
    }
    
    // Get the similarity of each molecule to the centroid
    // NOTE: column_sum adds all cols in an elementwise fashion per row. very confusing... 
    let pop_counts: Vec<f32> = centroids.column_sum().as_slice().to_vec();
    
    // If centroid is all zeros, pick the most central molecule via medoid calculation
    // as centroid. We should make this robust to all zero centroids.
    // Also, we need to make sure centroid is not all zeros.
    // Also, we need to make a medoid calculation function, which is shown below.
    // Remember, this only happens with your have a sparse voxel space.
    // A sparse voxel space could be generated by low resolution, 
    // and/or molecules spread out too thinly in a large box.
    let centroid_is_zero: bool = centroid.row(0).iter().all(|&x| x == 0.0);
    debug!("centroid.iter().any(|x| x.is_nan()) before: = {}", centroid.iter().any(|x| x.is_nan()));
    if centroid_is_zero {

        warn!("({}:{}) The centroid has all zeros, probably due to a sparse voxel space: {}", 
            file!(), 
            line!(),
            centroid_is_zero
        );

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
    let centroid_has_non_zero = centroid.row(0).iter().any(|&x| x != 0.0);

    let centroids_has_non_zero = centroids.iter().any(|&x| x != 0.0);

    //This needs to calculate dot products for each row of centroids
    // with centroid. centrpoid has a shape of n_mols, n_features
    debug!("centroids * &centroid.transpose(): = {}", centroids * &centroid.transpose());
    let a_centroid: Vec<f32>= (centroids * &centroid.transpose()).column(0).iter().copied().collect();

    debug!(
        "max_separation: 
        centroid_non_zero = {},
        centroid_nans = {},
        centroid_inf = {},
        centroids_non_zero = {},
        centroids_nans = {}, 
        centroids_inf = {},
        centroid.nrows(), ncols() = {},{},
        centroids.nrows(), ncols() = {},{},
        a_centroid.len() = {},
        a_centroid = {:?}",
        centroid_has_non_zero, 
        centroid.iter().any(|x| x.is_nan()),
        centroid.iter().any(|x| x.is_infinite()),
        centroids_has_non_zero,
        centroids.iter().any(|x| x.is_nan()),
        centroids.iter().any(|x| x.is_infinite()),
        centroid.nrows(), centroid.ncols(),
        centroids.nrows(), centroids.ncols(),
        a_centroid.len() ,
        a_centroid
    );

    let centroid_norms_sq = centroids.row_iter()
        .map(|row| row.dot(&row))
        .collect::<Vec<f32>>();

    let mut query_norm_sq = centroid.row(0).dot(&centroid.row(0));

    let mut sims_med : Vec<f32> = Vec::with_capacity(a_centroid.len());

    for i in 0..a_centroid.len() {
        let union = centroid_norms_sq[i] + query_norm_sq - a_centroid[i];
        let score = a_centroid[i] / (union);

        debug!(
        "\n--- centroid {} ---
        a_centroid:        {:?},
        centroid_norm_sq:  {:?},
        query_norm_sq:     {:?},
        union:             {:?},
        tanimoto:          {:?}",
        i,
        a_centroid[i],
        centroid_norms_sq[i],
        query_norm_sq,
        union,
        score
        );
        debug_assert!(score >= 0.0 && score <= 1.0 );
        sims_med.push(score);
    }

    let mol1 = sims_med
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap().0;

    let a_mol1 : Vec<f32> = (centroids * centroids.row(mol1).transpose()).as_slice().to_vec();

    debug!("a_mol1: = {:?}, centroids * centroids.row(mol1).transpose(): = {}", 
    a_mol1, 
    centroids * centroids.row(mol1).transpose());

    let mut sims_mol1 : Vec<f32> = Vec::with_capacity(a_mol1.len());

    query_norm_sq = centroids.row(mol1).dot(&centroids.row(mol1));

    // Get the similarity of each molecule to mol1 
    for i in 0..a_mol1.len() {
        let union = centroid_norms_sq[i] + query_norm_sq - a_mol1[i];
        let score = a_mol1[i] / (union);

        debug!(
        "\n--- centroid {} ---
        a_mol1:            {:?},
        centroid_norm_sq:  {:?},
        query_norm_sq:     {:?},
        union:             {:?},
        tanimoto:          {:?}",
        i,
        a_mol1[i],
        centroid_norms_sq[i],
        query_norm_sq,
        union,
        score
        );
        debug_assert!(score >= 0.0 && score <= 1.0 );
        sims_mol1.push(score);
    }


    // # Get the least similar molecule to mol1
    let mol2 = sims_mol1
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap().0;

    // # Get the similarity of each molecule to mol2
    let a_mol2 : Vec<f32> = (centroids * centroids.row(mol2).transpose()).as_slice().to_vec();
    debug!(
        "a_mol2: = {:?}, centroids * centroids.row(mol2).transpose(): = {}", 
        a_mol2, 
        centroids * centroids.row(mol2).transpose()
    );
    let mut sims_mol2 : Vec<f32> = Vec::with_capacity(a_mol2.len());

    query_norm_sq = centroids.row(mol2).dot(&centroids.row(mol2));

    for i in 0..a_mol2.len() {
        let union = centroid_norms_sq[i] + query_norm_sq - a_mol2[i];
        let score = a_mol2[i] / (union);

        debug!(
        "\n--- centroid {} ---
        a_mol2:            {:?},
        centroid_norm_sq:  {:?},
        query_norm_sq:     {:?},
        union:             {:?},
        tanimoto:          {:?}",
        i,
        a_mol2[i],
        centroid_norms_sq[i],
        query_norm_sq,
        union,
        score
        );
        debug_assert!(score >= 0.0 && score <= 1.0 );
        sims_mol2.push(score);
    }

    return (mol1, mol2, sims_mol1, sims_mol2)
}

// NOTE: Original did something different here. 
// This function doesn't take in max_branches or n_features.
// But because the original Code language can mutate any variable
// with any other data type. We need max_branches and n_featires 
// to retern a generated Dmatrix for self.centroid. 

// TODO: change this to a real-based centroid calculation
fn calc_centroid( ls: &Vec<f32>, nj: u32, max_branches: usize, n_features: usize) -> DMatrix<f32> {
    
    debug_assert_eq!(ls.len(), n_features);

    let mut centroid = DMatrix::<f32>::zeros(max_branches + 1, n_features);

    for (j, &x) in ls.iter().enumerate() {
        centroid[(0, j)] = x / nj as f32;
        // debug!("
        //     --calc_centroid-col-{j}--
        //     x: = {},
        //     nj: = {}",
        //     x,
        //     nj
        // );
        debug_assert!(!(x / nj as f32).is_nan());
    }

    centroid
}

fn split_node(
    node_child: &Option<Rc<RefCell<BFNode>>>,
    threshold: f32,
    max_branches: usize
)-> (BFSubcluster, BFSubcluster){


    let mut node = node_child.as_ref().unwrap().borrow_mut();
    debug!("Before split_node procress");
    debug!("node.subclusters.as_mut().unwrap().len(): {}", node.subclusters.as_mut().unwrap().len());
    debug!("node.subclusters.as_mut().unwrap().len(): {}", node.subclusters.as_mut().unwrap().len());
    for (idx, subcluster) in node.subclusters.as_mut().unwrap().iter_mut().enumerate() {
        debug!("idx_of_subcluster: {}", idx);
        debug!("subcluseter.ls: {:?}", subcluster.clone().ls.unwrap());
        debug!("subcluseter.nj: {:?}", subcluster.clone().nj);
        debug!("subcluseter.ss: {:?}", subcluster.clone().ss.unwrap());
        debug!("subcluseter.mols: {:?}", subcluster.clone().mols.unwrap());

        debug_assert!(!(subcluster.nj == 0));
    }
    let mut new_subcluster1 = BFSubcluster::new(
        None, 
        Vec::<String>::new(), 
        max_branches, 
        node.n_features
    );
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
                    node.d_type,
                )));

    let  new_node2 = Rc::new(RefCell::new(BFNode::new(
                    threshold,
                    max_branches,
                    node.is_leaf,
                    n_features,
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
        _,
        node1_dist,
        node2_dist
    ) = max_seperation(&node.centroids.clone().unwrap(), max_branches);


    let mut node1_closer: Vec<bool> = node1_dist.iter()
        .zip(node2_dist.iter())
        .map(|(&d1, &d2)| d1 > d2)
        .collect();

    // FROM BITBIRCH: "Make sure node1 is closest to itself even if all distances are equal.
    // This can only happen when all node.centroids_ are duplicates leading to all
    // distances between centroids being zero.""
    node1_closer[farthest_idx1] = true; 

    debug!("
    node1_closer.len() = {},
    node1_closer: = {:?},
    node1_closer all true = {}",
    node1_closer.len(),
    node1_closer,
    node1_closer.iter().all(|&x| x == true));


    if node1_closer.iter().all(|&x| x == true) {
        panic!("node1_closer is all true. new_subcluster2 would never be instantiated");
    }
    
    for (idx, subcluster) in node.subclusters.as_mut().unwrap().iter_mut().enumerate() {
        assert!(!(subcluster.nj == 0));
        if node1_closer[idx] {
            new_node1.borrow_mut().append_subcluster(subcluster.clone());
            new_subcluster1.update(subcluster, max_branches, n_features);

        // It could be possible that the second sublcuster was never updated due
        // to this if then logic. 
        } else {
            new_node2.borrow_mut().append_subcluster(subcluster.clone());
            new_subcluster2.update(subcluster, max_branches, n_features);
        }
    }

    (new_subcluster1, new_subcluster2)
} 

// Find the closest subcluster among all subclusters
// via index so we can insert our new subcluster there
// NOTE: this finds the closest subcluster via tani sim
// via real numbers
fn find_closest_subluster(
    centroids: &DMatrix<f32>, 
    centroid: &DMatrix<f32>
) -> usize {

    // perform dot product between two matrices. not inner dot product
    let a = centroids * &centroid.transpose();

    debug!(
        "a: = {},
        centroids non zero?: = {},
        centroid non zero?: = {}",
        a,
        centroids.iter().any(|&x| x != 0.0),
        centroid.iter().any(|&x| x != 0.0)
    );

    let centroid_norms_sq = centroids.row_iter()
        .map(|row| row.dot(&row))
        .collect::<Vec<f32>>();

    let query_norm_sq = centroid.row(0).dot(&centroid.row(0));

    let mut similarity = Vec::with_capacity(a.nrows());

    for i in 0..a.nrows() {
        let union = centroid_norms_sq[i] + query_norm_sq - a[(i,0)];
        let score = a[(i,0)] / (union);
        debug!(
        "\n--- centroid {} ---
        dot:               {:?},
        centroid_norm_sq:  {:?},
        query_norm_sq:     {:?},
        union:             {:?},
        tanimoto:          {:?}",
        i,
        a[i],
        centroid_norms_sq[i],
        query_norm_sq,
        union,
        score
        );
        debug_assert!(score >= 0.0);

        similarity.push(score);
    }

    let mut max_val = f32::MIN;
    let mut closest_index = 0;  

    for i in 0..similarity.len() {
        debug!("similarity[{}] = {}", i, similarity[i]);
        if similarity[i] > max_val {
            max_val = similarity[i];
            closest_index = i;
        }
        
    }

    debug!("find_closest_subluster: 
    centroids shape = {:?}, 
    centroid shape = {:?}, 
    similarity len = {:?}, 
    similarity = {:?}",
        centroids.shape(),
        centroid.shape(),
        similarity.len(),
        similarity
    );

    closest_index
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
       
        debug!("Before split_node procress");
        if !self.subclusters.is_none() {
            debug!("node.subclusters.as_mut().unwrap().len(): {}", self.subclusters.as_mut().unwrap().len());
            debug!("node.subclusters.as_mut().unwrap().len(): {}", self.subclusters.as_mut().unwrap().len());
            for (idx, sc) in self.subclusters.as_mut().unwrap().iter_mut().enumerate() {
                debug!("idx_of_subcluster: {}", idx);
                debug!("subcluseter.ls: {:?}", sc.clone().ls.unwrap());
                debug!("subcluseter.nj: {:?}", sc.clone().nj);
                debug!("subcluseter.ss: {:?}", sc.clone().ss.unwrap());
                debug!("subcluseter.mols: {:?}", sc.clone().mols.unwrap());

                debug_assert!(!(sc.nj == 0));
            }
        }

        debug!("###############");
        debug!("subcluseter.ls: {:?}", subcluster.clone().ls.unwrap());
        debug!("subcluseter.nj: {:?}", subcluster.clone().nj);
        debug!("subcluseter.ss: {:?}", subcluster.clone().ss.unwrap());
        debug!("subcluseter.mols: {:?}", subcluster.clone().mols.unwrap());

        debug_assert!(!(subcluster.nj == 0));
        // This also returns the index for the last subcluster via index. 
        // We need this to update the init_centroids and centroids matrices.
        // remember n_samples is zero-indexed. so 5 means the 6th index. 
        // This is important for copying the centroid row to the correct index.
        // n_samples is actually used to update centroid matrices for 
        // the newly appended subcluster.
        let n_samples: usize = match &self.subclusters {
            Some(subclusters) => subclusters.len(),
            None => 0,
        };

        // --- take the centroid BEFORE moving subcluster ---
        let centroid = subcluster.centroid.as_ref().unwrap().clone();

        debug!("append_sublcuster: centroid non zero?: = {},", centroid.iter().any(|&x| x != 0.0));
        // append the subcluster, which is gonna the n_samples index now.
        match &mut self.subclusters {
            Some(subclusters) => subclusters.push(subcluster),
            None => self.subclusters = Some(vec![subcluster]),
        };

        debug!("subclusters.len() = {}, after = {}", n_samples , &self.subclusters.as_ref().unwrap().len() );

        // copy the centroid into the init_centroid
        // n_samples is the last index BEFORE we pushed into subclusters. 
        let num_rows = self.init_centroids.as_ref().unwrap().nrows();
        if n_samples < num_rows {
            self.init_centroids.as_mut().unwrap().row_mut(n_samples).copy_from(&centroid.row(0));
        } else {
            error!("n_samples out of bounds. init_centroids has {} rows, but tried to access {}th row", num_rows, n_samples+1);
        }

        // I am gonna keep this here for now.
        // This code copies all centroids at once.
        // I wonder why would you do this? 
        // maybe there is a reason to overwrite all centroids at once.
        self.centroids.as_mut().unwrap().rows_mut(0, n_samples+1).copy_from(
            &self.init_centroids.as_ref().unwrap().view(
            (0,0), 
            (n_samples+1, self.init_centroids.as_ref().unwrap().ncols())
        ));
        
    }

    fn update_split_subclusters(
        & mut self,
        row_idx: usize , 
        new_subcluster1: BFSubcluster ,
        new_subcluster2: BFSubcluster
    ) {

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
        parent: BFSubcluster,
        write_out: &mut impl Write
    ) -> bool {

        debug!("
        subcluster insert_bf_subcluster,
        ls.len() = {:?},
        nj = {:?},
        mols.len() = {:?}",
        subcluster.clone().ls.unwrap().len(),
        subcluster.clone().nj,
        subcluster.clone().mols.unwrap().len()
        );

        debug_assert!(!(subcluster.nj == 0));

        if self.subclusters.is_none() {
            self.append_subcluster(subcluster.clone());
            return false;      
        }

        let threshold = self.threshold;
        let max_branches = self.max_branches;

        let row_idx = find_closest_subluster(
            self.centroids.as_ref().unwrap(),
            subcluster.centroid.as_ref().unwrap(),
        ); 

        // print statistics of length of various fields here
        debug!("
            Inserting BFSubcluster: 
            threshold = {},
            max_branches = {},
            subclusters.len() = {:?}, 
            closest_index/row_idx = {:?},
            subclusters.child[row_idx].sublcusters.len() = {:?},
            init_centroids shape = {:?}, 
            centroids shape = {:?},
            subclusters.child[row_idx].child is none? = {:?}",
            threshold,
            max_branches,
            self.subclusters.as_ref().map_or(0, |s| s.len()),
            row_idx, 
            self.subclusters.as_mut().unwrap()[row_idx].child.as_ref().map_or(
                0, 
                |c| c.borrow().subclusters.as_ref().map_or(0, |sc| sc.len())),
            self.init_centroids.as_ref().map_or((0,0), |c| (c.nrows(), c.ncols())
            ),
            self.centroids.as_ref().map_or((0,0), |c| (c.nrows(), c.ncols())),
            self.subclusters.as_mut().unwrap()[row_idx].child.is_none()
        );

        let closest_node_is_none = self.subclusters
            .as_mut()
            .unwrap()[row_idx]
            .child
                .is_none();

        if closest_node_is_none {
            let merged = self.subclusters.as_mut().unwrap()[row_idx].merge_subcluster(
                subcluster.clone(), max_branches, threshold, write_out
            );

            debug!(
                "Merged status: {}", merged
            );

            if !merged {
                if self.subclusters.as_ref().unwrap().len() >= self.max_branches + 1 {
                    println!(
                        "
                        The child node of the closest subcluster is None, and
                        merging criteria has not met. With this: VoxBirch reached
                        max branches. This typically means your molecules are  
                        very diverse, and requires more lenient clustering. 
                        Consider lowering threshold or more branches to cluster 
                        effectively.
                        "
                    );
                }
                self.append_subcluster(subcluster.clone());
                return self.subclusters.as_ref().unwrap().len() > self.max_branches
            }
            
            let closest_subcluster = self.subclusters.as_mut().unwrap(); // Unwrap Option to get a mutable reference
            let row = & closest_subcluster[row_idx].centroid;

            // NOTE: this saves the first row of centroids
            self.centroids
                .as_mut()
                .unwrap()
                .row_mut(row_idx)
                .copy_from(
                    &row.clone().unwrap().row(0)
            );
            self.init_centroids
                .as_mut()
                .unwrap()
                .row_mut(row_idx)
                .copy_from(
                    &row.clone().unwrap().row(0)
            );

            return false

        }


        let split_child = 
            self.subclusters
                .as_mut()
                .unwrap()[row_idx]
                    .child
                    .as_ref()
                    .unwrap()
                    .borrow_mut()
                    .insert_bf_subcluster(
                        subcluster.clone(),
                        parent.clone(),
                        write_out,
                );
        
        if split_child {
            let (
                new_subcluster1, new_subcluster2
            ) = split_node(
                &self.subclusters.as_mut().unwrap()[row_idx].child, // by reference. meaning this must be mutated
                threshold,
                max_branches
            );

            // this includes append_subcluster
            self.update_split_subclusters(
                row_idx,
                new_subcluster1, 
                new_subcluster2
            );

            return self.subclusters.as_ref().unwrap().len() > self.max_branches    
        }

        self.subclusters
            .as_mut()
            .unwrap()[row_idx]
            .update(
                &subcluster, 
                self.max_branches, 
                self.n_features 
            );

        // NOTE: this saves the first row of centroids
        self.init_centroids.as_mut()
            .unwrap()
            .row_mut(row_idx)
            .copy_from(
                &self.subclusters
                    .as_ref()
                    .unwrap()[row_idx]
                    .centroid
                    .clone()
                    .unwrap()
                    .row(0)
            );
        
        self.centroids.as_mut()
            .unwrap()
            .row_mut(row_idx)
            .copy_from(
                &self.subclusters
                    .as_ref()
                    .unwrap()[row_idx]
                    .centroid
                        .clone()
                        .unwrap()
                        .row(0)
            );

        return false
        
    }
}

#[derive(Debug, Clone)]
struct BFSubcluster {
    nj: u32,
    ls: Option<Vec<f32>>,
    ss: Option<Vec<f32>>,
    mols: Option<Vec<String>>,
    // cj: Option<Vec<f32>>,
    child: Option<Rc<RefCell<BFNode>>>,
    centroid: Option<DMatrix<f32>>,
}

impl BFSubcluster {
    pub fn new(linear_sum: Option<Vec<f32>>, mol_titles: Vec<String>, max_branches: usize, n_features: usize) -> Self {

        if linear_sum == None {

            BFSubcluster {
                nj: 0,
                ls: Some(vec![0.0; n_features]),
                ss: Some(vec![0.0; n_features]),
                mols: Some(Vec::<String>::new()),
                child: None,
                centroid:  Some(DMatrix::<f32>::zeros(max_branches+1, n_features)),
            }
            
        } else {

            let mut centroid_zeros = DMatrix::<f32>::zeros(max_branches + 1, n_features);

            // Convert linear_sum Vec<f32> into a dynamic RowVector
            let row_vec = RowVector::<f32, Dyn, VecStorage<f32, U1, Dyn>>::from_vec(linear_sum.clone().unwrap());

            // Set the first row
            centroid_zeros.row_mut(0).copy_from(&row_vec);

            debug!("BFSubcluster.centroid = {}.", centroid_zeros.iter().any(|&x| x != 0.0));

            BFSubcluster {
                nj: 1,
                ls: Some(linear_sum.clone().unwrap()),
                ss: Some(
                    linear_sum.clone().unwrap()
                    .iter()
                    .map(|&x| x * x)
                    .collect()
                ),
                mols: Some(mol_titles),
                child: None,
                centroid: Some(centroid_zeros),
            }

        }
        
    }

    pub fn update(& mut self, subcluster: &BFSubcluster, max_branches: usize, n_features: usize) {


        debug!("self.nj: = {} | subcluster.nj: = {}", self.nj,subcluster.nj);
        debug_assert!(self.nj + subcluster.nj > 0);

        self.nj += subcluster.nj;

        // This IS the same as elementwise increments between two ndarrays 
        if let (Some(a), Some(b)) = (self.ls.as_mut(), subcluster.ls.as_ref()) {
            assert_eq!(a.len(), b.len());

            for (x, y) in a.iter_mut().zip(b.iter()) {
                *x += *y;
            }
        }

        // update the squared sums by taking the self.ss and adding subcluster.ss elementwise
        self.ss = Some(
            self.ss.as_ref().unwrap().iter()
            .zip(subcluster.ss.as_ref().unwrap().iter())
            .map(|(&x, &y)| x + y)
            .collect()
        );

        // This is NOT the same as elementwise increments. This is an extension between two
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
        // NOR IS IT THE SAME LOGIC AS the resolution to deal with zero vectors
        // in the max_seperation FUNCTION 
        // Remember, this only happens with your have a sparse voxel space.
        // A sparse voxel space could be generated by low resolution, 
        // and/or molecules spread out too thinly in a large box.
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
        threshold: f32,
        write_out: & mut impl Write
    ) -> bool {

        let new_ls: Vec<f32> = self.ls.clone().unwrap().iter()
            .zip(nominee_cluster.ls.clone().unwrap().iter())
            .map(|(x, y)| *x + *y)
            .collect();

        let new_ss: Vec<f32> = self.ss.clone().unwrap().iter()
            .zip(nominee_cluster.ss.clone().unwrap().iter())
            .map(|(x, y)| *x + *y)
            .collect();

        let new_n = self.nj + nominee_cluster.nj;

        debug!("merge_subcluster:");
        debug_assert!(self.nj + nominee_cluster.nj > 0);

        // BUG: Here we know that calc_centroids can return zero.
        let new_centroid = calc_centroid(&new_ls, new_n, max_branches, new_ls.len() );


        // TODO: this needs to be changed where set_merge can called anywhere.
        let mut merge_accept = set_merge(MergeCriterion::Diameter, 0.05);

        // Collect the row into a Vec<f32> instead of trying to call as_slice()
        let row_as_vec: Vec<f32> = new_centroid.row(0).iter().cloned().collect();

        if merge_accept(
            threshold, 
            &new_ls, 
            &new_ss,
            &row_as_vec, 
            new_n as usize, 
            &self.ls.as_ref().unwrap(), 
            &self.ss.as_ref().unwrap(), 
            &nominee_cluster.ls.unwrap(), 
            self.nj as usize, 
            nominee_cluster.nj as usize,
            write_out
        ) {

            writeln!(
                write_out, 
                " Accepted between clusters with sizes {} and {}", self.nj, nominee_cluster.nj
            ).unwrap();

            self.nj = new_n;
            self.ls = Some(new_ls);
            self.ss = Some(new_ss);
            self.centroid = Some(new_centroid); 
            // self.mol_indices = self.mol_indices + nominee_cluster.mol_indices;
            if let (Some(a), Some(b)) = (self.mols.as_mut(), nominee_cluster.mols.as_ref()) {
                writeln!(write_out,"  Merging mols {:?} with {:?}", a, b).unwrap();
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
        stdout: & mut Box<dyn Write>
        ) -> &mut VoxBirch 
        {

        assert_eq!(grids.nrows(), titles.len());

        writeln!(stdout,"\n#############################\nFitting grids with the VoxBirch Clustering\n#############################\n\n").unwrap();
        let n_features = grids.ncols();

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


        for iter in 0..grids.nrows() {

            let mol_title = titles.get(iter).unwrap().to_string();

            writeln!(stdout, "\nInserting {}: {}/{}", 
                mol_title, 
                iter+1, grids.nrows()).unwrap();

            let grid: Option<Vec<f32>> = Some(
                grids.row(iter).iter().copied().collect()
            );
            let mol_indices: Vec<String> = vec![mol_title];
            let subcluster = BFSubcluster::new(
                grid.clone(),
                mol_indices, 
                self.max_branches, 
                grid.unwrap().len());

            // printout the inserting subcluster stats 
            debug!("
            subcluster --{}--
            ls.len() = {:?},
            nj = {:?},
            mols.len() = {:?}",
            iter+1,
            subcluster.clone().ls.unwrap().len(),
            subcluster.clone().nj,
            subcluster.clone().mols.unwrap().len()
            );
            debug_assert!(subcluster.nj == 1);

            let split = self.root
                .as_ref()
                .unwrap()
                .borrow_mut()
                .insert_bf_subcluster(
                    subcluster.clone(),
                    // Here, BitBirch feeds in subcluster.parent_
                    // But since Rust is a strongly typed language
                    // we must send in BFSubcluster rather than. BFNode.
                    subcluster.clone(),
                    stdout
            );

            if split {

                let (new_subcluster1, new_subcluster2) = split_node(
                    &self.root,
                    self.threshold,
                    self.max_branches
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



