# VoxBirch

VoxBirch is an voxel-based ultra-fast clustering algorithm 
that scales at large-scale for molecules in 3D. The 
VoxBirch takes advantage of the Balanced Iterative 
Reducing and Clustering using Hierarchies or BIRCH 
algorithm to efficently cluster 3D molecules at linear 
time. 

## How to build

### Install `rustup`
You must install `cargo` by install `rustup`to build voxbirch

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

or visit `https://rust-lang.org/learn/get-started/`

Then pull the voxbirch repository from github.

`git clone https://github.com/Pakman450/voxbirch`

`cd` into `voxfing`, then run 

```
make install
```

and you will see your `voxbirch` binary in `./bin`

To clean up your binary, run

```
make clean
```

### Example commands

If you want to cluster molecules with 10 max_branches:

```
voxbirch -p molecules.mol2 -m 10
```

Then you should get a list of clustered molecules via id, `clusters_mol_ids.txt`

If you must see the debug logs for the purposes of development:

```
voxbirch -p molecules.mol2 -vv > output.out 2>&1
```


