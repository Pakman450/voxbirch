# VoxBirch

VoxBirch is an voxel-based ultra-fast clustering algorithm 
that can categorizes molecules in 3D space. The 
VoxBirch takes advantage of the Balanced Iterative 
Reducing and Clustering using Hierarchies or BIRCH 
algorithm to efficently cluster 3D molecules at linear 
time. 

## How to build

### Install `rustup`
You must install `cargo` by installing `rustup` to 
build voxbirch

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

or visit `https://rust-lang.org/learn/get-started/`
for more information. 

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

If you want to make a fresh copy (` make clean && make install`), run

```
make fresh
```

### Example commands

If you want to cluster molecules with 10 max_branches:

```
voxbirch -p molecules.mol2 -m 10
```

Then you should get a list of clustered molecules via id, `clusters_mol_ids.txt`. And, corresponding mol2s files for 
each cluster index in `./molecules`, ranging from most 
populated to least populated (defaulted to write
10 top ranking clusters).

If you want to cluster molecules with atom typing:

```
voxbirch -p molecules.mol2 -m 10 --atom-typing
```

If you must see the debug logs for the purposes of development:

```
voxbirch -p molecules.mol2 -vv > output.out 2>&1
```

## Outputs

There are three main outputs
1. tmp/
2. molecules/
2. clustered_mol_ids.txt

`voxbirch` writes out grids and molecules into `tmp` 
temporarily so users can read in large molecules 
on-the-fly. This makes the program memory efficient. 

By default, `voxbirch` will write out the top ten 
clusters based on size into `molecules/`

If you require whole list of clusters, it is written
in `clustered_mol_ids.txt`

