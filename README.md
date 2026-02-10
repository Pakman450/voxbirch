# VoxBirch

VoxBirch is an voxel-based ultra-fast clustering algorithm 
that can categorizes molecules in 3D space. The 
VoxBirch takes advantage of the Balanced Iterative 
Reducing and Clustering using Hierarchies or BIRCH 
algorithm to efficently cluster 3D molecules at linear 
time. 

## Code is inspired by...
Please check out these works, which are the basis for this code base

- BitBIRCH: efficient clustering of large molecular libraries:
    https://doi.org/10.1039/D5DD00030K
- BitBIRCH Clustering Refinement Strategies:
    https://doi.org/10.1021/acs.jcim.5c00627
- BitBIRCH-Lean:
    (preprint) https://www.biorxiv.org/content/10.1101/2025.10.22.684015v1
- iSIM: instant similarity
    https://pubs.rsc.org/en/content/articlelanding/2024/dd/d4dd00041b


![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pakman450/voxbirch/voxbirch.yml?logo=github&label=ci)
![GitHub License](https://img.shields.io/github/license/pakman450/voxbirch)


## How to build the voxbirch binary

### Install `rustup`
You must install `cargo` by installing `rustup` to 
build voxbirch

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

or visit `https://rust-lang.org/learn/get-started/`
for more information. 

Then pull the voxbirch repository from github.

```
git clone https://github.com/Pakman450/voxbirch
```

`cd` into `voxbirch`, then run 

```
make install
```

and you will see your `voxbirch` binary in `./bin`

To clean up your binary, run

```
make clean
```

If you want to make a fresh copy (`make clean && make install`), run

```
make fresh
```

### How to build python wrapper functions

Voxbirch allows users to install a python package. To install: 

```
make bindings
```

Then, you have to activate the virtual python env:

```
source python/.venv/bin/activate
```

This will allow you to interact with the package
on your favorite python environment:

```
python
>>> import voxbirch
>>> voxbirch.read_mol2("path/to/mol.mol2")
>>> ...
>>> ...
>>> vb = voxbirch.Voxbirch( num_features = condensed_num_cols )
>>> vb.cluster()
```


### Example commands

If you want to cluster molecules with 10 max_branches,
provide path to mol file and specify the num of max_branches:

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


## Usage
```
Usage: voxbirch [OPTIONS] --path <PATH>

Options:
  -p, --path <PATH>
          Path to the MOL2 file (required) [default: mol.file]
  -d, --dims <DIMS>
          Dimensions of the voxel grid (x, y, z), comma separated [default: 20,20,20]
  -r, --resolution <RESOLUTION>
          Resolution of the voxel grid in Angstroms [default: 1.4142135]
  -o, --origin <ORIGIN>
          Target origin x0 y0 z0 via comma separated string [default: 0.0,0.0,0.0]
  -t, --threshold <THRESHOLD>
          Threshold of similarity [default: 0.65]
  -m, --max-branches <MAX_BRANCHES>
          Number of max branches [default: 50]
      --clustered-ids-path <CLUSTERED_IDS_PATH>
          Clustered mol ids output name
      --output-path <OUTPUT_PATH>
          Clustered mol ids output name
  -c, --cluster-write-limit <CLUSTER_WRITE_LIMIT>
          Number of clusters to write out [default: 10]
      --atom-typing
          Add atom typing
      --no-condense
          Do not condense voxel grids. Leaving this out condenses grids
  -v, --verbosity...
          Verbosity level. -v means level 1, -vvv means level 3
  -q
          Quiet mode. -q means nothign will printout to screen or to an file
      --rss
          Memory logging mode
  -h, --help
          Print help
  -V, --version
          Print version
```



