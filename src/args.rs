use clap::Parser;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct ArgsV {
    /// Path to the MOL2 file (required)
    #[arg(short, long, required = true, default_value = "mol.file")]
    pub path: String,

    /// Dimensions of the voxel grid (x, y, z), comma separated
    #[arg(short, long, value_delimiter = ',', default_value = "20,20,20")]
    pub dims: Vec<usize>,

    /// Resolution of the voxel grid in Angstroms
    #[arg(short, long, default_value_t = 2.0)]
    pub resolution: f32,

    /// Target origin x0 y0 z0 via comma separated string
    #[arg(short, long, value_delimiter = ',', default_value = "0.0,0.0,0.0")]
    pub origin: Vec<f32>,

    /// Threshold of similarity
    #[arg(short, long, default_value_t = 0.65)]
    pub threshold: f32,

    /// Number of max branches
    #[arg(short, long, default_value_t = 50)]
    pub max_branches: usize,

    /// Clustered mol ids output name
    #[arg(long, default_value = None)]
    pub clustered_ids_path: Option<std::string::String>,

    /// Clustered mol ids output name
    #[arg(long, default_value = None)]
    pub output_path: Option<std::string::String>,

    /// Number of clusters to write out
    #[arg(short, long, default_value_t = 10)]
    pub cluster_write_limit: usize,

    /// Add atom typing
    #[arg(long)]
    pub atom_typing: bool,

    /// Do not condense voxel grids. Leaving this out condenses grids.
    #[arg(long)]
    pub no_condense: bool,

    /// Verbosity level. -v means level 1, -vvv means level 3 
    #[arg(short, long, action = clap::ArgAction::Count, default_value_t = 0)]
    pub verbosity: u8,

    /// Quiet mode. -q means nothign will printout to screen or to an file
    #[arg(short)]
    pub quiet: bool,
}