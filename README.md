# scRNA-Seq
## Author: Brian Wray

### Overview

This repository holds the pipeline and reporting templates for single-cell RNA-seq analysis using Seurat v5. The project is designed to run within the Quest environment, using renv for reproducibility and SLURM for execution.

### Setup

1. Open a browser window and navigate to a quest analytics node
2. In Rstudio: 
* File → New Project → Version Control → Git
* Repository URL:
* https://github.com/brianFSM/scRNA-seq_seurat_V5.git

### Configure the Environment
Once the project exists locally, SSH into Quest and move into the cloned repository directory. Copy the shared renv.lock file:

`cp /projects/b1197/PROJECTS/Seurat_v5_renv/December_2025/renv.lock .`

Return to the RStudio session and initialize the environment:

`renv::init(bare=TRUE)`

`renv::restore()`

If you’ve installed these packages before, they’ll symlink from cache; otherwise, they’ll build fresh.

### PDF rendering requirements

These reports build to PDF, not HTML. If TinyTeX isn’t installed, do it now:

`install.packages("tinytex")`

`tinytex::install_tinytex()   # downloads + sets up ~/.TinyTeX`

`tinytex::is_tinytex()        # should return TRUE`

### Configure the Run

Edit the configuration YAML to match your data and paths.

###Running Templates on Quest

Submit a report via SLURM:

`$ sbatch run_templates.sh scRNA_part1_QC.Rmd`

In this example, the output will be scRNA_part1_QC.pdf

### Notes
The goal is reproducibility. Don’t install random packages unless you want future-you cursing present-you.
If something breaks, it’s usually the config YAML, missing modules, or someone messing with the lockfile.
