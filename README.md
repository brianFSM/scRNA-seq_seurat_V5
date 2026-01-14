# scRNA-Seq
## Author: Brian Wray

### Overview

This repository holds the pipeline and reporting templates for single-cell RNA-seq analysis using Seurat v5. The project is designed to run within the Quest environment, using renv for reproducibility and SLURM for execution.

### Setup

Option 1: analytics node
1. Open a browser window and navigate to a quest analytics node
2. In Rstudio: 
* File → New Project → Version Control → Git
* Repository URL:
* https://github.com/brianFSM/scRNA-seq_seurat_V5.git

Option 2: command line on quest
* git clone https://github.com/brianFSM/scRNA-seq_seurat_V5.git
* cd scRNA-seq_seurat_V5
* module purge all
* module load R/4.4.0

### Configure the Environment
Once the project exists locally, SSH into Quest and move into the cloned repository directory. Copy the shared renv.lock file:

``` bash
cp /projects/b1197/PROJECTS/Seurat_v5_renv/December_2025/renv.lock .
```

Return to the RStudio session in the analytics node, or start R on the command line, and initialize the environment:

```R
renv::init(bare=TRUE)
```

```R
renv::restore()
```

If you’ve installed these packages before, they’ll symlink from cache; otherwise, they’ll build fresh.

### PDF rendering requirements

These reports build to PDF, not HTML. If TinyTeX isn’t installed, do it now:

```R
install.packages("tinytex")

tinytex::install_tinytex()   # downloads + sets up ~/.TinyTeX

tinytex::is_tinytex()        # should return TRUE
```

### Configure the Run

```bash
cp config_template.yaml config.yaml
```

Edit the configuration.yaml to match your data and paths.

### Running Templates on Quest

Submit a report via SLURM:

```bash
 sbatch run_templates.sh scRNA_part1_QC.Rmd
```

In this example, the output will be scRNA_part1_QC.pdf

### Notes
* This pipeline is setup to minimize user input. Render part 1 and look at the resulting pdf to decide cutoffs for QC. Enter those in to the comfig.yaml and then run part2. Part 2 will do the rest! It picks the ideal clustering resolution, and even enters that resolution in to the config.yaml after rendering. 
* The goal is reproducibility. Don’t install random packages unless you want future-you cursing present-you.
* If something breaks, it’s usually the config YAML, missing modules, or someone messing with the lockfile.
