# scRNA-Seq
# Author: Brian Wray

### General instructions

- Open a window and navigate to a quest analytics node
- File > New Project... 
- (In new project wizard) Version Control > Git 
- Repository URL: https://github.com/brianFSM/scRNA-seq_seurat_V5.git

After project is created, go to quest on the command line, navigate to the newly created repository folder (which should be a subdirectory in the directory of whatever project you're working on), and copy over the renv.lock file we have on quest. This will tell the project which versions of which packges you'll use:
$ cp /projects/b1197/PROJECTS/Seurat_v5_renv/November_2025/renv.lock . 

Go back to your project in the analytics node, and at the conosole:
> renv::init(bare=TRUE)
> renv::restore()

If you've ever installed the libraries before, they will link from the cache in a flash. Otherwise they will be downloaded and compiled. 

If this is the first time you've ever run these reports, you will need to run the following from the console in the analytics node in order to get pdfs to render (the reports are now pdfs, not html files):

> install.packages("tinytex")
> tinytex::install_tinytex()   # downloads + sets up ~/.TinyTeX
> tinytex::is_tinytex()        # should return TRUE

Finally, fill out your config file. 

To run the 
