# scRNA-Seq

### General instructions

There are 2 sets of files, one set (*_PART1.Rmd, *_PART2.Rmd, etc) that are there for you. Start
by filling out as much of 'config_template.yaml' as you can, and changing the name to 'config.yaml'.
Then run scRNA_template_PART1.Rmd, check out the results, and fill out the cutoffs you need for part2.
Keep moving through those files until you get the results for scRNA_template_PART5.Rmd. As these files
progress they will save all of the relevant figures and information for making reports for the users.

Once you have completed scRNA_template_PART5.Rmd and you are satisfied with the results, you can generate
the QC and clustering reports for the users. These run automatically based off of your config.yaml. Send
the QC and clustering html reports to the users.

You can run the templates, such as number 2, by entering

`sbatch run_templates.sh 2`

as long as you have the default config (./config.yaml) 

By default, if you run the script from the command line, or use the 'knit' button in RMarkdown, the script
will use ./config.yaml as the config file. To pass it an alternate config file at the command line check out
execute_pipeline.R



### To run on quest
[Instructions for installing Seurat on Quest](https://kb.northwestern.edu/page.php?id=98203#Seurat)

Make sure to also install hdf5r in R. First

load module hdf5/1.8.19-serial

then open R/4.1.1 and install.packages("hdf5r")

NOTE: hdf5 is not available on quest analytics. No biggie, just indicate that in the config file.

When running as a submit job on quest, you need to load the following modules in your slurm script:

module purge all

module load R/4.1.1

module load geos/3.8.1

module load hdf5/1.8.19-serial

module load pandoc/2.2.1


### Setting up GitHub on Quest
Follow [instructions for using Git on Quest](https://kb.northwestern.edu/page.php?id=78598)
Clone scRNA-Seq repository from GitHub with HTTPS

To work with different branches and to use 'commit' & 'push', you'll need Personal Access Token
Create Personal Access Token by navigating to Settings > Developer Settings > Personal Access Tokens > Generate new token > Check All > Generate token
Make sure to note down the token as it won't be available in the future.
Whenever prompted while 'commit' or 'push', use the Personal Access Token as your password. This applies only for 'commit' or 'push' directly from Quest.
