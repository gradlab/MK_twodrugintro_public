#!/bin/bash
#SBATCH -c 10                # Number of cores (-c)
#SBATCH -t 0-36:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p sapphire --contiguous   # Partition to submit to
#SBATCH --mem=200G            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o logs/baselineruns_seq_10k_sapphire_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e logs/baselineruns_seq_10k_sapphire_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=END
#SBATCH --mail-user=madeleine_kline@hms.harvard.edu

# set R packages and rstudio server singularity image locations
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_19
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_19.sif"

# Run the R script
singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} Rscript baselineparms_10kruns_seq_sapphire.R
