#!/bin/bash
#SBATCH -c 8                # Number of cores (-c)
#SBATCH -t 0-10:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=16G            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o logs/threshold_analysis_output_plots_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e logs/threshold_analysis_output_plots_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=END
#SBATCH --mail-user=madeleine_kline@hms.harvard.edu

# set R packages and rstudio server singularity image locations
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_19
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_19.sif"

# Run the R script
singularity exec --cleanenv --env R_LIBS_USER=${my_packages} ${rstudio_singularity_image} Rscript threshold_analysis_output_plots.R
~
~
