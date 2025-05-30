#!/bin/bash
# submit_job.sh

#SBATCH --account=def-bolker

# Job parameters
#SBATCH --job-name=my_job
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10        # number of processes
#SBATCH --mem-per-cpu=100G      # memory; default unit is megaby


#SBATCH --array=1-2
#SBATCH --output=array-%j-%a.out               # log file
#SBATCH --error=array-%j-%a.err                # error file

#SBATCH --mail-user=agronahm@mcmaster.ca    # who to email
#SBATCH --mail-type=FAIL                  # when to email


# Load R module
module load r/4.3.1

# Run R script with array index as argument

Rscript nbmm_combine.R $SLURM_ARRAY_TASK_ID

