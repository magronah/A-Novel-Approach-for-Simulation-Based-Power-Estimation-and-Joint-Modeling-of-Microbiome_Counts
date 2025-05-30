#!/bin/bash
#SBATCH --job-name="Single job"     # job name
#SBATCH --nodes=1                # number of node MUST be 1
#SBATCH --cpus-per-task=40        # number of processes
#SBATCH --mem-per-cpu=50G      # memory; default unit is megabytes
#SBATCH --time=0-15:15           # time (DD-HH:MM)
#SBATCH --output=sim-%j.log               # log file
#SBATCH --error=sim-%j.err                # error file
#SBATCH --mail-user=agronahm@mcmaster.ca    # who to email
#SBATCH --mail-type=FAIL                  # when to email
#SBATCH --account=def-bolker
module load r/4.3.1

R CMD BATCH  wald_cov_read.R > wald_cov_read.Rout
