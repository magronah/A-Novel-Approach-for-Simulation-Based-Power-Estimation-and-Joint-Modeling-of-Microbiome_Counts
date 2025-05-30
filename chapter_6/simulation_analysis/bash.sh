#!/bin/bash

#SBATCH --nodes=1               # number of nodes to use
#SBATCH --time=00-00:30:00         # time (DD-HH:MM:SS)
#SBATCH --job-name="Single job"     # job name

#SBATCH --cpus-per-task=1         # Number of CPU cores to use

#SBATCH --mem=5G                       # memory per node
#SBATCH --output=sim-%j.log               # log file
#SBATCH --error=sim-%j.err                # error file
#SBATCH --mail-user=agronahm@mcmaster.ca    # who to email
#SBATCH --mail-type=FAIL                  # when to email
#SBATCH --account=def-bolker

module load r/4.4.0

R CMD BATCH zinbres.R   zinbres.Rout


