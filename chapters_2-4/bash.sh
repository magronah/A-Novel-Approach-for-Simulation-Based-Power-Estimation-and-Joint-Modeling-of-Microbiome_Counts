#!/bin/bash

#SBATCH --nodes=1               # number of nodes to use
#SBATCH --time=00-03:00:00         # time (DD-HH:MM:SS)
#SBATCH --job-name="Single job"     # job name

#SBATCH --ntasks-per-node=1              # tasks per node
#SBATCH --mem=150G                       # memory per node
#SBATCH --output=sim-%j.log               # log file
#SBATCH --error=sim-%j.err                # error file
#SBATCH --mail-user=agronahm@mcmaster.ca    # who to email
#SBATCH --mail-type=FAIL                  # when to email
#SBATCH --account=def-bolker
module load r/4.3.1

R CMD BATCH  mod_check.R > mod_check.Rout













