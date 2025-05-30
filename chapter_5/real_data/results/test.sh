#!/bin/bash
#SBATCH --job-name=SingleJob
#SBATCH --output=output.txt
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G

# Check OpenMP thread count
echo "OMP_NUM_THREADS is set to: $OMP_NUM_THREADS"

# Run your job here
