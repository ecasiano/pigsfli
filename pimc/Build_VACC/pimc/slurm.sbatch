#!/bin/bash
#
#SBATCH -J pigsl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=30:00:00
#SBATCH --array=0-999
#SBATCH -o out/pigsl_%A_%4a.out
line_number=$((1 + SLURM_ARRAY_TASK_ID))
LINE=$(sed -n "$line_number"p pimc_commands)
echo $LINE
time eval "$LINE"
