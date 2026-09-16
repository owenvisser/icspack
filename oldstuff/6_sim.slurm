#!/bin/bash
#SBATCH --job-name=sim_NHANES# Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -A somnath.datta
#SBATCH --qos=somnath.datta-b
#SBATCH --mail-user=owenvisser@ufl.edu # Where to send mail	
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=20:00:00
#SBATCH --output=sims_%j.log        # Standard output and error log

pwd; hostname; date

module load R

echo "Running simulation of NHANES Data"

R CMD BATCH 6_Simulator.R ## the file to run

date