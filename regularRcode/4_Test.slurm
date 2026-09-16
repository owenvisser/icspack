#!/bin/bash
#SBATCH --job-name=ISS_Test_NHANES# Job name
#SBATCH --mail-type=ALL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -A somnath.datta
#SBATCH --qos=somnath.datta-b
#SBATCH --mail-user=owenvisser@ufl.edu # Where to send mail	
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=60:00:00
#SBATCH --output=iss_test_%j.log        # Standard output and error log

pwd; hostname; date

module load R

echo "Running Sam's ISS Test on NHANES Data"

R CMD BATCH 4_Subgroup_Informativeness_Test.R ## the file to run

date