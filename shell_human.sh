#!/bin/bash
#BATCH --job-name=human_trial
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kaiyiwu0124@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --mem=250gb
#SBATCH --time=12:00:00
#SBATCH --output=../melissa/result.log

module load matlab/2020a
matlab -nosplash -nodesktop <run_human_tests.m> output.txt
