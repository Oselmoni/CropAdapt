#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --time=0-20:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=40000
#SBATCH --output=./slurmOutput/R-%x.%j.out
#SBATCH --wait

orthofinder -f DATA/RAW/PROTEINS/ -t 40

mv DATA/RAW/PROTEINS/*/* DATA/ORTHOFINDER/