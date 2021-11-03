#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=dedup
#SBATCH --error=dedup.err
#SBATCH --output=dedup.out
#SBATCH --time=0-02:00:00
#SBATCH --cpus-per-task=8

conda activate bgmp_py39
/usr/bin/time -v /projects/bgmp/wrosales/bioinfo/Bi624/Deduper-WesleyRosales/rosales_deduper.py -f ./C1_SE_uniqAlign.sorted.sam -u STL96.txt -o ./C1_SE_uniqAlign_deduped.sam