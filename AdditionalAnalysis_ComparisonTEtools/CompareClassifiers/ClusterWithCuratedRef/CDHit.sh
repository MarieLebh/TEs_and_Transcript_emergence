#!/bin/bash
#SBATCH --mem 10G
#SBATCH --cpus-per-task=1


/global/projects/programs/source/cd-hit-v4.6.8-2017-1208/cd-hit-est-2d -i ../ClassificationErrorCheck/MCTE_RechEtAl22.fasta -i2 ../ReClassifyTransposonUltimateLibs/FinalAnnotations_TransposonSequences.fasta -o ClusterTU_MC -c 0.8 -G 0 -aS 0.8 -n 5 -d 0
