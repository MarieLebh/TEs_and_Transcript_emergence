#!/bin/bash
#SBATCH --mem 10G
#SBATCH --cpus-per-task=1


RepeatClassifier -consensi /global/group/research/m_lebh01/TEpaperCorrection_Revisions/ClassificationErrorCheck/MCTE_RechEtAl22.fasta
