#!/bin/bash
#SBATCH --mem 10G
#SBATCH --cpus-per-task=1


transposon_classifier_RFSB -mode classify -fastaFile /global/group/research/m_lebh01/TEpaperCorrection_Revisions/ClassificationErrorCheck/MCTE_RechEtAl22.fasta -outputPredictionFile MCTE_RFSB_results.txt


