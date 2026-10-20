#!/bin/bash
#SBATCH --mem 10G
#SBATCH --cpus-per-task=1


RepeatClassifier -consensi FinalAnnotations_TransposonSequences.fasta
