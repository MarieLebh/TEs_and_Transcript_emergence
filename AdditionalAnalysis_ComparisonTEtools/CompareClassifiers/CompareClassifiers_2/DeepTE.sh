#!/bin/bash
#SBATCH --mem 50G
#SBATCH --cpus-per-task=1


python /global/students/homes/m_lebh01/Programs/DeepTE/DeepTE.py -d DeepTErun -i /global/group/research/m_lebh01/TransposonUltimateTE/DmelDK/FinalAnnotations_TransposonSequences.fasta -sp M -m_dir /global/students/homes/m_lebh01/Droso_GenomeAnnotations/Metazoans_model

