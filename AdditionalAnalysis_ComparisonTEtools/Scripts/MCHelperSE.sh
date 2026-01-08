#!/bin/bash

#MCHelper command: Run similarly for all other genomes
#First run MCHelper
python3  /global/students/homes/m_lebh01/Programs/MCHelper/MCHelper.py -l Data/DmelSE_TE.fa -o Output/DmelSE -g Data/DmelSE_Genome.fa --input_type fasta -b Data/references.hmm -a F -t 5 -r A

#Then (from a separate folder where you copy the genome and the MCHelper library) run RepeatMasker

RepeatMasker -xsmall -lib curated_sequences_NR.fa -s -gff  DmelSE_Genome.fa

