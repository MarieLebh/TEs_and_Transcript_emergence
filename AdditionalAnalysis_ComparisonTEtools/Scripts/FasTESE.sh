#!/bin/bash

#Run EDTA2 
perl EDTA//EDTA.pl --genome DmelSE_Genome.fa  --species others --sensitive 1 --threads 42 

#Reclassify the EDTA library
python DeepTE.py -d DeepTErun -i ../S1_EDTA/DmelSE_Genome.fa.mod.EDTA.TElib.fa -sp M -m_dir ../../Metazoans_model

#Clean the headers after running
sed -e 's/\(#\).*\(__\)/\1\2/' opt_DeepTE.fasta > DmelSEclean.fasta

#Run RepeatMasker
RepeatMasker ../S1_EDTA/DmelSE_Genome.fa -s -a -lib ../S2_DeepTE/DmelSEclean.fasta -xsmall

#For later analysis turn the output to a bedfile
RM2Bed.py DmelSE_Genome.fa.out

