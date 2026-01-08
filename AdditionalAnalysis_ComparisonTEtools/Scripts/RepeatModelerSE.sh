#!/bin/bash
#Run RepeatModeler (also see: Grandchamp et al, 2023)

#Repeat Modeler 
/global/projects/programs/source/RepeatModeler-open-1.0.11/BuildDatabase -name SW5 SW5finalGenome.fa
/global/projects/programs/source/RepeatModeler-open-1.0.11/RepeatModeler -database SW5 -pa 1

#RepeatMasker
/global/projects/programs/source/RepeatMasker/RepeatMasker -xsmall -lib consensi.fa.classified SW5finalGenome.fa -pa 1

