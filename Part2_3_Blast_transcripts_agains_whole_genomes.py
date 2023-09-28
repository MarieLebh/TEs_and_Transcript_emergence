#!/usr/bin/env python

"""
Find transcript on expressed homologous in the other genomes
"""

#import all necessary modules
import subprocess
import pandas as pd
from Bio import SeqIO
import os

"""
This function takes a genome fasta file as an input.
It returns a fasta file where the population is added to the seqname (=chromosome).
This will later be helpful to know from which population a corresponding blast hit is.
"""

def make_fasta_for_db(Fasta, population, filename)-> "modified fasta file":
      with open(filename, "w") as handle:
            records = SeqIO.parse(Fasta, "fasta")
            for record in records:
                  record.id = population + "_" + record.id
                  record.description = record.id
                  SeqIO.write(record, handle, "fasta")

"""
This function creates seven different blast databases for the given populations.
"""

def create_blast_databases(AK5, DK5, GI5, SW5, UM, YE, Zamb):
      subprocess.call("makeblastdb -in " + AK5 + " -dbtype nucl", shell=True)
      subprocess.call("makeblastdb -in " + DK5 + " -dbtype nucl", shell=True)
      subprocess.call("makeblastdb -in " + GI5 + " -dbtype nucl", shell=True)
      subprocess.call("makeblastdb -in " + SW5 + " -dbtype nucl", shell=True)
      subprocess.call("makeblastdb -in " + UM + " -dbtype nucl", shell=True)
      subprocess.call("makeblastdb -in " + YE + " -dbtype nucl", shell=True)
      subprocess.call("makeblastdb -in " + Zamb + " -dbtype nucl", shell=True)

"""
This function takes a seq_id and a fasta file as an input.
This function checks if a given fasta file contains a certain seq id.
If it does, it returns the id and its corresponding fasta sequence.
"""

def getSeqRec(seq_id, fasta):
        fasta=SeqIO.parse(fasta,"fasta")
        for record in fasta:
                if seq_id==record.id:
                        return ">" + record.id + "\n" + record.seq

"""
This function is the blast function. It takes two files as an input:

population file = a file with all transcript ids and a list of all population that it should be blasted against 
de_novo_fa = a multi-fasta file containing all de novo transcript sequences

This function then opens the population file and stores the corresponding fasta sequence in a temporary file.
This fasta sequence is then blasted against all genomes that were in the list of the population file.
So if a transcript has no population in the list, no blast will be run.

This part of the script should be run via slurm/the terminal. The results will then automatically be stored in the (slurm) output file as tabular output.
The output format will have the following columns:

qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp

"""
def blast_seqs(population_file, de_novo_fa):
      with open(population_file) as file_in:
          lines = file_in.readlines()[1:]
          for line in lines:
                l = line.split("\t")
                fasta = getSeqRec(l[0], de_novo_fa)
                with open("query_temp.fa", "w") as query_file:
                      fasta2 = str(fasta)
                      header, seq = fasta2.split("\n", 1)
                      query_file.write(header + "\n")
                      query_file.write(seq)
                      query_file.close()
                if "AK5" in l[1]:
                     subprocess.call("blastn -db AK5Genome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)
                if "DK5" in l[1]:
                      subprocess.call("blastn -db DK5Genome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)
                if "GI5" in l[1]:
                      subprocess.call("blastn -db GI5Genome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)
                if "SW5" in l[1]:
                      subprocess.call("blastn -db SW5Genome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)
                if "UM" in l[1]:
                      subprocess.call("blastn -db UMGenome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)
                if "YE" in l[1]:
                      subprocess.call("blastn -db YEGenome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)
                if "Zamb" in l[1]:
                      subprocess.call("blastn -db ZambGenome.fa -query query_temp.fa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp' -num_threads 5", shell = True)


#SPECIFY ALL INPUT FILES HERE
#path to de novo fasta files:
de_novo_fa = "/global/students/research/m_lebh01/Synteny_Analysis/Prep_for_blast/Fasta_files/AK5/AK5_de_novo.fa"

#STEP1: CREATE THE BLAST DATABASES FOR THE BLAST SEARCH
#path to all genomes:
AK5_genome = "AK5finalGenome.masked.fa"
DK5_genome = "DK5finalGenome.masked.fa"
GI5_genome = "GI5finalGenome.masked.fa"
SW5_genome = "SW5finalGenome.masked.fa"
UM_genome = "UMfinalGenome.masked.fa"
YE_genome = "YEfinalGenome.masked.fa"
Zamb_genome = "ZambfinalGenome.masked.fa"


AK5 = "AK5Genome.fa"
DK5 ="DK5Genome.fa"
GI5 = "GI5Genome.fa"
SW5 = "SW5Genome.fa"
UM = "UMGenome.fa"
YE = "YEGenome.fa"
Zamb = "ZambGenome.fa"

make_fasta_for_db(AK5_genome, "AK5", AK5)
make_fasta_for_db(DK5_genome, "DK5", DK5)
make_fasta_for_db(GI5_genome, "GI5", GI5)
make_fasta_for_db(SW5_genome, "SW5", SW5)
make_fasta_for_db(UM_genome, "UM", UM)
make_fasta_for_db(YE_genome, "YE", YE)
make_fasta_for_db(Zamb_genome, "Zamb", Zamb)

#create_blast_databases(AK5, DK5, GI5, SW5, UM, YE, Zamb)

#STEP2: Run the blasts:
#path to the population files:
population_file_AK5 = "AK5_file_for_blast.csv"
population_file_DK5 = "DK5_file_for_blast.csv"
population_file_GI5 = "GI5_file_for_blast.csv"
population_file_SW5 = "SW5_file_for_blast.csv"
population_file_UM = "UM_file_for_blast.csv"
population_file_YE = "YE_file_for_blast.csv"
population_file_Zamb = "Zamb_file_for_blast.csv"

#path to de novo fasta files:
de_novo_fa_AK5 = "AK5_de_novo.fa"
de_novo_fa_DK5 = "DK5_de_novo.fa"
de_novo_fa_GI5 = "GI5_de_novo.fa"
de_novo_fa_SW5 = "SW5_de_novo.fa"
de_novo_fa_UM = "UM_de_novo.fa"
de_novo_fa_YE = "YE_de_novo.fa"
de_novo_fa_Zamb = "Zamb_de_novo.fa"

#Run this 
blast_seqs(population_file_AK5, de_novo_fa_AK5)
blast_seqs(population_file_DK5, de_novo_fa_DK5)
blast_seqs(population_file_GI5, de_novo_fa_GI5)
blast_seqs(population_file_SW5, de_novo_fa_SW5)
blast_seqs(population_file_UM, de_novo_fa_UM)
blast_seqs(population_file_YE, de_novo_fa_YE)
blast_seqs(population_file_Zamb, de_novo_fa_Zamb)


