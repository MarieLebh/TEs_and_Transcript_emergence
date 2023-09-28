#!/usr/bin/env python

"""
Search for motifs in upstream regions
"""

#import the necessary packages and modules
from multiprocessing import Pool
from Bio import motifs
import csv
import pandas as pd
import random


"""
This function takes a fasta file and stores all sequences and their corresponding header as a dictionary.

fasta = fasta file containing a set of DNA sequences

"""
def create_seq_dict(fasta: "fasta_file")-> "dictionary":
        f=open(fasta)
        seqs={}
        for line in f:
                line=line.rstrip()
                if line[0]=='>':
                        words=line.split()
                        name=words[0][1:] #get transcript id 
                        seqs[name]=''
                else:
                        seqs[name]=seqs[name]+line #add sequence
        return seqs #output dictionary

"""
This function will be used for the parallelisation of the different files.
The output will be a list of tuples containing the key and value in a list paired with the db name (always the same for each key value pair)


seqs = sequence dictionary (key = name, value = sequence)
db = name/path to the database for the motif search

"""

def create_parameter_list(seqs, db):
    seqs_tuples = [[key] + [val] for key, val in seqs.items()]
    db_list = [db] * len(seqs_tuples)
    parameter_list = list(zip(seqs_tuples, db_list))
    return parameter_list

"""
This function calculates the occurence of a set of motifs in a set of dna seqs.

seqname = sequence id
fasta = a fasta sequence
file = a jaspar motifs file (txt format) containing one or more pfm matrices of TF motifs
        - databases to use: JASPAR2022_CORE_insects_non-redundant_pfms_jaspar or JASPAR2020_POLII

"""

def searchmotifs(seqname, fasta, file):
    with open(file) as f:
        sequence = fasta
        seq_name = seqname
        motif_id = []
        positions = []
        scores = []
        rel_score = []
        threshold = 0
        for m in motifs.parse(f, "jaspar"):
            pseudocounts = motifs.jaspar.calculate_pseudocounts(m)
            pwm = m.counts.normalize(pseudocounts)
            background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25} 
            pssm = pwm.log_odds(background) #*
            #mean = pssm.mean()
            #std = pssm.std()
            max_score = pssm.max
            min_score = pssm.min
            threshold = (max_score - min_score) * 0.7 + min_score  
            for position, score in pssm.search(sequence, threshold):
                    srel = (score - min_score) / (max_score - min_score)
                    rel_score.append(srel)
                    positions.append(position)
                    scores.append(score)
                    motif_id.append(m.name)
    df = pd.DataFrame()
    df["position"] = positions
    df["score"] = scores
    df["rel_score"] = rel_score
    df["motif_id"] = motif_id
    df["seq_name"] = seq_name
    return df

"""
This is the multiwrapper function.
It takes a tuple as a parameter and splits it into the parameters needed for the motif search.
It then applies the motif search function.
"""

def multi_wrapper(t): #t = tuple with parameter
        key_value, db = t
        key = key_value[0]
        value = key_value[1]
        x = searchmotifs(key, value, db) #Funktion hier ausführen
        return x



#Run like this

#AK5
fasta = "/global/students/research/m_lebh01/Bedtools/Motif_search_fasta_files/Fasta_files/AK5_intergenic_transcripts_upstream.fa" #fasta sequences
db = "JASPAR_merged_motifs_db.txt" #motif database 
seqs = create_seq_dict(fasta)
parameter_list = create_parameter_list(seqs, db)
p = Pool(20)  #number of cores #write the same number in the slurm script
results_variable = p.map(multi_wrapper, parameter_list)  #apply the multi wrapper function (result = pd.df)
p.close()
p.join()
results_df = pd.concat(results_variable) #join all pd.dfs for every sequence together to one big df with all sequences
results_df = results_df.to_csv("AK5_intergenic_transcripts_motifs.csv", index = False, sep = "\t")




