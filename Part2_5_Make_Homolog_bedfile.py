#!/usr/bin/env python

"""
Filter blast output:

-Input files: 

Blast output file
Population file for Blast (contains all transcripts and what populations need to be searched for in blast)
"""

#import the necessary packages and modules
import pandas as pd
import os
import numpy as np
import subprocess

#Exclude all blast matches that do not meet the wanted criteria

def filter_matches(blast_file):
        with open(blast_file) as file_in:
                qseqid = [] #query_id
                sseqid = [] #target seq id
                pident = [] #percent identity
                evalue = [] #E value
                qstart = [] #start of alignment in query
                qend = [] #end of alignment in query
                sstart = [] #start of alignment in target
                send = [] #end of alignment in target
                coverage = [] #coverage
                strand = [] #strand
                population_query = [] #query population
                population_hit = [] #target population
                for line in file_in:
                        l = line.split()
                        #Now we check that all conditions are fullfilled
                        #Filter for coverage (l14 and percent  identity l2)
                        if float(l[14]) >= 80 and float(l[2]) >= 80:
                              qseqid.append(l[0])
                              sseqid.append(l[1].split('_', 1)[1])
                              pident.append(l[2])
                              qstart.append(l[6])
                              qend.append(l[7])
                              evalue.append(l[10])
                              sstart.append(l[8])
                              send.append(l[9])
                              if l[12] == "plus":
                                    strand.append("+")
                              elif l[12] == "minus":
                                    strand.append("-")
                              coverage.append(l[14])
                              population_query.append(l[0].split('::')[-1])
                              population_hit.append(l[1].split('_', 1)[0])
                #All query-subject matches that fullfill the conditions are appended to a df
                df = pd.DataFrame()
                df["qseqid"] = qseqid
                df["sseqid"] = sseqid
                df["pident"] = pident
                df["qstart"] = qstart
                df["qend"] = qend
                df["sstart"] = sstart
                df["send"] = send
                df["cov"] = coverage
                df["strand"] = strand
                df["e_value"] = evalue
                df["population_query"] = population_query
                df["population_hit"] = population_hit
                #Now a second df is made where for each unique seq id it is stored in which unique population it has a match
                df.to_csv("filtered_blast_hits.txt", index = None, header = None, sep = "\t")

#Make a bedfile for the homologues 
#Homologues from the same population will be stored in one bedfile (to be able to run bedtools...)

def merge_dataframes(best_homologues, population_file):
      my_colnames = ["qseqid", "chromosome", "pident", "qstart", "qend", "sstart", "send", "cov", "strand2", "evalue", "pop_query", "pop_hit"]
      df1 = pd.read_csv(population_file, sep = "\t")
      df2 = pd.read_csv(best_homologues, names = my_colnames, sep = "\t")
      df = pd.merge(df1,df2,on = "qseqid", how = "right") #Merge the two dataframes (one containing the transcripts that were blasted and the other one the blast matches)
      df["sstart2"] = np.where(df["strand2"]== "-", df["send"] , df["sstart"])
      df["send2"] = np.where(df["strand2"]== "-", df["sstart"] , df["send"])
      df.to_csv("NE_homologues.bed", index = None, header = None, sep = "\t")
      df["ID"] = "NonExpressedHomolog" + "_" + df["qseqid"] +  "_" + df['pop_hit'] + "_" +  df['pop_query'] #Take id of the corresponding transcript 
      df["sstart"] = pd.to_numeric(df["sstart"])
      df["send"] = pd.to_numeric(df["send"])
      #Split the df by population and make a bedfile for the homologues
      #AK5
      AK5 = df[df['pop_hit'] == "AK5"]
      AK5 = AK5[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      AK5.to_csv("NE_homologuesAK51.bed", index = None, header = None, sep = "\t")
      #DK5
      DK5 = df[df['pop_hit'] == "DK5"]
      DK5 = DK5[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      DK5.to_csv("NE_homologuesDK51.bed", index = None, header = None, sep = "\t")
      #GI5
      GI5 = df[df['pop_hit'] == "GI5"]
      GI5 = GI5[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      GI5.to_csv("NE_homologuesGI51.bed", index = None, header = None, sep = "\t")
      #SW5
      SW5 = df[df['pop_hit'] == "SW5"]
      SW5 = SW5[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      SW5.to_csv("NE_homologuesSW51.bed", index = None, header = None, sep = "\t")
      #UM
      UM = df[df['pop_hit'] == "UM"]
      UM = UM[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      UM.to_csv("NE_homologuesUM1.bed", index = None, header = None, sep = "\t")
      #YE
      YE = df[df['pop_hit'] == "YE"]
      YE = YE[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      YE.to_csv("NE_homologuesYE1.bed", index = None, header = None, sep = "\t")
      #Zamb
      Zamb = df[df['pop_hit'] == "Zamb"]
      Zamb = Zamb[["chromosome", "sstart2", "send2", "ID", "score", "strand2"]]
      Zamb.to_csv("NE_homologuesZamb1.bed", index = None, header = None, sep = "\t")
      return df


#Run this in python

#Input files 
population_file = "all_de_novo_transcripts.csv" #File that was created for the blast search
blast_file = "blast_all_output.txt" #Blast output file

filter_matches(blast_file)

#After this run this line on the output:
#This will sort the hits by e-value, percent identity and coverage and (if there are multiple sufficient blast hit) take the best. 
#If two hits would be exactly equal with regards to these parameters, the first will be taken.

subprocess.call("sort -k1,1 -k12,12 -k10,10 -rk3,3  -rk8,8  filtered_blast_hits.txt  | sort --merge -u  -k1,1 -k12,12 > best_non_expressed_homologues.txt", shell=True)

best_homologues = "best_non_expressed_homologues.txt"#File containing the best homologs matching the criteria

merge_dataframes(best_homologues, population_file)


