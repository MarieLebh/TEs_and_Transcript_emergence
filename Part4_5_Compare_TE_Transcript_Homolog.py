#!/usr/bin/env python

"""
Compare transcript and homolog TEs 

"""

#Input files (specify here)
transcripts = "All_transcripts_TE_overlap.bed" #file with info on transcript TE overlap (bedtools intersect transcripts and TE)
homologs = "All_ne_homologues_TE_overlap.bed" #file with info on homolog TE overlap (bedtools intersect homologs and TE)

#import the necessary packages and modules
import pandas as pd
import os
import numpy as np

#Merge the transcript file with the homolog file

def merge_dataframes(transcripts, homologs):
      colnames=["chrom", "start", "end", "id", "score", "strand", "chrom2", "start2", "end2", "id2", "score2", "strand2", "overlap"]
      #Prepare transcript df
      tr = pd.read_csv(transcripts, sep = "\t", names = colnames)
      tr[["te","class", "description"]] = tr['id2'].str.split(';',expand=True)
      tr2 = pd.DataFrame({"id":tr.id.unique()})
      tr2["tes_transcript"] = [list(set(tr["class"].loc[tr["id"] == x["id"]])) #get list with all TE per transcript
            for _, x in tr2.iterrows()]
      #Prepare homolog df
      ho = pd.read_csv(homologs, sep = "\t", names = colnames)
      ho[["te","class", "description"]] = ho["id2"].str.split(';',expand=True)
      ho2 = pd.DataFrame({"id":ho.id.unique()})
      ho2["tes_homolog"] = [list(set(ho["class"].loc[ho["id"] == x["id"]])) #get list with all TE per homolog
            for _, x in ho2.iterrows()]
      ho2[["type","id", "pop_homolog", "pop_transcript"]] = ho2['id'].str.split('_',expand=True)
      df = pd.merge(tr2,ho2,on = "id", how = "right") #Merge the two dataframes 
      #Remove all None characters
      df["tes_transcript"] = df["tes_transcript"].apply(lambda el: [x for x in el if pd.notna(x)])
      df["tes_homolog"] = df["tes_homolog"].apply(lambda el: [x for x in el if pd.notna(x)])
      df.to_csv("temp_file.txt", index = None, header = None, sep = "\t") #save this df in a temporary file (will be deleted after use)
      #return df

#Makes the outfile and sorts the transcript - homolog pairs into categories 

def make_final_file():
      File = "temp_file.txt" #temporary file from before
      Out = open("TE_analysis_homologs.txt", "w") #open the output file
      Out.write("transcript_id,category,population_homolog,population_transcript" + "\n") #colnames of the output file
      category = 0
      with open(File) as File_in:
            for line in File_in:
                  l = line.split("\t")
                  temp = l[1] #Turn list column back to list (as it was recognized as string)
                  if temp == "[]":
                        l[1] = []
                  else:
                        l[1] = [i.strip() for i in l[1][1:-1].replace('"',"").split(',')]
                  temp = l[2] #Turn 2nd list column back to list (as it was recognized as string)
                  if temp == "[]":
                        l[2] = []
                  else:
                        l[2] = [i.strip() for i in l[2][1:-1].replace('"',"").split(',')]
                  if not l[1] and not l[2]:
                        category = "No TE in transcript and homolog"
                  elif not l[2]:
                        category = "TE(s) in transcript but not homolog"
                  elif not l[1]:
                        category = "TE(s) in homolog but not transcript"
                  elif sorted(l[1]) == sorted(l[2]):
                        category = "Same TE(s) in transcript and homolog"
                  elif any(i in l[1] for i in l[2]):
                        category = "Transcript and homolog share some but not all TEs"
                  else:
                        category = "Transcript and homolog have different TEs"
                  Out.write(l[0] + "," + category + ","  + l[4] + ","+ l[5]) #write the results in the outfile
                  category = 0 #reset the parameters
      Out.close()

#Merge all functions created before

def final_function(transcripts, homologs):
      merge_dataframes(transcripts, homologs)
      make_final_file()
      os.remove("temp_file.txt")

#Run this
final_function(transcripts, homologs)
