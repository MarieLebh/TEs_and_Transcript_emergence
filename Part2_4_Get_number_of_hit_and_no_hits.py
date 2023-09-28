#!/usr/bin/env python

#import the necessary packages and modules
import pandas as pd
import os

"""
Script: Search for non expressed homologes
> Checks how many transcripts (unspliced) have a blast hit with 80 % coverage and 80 % identity (both changeable) in each population (that it was searched against)

- Input: 
	1. Blast tabular output file 
	2. Population file(tab-delimited, contains the transcript id, population and list of popultions not sharing the transcript)
- Output: 
	1. Csv file checking for each transcript in each population: 
		- if they have a match in that population (= 1) 
		- or not (= 0) 
		- or if they don't need one cause the populations share a transcript (= s)

"""

#Input files 
population_file = "All_merged_blast_populations.csv"
blast_file = "blast_all_output.txt"


"""
FILTER MATCHES:
Blast output parameters: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qcovs qcovhsp
"""

def filter_matches(blast_file):
        with open(blast_file) as file_in:
                qseqid = [] #query_id
                sseqid = [] #target seq id
                pident = [] #percent identity
                alen = [] #alignment length
                mismatch = [] #number mismatches
                gapopen = [] #number of gap openings
                qstart = [] #start of alignment in query
                qend = [] #end of alignment in query
                sstart = [] #start of alignment in target
                send = [] #end of alignment in target
                evalue = [] #e(xpect) value
                bitscore = [] #bitscore
                population_query = [] #query population
                population_hit = [] #target population
                for line in file_in:
                        l = line.split()
			#Now we check that all conditions are fullfilled
 			#Filter for coverage (l14 and percent  identity l2)
                        if float(l[14]) >= 80 and float(l[2]) >= 80:
                              qseqid.append(l[0])
                              sseqid.append(l[1])
                              pident.append(l[2])
                              alen.append(l[3])
                              mismatch.append(l[4])
                              gapopen.append(l[5])
                              qstart.append(l[6])
                              qend.append(l[7])
                              sstart.append(l[8])
                              send.append(l[9])
                              evalue.append(l[10])
                              bitscore.append(l[11])
                              population_query.append(l[0].split('::')[-1])
                              population_hit.append(l[1].split('_', 1)[0])
                #All query-subject matches that fullfill the conditions are appended to a df
                df = pd.DataFrame()
                df["qseqid"] = qseqid
                df["sseqid"] = sseqid
                df["pident"] = pident
                df["alignment_length"] = alen
                df["mismatch"] = mismatch
                df["gapopen"] = gapopen
                df["qstart"] = qstart
                df["qend"] = qend
                df["sstart"] = sstart
                df["send"] = send
                df["evalue"] = evalue
                df["bitscore"] = bitscore
                df["population_query"] = population_query
                df["population_hit"] = population_hit
		#Now a second df is made where for each unique seq id it is stored in which unique population it has a match
                df2 = pd.DataFrame({"qseqid":df.qseqid.unique()})
                df2["non_expressed_homologue"] = [list(set(df["population_hit"].loc[df["qseqid"] == x["qseqid"]]))
                    for _, x in df2.iterrows()]
                print(df2.head())
        return df2 #This is then returned


def merge_dataframes(df2, population_file):
     #colnames=["qseqid", "not_shared", "number_shared", "number_not_shared", "population"]
      df1 = pd.read_csv(population_file, sep = "\t")
      print(df1)
      df = pd.merge(df1,df2,on = "qseqid", how = "left") #Merge the two dataframes (one containing the transcripts that were blasted and the other one the blast matches)
      #Turn all NAs to an empty list (= transcripts that have no blastp hit either because they have none or because they had a match in all populations)
      df.loc[df["non_expressed_homologue"].isnull(),["non_expressed_homologue"]] = df.loc[df["non_expressed_homologue"].isnull(),"non_expressed_homologue"].apply(lambda x: [])
      df = df.drop(columns =  ["Number_of_populations_sharing_transcript", "Number_of_populations_without_the_transcript"]) #drop unecessary columns
      df.to_csv("temp_file.txt", index = None, header = None, sep = "\t") #save this df in a temporary file (will be deleted after use)
      return df

def make_final_file():
      File = "temp_file.txt" #temporary file from before
      Out = open("Non_expressed_homologues.txt", "w") #open the output file
      Out.write("transcript_id,population,matchAK5,matchDK5,matchGI5,matchSW5,matchUM,matchYE,matchZamb" + "\n") #colnames of the output file
      AK5 = 0
      DK5 = 0
      GI5 = 0
      SW5 = 0
      UM = 0
      YE = 0
      Zamb = 0
      with open(File) as File_in:
            for line in File_in:
                  l = line.split("\t")
                  Out.write(l[0] + "," + l[2]+ ",") #write seq id and population
                  for item in l[1]:
		  #check for each population that had to be blasted if it had a hit (=1) or not (=0)
                        if "AK5" in l[1]:
                              if "AK5" in l[3]:
                                     AK5 = "1"
                              else:
                                     AK5 = "0"
                        if "DK5" in l[1]:
                              if "DK5" in l[3]:
                                     DK5 = "1"
                              else:
                                     DK5 = "0"
                        if "GI5" in l[1]:
                              if "GI5" in l[3]:
                                     GI5 = "1"
                              else:
                                     GI5 = "0"
                        if "SW5" in l[1]:
                              if "SW5" in l[3]:
                                     SW5 = "1"
                              else:
                                     SW5 = "0"
                        if "UM" in l[1]:
                              if "UM" in l[3]:
                                     UM = "1"
                              else:
                                     UM = "0"
                        if "YE" in l[1]:
                              if "YE" in l[3]:
                                     YE = "1"
                              else:
                                     YE = "0"
                        if "Zamb" in l[1]:
                              if "Zamb" in l[3]:
                                     Zamb = "1"
                              else:
                                     Zamb = "0"
			#if a population was not in the list to be blasted against it means that it shared a transcript with the query
                        if "AK5" not in l[1]:
                              AK5 = "s"
                        if "DK5" not in l[1]:
                              DK5 = "s"
                        if "GI5" not in l[1]:
                              GI5 = "s"
                        if "SW5" not in l[1]:
                              SW5 = "s"
                        if "UM" not in l[1]:
                              UM = "s"
                        if "YE" not in l[1]:
                              YE = "s"
                        if "Zamb" not in l[1]:
                              Zamb = "s"
                  Out.write(AK5 + "," + DK5 + ","+ GI5 + "," + SW5 + "," + UM + "," + YE + "," + Zamb + "\n") #write the results in the outfile
                  AK5 = 0 #reset the parameters
                  DK5 = 0
                  GI5 = 0
                  SW5 = 0
                  UM = 0
                  YE = 0
                  Zamb = 0
      Out.close()

#Merge all functions created before

def final_function(blast_file, population_file):
      df2 =  filter_matches(blast_file)
      merge_dataframes(df2, population_file)
      make_final_file()
      os.remove("temp_file.txt") #remove the temporary file

#run this in python

final_function(blast_file, population_file)

