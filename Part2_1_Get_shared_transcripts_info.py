#!/usr/bin/env python

"""
Filter the blast matches and get all populations to blast against
This script will check if a transcripts is shared with another populations. 
The output will be a file containing a list for each transcripts containing all populations that do NOT share the transcripts.
The second output will be a file containing for each transcript in how many populations it has a hit matching the criteria from before.
"""

#import the necessary packages and modules
import csv
import pandas as pd

#For each population two inputs are needed:
#1: Blast output (of blasting spliced transcripts against all populations)
#2: Text file with all de novo transcript ids

"""
These are the parameters of the blast output files (in this order):

1. qseqid: id of query
2. qlen: length of query
3. sseqid: id of subject (=hit)
4. slen: length of the subject
5.pident: percent identity
6.length: length of alignment
7. mismatch: number of mismatches
8. gapopen: number of gaps
9. qstart: start alignment in query
10. qend: end alignment in query
11. sstart: start alignment in subject
12. send: end alignment in subject
13. evalue: E-value
14. bitscore: Bitscore

"""

"""
FILTER MATCHES:

Input: blast tabular output file

The following criteria need to be fullfilled:
> When filtering only for a shared end exclude condition 1 and 2
> When filtering for a shared start exclude condition 3-7

1. Start alignment query <= 200
2. Start alignment target <= 200
3. Length query - 200 <= End alignment query <= Length query +200
4. Length query - 200 <= End alignment target <= Length query +200
5. Length target - 200 <= End alignment query <= Length target +200
6. Length target - 200 <= End alignment target <= Length target +200
7. Abs(Length target - length query) < 200

The number 200 can be replaced by a different (smaller/higher) number if required.

Output: pandas df (containing only filtered matches), list of all ids with a hit
"""

def filter_matches(blast: "blast tabular output file")-> "df":
        with open(blast) as file_in:
                qseqid = [] #query_id
                qlen = [] #query length (= end of query)
                sseqid = [] #target seq id
                slen = [] #target seq length
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
                        if l[0] != l[2]  and abs(int(l[1]) - int(l[3])) < 200 and float(l[4]) >= 80 and int(l[1])-200 <= int(l[11]) <= int(l[1])+200 and int(l[1])-200 <= int(l[9]) <= int(l[1])+200 and int(l[10]) <= 200 and int(l[3])-200 <= int(l[11]) <= int(l[3])+200 and int(l[8]) <= 200 and int(l[3])-200 <= int(l[9]) <= int(l[3])+200:
                              qseqid.append(l[0])
                              qlen.append(l[1])
                              sseqid.append(l[2])
                              slen.append(l[3])
                              pident.append(l[4])
                              alen.append(l[5])
                              mismatch.append(l[6])
                              gapopen.append(l[7])
                              qstart.append(l[8])
                              qend.append(l[9])
                              sstart.append(l[10])
                              send.append(l[11])
                              evalue.append(l[12])
                              bitscore.append(l[13])
                              population_query.append(l[0].split('::')[-1])
                              population_hit.append(l[2].split('::')[-1])
                #All query-subject matches that fullfill the conditions are appended to a df
                df = pd.DataFrame()
                df["qseqid"] = qseqid
                df["qlen"] = qlen
                df["sseqid"] = sseqid
                df["slen"] = slen
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
        return qseqid, df #returns a list of all hits as well as a df with all hits

"""
LIST DE NOVO

This function returns a list of all de novo transcripts from a text file containing the respective ids.

Input: file containing ids of all de novo transcripts, population

Output: list containing all de novo transcripts (with added population)
"""

def list_de_novo(de_novo, population)-> "list":
        with open(de_novo) as file_in:
                lines = [line.rstrip() for line in file_in]
                lines2 = {i+"::" + population for i in lines}
        return lines2

"""
ADD NO MATCHES

This function compares the two lists (ids with blast match vs de novo transcripts)
with each other and return a new list with only the ones that are not in the list with all hits.

It then merges the dataframe with all hits with a second df with all transcript ids without any hit.

This df will then be used in the following function.

Input: list hits, list with all de novo transcripts, filtered df, population

Output: df where all de novo ids with no match are added (with empty columns except for id and population)
"""

def add_no_matches(a, b, df, population):
    a = set(a) #list of all transcripts
    b = set(b) # list of all blast hits
    c = list(a - b) #list of all transcript without blast hits
    population = len(c) *  [population]
    df2 = pd.DataFrame()
    df2["qseqid"] = c
    df2["population_query"] = population
    final_df = pd.concat([df, df2])
    return final_df

"""
GET ALL HIT POPULATIONS

Input:

-df = a dataframe containing all transcripts and their filtered hits (if they have any)
- population = the population (this population will be substracted from the list of populations to blast against)

The function  checks for each transcript in which population it has at least one hit and saves these in a list.
If it doesn't have any hits this list will be empty.

In a second step the program takes a list of all populations and substracts the hit transcripts (and the transcript population).

Output:

-A df (final) containing the query id and a list of populations with no blast hit.
- A second df (df) containing the number of shared and non shared populations (can be analyzed in R).

"""

def get_all_hit_populations(df, population):
    df2 = pd.DataFrame({"qseqid":df.qseqid.unique()})
    df2['population_hit'] = [list(set(df['population_hit'].loc[df['qseqid'] == x['qseqid']]))
        for _, x in df2.iterrows()]
    x = ['AK5', 'DK5', 'GI5', 'SW5', 'UM', 'YE', 'Zamb']
    x.remove(population)
    l = len(df2) *[x]
    df2['control'] = l
    df2['population_no_hit'] = df2['control'].map(set) - df2['population_hit'].map(set)
    df2['populations_to_blast_against'] = df2['population_no_hit'].apply(list)
    final = df2.drop(columns=['control', 'population_hit', 'population_no_hit'])
    df = final
    df['Number_of_populations_sharing_transcript'] = 6 - df['populations_to_blast_against'].str.len()
    df['Number_of_populations_without_the_transcript'] = df['populations_to_blast_against'].str.len()
    df['Population'] = population
    return final, df


#Run this in python

#AK5
blast = "matches_AK5.txt"

b, df  = filter_matches(blast)

de_novo = "AK5_de_novo_transcripts.txt"
population = "AK5"

a = list_de_novo(de_novo, population)

final_dfa = add_no_matches(a, b, df, population)

final_df2, p1  = get_all_hit_populations(final_dfa, population)

final_df2.to_csv("AK5_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes

#DK5
blast = "matches_DK5.txt"
de_novo = "DK5_de_novo_transcripts.txt"
population = "DK5"

b, df  = filter_matches(blast)

a = list_de_novo(de_novo, population)

final_dfb = add_no_matches(a, b, df, population)

final_df2, p2 = get_all_hit_populations(final_dfb, population)

final_df2.to_csv("DK5_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes

#GI5
blast = "matches_GI5.txt"
de_novo = "GI5_de_novo_transcripts.txt"
population = "GI5"

b, df  = filter_matches(blast)

a = list_de_novo(de_novo, population)

final_dfc = add_no_matches(a, b, df, population)

final_df2, p3 = get_all_hit_populations(final_dfc, population)

final_df2.to_csv("GI5_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes

#SW5
blast = "matches_SW5.txt"
de_novo = "SW5_de_novo_transcripts.txt"
population = "SW5"

b, df  = filter_matches(blast)

a = list_de_novo(de_novo, population)

final_dfd = add_no_matches(a, b, df, population)

final_df2, p4 = get_all_hit_populations(final_dfd, population)

final_df2.to_csv("SW5_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes

#UM
blast = "matches_UM.txt"
de_novo = "UM_de_novo_transcripts.txt"
population = "UM"

b, df  = filter_matches(blast)

a = list_de_novo(de_novo, population)

final_dfe = add_no_matches(a, b, df, population)

final_df2, p5 = get_all_hit_populations(final_dfe, population)

final_df2.to_csv("UM_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes

#YE
blast = "matches_YE.txt"
de_novo = "YE_de_novo_transcripts.txt"
population = "YE"

b, df  = filter_matches(blast)

a = list_de_novo(de_novo, population)

final_dff = add_no_matches(a, b, df, population)

final_df2, p6 = get_all_hit_populations(final_dff, population)

final_df2.to_csv("YE_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes

#Zamb
blast = "matches_Zamb.txt"
de_novo = "Zamb_de_novo_transcripts.txt"
population = "Zamb"

b, df  = filter_matches(blast)

a = list_de_novo(de_novo, population)

final_dfg = add_no_matches(a, b, df, population)

final_df2, p7 = get_all_hit_populations(final_dfg, population)

final_df2.to_csv("Zamb_file_for_blast.csv", index = False, sep = "\t") #this will be used for the blast against the other genomes


#All dataframes are merged together.
#This file then shows for each transcript in how many populations it has a hit matching the criteria from before.
all = pd.concat([p1, p2, p3, p4, p5, p6, p7])
all.to_csv("Orthogroup_file1.csv", index = False, sep = "\t") 

#Merge the dataframes containing all filtered hits together to build a final one with all populations
#This script is necessary for the second script "build_orthogroups.py"
i =  pd.concat([final_dfa, final_dfb, final_dfc, final_dfd, final_dfe, final_dff, final_dfg])
i.to_csv("transcripts_with_blast.csv", index = False, sep = "\t")
