#!/usr/bin/env python

"""
Count the number of all individual motifs in the sequence upstream regions (-1000/+100)

Input: Result csv file from motif search
"""

#import the necessary packages and modules
import csv
import pandas as pd


def get_motif_per_transcript(df, population, t):
    df = pd.read_csv(df, sep = "\t")
    df = df[df["position"] >= 0]
    df["population"] = population
    df["type"] = t
    df["seq_name"] = df["seq_name"].apply(str)
    df["unique_id"] = df["seq_name"] + df["population"]
    df["unique_id2"] = df["unique_id"] + df["motif_id"]
    df["number"] = 0
    df["number_high"] = 0
    df.loc[df["rel_score"] >= 0.8 ,  "number"] = 1
    df.loc[df["rel_score"] >= 0.95 , "number_high"] = 1
    df['number_per_transcript'] = df.groupby(df['unique_id2'].ne(df['unique_id2'].shift()).cumsum())['number'].transform('sum')
    df['number_per_transcript_high'] = df.groupby(df['unique_id2'].ne(df['unique_id2'].shift()).cumsum())['number_high'].transform('sum')
    df = df.drop_duplicates(subset=['unique_id2'])
    final = df[["motif_id", "seq_name", "type", "unique_id", "population", "number_per_transcript", "number_per_transcript_high"]]
    return final

#Run like this

name = "Motif_number_of_intergenic_de_novo_transcripts.csv"
t = "intergenic_de_novo_transcript"

AK5 ="AK5_intergenic_transcripts_motifs.csv"
population = "AK5"
ak5 = get_motif_per_transcript(AK5, population, t)

DK5= "DK5_intergenic_transcripts_motifs.csv"
population = "DK5"
dk5= get_motif_per_transcript(DK5, population, t)

GI5=  "GI5_intergenic_transcripts_motifs.csv"
population = "GI5"
gi5 = get_motif_per_transcript(GI5, population, t)

SW5 = "SW5_intergenic_transcripts_motifs.csv"
population = "SW5"
sw5 = get_motif_per_transcript(SW5, population, t)

UM = "UM_intergenic_transcripts_motifs.csv"
population = "UM"
um = get_motif_per_transcript(UM, population, t)

YE = "YE_intergenic_transcripts_motifs.csv"
population = "YE"
ye = get_motif_per_transcript(YE, population, t)

Zamb =  "Zamb_intergenic_transcripts_motifs.csv"
population = "Zamb"
zamb = get_motif_per_transcript(Zamb, population, t)

final = pd.concat([ak5, dk5, gi5, sw5, um, ye, zamb])
final.to_csv(name , index = False, sep = "\t")


#Now correct the motif number (for cases where a motif was not at all detected even with a 0.7 score)

df = pd.read_csv("Motif_number_of_intergenic_de_novo_transcripts.csv",delim_whitespace=True)
df1 = df.set_index(['unique_id','motif_id']).unstack(fill_value=0).stack().reset_index()
df = df1
df[['id','position']] = df['unique_id'].str.split('::',expand=True)
df[['chrom','position']] = df['position'].str.split(':',expand=True)
df[['position','population']] = df['position'].str.split('(',expand=True)
df[['start','end']] = df['position'].str.split('-',expand=True)
df['population'] = df['population'].str[2:]
df['size'] = pd.to_numeric(df['end']) - pd.to_numeric(df['start'])
df['type'] = 'de_novo_transcript'

df.to_csv("Motif_number_de_novo_corrected.csv")


