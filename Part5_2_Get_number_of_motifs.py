#!/usr/bin/env python

"""
Count the number of motifs in the sequence upstream regions (-1000/+100)

Input: Result csv file from motif search
"""

#import the necessary packages and modules
import csv
import pandas as pd

def merge_dfs(df, population, t):
    df = pd.read_csv(df, sep = "\t")
    df = df[df["position"] >= 0]
    df["population"] = population
    df["type"] = t
    df["seq_name"] = df["seq_name"].apply(str)
    df["unique_id"] = df["seq_name"] + df["population"]
    df["number_of_motifs"] = 0
    df["number_of_motifs_high"] = 0
    df["number_of_core"] = 0
    df["number_of_core_high"] = 0
    df["core_id"] = df["motif_id"]
    df.loc[~df["core_id"].isin(motifs_core), "core_id"] = "Not_core"
    df.loc[df["core_id"].isin(motifs_core), "core_id"] = "core"
    df.loc[(df["rel_score"] > 0.8) & (df["core_id"] == "Not_core"), "number_of_motifs"] = 1
    df.loc[(df["rel_score"] > 0.95) & (df["core_id"] == "Not_core"), "number_of_motifs_high"] = 1
    df.loc[(df["rel_score"] > 0.8) & (df["core_id"] == "core"), "number_of_core"] = 1
    df.loc[(df["rel_score"] > 0.95) & (df["core_id"] == "core"), "number_of_core_high"] = 1
    df['number_motifs_per_transcript'] = df.groupby(df['unique_id'].ne(df['unique_id'].shift()).cumsum())['number_of_motifs'].transform('sum')
    df['number_motifs_high_per_transcript'] = df.groupby(df['unique_id'].ne(df['unique_id'].shift()).cumsum())['number_of_motifs_high'].transform('sum')
    df['number_core_per_transcript'] = df.groupby(df['unique_id'].ne(df['unique_id'].shift()).cumsum())['number_of_core'].transform('sum')
    df['number_core_high_per_transcript'] = df.groupby(df['unique_id'].ne(df['unique_id'].shift()).cumsum())['number_of_core_high'].transform('sum')
    df = df.drop_duplicates(subset=['unique_id'])
    final = df[["seq_name","population", "unique_id", "type", "number_motifs_per_transcript", "number_motifs_high_per_transcript", "number_core_per_transcript", "number_core_high_per_transcript"]]
    return final

motifs_core = ["POL001.1.MTE", "POL002.1.INR", "POL003.1.GC-box", "POL004.1.CCAAT-box", "POL005.1.DPE", "POL006.1.BREu", "POL007.1.BREd", "POL008.1.DCE_S_I", "POL009.1.DCE_S_II", "POL010.1.DCE_S_III", "POL011.1.XCPE1", "POL012.1.TATA-Box", "POL013.1.MED-1"]


#Run like this
name = "Intergenic_de_novo_transcripts.csv"
t = "intergenic_de_novo_transcript"

AK5 ="/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/AK5_intergenic_transcripts_motifs.csv"
population = "AK5"
ak5 = merge_dfs(AK5, population, t)

DK5= "/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/DK5_intergenic_transcripts_motifs.csv"
population = "DK5"
dk5= merge_dfs(DK5, population, t)

GI5=  "/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/GI5_intergenic_transcripts_motifs.csv"
population = "GI5"
gi5 = merge_dfs(GI5, population, t)

SW5 = "/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/SW5_intergenic_transcripts_motifs.csv"
population = "SW5"
sw5 = merge_dfs(SW5, population, t)

UM = "/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/UM_intergenic_transcripts_motifs.csv"
population = "UM"
um = merge_dfs(UM, population, t)

YE = "/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/YE_intergenic_transcripts_motifs.csv"
population = "YE"
ye = merge_dfs(YE, population, t)

Zamb =  "/global/students/research/m_lebh01/MotifAnalysis/Transcript_motifs/Zamb_intergenic_transcripts_motifs.csv"
population = "Zamb"
zamb = merge_dfs(Zamb, population, t)

final = pd.concat([ak5, dk5, gi5, sw5, um, ye, zamb])
final.to_csv(name , index = False, sep = "\t")
