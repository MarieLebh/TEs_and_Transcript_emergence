#!/usr/bin/env python

"""
Merge the individual output files from bedtools intersect for each population
"""

import pandas as pd
import csv
import random

def make_bedfile(file1, population):
      colnames = ["chrom_seq", "start_seq", "end_seq", "id_seq", "score" , "strand_seq"  , "chrom_te", "start_te", "end_te", "id_te", "score_te", "strand_te" , "overlap"]
      df = pd.read_csv(file1, sep='\t', names = colnames, header = None)
      df["population"] = population
      #df["id"] = df["id_seq"] + "_" + type + "_" + population
      #df.to_csv(outfile, sep='\t', index = False)
      return df

#Run like this 


ak5 = make_bedfile("AK5_intergenicTE.bed", "AK5")
dk5 = make_bedfile("DK5_intergenicTE.bed", "DK5")
gi5 = make_bedfile("GI5_intergenicTE.bed", "GI5")
sw5 = make_bedfile("SW5_intergenicTE.bed", "SW5")
um = make_bedfile("UM_intergenicTE.bed", "UM")
ye = make_bedfile("YE_intergenicTE.bed", "YE")
zamb = make_bedfile("Zamb_intergenicTE.bed", "Zamb")

all = pd.concat([ak5, dk5, gi5, sw5, um, ye, zamb])


all.to_csv("All_noncoding.csv", sep='\t', index = False) 