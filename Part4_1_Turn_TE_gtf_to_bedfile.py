#!/usr/bin/env python

"""
Turn the TE annotation gtf into a bedfile

Input:

- text file that contains info which chromosome corresponds to which seqname in the gtf
- TE Annotation gtf file (TransposonUltimate pipeline output)

"""

import pandas as pd
import csv
import random

#Read in the seq head file and safe the name as id and the corresponding chrom as value

def make_chrom_dictionary(file_chrom):
      colnames = ['id', 'chrom']
      df = pd.read_csv(file_chrom, sep='\t', names = colnames)
      df['id'] = df['id'].str.replace('>', '')
      df['chrom'] = df['chrom'].str.replace('>', '')
      chrom_dict = dict(zip(df['id'], df['chrom']))
      return chrom_dict

#Make a bedfile (replace the id (seq1 etc) with the corresponding chromosome

def make_bedfile(file1, dic, outfile):
      colnames = ['chrom', 'program', 'type', 'start', 'end', 'sth', 'strand', 'score', 'id']
      df = pd.read_csv(file1, sep='\t', names = colnames)
      df['chrom'] = df['chrom'].replace(dic)
      df = df[["chrom", "start", "end", "id", "score", "strand"]]
      df.to_csv(outfile, index = False, header = False, sep = "\t")



#Run like this

seq_file = "sequence_heads.txt"
te_annotation = "FinalAnnotations_Transposons.gff3"
outfile = "AK5_te_annotation.bed"

chrom_dict = make_chrom_dictionary(seq_file)
make_bedfile(te_annotation, chrom_dict, outfile)


