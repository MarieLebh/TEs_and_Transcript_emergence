#!/usr/bin/env python

"""
A)
Create a bedfile for up- and downstream regions.

Input: Bedfile of Transcripts (or genes or homologs)

B)
Create a bedfile for the randomly sampled intergenic regions from bedtools sample.

Input: Bedoutfile from bedtools sample (contains intergenic regions), population

Run bedtools like this:

-w = 1100 (same size as gene/transcript upstream)
- n = 7910 (10 x the size of de novo transcripts AK5 (with TE transcripts))

bedtools makewindows -b All_intergenic_regions.bed -w 1100 \
     | awk '($3-$2) == 1100' \
     | bedtools sample -i - -n 7910  > AK5_random.bed


"""

import pandas as pd
import csv
import random

#Get the upstream regions

def get_promotor_bed(bedfile, outfile):
    with open(bedfile) as file_in:
        seqname = [] #Chrom
        start = [] #start of transcript
        end = [] #end of transcript
        start2 = [] #start of closest gene
        end2 = [] #end of closest gene
        score = [] #score
        strand = [] #strand
        trans_id = [] #final transcript id or gene id
        new_start = [] #promotor region start
        new_end = []  #promotor region end
        for line in file_in:
                l = line.split()
                seqname.append(l[0])
                start.append(l[1])
                end.append(l[2])
                if l[5] == "-": #For negative strand
                    new_start.append(int(l[2])- 100) 
                    new_end.append(int(l[2])+ 1000)
                else: #For positive strand
                     new_start.append(int(l[1])- 1000)
                     new_end.append(int(l[1]) + 100)
                score.append(l[4])
                strand.append(l[5])
                trans_id.append(l[3])
    new_list = [] #for all transcripts at the start of the chrom: if new_start >0 just take 0 as new start
    for num in new_start:
        if int(num) < 0:
             new_list.append(0)
        else:
            new_list.append(num)
    df = pd.DataFrame()
    df["seqname"] = seqname
    df["start"] = new_list
    df["end"] = new_end
    df["transcript_id"] = trans_id #transcript or gene id
    df["score"] = score
    df["strand"] = strand
    df.rename_axis(None, axis = 1)
    df.to_csv(outfile, header = False, index = False, sep = "\t") #create the outfile
    return outfile

#Get the downstream regions

def get_end_bed(bedfile, outfile):
    with open(bedfile) as file_in:
        seqname = [] #Chrom
        start = [] #start of transcript
        end = [] #end of transcript
        start2 = [] #start of closest gene
        end2 = [] #end of closest gene
        score = [] #score
        strand = [] #strand
        trans_id = [] #final transcript id 
        new_start = [] #start
        new_end = []  #end
        for line in file_in:
                l = line.split()
                seqname.append(l[0])
                start.append(l[1])
                end.append(l[2])
                if l[5] == "-":
                    new_start.append(int(l[1])- 1000) 
                    new_end.append(int(l[1]) + 100)
                else:
                     new_start.append(int(l[2]) - 100)
                     new_end.append(int(l[2]) + 1000)
                score.append(l[4])
                strand.append(l[5])
                trans_id.append(l[3])
    new_list = [] #for all transcripts at the start of the chrom: if new_start >0 just take 0 as new start
    for num in new_start:
        if int(num) < 0:
             new_list.append(0)
        else:
            new_list.append(num)
    df = pd.DataFrame()
    df["seqname"] = seqname
    df["start"] = new_list
    df["end"] = new_end
    df["transcript_id"] = trans_id #transcript or gene id
    df["score"] = score
    df["strand"] = strand
    df.rename_axis(None, axis = 1)
    df.to_csv(outfile, header = False, index = False, sep = "\t") #create the outfile
    return outfile

#Make a bedfile for the intergenics and randomly assign a strand

def make_intergenic_bedfile(bedfile, pop, outfile): 
    with open(bedfile) as file_in:
        seqname = []
        start = []
        end = []
        score = []
        name = []
        n = 0
        strand = ("+", "-")
        strand2 = []
        population = []
        for line in file_in:
            l = line.split()
            seqname.append(l[0])
            start.append(l[1])
            score.append(".")
            end.append(l[2])
            name.append("Noncoding_" + str(n+1))
            n = n +1
            strand2.append(random.choice(strand))
    df = pd.DataFrame()
    df["seqname"] = seqname
    df["start"] = start
    df["end"] = end
    df["transcript_id"] = name
    df["score"] = score
    df["strand"] = strand2
    df.rename_axis(None, axis = 1)
    df.to_csv(outfile, header = False, index = False, sep = "\t") #create the outfile
    return outfile

################################################
#Run like this

#AK5
outfile = "AK5_genes_upstream.bed"
bedfile = "AK5_overlap_genes.bed"
x = get_promotor_bed(bedfile, outfile)

outfile = "AK5_intergenic_transcripts_upstream.bed"
bedfile = "AK5_overlap_intergenic.bed"
x = get_promotor_bed(bedfile, outfile)

outfile = "AK5_intergenic_transcripts_end.bed"
bedfile = "AK5_overlap_intergenic2.bed"
x = get_end_bed(bedfile, outfile)

outfile = "AK5_intergenic_random.bed"
bedfile = "AK5_random.bed"
pop = "AK5"
make_intergenic_bedfile(bedfile, pop, outfile)
