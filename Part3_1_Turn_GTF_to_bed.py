#!/usr/bin/env python

import pandas as pd
import csv

"""
Turn a gene/transcript gtf file to a bedfile.

Run this script before running bedtools closest.

These two functions take a gtf/gff file as an input file and return a bedfile (tab separated file).

- The function make_transcript_bedfile() takes a gtf file with transcripts.
- The function make_gene_bedfile() takes a gtf file with genes.

The following columns of the input gtf will be taken:

name: name of the chromosome (e.g. 2R)
start:start position of the sequence
end: end position of the sequence
id: sequence id
score: score from 0 to 1000
strand: which strand the sequence is on

The output file can be used to run bedtools closest to get the closest up/downstream gene.
With this data the corresponding up/downstream region can be motified to exclude the gene overlapping parts.

"""

def make_transcript_bedfile(GTF, outfile):
    with open(GTF) as file_in:
        seqname = [] #Chrom
        start = [] #start of transcript
        end = [] #end of transcript
        score = [] #score
        strand = [] #strand
        x = [] #transcript id not modified
        trans_id = [] #final transcript id
        for line in file_in:
                l = line.split()
                if l[0]!="#" and l[2] == "transcript": #append all values, only for transcripts NOT exons
                        seqname.append(l[0])
                        start.append(l[3])
                        end.append(l[4])
                        score.append(l[5])
                        strand.append(l[6])
                        x.append(l[11])
                        trans_id = [i[1:-2] for i in x] #final transcript id
    df = pd.DataFrame()
    df["seqname"] = seqname
    df["start"] = start
    df["end"] = end
    df["transcript_id"] = trans_id
    df["score"] = score
    df["strand"] = strand
    df.rename_axis(None, axis = 1)
    df.to_csv(outfile, header = False, index = False, sep = "\t") #create the outfile
    return outfile


def make_gene_bedfile(GTF, outfile):
    with open(GTF) as file_in:
        seqname = [] #Chrom
        start = [] #start of gene
        end = [] #end of gene
        score = [] #score
        strand = [] #strand
        x = [] #originally meant for id but for simplicity it will just print gene
        for line in file_in:
            if line[0] != "#":
                l = line.split()
                if l[2] == "gene": #append all values, only for transcripts NOT exons
                            seqname.append(l[0])
                            start.append(l[3])
                            end.append(l[4])
                            score.append(l[5])
                            strand.append(l[6])
                            x.append("gene")
        df = pd.DataFrame()
        df["seqname"] = seqname
        df["start"] = start
        df["end"] = end
        df["type"] = x
        df["score"] = score
        df["strand"] = strand
        df.rename_axis(None, axis = 1)
        df.to_csv(outfile, header = False, index = False, sep = "\t")
        return outfile

#AK5

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/AK5_intergenic_transcripts.gtf"
outfile = "AK5_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/AK5/Step3GeMoMa/final_annotation.gff"
outfile =  "AK5_genes.bed"
x = make_gene_bedfile(GTF, outfile)

#DK5

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/DK5_intergenic_transcripts.gtf"
outfile = "DK5_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/DK5/Step3GeMoMa/final_annotation.gff"
outfile =  "DK5_genes.bed"
x = make_gene_bedfile(GTF, outfile)

#GI5

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/GI5_intergenic_transcripts.gtf"
outfile = "GI5_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/GI5/Step3GeMoMa/final_annotation.gff"
outfile =  "GI5_genes.bed"
x = make_gene_bedfile(GTF, outfile)

#SW5

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/SW5_intergenic_transcripts.gtf"
outfile = "SW5_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/SW5/Step3GeMoMa/final_annotation.gff"
outfile =  "SW5_genes.bed"
x = make_gene_bedfile(GTF, outfile)

#UM

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/UM_intergenic_transcripts.gtf"
outfile = "UM_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/UM/Step3GeMoMa/final_annotation.gff"
outfile =  "UM_genes.bed"
x = make_gene_bedfile(GTF, outfile)

#YE

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/YE_intergenic_transcripts.gtf"
outfile = "YE_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/YE/Step3GeMoMa/final_annotation.gff"
outfile =  "YE_genes.bed"
x = make_gene_bedfile(GTF, outfile)

#Zamb

GTF = "/global/students/research/m_lebh01/Bedtools/De_novo_transcript_final_bedfiles/Zamb_intergenic_transcripts.gtf"
outfile = "Zamb_intergenic_transcripts.bed"
x = make_transcript_bedfile(GTF, outfile)
GTF = "/global/students/research/m_lebh01/RealWork/Zamb/Step3GeMoMa/final_annotation.gff"
outfile =  "Zamb_genes.bed"
x = make_gene_bedfile(GTF, outfile)

