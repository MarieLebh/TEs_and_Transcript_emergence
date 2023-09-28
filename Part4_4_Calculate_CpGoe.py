#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 19:35:41 2021

@author: root
"""
from Bio import SeqIO
import os


def getCpGoe(seqsfile, spp, outfile1):#, outfile2):
    print(seqsfile)
    for record in SeqIO.parse(seqsfile, "fasta"):
        Len=len(record.seq)
        Ccount=record.seq.upper().count("C")
        Gcount=record.seq.upper().count("G")
        CGcount=record.seq.upper().count("CG")

        Cfreq=Ccount/Len
        Gfreq=Gcount/Len
        CGfreq=CGcount/(Len-1) #-1  becuase number of pairs=len-1

        if ((Cfreq*Gfreq) == 0):
            print(record.id, CGcount, record.seq, "ZEROOO")# can't devide by zero
            CGoe="NA"
        else:
            CGoe=round(CGfreq/(Cfreq*Gfreq),  7)

        print(spp, record.id,Len,round(Gfreq, 7),round(Cfreq, 7),round(CGfreq, 7), CGoe )
        name = str(record.id).replace('.model.','.TU.')
        outfile1.write(spp+"\t"+name+"\t"+str(Len)+"\t"+str(round(Gfreq, 7))+"\t"+str(round(Cfreq, 7))+"\t"+str(round(CGfreq, 7))+"\t"+str(CGoe)+"\n")
#        outfile2.write(spp+"\t"+record.id+"\t"+str(Len)+"\t"+str(round(Gfreq, 7))+"\t"+str(round(Cfreq, 7))+"\t"+str(round(CGfreq, 7))+"\t"+str(CGoe)+"\n")

#Run like this
path = '/media/bebs/DATA/dpun_clean/Dpun_genome/'
out = open('%sCpGoe_dpun.tsv' % path,'w+')
out.write("Spp"+"\t"+"ID"+"\t"+"Len"+"\t"+"Gfreq"+"\t"+"Cfreq"+"\t"+"CGfreq"+"\t"+"CpGoe"+"\n")
getCpGoe('%sEVM_OGSprelim3_cds.fa' % path,'Dpun',out)
out.close()
