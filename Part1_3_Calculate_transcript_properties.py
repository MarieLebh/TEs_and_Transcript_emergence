# ! / usr / bin / python3

"""
Calculate transcript characteristics:

Input:

- Fasta file of the transcriptome assembly
- GTF file (with either de novo or annotated transcripts)
- Bedtools intersect output gtf file (with info on gene overlap)

"""
import pandas as pd
import os
import csv

#This function adds the population name to a file

def add_population(GTF, population: "gtf_file")-> "pandas_file":
        with open(GTF) as file_in:
                x = [] #transcipt id
                p = [] #population
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                x.append(l[11])
                                p.append(population)
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        p_dict = {newx[i]: p[i] for i in range(len(x))} #make a dictionary out of the two lists
        p_df = pd.DataFrame(p_dict.items(), columns=['Transcript_ID', 'population'])
        return p_df

#Add the type (annotated or de novo)

def add_type(GTF, type: "gtf_file")-> "pandas_file":
        with open(GTF) as file_in:
                x = [] #transcipt id
                t = [] #population
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                x.append(l[11])
                                t.append(type)
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        t_dict = {newx[i]: t[i] for i in range(len(x))} #make a dictionary out of the two lists
        t_df = pd.DataFrame(t_dict.items(), columns=['Transcript_ID', 'type'])
        return t_df

#Add the chromosome

def add_chrom(GTF: "gtf_file")-> "pandas_file":
        with open(GTF) as file_in:
                x = [] #transcipt id
                c = [] #chromosome
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                x.append(l[11])
                                c.append(l[0])
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        c_dict = {newx[i]: c[i] for i in range(len(x))} #make a dictionary out of the two lists
        c_df = pd.DataFrame(c_dict.items(), columns=['Transcript_ID', 'chromosome'])
        return c_df

#This function opens the gtf file and returns a pandas_df with the gene overlap

def get_overlap(GTF: "gtf_file")-> "pandas_file":
        with open(GTF) as file_in:
                o = [] #overlap yes or no
                x = [] #transcipt id
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                x.append(l[11])
                                if l[27] != "0":
                                        o.append("gene_overlap")
                                elif l[27] == "0":
                                        o.append("intergenic")
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        overlap = {newx[i]: o[i] for i in range(len(x))} #make a dictionary out of the two lists
        o_df = pd.DataFrame(overlap.items(), columns=['Transcript_ID', 'overlap'])
        return o_df

#This function opens the gtf file and returns a csv file with the transcript id its corresponding gene_id

def get_gene_id(GTF: "gtf_file")-> "pandas_df":
        with open(GTF) as file_in:
                g = [] #gene id
                x = [] #transcipt id
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                g.append(l[9])
                                x.append(l[11])
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        newg = [i[1:-2] for i in g] #same just for gene_id
        gene_id = {newx[i]: newg[i] for i in range(len(x))} #make a dictionary out of the two lists
        gene_df = pd.DataFrame(gene_id.items(), columns=['Transcript_ID', 'Gene_ID'])
        return gene_df


#This function takes a gtf file as the input and returns a pandas df with transcript id and TPM value

def get_tpm(GTF: "gtf_file")-> "pandas df":
        with open(GTF) as file_in:
            tpm = [] #tpm of the transcript
            t = [] #transcipt id
            for line in file_in:
                    l = line.split()
                    if l[0]!="#" and l[2] == "transcript":
                            tpm.append(l[17])
                            t.append(l[11])
        newt = [i[1:-2] for i in t] #remove the first and the last two letters to get the same transcript id as in the assembly
        newtpm = [i[1:-2] for i in tpm] #remove the first and the last two letters
        TPM = {newt[i]: float(newtpm[i]) for i in range(len(t))} #make a dictionary out of the two lists
        TPM_df = pd.DataFrame(TPM.items(), columns=['Transcript_ID', 'TPM'])
        return TPM_df

#This function opens the gtf file and returns a csv file with the transcript id and the size of the unspliced transcript

def get_unspliced_transcript(GTF:"gtf_file")-> "csv_file":
        with open(GTF) as file_in:
                start = [] #start of the transcript
                end = [] #end of the transcript
                size = [] #size of the transcript  (end - start)
                x = [] #transcipt id
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                start.append(l[3])
                                end.append(l[4])
                                size.append(int(l[4])-int(l[3]))
                                x.append(l[11])
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        unspliced_length = {newx[i]: size[i] for i in range(len(x))} #make a dictionary out of the two lists
        us_df = pd.DataFrame(unspliced_length.items(), columns=['Transcript_ID', 'unspliced_length'])
        return us_df

#This function opens a gtf file and returns a csv file with the transcript ID and the number of exons of that transcript

def get_exon_number(GTF: "gtf_file")-> "csv_file":
        with open(GTF) as file_in:
                x = [] #transcipt id
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "exon":
                                x.append(l[11])
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly file
        exon_number = {}
        for item in newx:
                if (item in exon_number):
                        exon_number[item] += 1
                else:
                        exon_number[item] = 1
        exon_df = pd.DataFrame(exon_number.items(), columns=['Transcript_ID', 'exon_number'])
        return exon_df

#This function opens a gtf file and returns a csv file with the transcript ID and the number of introns of that transcript

def get_intron_number(GTF: "gtf_file")-> "csv_file":
        with open(GTF) as file_in:
                x = [] #transcipt id
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "exon":
                                x.append(l[11])
        newx = [i[1:-2] for i in x]
        exon_number = {}
        for item in newx: #creates a dictionay out of the list; whenerver an id appears more than once its value is summed up with the value of the same id
                if (item in exon_number):
                        exon_number[item] += 1
                else:
                        exon_number[item] = 1
        intron_number = {k: v - 1 for k, v in exon_number.items()} #intron number = exon number -1
        intron_df = pd.DataFrame(intron_number.items(), columns=['Transcript_ID', 'intron_number'])
        return intron_df

#This function reads in a fasta file and stores it in a dictionary. The key is transcript id and the value contain the fasta sequence

def create_seq_dict(Fasta: "fasta_file")-> "dictionary":
        f=open(Fasta)
        seqs={}
        for line in f:
                line=line.rstrip()
                if line[0]=='>':
                        words=line.split()
                        name=words[0][1:] #get transcript id
                        seqs[name]=''
                else:
                        seqs[name]=seqs[name]+line #add sequence
        return seqs


#This function can calculate the GC count of a DNA sequence

def gc(seq: "fasta_sequence")-> "gc content as float":
        seq = seq.upper() #in case there are lower case transcripts
        A = seq.count("A")
        T = seq.count("T")
        G = seq.count("G")
        C = seq.count("C")
        return (G+C)/(G+C+T+A)

#This function applies the GC count function to all values in the dictionary and returns a csv file whith the transcript id  and the GC content of that sequence

def get_gc(seqs: "dictionary containing fasta sequences and their IDs")-> "csv file":
        GC = {}
        for k, v in seqs.items():
                GC[k] = gc(v) #apply the gc count functions to all values in the sequence dictionary and  store them in a new dictionary
        gc_df = pd.DataFrame(GC.items(), columns=['Transcript_ID', 'gc_content'])
        return gc_df

#This function takes the sequence dictionary and returns the spliced length of the sequences as a df.

def get_spliced_transcript(seqs: "dictionary containing fasta sequences and their IDs")-> "csv file":
        sl = {}
        for k, v in seqs.items():
                sl[k] = len(v) -1  #apply the gc count functions to all values in the sequence dictionary and  store them in a new dictionary
        s_df = pd.DataFrame(sl.items(), columns=['Transcript_ID', 'spliced_length'])
        return s_df

#This function converts each of the created csv files into a pandas dataframe.
#It then uses pandas to merge all files together based on their transcript id (or in one case gene id).
#The output is a merged csv file containing all the data.

def merge_all_csv(GTF, GTF_Overlap, Fasta, population, type, outfile)-> "csv file":
        gene_df = get_gene_id(GTF)
        TPM_df = get_tpm(GTF)
        us_df = get_unspliced_transcript(GTF)
        exon_df = get_exon_number(GTF)
        intron_df = get_intron_number(GTF)
        seqs = create_seq_dict(Fasta)
        gc_df = get_gc(seqs)
        s_df = get_spliced_transcript(seqs)
        o_df = get_overlap(GTF_Overlap)
        p_df = add_population(GTF, population)
        t_df = add_type(GTF, type)
        c_df = add_chrom(GTF)
        data_all2 = gene_df.merge(TPM_df, on = 'Transcript_ID').merge(us_df, on = 'Transcript_ID').merge(exon_df, on = 'Transcript_ID', how = "left").merge(intron_df, on = 'Transcript_ID', how = "left") #merge with gene id file
        data_all = pd.merge(data_all2, gc_df, on="Transcript_ID", how="left").merge(s_df, on="Transcript_ID", how="left").merge(o_df, on = "Transcript_ID", how = "left").merge(t_df, on = "Transcript_ID", how = "left").merge(p_df, on = "Transcript_ID", how = "left").merge(c_df, on = "Transcript_ID", how = "left")
        data_all.fillna("NA", inplace = True)
        final_file = data_all.to_csv(outfile) 
        return data_all2



"""
Run like this
"""

#AK5
GTF = "AK5_transcripts.gtf"  #gtf file
GTF_Overlap = "AK5_transcripts_bed.gtf" #overlap gtf
Fasta = "AK5transcriptomeAssembly.fa"  #Assembly_fasta file
outfile = "AK5Data_denovo.csv" #output file here
population = "AK5"
type = "de_novo"

all = merge_all_csv(GTF, GTF_Overlap, Fasta, population, type, outfile)
print(len(all))


