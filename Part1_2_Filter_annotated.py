#!/usr/bin/env python

"""
Filter annotated transcripts for TPM and splicing
"""


import pandas as pd
pd.options.mode.chained_assignment = None

# Step 1: Extract the gene id and the transcript id from the gtf file and stores them in a pandas df

def get_gene_id(GTF: "gtf_file")-> "pandas df":
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
        TPM = {newx[i]: float(newtpm[i]) for i in range(len(t))} #make a dictionary out of the two lists
        TPM_df = pd.DataFrame(TPM.items(), columns=['Transcript_ID', 'TPM'])
        Gene_df = pd.DataFrame(gene_id.items(), columns=['Transcript_ID', 'Gene_ID'])
        df = Gene_df.merge(TPM_df, on = 'Transcript_ID')
        return df

#counts the number of times each gene id occurs in the dataset
def return_occurence(df: "pandas df")-> "pandas df with additional column":
        df['count_all'] = df.groupby('Gene_ID')['Gene_ID'].transform('count')
        return df

# Take the text file with the de novo transcript ids and stores them in a list

def get_de_novo_ids(de_novo_ids):
    de_novo_list = open(de_novo_ids).read().splitlines()
    return de_novo_list

# compares the df to a list of de novo transcript ids
# returns a dataset containing only the ids that match the list

def get_de_novo_positions(df, ids):
    df1 = df.query('Transcript_ID not in @ids')
    return df1

# Calculates the occurence per gene id again and stores it in the df

def return_occurence2(df1)-> "pd_df":
    df1['count_de_novo'] = df1.groupby('Gene_ID')['Gene_ID'].transform('count')
    return df1

# calculates the difference between the number of occurence of all transcripts
# and the number of occurence after all non de novo transcripts have been
# removed
# if its not 0 then there are some spliceforms that had a blast hit
# while others did not

def get_difference(df1):
    df1['Score_diff'] = df1['count_all'] - df1['count_de_novo']
    return df1

# remove all transcripts where the column score diff is not zero
# Thus keep all the ones where it is 0

def remove_transcripts(df):
    df= df[df['Score_diff'] == 0]
    return df

def remove_TPM(df):
    df2 = df[df["TPM"] >= 0.5]
    with open(outfile, "w") as f_out:
        dfAsList = df2["Transcript_ID"].to_list()
        for item in dfAsList:
            f_out.write(item +  "\n")
    return f_out

#This function is just there to get the number of transcripts after TPM removal
def remove_TPM2(df):
    df2 = df[df["TPM"] >= 0.5]
    return df2

# This df is then the final one where all alternative spliceforms (if any exist)
# show no blast hit and are therefore used as de novo transcripts


#Run like this

Population = "AK5"
GTF = "AK5transcriptome.gtf" #original transcriptome assembly file with all transcripts
de_novo_ids = "AK5_transcripts_no_blast_hit_all.txt" #a list of all transcript ids that had no blast hits
outfile = "AK5_not_de_novo.txt" #name of the outfile containing the final annotated ids to use

x = get_gene_id(GTF)
print("Number of transcripts original gtf: ", len(x))
df = return_occurence(x)
ids = get_de_novo_ids(de_novo_ids)
df1 = get_de_novo_positions(df, ids)
df1 = return_occurence2(df1)
df = get_difference(df1)
print("Number of annotated transcripts: ", len(df))
df = remove_transcripts(df)
print("Number of annotated transcripts where all alternative splice forms have a hit: " ,len(df))
d = remove_TPM(df)
s = remove_TPM2(df)
print("Number of de novo transcripts where all alternative spliceforms have no hit and a TPM >= 0.5: ", len(s))
