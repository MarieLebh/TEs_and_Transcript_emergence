# ! / usr / bin / python3

"""
This script extracts the ids of all hits from a blast (tabular) output file.
It then extracts all transcript ids from a gtf file.
As a last step it returns all transcript ids that were not in the list of ids that
had one or multiple hits during the blast search.
"""

#Step1: Extract the query id stored in the first column of the output file and store them in a list

def list_matches(blast: "blast tabular output file")-> "list":
        with open(blast) as file_in:
                x = [] #matches
                for line in file_in:
                        l = line.split()
                        x.append(l[0])
        return x

#Step 2: Get a list of all transcript ids from the gtf file of the transcriptome assembly

def get_gene_id_csv(GTF: "gtf_file")-> "list":
        with open(GTF) as file_in:
                x = [] #transcipt id
                for line in file_in:
                        l = line.split()
                        if l[0]!="#" and l[2] == "transcript":
                                x.append(l[11])
        newx = [i[1:-2] for i in x] #remove the first and the last two letters to get the same transcript id as in the assembly
        return newx


#Step 3: Compare the two lists with each other and return a new list with only the ones that are not in the list with all hits

def return_not_matches(a, b):
    a = set(a) #list of all transcripts
    b = set(b) # list of all blast hits
    c = list(a - b) #list of all transcript without blast hits
    return c

#Step 4: Create a text file where each transcript id without a hit is stored as a new line

def file_no_hit(l, Filename)-> "txt_file":
        with open(Filename, 'w') as outfile:
                for item in l:
                        outfile.write(item + "\n")
        return outfile


#Run like this

Blast =  "blast_all_AK5.txt" #add path to the blast output file
GTF ="AK5transcriptome.gtf" #add path to the transcriptome assembly file

#Get a list with all transcript ids from the transcriptome assembly gtf file

all_ids = get_gene_id_csv(GTF)

#Get a list with the transcript ids that had blast hits

ids = list_matches(Blast)

#Compare each of the lists with the blast matches with the list of all transcripts
#Return the ids that had no match in the blast search

no_match = return_not_matches(all_ids, ids)

#Create a textfile from the list. It contains each transcript id without a match as a new line

name = "AK5_transcripts_no_blast_hit_all.txt"

file_plus = file_no_hit(no_match, name)
