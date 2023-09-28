#!/usr/bin/env python

"""
Build orthogroups based on the shared transcript definition
"""

"""
MAKE LIST OF PAIRS

-Input:
	A csv file containing all filtered matches from the Blast and the no hit transcripts.

- Output:
	A list where each query-subject hit pair is saved as a sublist.
	Transcripts with no hits have an empty spot in their sublist.
"""

def make_list_of_pairs(File):
      F=open(File, "r")
      L=F.readlines()
      final_list = []
      names = []
      names2 = []
      for line in  L:
            line = line.split("\t")
            names.append(line[0])
            names2.append(line[2])
            final_list = [list(pair) for pair in zip(names, names2)]
      return final_list


"""
BUILD ORTHOGROUPS

- Input:
	List with all blast query-subject pairs (filtered) and no hit - empty pairs 

- Output:
	List where all transcripts are grouped into orthogroups 
"""

def build_orthogroups(list_matches):
    orthogroups = []
    for pair in list_matches:
        pair = set(pair)
        equals_found = 0
        for idx, group in enumerate(orthogroups):
            if group.intersection(pair) and '' not in pair:
	    #make sure here that empty values do not count (else all non matches are grouped together)
                equals_found += 1
                if equals_found == 1: #first transcript in pair matches transcript in group
                    group.update(pair)
                    first_group = group
                elif equals_found == 2: #second transcript in pair matches transcript in group
                    first_group.update(group)
                    del orthogroups[idx]
                    break
        if not equals_found: #If none of the values is found, a new orthogroup is created
            orthogroups.append(pair)
    x = [list(sorted(group)) for group in orthogroups]
    return x #a list containing all transcripts sorted into orthogroups; each orthogroup is saved in a separate set


"""
BUILD FINAL FILE

- Input:
	List of orthogroups

- Output:
	Text file with the orthogroups; number of seqs per orthogroup and number of populations per orthogroup
"""
def build_final_file(x):
            F = open("Orthogroups.txt", "w")
            F.write("Number_Orthogroup; IDs_in_orthogroup; Number_seqs_in_orthogroup; Number_populations_in_orthogroup")
            F.write("\n")
            count = 0
            count_pops = 0
            for i in x:
                  count = count + 1
                  if '' in i: #Check if an empty value is in the orthogroup (needs to be substracted from sequence count)
                        if any("AK5" in s for s in i): #Count if each population occurs
                              count_pops = count_pops + 1
                        if any("DK5" in s for s in i):
                              count_pops = count_pops + 1
                        if any("GI5" in s for s in i):
                              count_pops = count_pops + 1
                        if any("SW5" in s for s in i):
                              count_pops = count_pops + 1
                        if any("UM" in s for s in i):
                              count_pops = count_pops + 1
                        if any("YE" in s for s in i):
                              count_pops = count_pops + 1
                        if any("Zamb" in s for s in i):
                              count_pops = count_pops + 1
                        F.write("Orthogroup nr. " + str(count) + ";" + ",".join(i) + ";" + str(len(i)-1)+ ";" + str(count_pops))
                        F.write("\n")
                        count_pops = 0
                  else:
                        if any("AK5" in s for s in i): #Count if each population occurs
                              count_pops = count_pops + 1
                        if any("DK5" in s for s in i):
                              count_pops = count_pops + 1
                        if any("GI5" in s for s in i):
                              count_pops = count_pops + 1
                        if any("SW5" in s for s in i):
                              count_pops = count_pops + 1
                        if any("UM" in s for s in i):
                              count_pops = count_pops + 1
                        if any("YE" in s for s in i):
                              count_pops = count_pops + 1
                        if any("Zamb" in s for s in i):
                              count_pops = count_pops + 1
                        F.write("Orthogroup nr. " + str(count) + ";" + ",".join(i) + ";" + str(len(i))+ ";" + str(count_pops))
                        F.write("\n")
                        count_pops = 0
            F.close()



#Run this in python
x = make_list_of_pairs("transcripts_with_blast.csv")
y = build_orthogroups(x)
final =  build_final_file(y)
