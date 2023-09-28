#Prepare files for bedtools merge (works similarly for the homolog files)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library( data.table )
library( intervals )

#Load the input files (= bedtools intersect output files)

noncoding <- read_delim("all_noncoding_te_overlap.csv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

upstream <- read_delim("all_upstream_te_overlap.csv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

end <- read_delim("all_end_regions.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

transcript <- read_delim("all_transcripts.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

end$type = "end"
noncoding$type = "noncoding"
transcript$type = "transcript"
upstream$type = "upstream"

df4 = rbind(end, noncoding, transcript, upstream)

df4$id_seq = paste(df4$type, df4$id_seq, sep = "_")
df4$id_seq = paste(df4$population, df4$id_seq, sep = "_")


#Calculate size of the region 
df4$start_seq <- as.numeric(df4$start_seq)
df4$end_seq <- as.numeric(df4$end_seq)
df4$size <- df4$end_seq - df4$start_seq

#Prepare the files for bedtools analysis

#First deal with TE within a transcript that overlap each other
#Get the overlap of each TE with the actual transcript

#Get the start of the overlap
for(i in 1:nrow(df4)){
  if(df4$start_seq[i] >= df4$start_te[i]){
    df4$start_overlap[i]= df4$start_seq[i]
  } else if(df4$start_seq[i] < df4$start_te[i]){
    df4$start_overlap[i]= df4$start_te[i]
  } else {df4$start_overlap[i]= df4$start_te[i]}
}

#Get the end of the overlap
for(i in 1:nrow(df4)){
  if(df4$end_seq[i] <= df4$end_te[i]){
    df4$end_overlap[i]= df4$end_seq[i]
  } else if(df4$end_seq[i] > df4$end_te[i]){
    df4$end_overlap[i]= df4$end_te[i]
  } else {df4$end_seq[i]= df4$end_te[i]}
}

df4$start_overlap = as.numeric(df4$start_overlap)
df4$end_overlap = as.numeric(df4$end_overlap)

#short control (check that the overlap from the newly calculated coordinates is the same as the bedtools overlap)
df = df4
df$overlap2 = df$end_overlap - df$start_overlap
df[which(df$overlap == 0), "overlap2"] <- 0
df$difference = df$overlap -df$overlap2
control = subset(df, difference != 0) #This dataframe should be empty if correct
rm(control)

#Only choose the transcripts with overlap
df4[which(df4$overlap== 0), "start_overlap"] <- 0
df4[which(df4$overlap== 0), "end_overlap"] <- 0

#Exclude transcripts not overlapping on the same strand
df4$same_strand = paste(df4$strand_seq ,df4$strand_te,sep="_")
df4[which(df4$same_strand == "-_+"), "overlap"] <- 0
df4[which(df4$same_strand == "+_-"), "overlap"] <- 0
df4[which(df4$same_strand == "._-"), "overlap"] <- 0
df4[which(df4$same_strand == "-_."), "overlap"] <- 0

x = subset(df4, overlap != 0) #make a subset of all data where the overlap is not 0
df4 = x
df4$id2 = paste(df4$id_seq, df4$id_te, sep = "##")
data = data.frame(df4$chrom_seq, df4$start_overlap, df4$end_overlap, df4$id2, ".", df4$strand_te, df4$population)
class(data$df4.seqname)
data$df4.chrom_seq = as.factor(data$df4.chrom_seq)

#Make a bedfile for each population that includes the overlapping positions for all TE
#Use these files to run bedtools merge. This will merge the overlapping intervals together. 

AK5 = subset(data, df4.population == "AK5")
AK5 = subset(AK5, select = -df4.population)
write.table(AK5, sep = "\t", file = "AK5.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

DK5 = subset(data, df4.population == "DK5")
DK5 = subset(DK5, select = -df4.population)
write.table(DK5, sep = "\t", file = "DK5.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

GI5 = subset(data, df4.population == "GI5")
GI5 = subset(GI5, select = -df4.population)
write.table(GI5, sep = "\t", file = "GI5.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

SW5 = subset(data, df4.population == "SW5")
SW5 = subset(SW5, select = -df4.population)
write.table(SW5, sep = "\t", file = "SW5.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

UM = subset(data, df4.population == "UM")
UM = subset(UM, select = -df4.population)
write.table(UM, sep = "\t", file = "UM.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

YE = subset(data, df4.population == "YE")
YE = subset(YE, select = -df4.population)
write.table(YE, sep = "\t", file = "YE.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)

Zamb = subset(data, df4.population == "Zamb")
Zamb = subset(Zamb, select = -df4.population)
write.table(Zamb, sep = "\t", file = "Zamb.bed", col.names = FALSE, row.names = FALSE, quote = FALSE)


