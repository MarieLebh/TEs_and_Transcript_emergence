#Get the relative TE overlap 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library( data.table )
library( intervals )
library(ggpattern)

#Merge all population files (output of bedtools merge) to one big file
new_colnames = c("chromosome", "transcript_start", "transcript_end","chrom", "score", "strand_transcript", "chromosome", "overlap_start", "overlap_end","name", "overlap")

AK5 <- read.csv("AK5_final_file.bed", header = FALSE, sep="\t")
colnames(AK5) <- new_colnames
AK5$Population = "AK5"

DK5 <- read.csv("DK5_final_file.bed", header = FALSE, sep="\t")
colnames(DK5) <- new_colnames
DK5$Population = "DK5"

GI5 <- read.csv("GI5_final_file.bed", header = FALSE, sep="\t")
colnames(GI5) <- new_colnames
GI5$Population = "GI5"

SW5 <- read.csv("SW5_final_file.bed", header = FALSE, sep="\t")
colnames(SW5) <- new_colnames
SW5$Population = "SW5"

UM <- read.csv("UM_final_file.bed", header = FALSE, sep="\t")
colnames(UM) <- new_colnames
UM$Population = "UM"

YE <- read.csv("YE_final_file.bed", header = FALSE, sep="\t")
colnames(YE) <- new_colnames
YE$Population = "YE"

Zamb <- read.csv("Zamb_final_file.bed", header = FALSE, sep="\t")
colnames(Zamb) <- new_colnames
Zamb$Population = "Zamb"

All_merged <- rbind(AK5, DK5, GI5, SW5, UM, YE, Zamb)

#Restructure the data and split the columns
All_merged$index <- 1:nrow(All_merged)
All_merged = All_merged[,c(ncol(All_merged),1:(ncol(All_merged)-1))]
All_merged$seq_id = All_merged$chrom
All_merged = All_merged %>% 
  separate("chrom", c("Id", "Type", "population", "x"), sep = "_")

All_merged[which(All_merged$Id == "Noncoding"), "Type"] <- "Intergenic"

All_merged[which(All_merged$Type == "End"), "Type"] <- "Downstream"

df = All_merged
rm(AK5, DK5, GI5, SW5, UM, YE, Zamb, All_merged) #rm unecessary

#Exclude seqs where overlap is not on the same strand
df[which(df$strand_transcript == "."), "strand_transcript"] <- "+"
df$same_strand = paste(df$strand_transcript ,df$name,sep="_")
df[which(df$same_strand == "-_+"), "overlap"] <- 0
df[which(df$same_strand == "+_-"), "overlap"] <- 0

#Get the overlap
df = df %>% group_by(df$seq_id) %>% mutate(total_overlap = sum(overlap)) 
df$rel_overlap = df$total_overlap/(df$transcript_end-df$transcript_start) #Get the relativeoverlap

#Make these categories for final plot
df$seq_id2 = paste(df$seq_id, df$strand_transcript)
merged = distinct(df, seq_id2, .keep_all = TRUE) #Keep distinct ids

#Rename the columns
merged[which(merged$Type == "Noncoding"), "group"] <- "Noncoding"
merged[which(merged$Type == "Transcript"), "group"] <- "De novo"
merged[which(merged$Type == "Upstream"), "group"] <- "De novo"
merged[which(merged$Type == "End"), "group"] <- "De novo"
merged$ID = paste(merged$Id, merged$population, sep = "::")

merged[which(0 < merged$rel_overlap & merged$rel_overlap  <= 0.2), "overlap_te"] <- "20 % and less overlap"
merged[which(0.2 < merged$rel_overlap & merged$rel_overlap  <= 0.4), "overlap_te"] <- "20% - 40% overlap"
merged[which(0.4 < merged$rel_overlap & merged$rel_overlap  <= 0.6), "overlap_te"] <- "40% - 60% overlap"
merged[which(0.6 < merged$rel_overlap & merged$rel_overlap  <= 0.8), "overlap_te"] <- "60% - 80% overlap"
merged[which(0.8 < merged$rel_overlap & merged$rel_overlap  < 1), "overlap_te"] <- "80% and above overlap"
merged[which(merged$rel_overlap == 0), "overlap_te"] <- "0 % overlap"
merged[which(merged$rel_overlap == 1), "overlap_te"] <- "Inside TE"

to_save = select(merged, Id, Type, Population, rel_overlap, overlap_te)

to_save$`df$seq_id` = paste(to_save$Type, to_save$Id, to_save$Population, sep = "_")
to_save$Unique_ID = to_save$`df$seq_id`

write.csv(to_save, "Data_rel_TE_overlap.csv")



