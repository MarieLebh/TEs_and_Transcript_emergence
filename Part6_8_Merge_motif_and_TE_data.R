#Merge motif and TE data

library(dplyr)
library(readr)
library(tidyverse)
library(data.table )
library(car)
library(stringr)
library(ggforce)
library(rstatix)
library(ggpattern)

#Get the TE data
TEs = read.csv("TE.csv") #Contains all transcripts with >= 0.8 TE overlap (=TE)
id_TE = TEs$x

OverlapTE = read.csv("TE_overlap.csv") #Transcript ids with > 0.8 overlap
id_TE_overlap = OverlapTE$x

NoTE = read.csv("no_TE.csv") #No TE overlap
id_TE_no = NoTE$x

df$type_te = df$type
df$Id = paste(df$id, df$population, sep = "_")

#Load  and prepare the motif data
Genes <- read.delim("All_genes.csv")
Noncoding <- read.delim("All_noncoding.csv")
Transcripts <- read.delim("All_intergenic_de_novo_transcripts.csv")
Final_df <- rbind(Genes, Noncoding, Transcripts)
Final_df$type2 = Final_df$type

#Separate the columns
Final_df = Final_df %>% 
  separate("seq_name", c("type", "chrom", "position"), sep = "::")
Final_df = Final_df %>% 
  separate("chrom", c("chrom", "position"), sep = ":")
Final_df$final = str_sub(Final_df$position, end = -4)
Final_df = Final_df %>% 
  separate("final", c("start", "end"), sep = "-")
Final_df$size = as.numeric(Final_df$end) - as.numeric(Final_df$start)
Final_df[which(Final_df$type2 == "noncoding"), "size"] <- 1100 #As all noncodings have this size
Final_df$type = Final_df$type2
Final_df$id = Final_df$unique_id
Final_df = Final_df %>% 
  separate("id", c("seq_id", "pos"), sep = "::")
Final_df$id = paste(Final_df$seq_id, Final_df$population, sep = "_")

#Only keep the ones with 1100 size
Final_df = subset(Final_df, size == "1100")

#Join with TE df
Final_df$type_te = Final_df$type

Final_df[which(Final_df$id %in% id_TE), "type_te"] <- "dn_TE"
Final_df[which(Final_df$id %in% id_TE_no), "type_te"] <- "dn_no_TE"
Final_df[which(Final_df$id %in% id_TE_overlap), "type_te"] <- "dn_overlap_TE"

Final_df$type_te2 = Final_df$type_te 

#Shorten the name
Final_df[which(Final_df$type2== "intergenic_de_novo_transcript"), "type"] <- "de_novo"
Final_df$population_type = paste(Final_df$population, Final_df$type, sep = "_")

#Save final df
to_save = select(Final_df, unique_id, population, type, type_te, number_motifs_high_per_transcript, number_motifs_per_transcript)
write.csv(to_save, "Data_motifs.csv")


#Load the data for promotor regions to get the number of core
Genes <- read.delim("Promotor_All_genes.csv")
Noncoding <- read.delim("Promotor_All_noncoding.csv")
Transcripts <- read.delim("Promotor_Intergenic_de_novo_transcripts.csv")
Final_df2 <- rbind(Genes, Noncoding, Transcripts)
Final_df2$type2 = Final_df2$type

Final_df2 = Final_df2 %>% 
  separate("seq_name", c("type", "chrom", "position"), sep = "::")
Final_df2 = Final_df2 %>% 
  separate("chrom", c("chrom", "position"), sep = ":")
Final_df2$final = str_sub(Final_df2$position, end = -4)
Final_df2 = Final_df2 %>% 
  separate("final", c("start", "end"), sep = "-")
Final_df2$size = as.numeric(Final_df2$end) - as.numeric(Final_df2$start)
Final_df2[which(Final_df2$type2 == "noncoding"), "size"] <- 1100 #As I know all noncoding ones have this size
Final_df2$type = Final_df2$type2
Final_df2$id2 = Final_df2$unique_id 
Final_df2 = Final_df2 %>% 
  separate("id2", c("seq_id", "pos"), sep = "::")

Final_df2$promotor_size = as.numeric(Final_df2$size)
Final_df2[which(Final_df2$size >= 300), "promotor_size"] <- 300

Final_df2$id = paste(Final_df2$seq_id, Final_df2$population, sep = "_")
Final_df2$type_te = Final_df2$type


Final_df2 = subset(Final_df2, promotor_size >= 300)

#Again merge with TE data
Final_df2[which(Final_df2$id %in% id_TE), "type_te"] <- "dn_TE"
Final_df2[which(Final_df2$id %in% id_TE_no), "type_te"] <- "dn_no_TE"
Final_df2[which(Final_df2$id %in% id_TE_overlap), "type_te"] <- "dn_overlap_TE"
Final_df2$type_te2 = Final_df2$type_te

to_save = select(Final_df2, unique_id, population, type, type_te, number_core_high_per_transcript, number_core_per_transcript)
y = subset(to_save, type == "transcript")

write.csv(to_save, "Data_core.csv")
