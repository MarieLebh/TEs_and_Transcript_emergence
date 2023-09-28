#########################################
#Make dataframes to merge Homolog + transcripts
#Will be later used for plots & stats
#########################################

library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library( data.table )
library( intervals )

#########################################
#Get dataframe with TE family and number
#########################################

upstream <- read_delim("all_upstream_te_overlap.csv", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

end <- read_delim("all_end_regions.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

transcript <- read_delim("all_transcripts.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

upstream_h <- read_delim("all_homologs_TE_us.csv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)

end_h <- read_delim("all_homologs_TE_ds.csv", 
                  delim = "\t", escape_double = FALSE, 
                  trim_ws = TRUE)

transcript_h <- read_delim("all_homologs_TE.csv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)


end$type = "Downstream(t)"
transcript$type = "Transcript(= t)"
upstream$type = "Upstream(t)"
end_h$type = "Downstream(h)"
transcript_h$type = "Homolog(= h)"
upstream_h$type = "Upstream(h)"

df4 = rbind(end, transcript, upstream, end_h, upstream_h, transcript_h)

df4$id_seq = paste(df4$type, df4$id_seq, sep = "_")
df4$id_seq = paste(df4$population, df4$id_seq, sep = "_")

#Only choose the transcripts with overlap
df4[which(df4$overlap== 0), "start_overlap"] <- 0
df4[which(df4$overlap== 0), "end_overlap"] <- 0

#Make sure its the same strand
df4$same_strand = paste(df4$strand_seq ,df4$strand_te,sep="_")
df4[which(df4$same_strand == "-_+"), "overlap"] <- 0
df4[which(df4$same_strand == "+_-"), "overlap"] <- 0
df4[which(df4$same_strand == "._-"), "overlap"] <- 0
df4[which(df4$same_strand == "-_."), "overlap"] <- 0

df4$number_repeats = 1
df4[which(df4$overlap == 0), "number_repeats"] <- 0
df4 = df4 %>% group_by(id_seq) %>% mutate(number_repeats2 = sum(number_repeats))

#Split df column to get class, subclass and family
df4$TE <- sub("(.*\\()(.*)(\\))", "\\2", df4$id_te)
df4 = df4 %>% 
  separate("TE", c("family", "type2"), sep = ";")
df4 = df4 %>% 
  separate("family", c("family", "subclass", "class"), sep = ",")
df4 = df4 %>% 
  separate("type2", c("useless", "completeness"), sep = "=")

df4$unique = df4$id_seq 

df4 = df4 %>% 
  separate("id_seq", c("Population", "Type", "ID"), sep = "_")

df4$ID = paste(df4$ID, df4$population, sep = "::")

df4[which(df4$overlap == 0), "family"] <- NA

#Group by unique ids, collapse all columns together that are doubled (because multiple TE overlap)
df4 = df4 %>% group_by(unique) %>% mutate(family2 = n_distinct(family, na.rm = TRUE))

#Create multi TE column (if different families overlap)
df4$family_multi = df4$family
df4[which(df4$family2 > 1), "family_multi"] <- "Multiple_TE"

temp = df4 %>%
  group_by(unique) %>%
  summarise_all(funs(list(na.omit(.))))

temp$multi2_2 = sapply(temp$family_multi, `[`,1)

#Select final columns from df
Final = select(temp, unique, family, number_repeats2, multi2_2)

Final[which(Final$family == "character(0)"), "family"] <- NA

#Extract number of repeats
Final$number_repeats2 = as.vector(Final$number_repeats2)

Final <- Final %>%
  mutate(number_repeats2 = sapply(number_repeats2, `[`, 1))

Final = Final %>% 
  separate("unique", c("Population", "Type", "ID", "pop1", "pop2", "pop3"), sep = "_")

#Make the dataframe for merged data of transcripts and homologs
#Transcripts
Transcripts = subset(Final, Type == "Transcript(= t)")
Transcripts$Id = paste(Transcripts$ID, Transcripts$Population, sep = "::")
Transcripts$number_TE_transcript = paste(Transcripts$number_repeats2)
Transcripts$TE_family_transcript = paste(Transcripts$family)
Transcripts$TE_family_multi = Transcripts$multi2_2
Tr = select(Transcripts, Id, Population, family,TE_family_multi, number_repeats2)

Upstream = subset(Final, Type == "Upstream(t)")
Upstream$family_upstream_t = paste(Upstream$family)
Upstream$number_TE_upstream_t = paste(Upstream$number_repeats2)
Upstream$Id = paste(Upstream$ID, Upstream$Population, sep = "::")
Upstream$TE_family_multi_upstream = Upstream$multi2_2
U = select(Upstream, Id, number_TE_upstream_t, family_upstream_t, TE_family_multi_upstream)

Downstream = subset(Final, Type == "Downstream(t)")
Downstream$family_downstream_t = paste(Downstream$family)
Downstream$number_TE_downstream_t = paste(Downstream$number_repeats2)
Downstream$Id = paste(Downstream$ID, Downstream$Population, sep = "::")
Downstream$TE_family_multi_downstream = Downstream$multi2_2
D = select(Downstream, Id, number_TE_downstream_t, family_downstream_t, TE_family_multi_downstream)

temp1 = merge(Tr, U, by = "Id")
Merged_Transcript = merge(temp1, D, by = "Id")

Merged_Transcript = apply(Merged_Transcript, 2, as.character)

write.csv(Merged_Transcript, "TE_fam_number_transcript.csv") 


#Homologs

Transcripts = subset(Final, Type == "Homolog(= h)")
Transcripts$Id = Transcripts$pop1
Transcripts$number_TE_homolog = paste(Transcripts$number_repeats2)
Transcripts$TE_family_homolog = paste(Transcripts$family)
Transcripts$population_hom = Transcripts$Population
Transcripts$Population = Transcripts$pop2
Transcripts$TE_family_multi = Transcripts$multi2_2
Tr = select(Transcripts, Id, Population, population_hom, TE_family_homolog, TE_family_multi, number_TE_homolog)

Upstream = subset(Final, Type == "Upstream(h)")
Upstream$family_upstream_h = paste(Upstream$family)
Upstream$number_TE_upstream_h = paste(Upstream$number_repeats2)
Upstream$Id = Upstream$pop1
Upstream$population_hom = Upstream$Population
Upstream$Population = Upstream$pop2
Upstream$TE_family_multi_upstream = Upstream$multi2_2
U = select(Upstream, Id, population_hom, number_TE_upstream_h, family_upstream_h, TE_family_multi_upstream)

Downstream = subset(Final, Type == "Downstream(h)")
Downstream$family_upstream_d = paste(Downstream$family)
Downstream$number_TE_upstream_d = paste(Downstream$number_repeats2)
Downstream$Id = Downstream$pop1
Downstream$population_hom = Downstream$Population
Downstream$Population = Downstream$pop2
Downstream$TE_family_multi_downstream = Downstream$multi2_2
D = select(Downstream, Id, population_hom, number_TE_upstream_d, family_upstream_d, TE_family_multi_downstream)

temp1 = merge(Tr, U, by = c("Id", "population_hom"))
Merged_Transcript = merge(temp1, D, by = c("Id", "population_hom"))

Merged_Transcript = apply(Merged_Transcript, 2, as.character)
write.csv(Merged_Transcript, "TE_fam_number_homolog.csv")

#########################################
#Get dataframe with rel TE overlap
#########################################
#Load data (= output of bedtools merge)
new_colnames = c("chromosome", "transcript_start", "transcript_end","chrom", "score", "strand_transcript", "chromosome", "overlap_start", "overlap_end","name", "overlap")

#Transcripts
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

#Now homologs
new_colnames = c("chromosome", "transcript_start", "transcript_end","chrom", "score", "strand_transcript", "Population", "chromosome", "overlap_start", "overlap_end","name", "overlap")


AK52 <- read.csv("AK5_final_file_hom.bed", header = FALSE, sep="\t")
colnames(AK52) <- new_colnames
AK52$Population = "AK5"

DK52 <- read.csv("DK5_final_file_hom.bed", header = FALSE, sep="\t")
colnames(DK52) <- new_colnames
DK52$Population = "DK5"

GI52 <- read.csv("GI5_final_file_hom.bed", header = FALSE, sep="\t")
colnames(GI52) <- new_colnames
GI52$Population = "GI5"

SW52 <- read.csv("SW5_final_file_hom.bed", header = FALSE, sep="\t")
colnames(SW52) <- new_colnames
SW52$Population = "SW5"

UM2 <- read.csv("UM_final_file_hom.bed", header = FALSE, sep="\t")
colnames(UM2) <- new_colnames
UM2$Population = "UM"

YE2 <- read.csv("YE_final_file_hom.bed", header = FALSE, sep="\t")
colnames(YE2) <- new_colnames
YE2$Population = "YE"

Zamb2 <- read.csv("Zamb_final_file_hom.bed", header = FALSE, sep="\t")
colnames(Zamb2) <- new_colnames
Zamb2$Population = "Zamb"

#Merge transcripts
All_merged <- rbind(AK5, DK5, GI5, SW5, UM, YE, Zamb)
All_merged$index <- 1:nrow(All_merged)
All_merged = All_merged[,c(ncol(All_merged),1:(ncol(All_merged)-1))]
All_merged$seq_id = All_merged$chrom
All_merged = All_merged %>% 
  separate("chrom", c("Id", "Type", "Population", "x"), sep = "_")
All_merged$population_hom = "NA"
All_merged$x = NULL
All_merged$Category = "De_novo_Transcript"

#Merge homologs
All_merged2 <- rbind(AK52, DK52, GI52, SW52, UM2, YE2, Zamb2)
All_merged2$index <- 1:nrow(All_merged2)
All_merged2 = All_merged2[,c(ncol(All_merged2),1:(ncol(All_merged2)-1))]
All_merged2$seq_id = All_merged2$chrom
All_merged2 = All_merged2 %>% 
  separate("chrom", c("Category", "Id", "population_hom", "x"), sep = "_")
All_merged2 = All_merged2 %>% 
  separate("x", c("Population", "Type", "population_hom"), sep = "XX")

#Bind together and remove old data
All_merged_new = rbind(All_merged, All_merged2)
All_merged = All_merged_new

rm(AK5, DK5, GI5, SW5, UM, YE, Zamb,AK52, DK52, GI52, SW52, UM2, YE2, Zamb2, All_merged2, All_merged_new)

#Rename some columns
All_merged[which(All_merged$Id == "Noncoding"), "Category"] <- "Intergenic"
All_merged[which(All_merged$Id == "Noncoding"), "Type"] <- "Intergenic"
All_merged[which(All_merged$Type == "End"), "Type"] <- "Downstream"

#Get the overlap
df = All_merged
rm(All_merged) #rm unecessary

df[which(df$strand_transcript == "."), "strand_transcript"] <- "+"
df$same_strand = paste(df$strand_transcript ,df$name,sep="_")

#Exclude seqs where overlap is not on the same strand 
df[which(df$same_strand == "-_+"), "overlap"] <- 0
df[which(df$same_strand == "+_-"), "overlap"] <- 0

df = df %>% group_by(df$seq_id) %>% mutate(total_overlap = sum(overlap)) #Get the overlap
df$rel_overlap = df$total_overlap/(df$transcript_end-df$transcript_start) #Get the relativeoverlap
control = subset(df, rel_overlap > 1) #control that it worked 

df$seq_id2 = paste(df$seq_id, df$strand_transcript)

#Make these categories for final plot
merged = distinct(df, seq_id2, .keep_all = TRUE) #Only keep unique id

merged2 = select(merged, Id, chromosome, Type, Category, Population, rel_overlap)

merged2[which(merged2$Type == "Intergenic"), "Category"] <- "Intergenic"
merged2[which(merged2$Type == "us"), "Type"] <- "Upstream(h)"
merged2[which(merged2$Type == "ds"), "Type"] <- "Downstream(h)"
merged2[which(merged2$Type == "hom"), "Type"] <- "Homolog(= h)"
merged2[which(merged2$Type == "Upstream"), "Type"] <- "Upstream(t)"
merged2[which(merged2$Type == "Downstream"), "Type"] <- "Downstream(t)"
merged2[which(merged2$Type == "Transcript"), "Type"] <- "Transcript(= t)"

#Prepare subdata for big merged data
data = select(merged, Id, transcript_start, transcript_end, Type, Population)
Transcripts = subset(data, Type == "Transcript")
Transcripts$Id = paste(Transcripts$Id, Transcripts$Population, sep = "::")
Homologs = subset(data, Type == "hom")
Homologs = Homologs %>%
  separate(`df$seq_id`, c("seq_id", "bar","Population"), "XX")
Homologs = select(Homologs, Id, transcript_start, transcript_end, Type, Population)
Transcripts$`df$seq_id` = NULL
df = rbind(Transcripts, Homologs)
df = select(df, Id, transcript_start, transcript_end, Population)
write.csv(df, "Positions_trans_hom.csv")


Transcripts = subset(merged2, Type == "Transcript(= t)")
Transcripts$Id = paste(Transcripts$Id, Transcripts$Population, sep = "::")
Tr = select(Transcripts, chromosome , Id, Population, rel_overlap)

Upstream = subset(merged2, Type == "Upstream(t)")
Upstream$rel_overlap_upstream_t = Upstream$rel_overlap
Upstream$Id = paste(Upstream$Id, Upstream$Population, sep = "::")
U = select(Upstream,Id, rel_overlap_upstream_t)

Downstream = subset(merged2, Type == "Downstream(t)")
Downstream$rel_overlap_downstream_t = Downstream$rel_overlap
Downstream$Id = paste(Downstream$Id, Downstream$Population, sep = "::")
D = select(Downstream,Id, rel_overlap_downstream_t)

temp1 = merge(Tr, U, by = "Id")
Merged_Transcript = merge(temp1, D, by = "Id")
Merged_Transcript$`df$seq_id`= NULL
Merged_Transcript$`df$seq_id.x`= NULL
Merged_Transcript$`df$seq_id.y`= NULL

Merged_Transcript$Transcription_status = "Transcript"
Merged_Transcript[which(Merged_Transcript$rel_overlap >= 0.8), "Transcription_status"] <- "Transcribed_TE"

TEs = subset(Merged_Transcript, Transcription_status == "Transcribed_TE")
ID_TE = TEs$Id

write.csv(Merged_Transcript, "Rel_TE_transcript.csv")

#Now homologs
Transcripts = subset(merged2, Type == "Homolog(= h)")
Transcripts = Transcripts %>% 
  separate("df$seq_id", c("population", "Type", "population_hom"), sep = "XX")
Transcripts$rel_overlap_homolog = Transcripts$rel_overlap
Transcripts$Id2 = paste(Transcripts$Id, Transcripts$population_hom)
Tr = select(Transcripts, chromosome, Id, Id2, Population, population_hom, rel_overlap_homolog)

Upstream = subset(merged2, Type == "Upstream(h)")
Upstream = Upstream %>% 
  separate("df$seq_id", c("population", "Type", "population_hom"), sep = "XX")
Upstream$rel_overlap_upstream_h = Upstream$rel_overlap
Upstream$Id2 = paste(Upstream$Id, Upstream$population_hom)
U = select(Upstream, Id, Id2, rel_overlap_upstream_h)

Downstream = subset(merged2, Type == "Downstream(h)")
Downstream = Downstream %>% 
  separate("df$seq_id", c("population", "Type", "population_hom"), sep = "XX")
Downstream$rel_overlap_downstream_h = Downstream$rel_overlap
Downstream$Id2 = paste(Downstream$Id, Downstream$population_hom)
D = select(Downstream, Id, Id2, rel_overlap_downstream_h)

temp1 = merge(Tr, U, by = "Id2")
Merged_Homologs = merge(temp1, D, by = "Id2")
Merged_Homologs$`Id.x`= NULL
Merged_Homologs$`Id.y`= NULL
Merged_Homologs$`Id2`= NULL

Merged_Homologs = select(Merged_Homologs, Id, chromosome, population_hom,  rel_overlap_homolog, rel_overlap_upstream_h, rel_overlap_downstream_h)
Merged_Homologs$Transcription_status = "Non_transcribed_Homolog"
Merged_Homologs[which(Merged_Homologs$Id %in% ID_TE), "Transcription_status"] <- "Non_transcribed_Homolog_TE"

write.csv(Merged_Homologs, "Rel_TE_homologs.csv")

#########################################
#Get dataframe with Number of motifs
#########################################
Transcripts <- read.delim("All_intergenic_de_novo_transcripts.csv")

Transcripts$type2 = Transcripts$type

Transcripts = Transcripts %>% 
  separate("seq_name", c("type", "chrom", "position"), sep = "::")
Transcripts = Transcripts %>% 
  separate("chrom", c("chrom", "position"), sep = ":")
Transcripts$final = str_sub(Transcripts$position, end = -4)
Transcripts = Transcripts %>% 
  separate("final", c("start", "end"), sep = "-")
Transcripts$size = as.numeric(Transcripts$end) - as.numeric(Transcripts$start)
Transcripts[which(Transcripts$type2 == "noncoding"), "size"] <- 1100 #As I know all noncoding ones have this size
Transcripts$type = Transcripts$type2
Transcripts$id = Transcripts$unique_id
Transcripts = Transcripts %>% 
  separate("id", c("seq_id", "pos"), sep = "::")
Transcripts$id = paste(Transcripts$seq_id, Transcripts$population, sep = "_")

Transcripts = subset(Transcripts, size == "1100")

Transcripts$unique_id = paste(Transcripts$seq_id, Transcripts$population, sep = "::")
Transcripts$Population = Transcripts$population


#Now the homologs
Homologs <- read.delim("Motifs_homolog_upstream.csv")
Homologs = Homologs %>% 
  separate("seq_name", c("Id",  "Population", "position"), sep = "::")
Homologs = Homologs %>% 
  separate("Id", c("type", "Id"), sep = "_")
Homologs = Homologs %>% 
  separate("position", c("chrom", "position"), sep = ":")
Homologs$final = str_sub(Homologs$position, end = -4)
Homologs = Homologs %>% 
  separate("final", c("start", "end"), sep = "-")
Homologs$size = as.numeric(Homologs$end) - as.numeric(Homologs$start)
Homologs$type = "Homolog"
Homologs$type2 = "Homolog"
Homologs = Homologs %>% 
  separate("Population", c("Population", "population_hom", "population2"), sep = "_")
Homologs$unique_id = paste(Homologs$Id, Homologs$Population, sep = "::")
Homologs$Id = paste(Homologs$Id, Homologs$Population, sep = "::")
Homologs = subset(Homologs, size == 1100)

#Make df 
Homologs$number_motifs_per_transcript_h = Homologs$number_motifs_per_transcript
Homologs$number_motifs_high_per_transcript_h = Homologs$number_motifs_high_per_transcript
Hom = select(Homologs, Id, population_hom, number_motifs_per_transcript, number_motifs_high_per_transcript)
write.csv(Hom, "Motifs_upstream_homologs.csv")

Transcripts$Id = paste(Transcripts$seq_id, Transcripts$population, sep = "::")
Transcripts$Population = Transcripts$population
Trans = select(Transcripts, Id, Population,number_motifs_per_transcript, number_motifs_high_per_transcript)
write.csv(Trans, "Motifs_upstream_transcripts.csv")

#########################################
#Get dataframe with Number of core
#########################################
Transcripts <- read.delim("Promotor_Intergenic_de_novo_transcripts.csv")
Transcripts$type2 = Transcripts$type

Transcripts = Transcripts %>% 
  separate("seq_name", c("type", "chrom", "position"), sep = "::")
Transcripts = Transcripts %>% 
  separate("chrom", c("chrom", "position"), sep = ":")
Transcripts$final = str_sub(Transcripts$position, end = -4)
Transcripts = Transcripts %>% 
  separate("final", c("start", "end"), sep = "-")
Transcripts$size = as.numeric(Transcripts$end) - as.numeric(Transcripts$start)
Transcripts[which(Transcripts$type2 == "noncoding"), "size"] <- 1100 #As I know all noncoding ones have this size
Transcripts$type = Transcripts$type2
Transcripts = Transcripts %>% 
  separate("unique_id", c("seq_id", "pos"), sep = "::")
Transcripts$unique_id = paste(Final_df2$seq_id, Final_df2$population, sep = "_")
Transcripts$promotor_size = as.numeric(Transcripts$size)
Transcripts[which(Transcripts$size >= 300), "promotor_size"] <- 300
Transcripts = subset(Transcripts, size == "300")

Transcripts$id = paste(Transcripts$seq_id, Transcripts$population, sep = "_")
Transcripts$unique_id = paste(Transcripts$seq_id, Transcripts$population, sep = "::")
Transcripts$Population = Transcripts$population

Homologs <- read.delim("Promotor_homolog_upstream.csv")
Homologs = Homologs %>% 
  separate("seq_name", c("Id",  "Population", "position"), sep = "::")
Homologs = Homologs %>% 
  separate("Id", c("type", "Id"), sep = "_")
Homologs = Homologs %>% 
  separate("position", c("chrom", "position"), sep = ":")
Homologs$final = str_sub(Homologs$position, end = -4)
Homologs = Homologs %>% 
  separate("final", c("start", "end"), sep = "-")
Homologs$size = as.numeric(Homologs$end) - as.numeric(Homologs$start)
Homologs$type = "Homolog"
Homologs$type2 = "Homolog"
Homologs = Homologs %>% 
  separate("Population", c("Population", "population_hom", "population2"), sep = "_")
Homologs$unique_id = paste(Homologs$Id, Homologs$Population, sep = "::")
Homologs$Id = paste(Homologs$Id, Homologs$Population, sep = "::")
Homologs = subset(Homologs, size >= 300)

#Make df for both
#Transcripts
Transcripts$Id = paste(Transcripts$seq_id, Transcripts$population, sep = "::")
Transcripts$Population = Transcripts$population
Trans = select(Transcripts, Id, Population,number_core_per_transcript, number_core_high_per_transcript)
write.csv(Trans, "Core_upstream_transcripts.csv")

#Homologs
Homologs$number_core_per_transcript_h = Homologs$number_core_per_transcript
Homologs$number_core_high_per_transcript_h = Homologs$number_core_high_per_transcript
Hom = select(Homologs, Id, population_hom,number_core_per_transcript, number_core_high_per_transcript)
write.csv(Hom, "Core_upstream_homologs.csv")


#########################################
#Merge all data
#########################################

Core_homologs <- read_csv("F:/MasterThesisMarie_AllData/Make big data merged/Core_upstream_homologs.csv")
colnames(Core_homologs) = c("Index", "Seq_id", "Population", "Number_core", "Number_core_high")

Core_transcripts <- read_csv("F:/MasterThesisMarie_AllData/Make big data merged/Core_upstream_transcripts.csv")
colnames(Core_transcripts) = c("Index", "Seq_id", "Population", "Number_core", "Number_core_high")

Core = rbind(Core_homologs, Core_transcripts)

Motif_homologs <- read_csv("Motifs_upstream_homologs.csv")
colnames(Motif_homologs) = c("Index", "Seq_id", "Population", "Number_motifs", "Number_motifs_high")

Motif_transcripts <- read_csv("Motifs_upstream_transcripts.csv")
colnames(Motif_transcripts) = c("Index", "Seq_id", "Population", "Number_motifs", "Number_motifs_high")

Motifs = rbind(Motif_homologs, Motif_transcripts)

TE_overlap_homologs <- read_csv("Rel_TE_homologs.csv")
colnames(TE_overlap_homologs) = c("Index", "Seq_id", "Chromosome","Population", "Relative_TE_overlap", "Relative_TE_overlap_upstream", "Relative_TE_overlap_downstream", "Transcription_status")

TE_overlap_transcripts <- read_csv("Rel_TE_transcript.csv")
colnames(TE_overlap_transcripts) = c("Index", "Seq_id", "Chromosome","Population", "Relative_TE_overlap", "Relative_TE_overlap_upstream", "Relative_TE_overlap_downstream", "Transcription_status")

Overlap = rbind(TE_overlap_homologs, TE_overlap_transcripts)

TE_fam_number_homologs <- read_csv("TE_fam_number_homolog.csv", na = c("NA", "NULL", ""))
TE_fam_number_homologs$population_hom = NULL

colnames(TE_fam_number_homologs) = c("Index", "Seq_id", "Population", "TE_family", "TE_family2", "Number_TE", "Number_TE_upstream", "TE_family_upstream", "TE_family_upstream2", "Number_TE_downstream", "TE_family_downstream", "TE_family_downstream2")

TE_fam_number_transcript <- read_csv("TE_fam_number_transcript.csv", na = c("NA", "NULL", ""))
colnames(TE_fam_number_transcript) = c("Index", "Seq_id", "Population", "TE_family", "TE_family2", "Number_TE", "Number_TE_upstream", "TE_family_upstream2", "TE_family_upstream", "Number_TE_downstream", "TE_family_downstream", "TE_family_downstream2")
Family = rbind(TE_fam_number_homologs, TE_fam_number_transcript)

Temp = merge(Overlap, Family, by = c("Seq_id", "Population"), all = TRUE)
Temp2 = merge(Temp, Motifs, by = c("Seq_id", "Population"), all = TRUE)
Temp3 = merge(Temp2, Core, by = c("Seq_id", "Population"), all = TRUE)

Temp3$Index.x = NULL
Temp3$Index.y = NULL

#Now load the orthogroup file 1
Orthogroup_file1 <- read_delim("Orthogroup_file1.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

Orthogroup_file1$Seq_id = Orthogroup_file1$qseqid

DF = select(Orthogroup_file1, Population, Seq_id ,Number_of_populations_sharing_transcript)

Temp4 = merge(Temp3, DF, by = c("Seq_id", "Population"), all = TRUE)


#Now load the orthogroup file 2

df2 <- read.table("Orthogroups.txt", sep = ";", header = TRUE, stringsAsFactors = FALSE)

#Split the ID
df2$IDs_in_orthogroup <- strsplit(df2$IDs_in_orthogroup, ",")

# Explode the df to be able to merge it

df3 <- df2 %>%
  tidyr::unnest(IDs_in_orthogroup) %>%
  mutate(Seq_id = IDs_in_orthogroup) %>%
  select(-IDs_in_orthogroup)


Merged_df = merge(Temp4, df3, by = "Seq_id")

#Get rid of some unecessary columns
Merged_df$Number_of_populations_sharing_transcript = NULL
Merged_df$Index.x = NULL
Merged_df$Index.y = NULL


#Name some columns
Merged_df[which(Merged_df$Transcription_status %in%  c("Non_transcribed_Homolog_TE", "Non_transcribed_Homolog")), "To_name"] <- "Homolog of"

Merged_df[which(Merged_df$Transcription_status %in%  c("Transcribed_TE", "Transcript")), "To_name"] <- "Transcript"

Merged_df$Transcript_Id = Merged_df$Seq_id 

Merged_df$Seq_id = paste(Merged_df$To_name, Merged_df$Number_Orthogroup)

#Choose which columns to keep
Merged_df2 = Merged_df[,c(1,21, 25, 2,3,22,23,7, 4,5,6,8,9,10, 12, 13, 11,15,16,14,17,18,19,20)]

#Final df for statistics
write.csv(Merged_df2, "Merged_data_de_novo_transcript_homolog_correct.csv") 