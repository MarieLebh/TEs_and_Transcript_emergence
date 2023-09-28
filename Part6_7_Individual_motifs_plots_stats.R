#Get the number of individual motifs of de novo transcripts

#Load all packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyverse)
library(RColorBrewer)
library(data.table )
library(car)
library(stringr)
library(xtable)
library(readxl)
library(rstatix)

#Load the data and exclude core motifs from the data
#Three input files one for genes, noncoding (= intergenic) and de novo transcripts

Genes <- read.csv("Motif_number_genes_corrected.csv")
Genes[which(str_detect(Genes$motif_id, "POL") == TRUE), "Core"] <- "Core"
Genes[which(str_detect(Genes$motif_id, "POL") == FALSE), "Core"] <- "Other"
Genes = subset(Genes, Core == "Other")

Noncoding <- read.csv("Motif_number_noncoding_corrected.csv")
Noncoding[which(str_detect(Noncoding$motif_id, "POL") == TRUE), "Core"] <- "Core"
Noncoding[which(str_detect(Noncoding$motif_id, "POL") == FALSE), "Core"] <- "Other"
Noncoding = subset(Noncoding, Core == "Other")

Transcripts <- read.csv("Motif_number_de_novo_corrected.csv") #should include TE transcripts
Transcripts[which(str_detect(Transcripts$motif_id, "POL") == TRUE), "Core"] <- "Core"
Transcripts[which(str_detect(Transcripts$motif_id, "POL") == FALSE), "Core"] <- "Other"
Transcripts = subset(Transcripts, Core == "Other")

#Merge the data together
df <- rbind(Genes, Noncoding, Transcripts)
df$type2 = df$type #this is just for colouring

#discard all non 1100 seqs
df = subset(df, size == "1100")

rm(Genes, Noncoding, Transcripts)

#########################################################################################
#Load the TE infos from info_files
TEs = read.csv("TE.csv") #Contains all transcripts with >= 0.8 TE overlap (=TE)
id_TE = TEs$x

OverlapTE = read.csv("TE_overlap.csv") #Transcript ids with > 0.8 overlap
id_TE_overlap = OverlapTE$x

NoTE = read.csv("no_TE.csv") #No TE overlap
id_TE_no = NoTE$x

################################################
df$type_te = df$type
df$Id = paste(df$id, df$population, sep = "_")

#Merge imported ids with df
df[which(df$Id %in% id_TE), "type_te"] <- "dn_TE"
df[which(df$Id %in% id_TE_no), "type_te"] <- "dn_no_TE"
df[which(df$Id %in% id_TE_overlap), "type_te"] <- "dn_overlap_TE"

df2 <- transform(df, id = match(motif_id, unique(motif_id)))
df = df2

df$number_per_transcript_high = as.numeric(df$number_per_transcript_high)
df$number_per_transcript = as.numeric(df$number_per_transcript)
df$Legend = df$type
df[which(df$type == "de_novo_transcript"), "Legend"] <- "de_novo"

#Make the final plots

plot1 = ggplot(mean_se1, mapping = aes(x = motif_id, y = mean, col = type_te))+ 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position = position_dodge(width = 0.7))+
  geom_point(size = 3, position = position_dodge(width = 0.7))+
  scale_color_manual(values = c("darkgreen","lightgreen", "green", "#E57373", "#1E88E5"))+
  ylab(label = "Number of motifs per sequence (Rel_score >= 0.8)")+
  theme_pubclean()+
  coord_flip()+
  labs(col='Legend') +
  theme(legend.position = "top",  
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave("All_08.pdf", width = 16, height = 20)

plot2 = ggplot(mean_se2, mapping = aes(x = motif_id, y = mean, col = type_te))+ 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position = position_dodge(width = 0.7))+
  geom_point(size = 3, position = position_dodge(width = 0.7))+
  labs(fill='Legend') +
  scale_color_manual(values = c("darkgreen","lightgreen", "green", "#E57373", "#1E88E5"))+
  ylab(label = "Number of motifs per sequence (Rel_score >= 0.95)")+
  theme_pubclean()+
  coord_flip()+
  labs(col='Legend') +
  theme(legend.position = "top",  
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave("All_095.pdf", width = 16, height = 20)

#################################################################

#Do the stats

#Anova
stats1 <- df %>% 
  group_by(motif_id) %>%
  tukey_hsd(number_per_transcript ~ type_te)#treshold 0.8

stats2 <- df %>% 
  group_by(motif_id) %>%
  tukey_hsd(number_per_transcript_high ~ type_te) #treshold 0.95

#Summary stats

mean_se1 = df %>%
  group_by(motif_id, type_te) %>%
  get_summary_stats(number_per_transcript) #treshold 0.8


mean_se2 = df %>%
  group_by(motif_id, type_te) %>%
  get_summary_stats(number_per_transcript_high)#treshold 0.95


#Read in the info file about the motifs (Motif_id, family, class)

Motif_info <- read_excel("Motif_info.xlsx")

#Make results table for rel_score >= 0.8
other_means <- mean_se1 %>%
  filter(type_te != "noncoding")

noncoding_means <- mean_se1 %>%
  filter(type_te == "noncoding") %>%
  select(motif_id, mean_noncoding = mean)

df1_with_noncoding <- other_means %>%
  left_join(noncoding_means, by = "motif_id") %>%
  mutate(mean_diff = mean - mean_noncoding)

stats_non = subset(stats1, group2 == "noncoding")
stats_non$type_te = stats_non$group1

result <- df1_with_noncoding %>%
  inner_join(stats_non, by = c("motif_id", "type_te")) 


# Filter the dataframe 
filtered_df <- result %>%
  filter(mean_diff > 0 & p.adj < 0.05) %>%
  mutate(Significance = "X") %>%
  select(motif_id, type_te, Significance)

#Make format for table
result_df <- filtered_df %>%
  pivot_wider(names_from = type_te, values_from = Significance, values_fill = "-")

sub_result_df = subset(result_df, dn_overlap_TE == "X" |  dn_TE == "X" |  dn_no_TE == "X")
sub_result_df = sub_result_df[,c(1,3,2,4,5)]

Motif_table1 = merge(Motif_info, sub_result_df, by = "motif_id")


#Make results table for rel_score >= 0.95
other_means2 <- mean_se2 %>%
  filter(type_te != "noncoding")

noncoding_means2 <- mean_se2 %>%
  filter(type_te == "noncoding") %>%
  select(motif_id, mean_noncoding = mean)

df2_with_noncoding <- other_means2 %>%
  left_join(noncoding_means2, by = "motif_id") %>%
  mutate(mean_diff = mean - mean_noncoding)

stats_non2 = subset(stats2, group2 == "noncoding")
stats_non2$type_te = stats_non2$group1

result2 <- df2_with_noncoding %>%
  inner_join(stats_non2, by = c("motif_id", "type_te")) 

#Filter the dataframe 
filtered_df2 <- result2 %>%
  filter(mean_diff > 0 & p.adj < 0.05) %>%
  mutate(Significance = "X") %>%
  select(motif_id, type_te, Significance)

#Make the format right for table
result_df2 <- filtered_df2 %>%
  pivot_wider(names_from = type_te, values_from = Significance, values_fill = "-")

sub_result_df2 = subset(result_df2, dn_overlap_TE == "X" |  dn_TE == "X" |  dn_no_TE == "X")
sub_result_df2 = sub_result_df2[,c(1,2,3,5,4)]

#Merge the two tables
Motif_table1$rel_score = "0.8"
Motif_table2$rel_score = "0.95"
Motif_table_final = rbind(Motif_table1, Motif_table2) #Final table

write.csv(Motif_table_final, file = "Motif_table_supplement.csv")
