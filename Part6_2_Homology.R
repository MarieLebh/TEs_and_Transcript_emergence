#Make plots for homology

library(readr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(scales)

#Part I: Shared transcripts

#Load the data
Orthogroup<- read_delim("Orthogroup_file1.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

#Get number of blast searches
sum(Orthogroup$Number_of_populations_without_the_transcript) #with TE

#Exclude overlapping TE transcripts
TE <- read_csv("Data_rel_TE_overlap_new.csv") #Load the data with TE information

TE$qseqid = paste(TE$Id, TE$Population, sep = "::")
TE = subset(TE, Type == "Transcript") #Only keep transcripts

df = merge(TE, Orthogroup, by = "qseqid")

df_TE = subset(df, rel_overlap < 0.8) #Exclude TE transcripts
sum(df_TE$Number_of_populations_without_the_transcript) #get number of blast searches without TE

Orthogroup = df_TE #Change df for plot

#Count the percentage
shared = Orthogroup %>% group_by(Number_of_populations_sharing_transcript) %>%
  summarize(n = n_distinct(Number_of_populations_sharing_transcript), Count = n())

shared$Percentage = (shared$Count /3528) * 100 #3528 (= number de novo transripts)

#Plot Fig1_E
Orthogroup$Number_of_populations_sharing_transcript = as.numeric(Orthogroup$Number_of_populations_sharing_transcript)

ggplot()+
  geom_histogram(data = Orthogroup, aes(x = Number_of_populations_sharing_transcript), fill = "grey")+
  labs(x = "Number of Lines shared", y = "Number of de novo transcripts", size = 25)+
  theme_pubclean()
  
ggsave("Number_of_shared_transcripts.pdf", width = 11)


#Part II: Non expressed homologs

#Load the data 
df1 <- read_csv("Non_expressed_homologues_80cov.txt")
data <- data.frame(lapply(df1, function(x) {gsub("s", "shared_transcript", x)}))
data[2:9] <- data.frame(lapply(data[2:9], function(x) {gsub("0", "no_blast_hit", x)}))
data[2:9] <- data.frame(lapply(data[2:9], function(x) {gsub("1", "blast hit(s)", x)}))
df1 = data
rm(data)

#Load the TE data and exclude TE transcripts
TE <- read_csv("Data_rel_TE_overlap_new.csv")

TE$transcript_id = paste(TE$Id, TE$Population, sep = "::")
TE = subset(TE, Type == "Transcript")

df = merge(df1, TE, by.x = "transcript_id")

df_TE = subset(df, rel_overlap < 0.8)

df_TE = df_TE %>%
  separate(transcript_id, c("transcript_id", "population"), "::")

df_TE = df_TE[, c(1,2,3,4,5,6,7,8,9)] #Choose the columns to keep

df_TE = df_TE %>%
  arrange(population)  #Sort by population

#Change this for the final plot
df_TE [which(df_TE $population == "AK5"), "matchAK5"] <- "-"
df_TE [which(df_TE $population == "DK5"), "matchDK5"] <- "-"
df_TE [which(df_TE $population == "GI5"), "matchGI5"] <- "-"
df_TE [which(df_TE $population == "SW5"), "matchSW5"] <- "-"
df_TE [which(df_TE $population == "UM"), "matchUM"] <- "-"
df_TE [which(df_TE $population == "YE"), "matchYE"] <- "-"
df_TE [which(df_TE $population == "Zamb"), "matchZamb"] <- "-"

write.csv(x,file="Transcripts_Blast_Hits.csv", row.names=FALSE) 