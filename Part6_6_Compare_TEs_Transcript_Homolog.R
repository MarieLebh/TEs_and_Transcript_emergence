#Compare TEs in transcripts and homologs

library(readr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(scales)

#Load the data

df <- read_csv("TE_analysis_homologs.txt") #Load the python output

TE <- read_csv("Data_rel_TE_overlap_new.csv") #Load the TE overlap data

TE$transcript_id = paste(TE$Id, TE$Population, sep = "::")
TE = subset(TE, Type == "Transcript")
z = merge(df, TE, by.x = "transcript_id") #merge both dataframes 
sub = subset(z, rel_overlap < 0.8) #Exclude TE transcripts

#Count the number for each category
stats = sub %>% 
  group_by(category) %>%
  summarize(n = n_distinct(category), Count = n())

stats$percent = (stats$Count/14058) * 100 #Get the percentages of each category

#Plot this
ggplot()+
  geom_histogram(data = sub, aes(x = population_homolog, fill = category),stat = "count")+
  labs(x = "Line of non expressed homolog", y = "Number")+
  theme_pubclean()+
  labs(fill='Legend') +
  scale_fill_manual(values = c("grey", "chartreuse", "azure4", "orange", "red", "yellow" ))+
  facet_wrap(~population_transcript)+
  theme(legend.position = "top", nrow(3))+
  theme(axis.title = element_text(size = 15), axis.text = element_text(size = 13), legend.text = element_text(size = 11))+
  ggtitle(label = "Presence of TEs in nonexpressed homolog", subtitle = "Only homologs with percent identity > 80%, Coverage >= 80%")

ggsave("TEs_present_in_homolog.pdf", height = 10, width = 11.5) #Save the plot
