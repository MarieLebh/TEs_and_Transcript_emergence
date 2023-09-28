#Get number TE and fraction of overlapping TEs 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library( data.table )
library( intervals )

###################################
#Load the data
###################################

noncoding <- read_delim("All_noncoding.csv",
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


end$type = "Downstream"
noncoding$type = "Intergenic"
transcript$type = "Transcript"
upstream$type = "Upstream"

#Merge all dfs together
df4 = rbind(end, noncoding, transcript, upstream)

#Make seq id column
df4$id_seq = paste(df4$type, df4$id_seq, sep = "_")
df4$id_seq = paste(df4$population, df4$id_seq, sep = "_")

#Calculate size of the upstream region 
df4$start_seq <- as.numeric(df4$start_seq)
df4$end_seq <- as.numeric(df4$end_seq)
df4$size <- df4$end_seq - df4$start_seq

#Exluce sequences not overlapping the TE on the same strand
df4$same_strand = paste(df4$strand_seq ,df4$strand_te,sep="_")
df4[which(df4$same_strand == "-_+"), "overlap"] <- 0
df4[which(df4$same_strand == "+_-"), "overlap"] <- 0
df4[which(df4$same_strand == "._-"), "overlap"] <- 0
df4[which(df4$same_strand == "-_."), "overlap"] <- 0

#Get the number of TEs
df4$number_repeats = 1
df4[which(df4$overlap == 0), "number_repeats"] <- 0
df4 = df4 %>% group_by(id_seq) %>% mutate(number_repeats2 = sum(number_repeats))
df4$number_repeats2_size = df4$number_repeats2 / df4$size

#Split df column to get class, subclass and family
df4$TE <- sub("(.*\\()(.*)(\\))", "\\2", df4$id_te)
df4 = df4 %>% 
  separate("TE", c("family", "type2"), sep = ";")
df4 = df4 %>% 
  separate("family", c("family", "subclass", "class"), sep = ",")
df4 = df4 %>% 
  separate("type2", c("useless", "completeness"), sep = "=")

df4$unique = df4$id_seq #save this column for later

df4 = df4 %>% 
  separate("id_seq", c("Population", "Type", "ID"), sep = "_")

df4$ID = paste(df4$ID, df4$population, sep = "::")

df4[which(df4$overlap == 0), "family"] <- NA

#Make a subset with all TEs that do overlap
df4_2 = subset(df4, df4$overlap != 0)

#Split the unique column
df4_2 = df4_2 %>% 
  separate("unique", c("Population", "Type", "ID", "pop1", "pop2", "pop3"), sep = "_")

#Split the data (because the unique column is named a bit differently dor intergenic)
Noncoding = subset(df4_2, Type == "Intergenic")
Noncoding$ID = Noncoding$pop1
Else = subset(df4_2, Type != "Intergenic")

df4_2 = rbind(Noncoding, Else) #bind them together again

df4_2[which(df4_2$Type == "Intergenic"), "Type"] <- "Noncoding" #Change this so it matches with the TE data

df4_2$Id = paste(df4_2$Type, df4_2$ID, df4_2$Population, sep = "_")


###################################
#Load the TE info data
###################################

TE <- read_csv("Data_rel_TE_overlap_new.csv")
TE1 = subset(TE, Type == "Intergenic")
TE2 = subset(TE, Type != "Intergenic")
TE2$Id = paste(TE2$Type, TE2$Id, TE2$Population, sep = "_")
TE = rbind(TE1, TE2) #name the columns appropriately and merge back together

test = merge(df4_2, TE, by = "Id", all.y = FALSE, all.x = TRUE)

test2 = subset(test, Type.x == "Noncoding" |Type.x == "Transcript")
test2 = subset(test2, rel_overlap >= 0.8) #Make a subset with only TE seqs

test2$ID = paste(test2$ID, test2$Population.x, sep = "_") 
id = test2$ID #save ids of the TE seqs as list

test$ID = paste(test$ID, test$Population.x, sep = "_")

df_figure = subset(test, !ID %in% id) #Exclude TE seqs from df (use this for plot)

df_figure[which(df_figure$Type.x == "Noncoding"), "Type.x"] <- "Intergenic" #Rename the column again

###################################
#Make Figure 3B
###################################

my_colours = c(brewer.pal(7, "Set2"), brewer.pal(8, "Dark2"))

axis_order = c("Intergenic", "Upstream", "Transcript", "Downstream")

Bar = ggbarplot(df_figure, y = "overlap",
                position = position_fill(),
                col = "family", fill = "family")+
  facet_wrap(~factor(Type.x, levels = axis_order), ncol = 2, nrow = 2)+
  labs(fill='TE Family') +
  labs(col='TE Family') +
  theme_pubclean()+
  scale_color_manual(values = my_colours)+
  scale_fill_manual(values = my_colours)+
  labs(x= " ", y = " ")+
  theme(legend.position = "top",  
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  theme(strip.text = element_text(size = 20),
        axis.title = element_text(size = 20), legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),  
        legend.text = element_text(size = 20), axis.text.x = element_blank(), axis.ticks = element_blank())+
  ggtitle("Overlapping TE families")+
  theme(title = element_text(size = 20))

a = Bar + 
  coord_polar("y", start=0) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())

ggsave("TE_figure_new.pdf",a, height = 15, width = 15)

###################################
#Get the number of TE 
###################################

#Now restructure the dataset so that each unique sequence has only one column
df4 = df4 %>% group_by(unique) %>% mutate(family2 = n_distinct(family, na.rm = TRUE))

#If multiple TEs overlap set this column to Multiple TE
df4$family_multi = df4$family
df4[which(df4$family2 > 1), "family_multi"] <- "Multiple_TE"

temp = df4 %>%
  group_by(unique) %>%
  summarise_all(funs(list(na.omit(.))))

temp%>%
  mutate(family_multi = map_dbl(multi2, first))

temp$multi2_2 = sapply(temp$family_multi, `[`,1)

Final = select(temp, unique, family, number_repeats2, number_repeats2_size, multi2_2) #Select all columns for final df

Final[which(Final$family == "character(0)"), "family"] <- NA

Final$number_repeats2 = as.vector(Final$number_repeats2)
Final$number_repeats2_size = as.vector(Final$number_repeats2_size)

Final <- Final %>%
  mutate(number_repeats2 = sapply(number_repeats2, `[`, 1))#Extract number of TE

Final <- Final %>%
  mutate(number_repeats2_size = sapply(number_repeats2_size, `[`, 1)) #Extract number of TE divided by size

Final = Final %>% 
  separate("unique", c("Population", "Type", "ID", "pop1", "pop2", "pop3"), sep = "_")

#Make dfs fit
Noncoding = subset(Final, Type == "Intergenic")
Noncoding$ID = Noncoding$pop1
Else = subset(Final, Type != "Intergenic")

Final = rbind(Noncoding, Else) 

Final[which(Final$Type == "Intergenic"), "Type"] <- "Noncoding"

Final$Id = paste(Final$Type, Final$ID, Final$Population, sep = "_")

#Read in TE data (again)
TE <- read_csv("Data_rel_TE_overlap_new.csv")
TE1 = subset(TE, Type == "Intergenic")
TE2 = subset(TE, Type != "Intergenic")
TE2$Id = paste(TE2$Type, TE2$Id, TE2$Population, sep = "_")
TE = rbind(TE1, TE2)

merge = merge(Final, TE, by = "Id")

merge2 = subset(merge, Type.x == "Noncoding" |Type.x == "Transcript")
merge2 = subset(merge2, rel_overlap >= 0.8)
merge2$ID = paste(merge2$ID, merge2$Population.x, sep = "_")
id = merge2$ID

merge$ID = paste(merge$ID, merge$Population.x, sep = "_")

final = subset(merge, !ID %in% id) #Again exclude TE seqs (intergenic + transcript)

#Do the stats (here with number/size)
x = final%>%
  dunn_test(number_repeats2_size ~ Type.x)

#Calculate the mean number of TE
x = subset(final, Type.x == "Transcript")
mean(as.numeric(lol$number_repeats2))

x = subset(final, Type.x == "Upstream")
mean(as.numeric(lol$number_repeats2))

x = subset(final, Type.x == "Downstream")
mean(as.numeric(lol$number_repeats2))

x = subset(final, Type.x == "Noncoding")
mean(as.numeric(lol$number_repeats2))

