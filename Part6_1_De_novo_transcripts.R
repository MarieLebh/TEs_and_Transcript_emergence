#Plot the Blast results and distribution on chromosomes

#Import packages
library(readr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggpubr)
library(ggforce)
library(ggpattern)
library(dunn.test)
library(rstatix)

#Plot Blast results (Input is an excel sheet with the numbers of each filtering step)
novo <- read_excel("Blast.xlsx", sheet = "De_novo")#For de novo transcripts
anno <- read_excel("Blast.xlsx", sheet = "Annotated") #annotated data

#Plot de novo transcripts
blast_novo <- ggplot(novo, aes(y=Number, x=Population, fill = Dataset)) + 
  geom_bar(position="dodge", stat="identity", orientation = "x")+
  theme_pubclean()+
  labs(x = "Line")+
  scale_fill_manual(values = c("salmon2", "skyblue", "olivedrab", "blue"))+
  theme(legend.position = "top",  
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave("Fig1_A_new.pdf", height = 13, width = 14)

#Plot annotated transcripts
blast_a <- ggplot(anno, aes(y=Number, x=Population, fill = Dataset)) + 
  geom_bar(position="dodge", stat="identity", orientation = "x")+
  labs(x = "Line")+
  theme_pubclean()+
  scale_fill_manual(values = c("salmon2", "skyblue", "olivedrab"))+
  theme(legend.position = "top",  
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))



#Part 2: Plot the distribution on chromosomes (Input is the file containing all transcript properties)
Data <- read_csv("Transcript_features_corrected.csv") #Contains gene overlap, chrom and features

Data$transcript_id = paste(Data$Transcript_ID, Data$population, sep = "::")

#Join with TE info
TE <- read_csv("F:/MasterThesisMarie_AllData/Data/Data_rel_TE_overlap_new.csv")

TE$transcript_id = paste(TE$Id, TE$Population, sep = "::")
TE = subset(TE, Type == "Transcript")

Data = merge(Data, TE, by = "transcript_id", all.x = TRUE)

Data = subset(Data, rel_overlap < 0.8  | is.na(rel_overlap) == TRUE) #remove TE (ONLY intergenic de novo)

#Shorten group names
Data$Group <- paste(Data$overlap,Data$type,sep="_")
Data[which(Data$Group == "gene_overlap_de_novo"), "Group"] <- "Dn_g"
Data[which(Data$Group == "gene_overlap_annotated"), "Group"] <- "An_g"
Data[which(Data$Group == "intergenic_de_novo"), "Group"] <- "Dn_i"
Data[which(Data$Group == "intergenic_annotated"), "Group"] <- "An_i"

#Change population names
Data[which(Data$population == "AK5"), "population"] <- "FI"
Data[which(Data$population == "DK5"), "population"] <- "DK"
Data[which(Data$population == "GI5"), "population"] <- "ES"
Data[which(Data$population == "SW5"), "population"] <- "SE"
Data[which(Data$population == "UM"), "population"] <- "UA"
Data[which(Data$population == "YE"), "population"] <- "TR"
Data[which(Data$population == "Zamb"), "population"] <- "ZI"

#Split annotated an de novo
annotated <- Data[Data$type=="annotated",]
de_novo <- Data[Data$type=="de_novo",]

#Plot the distribution of the chromosomes:

de_novo2 <- de_novo[de_novo$overlap=="intergenic",] 

chrom <- ggplot(data = de_novo2, aes(y = population, fill = chromosome))+
  coord_flip()+
  scale_fill_manual(values = blue)+
  geom_bar(position = "dodge", width = 0.6)+
  #facet_wrap(~chromosome, scales = "free")+
  labs(x= "Number of transcripts", y = "Line")+
  theme_pubclean()+
  ggtitle(label = "Intergenic de novo transcripts")+
  theme(legend.position = "top",  
        panel.border = element_rect(colour = "black", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

ggsave("Supp_figure_chrom.pdf", height = 8.5, width = 8.5)

