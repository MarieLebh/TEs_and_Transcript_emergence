setwd("D:/Msc Data/Make circular plots")

library(rtracklayer)
library("circlize")
library(ggtree)
library("viridis")
library('dplyr')
library(ggplot2)
library(ggpubr)
library(readr)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library( data.table )
library( intervals )
library(ggpattern)

#Explanation what files I used:
# AK5finalGenome.masked.fa.fai = Genome index file
# AK5_intergenic_transcripts.gtf = Gtf file of all de novo transcrips
# AK5_chrom.txt = a bedfile with the chromosome positions
# AK5_te_annotation.gff = A gtf file of the ResonaTE annotation

## Read scaffold lengths from file
Scaff_len<-read.delim("AK5finalGenome.masked.fa.fai", sep="\t",stringsAsFactors = FALSE,header = FALSE)
Scaff_len_sorted<-Scaff_len[order(Scaff_len[,2], decreasing = TRUE),]
seqlengths<-Scaff_len[,2]
names(seqlengths)<-Scaff_len[,1]

## Set Genome window size
genomewidows = tileGenome(seqlengths, tilewidth=100000, cut.last.tile.in.chrom=T)# split genome 10 000

############ Genes
## read gff and format names
GbiGffgenes0<-read.delim("AK5_intergenic_transcripts.gtf", header = FALSE, sep="\t", stringsAsFactors = FALSE)
GbiGffgenes_2<-GbiGffgenes0[GbiGffgenes0$V3=="transcript",c(1,4,5,9)]
GbiGffgenes_2$name0<-sapply(strsplit(GbiGffgenes_2$V9, ";"), "[[", 2) 
GbiGffgenes_2$name<-sapply(strsplit(GbiGffgenes_2$name0, " "), "[[", 3)

GbiGffgenes<-GbiGffgenes_2[,c(1,2,3,6)]
#rename columns
colnames(GbiGffgenes)<-c("scaff", "start" ,"end", "name")


gff3genes<-import.gff3("AK5_intergenic_transcripts.gtf")
#select genes
gff3genes<-subset(gff3genes, type =="transcript" )## only genes with mRNA annotation
gff3genes


seqlengths<-Scaff_len[,2]
names(seqlengths)<-Scaff_len[,1]
seqlengths

# count genes inside genome windows
Overlapsgene<-findOverlaps(genomewidows, gff3genes)
UniqOverlapsgene<-Overlapsgene[!duplicated(subjectHits(Overlapsgene))]

# calcualte gene density within windows
genedensity=genomewidows
genedensity$totgenes<-0
genedensity[unique(queryHits(UniqOverlapsgene))]$totgenes<-table(queryHits(UniqOverlapsgene))
genedensity_df<-as.data.frame(genedensity)

### corrected follow bed genome file
Scaff_lenLOTS <- read.csv('AK5_chrom.txt',sep='\t',header = FALSE, col.names = c("name", "scafstart", "scafend"))

genedensity_df$scafstart <- Scaff_lenLOTS[match(genedensity_df$seqnames,Scaff_lenLOTS$name),2]

genedensity_df$newstart <- genedensity_df$start + genedensity_df$scafstart
genedensity_df$newend <- genedensity_df$end + genedensity_df$scafstart
ngenedens <- genedensity_df[,c(1,8,9,4,5,6)]
colnames(ngenedens) <- c('seqnames','start','end','width','strand','totgenes')

############ CpGoe
## Get CpGoe per CDS
CpGoe_Gbi<-read.csv("AK5_dn_CpGoe.tsv",sep='\t')
CpGoe_Gbi_list<-CpGoe_Gbi[,c("ID","CpGoe")]
head(CpGoe_Gbi_list)
dim(CpGoe_Gbi_list)
names(CpGoe_Gbi_list) <- c('name','CpGoe')
CpGoe_Gbi_coordinates<-merge(GbiGffgenes, CpGoe_Gbi_list, all.y=TRUE, by.x="name", by.y="name")
dim(CpGoe_Gbi_coordinates)
head(CpGoe_Gbi_coordinates)

Scaff_lenLOTS <- read.csv('AK5_chrom.txt',sep='\t', header = FALSE, col.names = c("name", "scafstart", "scafend"))
CpGoe_Gbi_coordinates$scafstart <- Scaff_lenLOTS[match(CpGoe_Gbi_coordinates$scaff,Scaff_lenLOTS$name),2]
CpGoe_Gbi_coordinates$newstart <-CpGoe_Gbi_coordinates$start + CpGoe_Gbi_coordinates$scafstart
CpGoe_Gbi_coordinates$newend <- CpGoe_Gbi_coordinates$end + CpGoe_Gbi_coordinates$scafstart
nCpGoe <- CpGoe_Gbi_coordinates[,c(1,2,7,8,5)]
colnames(nCpGoe) <- c('name','scaff','start','end','CpGoe')


############ TEs
gff3ALL<-import.gff3("AK5_te_annotation.gff")
gff3ALL
#select genes
gff3Rep<-subset(gff3ALL, source =="reasonaTE" )

# calculate overlaps withe genome windows
OverlapsRep<-findOverlaps(genomewidows, gff3Rep)
UniqOverlapsRep<-OverlapsRep[!duplicated(subjectHits(OverlapsRep))]

RepDensity=genomewidows
RepDensity$RepCounts<-0
RepDensity[unique(queryHits(UniqOverlapsRep))]$RepCounts<-table(queryHits(UniqOverlapsRep))
RepDensity_df<-as.data.frame(RepDensity)

head(RepDensity_df[order(RepDensity_df$RepCounts, decreasing = TRUE),])

RepMeandensity_chr<-aggregate(RepDensity_df$RepCounts, by=list(Category=RepDensity_df$seqnames), FUN=mean)
head(RepMeandensity_chr[order(RepMeandensity_chr$x, decreasing = TRUE),])

### corrected follow bed genome file
Scaff_lenLOTS <- read.csv('AK5_chrom.txt',sep='\t', header = FALSE,  col.names = c("name", "scafstart", "scafend"))
RepDensity_df$scafstart <- Scaff_lenLOTS[match(RepDensity_df$seqnames,Scaff_lenLOTS$name),2]
RepDensity_df$newstart <- RepDensity_df$start + RepDensity_df$scafstart
RepDensity_df$newend <- RepDensity_df$end + RepDensity_df$scafstart
nrepdens <- RepDensity_df[,c(1,8,9,4,5,6)]
colnames(nrepdens) <- c('seqnames','start','end','width','strand','RepCounts')

Scaff_lenLOTS <- read.csv('AK5_chrom.txt',sep='\t', header = FALSE, col.names = c("name", "scafstart", "scafend"))


library(RColorBrewer)
Selectedcolors<- c("orange","springgreen4")



##########################
## Circos Plot
###########################

pdf("Circos_Plot_AK52.pdf", bg = "white",height = 6, width = 6)

# Initialize Circus
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0, gap.after=0)
circos.genomicInitialize(Scaff_lenLOTS, major.by=5000000, tickLabelsStartFromZero=FALSE, labels.cex=0.45,plotType =  c("axis", "labels"), sector.names = c("2L", "2R", "3L", "3R", "4", "Mito", "X", "Y")) 

# Gene density track
genedensity_df_region<-ngenedens[ngenedens$seqnames %in% Scaff_lenLOTS$name,c(1,2,3,6)]

blue = brewer.pal(8, name ="Paired")
circos.track(ylim = c(0, 1), 
             bg.col = blue, 
             bg.border = NA, track.height = 0.05)

circos.genomicTrackPlotRegion(genedensity_df_region , numeric.column = 4,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "l", col =Selectedcolors[2], ...)
                              }, track.height = 0.2,bg.lty=0,bg.lwd=0.00001, bg.border = NA)


#Repeat overlap with transcripts
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
All_merged$index <- 1:nrow(All_merged)
All_merged = All_merged[,c(ncol(All_merged),1:(ncol(All_merged)-1))]
All_merged$seq_id = All_merged$chrom
All_merged = All_merged %>% 
  separate("chrom", c("Id", "Type", "population", "x"), sep = "_")

All_merged[which(All_merged$Id == "Noncoding"), "Type"] <- "Noncoding"

All_merged[which(All_merged$Type == "End"), "Type"] <- "Downstream"


#Get the overlap
df = All_merged
rm(AK5, DK5, GI5, SW5, UM, YE, Zamb, All_merged) #rm unecessary

df[which(df$strand_transcript == "."), "strand_transcript"] <- "+"
df$same_strand = paste(df$strand_transcript ,df$name,sep="_")
df[which(df$same_strand == "-_+"), "overlap"] <- 0
df[which(df$same_strand == "+_-"), "overlap"] <- 0


df = df %>% group_by(df$seq_id) %>% mutate(total_overlap = sum(overlap)) #Get the overlap
df$rel_overlap = df$total_overlap/(df$transcript_end-df$transcript_start) #Get the relativeoverlap

#Make these categories for final plot
df$seq_id2 = paste(df$seq_id, df$strand_transcript)
merged = distinct(df, seq_id2, .keep_all = TRUE)

merged[which(0 < merged$rel_overlap & merged$rel_overlap  <= 0.2), "overlap_te"] <- "20 % and less overlap"
merged[which(0.2 < merged$rel_overlap & merged$rel_overlap  <= 0.4), "overlap_te"] <- "20% - 40% overlap"
merged[which(0.4 < merged$rel_overlap & merged$rel_overlap  <= 0.6), "overlap_te"] <- "40% - 60% overlap"
merged[which(0.6 < merged$rel_overlap & merged$rel_overlap  <= 0.8), "overlap_te"] <- "60% - 80% overlap"
merged[which(0.8 < merged$rel_overlap & merged$rel_overlap  < 1), "overlap_te"] <- "80% and above overlap"
merged[which(merged$rel_overlap == 0), "overlap_te"] <- "0 % overlap"
merged[which(merged$rel_overlap == 1), "overlap_te"] <- "Inside TE"

merged[which(merged$Type == "Noncoding"), "group"] <- "Noncoding"
merged[which(merged$Type == "Transcript"), "group"] <- "De novo"
merged[which(merged$Type == "Upstream"), "group"] <- "De novo"
merged[which(merged$Type == "End"), "group"] <- "De novo"

merged <- ungroup(merged)

transcripts = subset(merged, Type == "Transcript")
sub2 = subset(transcripts, population== "AK5")
sub2$name = sub2$Id


TE_overlap = merge(nCpGoe, sub2, by = "name")

TE_100 <-TE_overlap[TE_overlap$overlap_te == "Inside TE", 2:5]
TE_80p <-TE_overlap[TE_overlap$overlap_te == "80% and above overlap", 2:5]
TE_60 <-TE_overlap[TE_overlap$overlap_te == "60% - 80% overlap", 2:5]
TE_40 <-TE_overlap[TE_overlap$overlap_te == "40% - 60% overlap", 2:5]
TE_20 <-TE_overlap[TE_overlap$overlap_te == "20% - 40% overlap", 2:5]
TE_20m <-TE_overlap[TE_overlap$overlap_te == "20 % and less overlap", 2:5]
TE_0 <-TE_overlap[TE_overlap$overlap_te == "0 % overlap", 2:5]

circos.genomicRainfall(list(TE_0, TE_20m, TE_20, TE_40, TE_60, TE_80p, TE_100), pch =16,stack = FALSE,
                       cex =0.5, col = c("darkgrey", "#FEFF89", "yellow2", "goldenrod2", "darkorange1", "red2", "darkred"), track.height = 0.2,
                       bg.lty=0,bg.lwd=0.00001, bg.border = NA)




# rep density
RepDensity_df_region<-nrepdens[nrepdens$seqnames %in%Scaff_lenLOTS$name,c(1,2,3,6)]


circos.genomicTrackPlotRegion(RepDensity_df_region , numeric.column = 4,
                              panel.fun = function(region, value, ...) {
                                #i = getI(...)
                                circos.genomicLines(region, value, type = "l", col =Selectedcolors[1], ...)
                              } , track.height = 0.1,bg.lty=0,bg.lwd=0.00001, bg.border = NA)



dev.off() 



