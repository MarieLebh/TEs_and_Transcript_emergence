#Make circular plots (here for AK5)--> For other populations just switch AK5 with the other name

library(rtracklayer)
library("circlize")
library(ggtree)
library("viridis")
library('dplyr')

#Explanation what files used:
# AK5finalGenome.masked.fa.fai = Genome index file
# AK5_intergenic_transcripts.gtf = Gtf file of all de novo transcrips
# AK5_chrom.txt = a bedfile with the chromosome positions
# AK5_te_annotation.gff = A gtf file of the ResonaTE annotation

## Read chrom lengths from file
Scaff_len<-read.delim("AK5finalGenome.masked.fa.fai", sep="\t",stringsAsFactors = FALSE,header = FALSE)
Scaff_len_sorted<-Scaff_len[order(Scaff_len[,2], decreasing = TRUE),]
seqlengths<-Scaff_len[,2]
names(seqlengths)<-Scaff_len[,1]

## Set Genome window size
genomewidows = tileGenome(seqlengths, tilewidth=100000, cut.last.tile.in.chrom=T)# split genome 10 000

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
gff3genes<-subset(gff3genes, type =="transcript" )#only transcripts

seqlengths<-Scaff_len[,2]
names(seqlengths)<-Scaff_len[,1]
seqlengths

# count genes inside genome windows
Overlapsgene<-findOverlaps(genomewidows, gff3genes)
UniqOverlapsgene<-Overlapsgene[!duplicated(subjectHits(Overlapsgene))]

#calcualte gene density within windows
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

write.table(ngenedens,"AK5_gene_density.txt",sep="\t",row.names=FALSE)

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

write.table(ngenedens,"AK5_CpGoe.txt",sep="\t",row.names=FALSE)

############ TEs
gff3ALL<-import.gff3("AK5_te_annotation.gff")
gff3ALL
#select genes
gff3Rep<-subset(gff3ALL, source =="reasonaTE" )
#gff3Rep$Name2<-sapply(strsplit(gff3Rep$Name, "genus:" ), `[`, 2)
#gff3Rep$Class<-sapply(strsplit(gff3Rep$Name2, '/' ), `[`, 1)
#table(gff3Rep$Class)

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

write.table(nrepdens,"AK5_Repeat_Counts.txt",sep="\t",row.names=FALSE)


Scaff_lenLOTS <- read.csv('AK5_chrom.txt',sep='\t', header = FALSE, col.names = c("name", "scafstart", "scafend"))


library(RColorBrewer)
Selectedcolors<- c("orange","springgreen4")

##########################
## Circos Plot
###########################

png("Circos_Plot_AK5.png", bg = "white", res=400, height = 6, width = 6, units = "in")

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

## CpGoe
Gbi_cpg_threshold=0.654
CpGoe_High<-nCpGoe[nCpGoe$CpGoe>Gbi_cpg_threshold ,2:5]
CpGoe_Low<-nCpGoe[nCpGoe$CpGoe<=Gbi_cpg_threshold, 2:5]

circos.genomicRainfall(list(CpGoe_High,CpGoe_Low ), pch =16,stack = FALSE,
                       cex =0.5, col = c("#E69F00", "#56B4E9"), track.height = 0.2,
                       bg.lty=0,bg.lwd=0.00001, bg.border = NA)

# rep density
RepDensity_df_region<-nrepdens[nrepdens$seqnames %in%Scaff_lenLOTS$name,c(1,2,3,6)]


circos.genomicTrackPlotRegion(RepDensity_df_region , numeric.column = 4,
                              panel.fun = function(region, value, ...) {
                                #i = getI(...)
                                circos.genomicLines(region, value, type = "l", col =Selectedcolors[1], ...)
                              } , track.height = 0.1,bg.lty=0,bg.lwd=0.00001, bg.border = NA)

dev.off() 



