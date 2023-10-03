library(ggplot2)
library(lme4)
library(cowplot)
library(tidyr)
library(dplyr)
library(car)

lPop <- c('AK5','DK5','SW5','GI5','UM','YE','Zamb')
df <- data.frame()
for(i in 1:7){
	filename <- paste(lPop[i],"Repeat_Counts.txt", sep="_")
	fd <- read.csv(filename,sep='\t',header=T)
	fd$Population <- rep(lPop[i],length(fd[,1]))
	df <- rbind(df,fd)
	} 
	
df2 <- data.frame()
for(i in 1:7){
	filename <- paste(lPop[i],"transcript_density.txt", sep="_")
	fd <- read.csv(filename,sep='\t',header=T)
	fd$Population <- rep(lPop[i],length(fd[,1]))
	df2 <- rbind(df2,fd)
	} 

ndf3 <- cbind(df,df2$totgenes)
colnames(ndf3) <- c(colnames(df),'totgenes')
df3 <- ndf3#[!ndf3$seqnames == 'mitochondrion_Chromosome',]
df3$chr_pop <- paste(df3$seqnames,df3$Population,sep="::")
sumchr <- tapply(df3$width,df3$chr_pop, sum)
df3$distcent <- ifelse(df3$start < round(sumchr[df3$chr_pop],-5)/2, ((round(sumchr[df3$chr_pop],-5)/2)-df3$start), (df3$start- (round(sumchr[df3$chr_pop],-5)/2)))

df3$densG <- df3$totgenes/df3$width
df3$densT <- df3$RepCounts/df3$width
df3$denst <- scale(df3$densT)

df3$distance <- scale(df3$distcent)
df3$repcount <- scale(df3$RepCounts)

md <- glmer(totgenes ~ repcount + (1|seqnames) + (1|Population),family=poisson,data = df3, control = glmerControl(optimizer ='bobyqa'))

mdD <- glmer(totgenes ~ distance + (1|seqnames) + (1|Population),family=poisson,data = df3, control = glmerControl(optimizer ='bobyqa'))

mdT <- glmer(RepCounts ~ distance + (1|seqnames) + (1|Population),family=poisson,data = df3, control = glmerControl(optimizer ='bobyqa'))

GG <- ggplot(df3,aes(x=totgenes,y=RepCounts, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm)+ facet_wrap(.~ seqnames,scales='free',ncol=2) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Number of de novo transcripts per 100 kb", color = "Population",y= "Number of TEs per 100 kb")

GGa <- ggplot(df3,aes(x=totgenes,y=RepCounts, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_blank(),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Number of de novo transcripts per 100 kb", color = "Population",y= "Number of TEs per 100 kb")
GGb <- ggplot(df3,aes(x=totgenes,y=RepCounts, color=seqnames, shape=seqnames,  linetype=seqnames, fill=seqnames)) + geom_point() +
	geom_smooth(method=lm) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Number of de novo transcripts per 100 kb",y= "Number of TEs per 100 kb")
	
G0 <- ggplot(df3,aes(y=totgenes,x=distcent, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) + facet_wrap(.~ seqnames,scales='free_x') +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Distance from center of chromosome",y= "Number of de novo transcripts per 100 kb")
G0a <- ggplot(df3,aes(y=totgenes,x=distcent, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) + 
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_blank(),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Distance from center of chromosome",y= "Number of de novo transcripts per 100 kb")
G0bb <- ggplot(df3,aes(y=totgenes,x=distcent, color=seqnames, shape=seqnames,  linetype=seqnames, fill=seqnames)) + geom_point() +
	geom_smooth(method=lm) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Distance from center of chromosome",y= "Number of de novo transcripts per 100 kb")
	
G1 <- ggplot(df3,aes(y=RepCounts,x=distcent, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) + facet_wrap(.~ seqnames,scales='free') +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Distance from center of chromosome",y= "Number of de novo transcripts per 100 kb")
G1a <- ggplot(df3,aes(y=RepCounts,x=distcent, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) + 
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_blank(),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Distance from center of chromosome",y= "Number of TEs per 100 kb")
G1b <- ggplot(df3,aes(y=RepCounts,x=distcent, color=seqnames, shape=seqnames,  linetype=seqnames, fill=seqnames)) + geom_point() +
	geom_smooth(method=lm) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "Distance from center of chromosome",y= "Number of TEs per 100 kb")
	
GraphGG <- plot_grid(GGa, GGb, labels = c('A', 'B'), nrow=2, label_size = 24)
GraphG0 <- plot_grid(G0a, G0b, labels = c('A', 'B'), nrow=2, label_size = 24)
GraphG1 <- plot_grid(G1a, G1b, labels = c('A', 'B'), nrow=2, label_size = 24)

ggsave('dnG_TE_chrom_scale.pdf',GraphGG,width=800,height=400,units="mm",dpi=300)
ggsave('dnG_TE_chrom_scale_detailed.pdf',GG,width=800,height=400,units="mm",dpi=300)
ggsave('Distribution of de novo genes from chromo center.pdf',GraphG0,width=800,height=400,units="mm",dpi=300)
ggsave('Distribution of de novo genes detailed from chromo center.pdf',G0,width=800,height=400,units="mm",dpi=300)
ggsave('Distribution of TEs from chromo center.pdf',GraphG0,width=800,height=400,units="mm",dpi=300)
ggsave('Distribution of TEs detailed from chromo center.pdf',G0,width=800,height=400,units="mm",dpi=300)

################################
## CpGoe with TE ovlp and TPM ##
################################

dc <- data.frame()
for(i in 1:7){
	filename <- paste(lPop[i],"CpGoe.txt", sep="_")
	fd <- read.csv(filename,sep='\t',header=T)
	fd$Population <- rep(lPop[i],length(fd[,1]))
	dc <- rbind(dc,fd)
	} 
dc$id_pop <- paste(dc$name,dc$Population,sep="::")
dc <- dc[!is.na(dc$CpGoe),]

dt <- read.csv('Data_rel_TE_overlap_new.csv')
dt$id_pop <- paste(dt$Id,dt$Population,sep="::")
dt[match(dc$id_pop,dt$id_pop),'rel_overlap']

dm <- read.csv('Data_motifs_new.csv')
dm$Id <- strsplit(dm$unique_id,"::")[[1]][1]
dm <- dm %>% mutate(Id = stringr::str_remove(unique_id, "::.*"))
dm$id_pop <- paste(dm$Id,dm$population,sep="::")

fd <- read.csv('Transcript_Feature_Comparison.csv')
fd$Population <- ifelse(fd$population == 'FI','AK5',ifelse(fd$population == 'DK','DK5',ifelse(fd$population == 'TR','YE',ifelse(fd$population == 'UA','UM',ifelse(fd$population == 'ZI','Zamb',ifelse(fd$population == 'ES','GI5','SW5'))))))
fd$id_pop <- paste(fd$Transcript_ID,fd$Population,sep="::")

dc$TE_ovlp <- dt[match(dc$id_pop,dt$id_pop),'rel_overlap']
dc$Motifs <- dm[match(dc$id_pop,dm$id_pop),'number_motifs_per_transcript']
dc$Motifs_high <- dm[match(dc$id_pop,dm$id_pop),'number_motifs_high_per_transcript']
dc$TPM <- fd[match(dc$id_pop,fd$id_pop),'TPM']
dc$GC_content <- fd[match(dc$id_pop,fd$id_pop),'gc_content']

dc$Meth <- scale(dc$CpGoe)
dc$tpm <- scale(dc$TPM)

mc <- lmer(CpGoe ~ TE_ovlp + (1|Population) + (1|scaff),data = dc, control = lmerControl(optimizer ="nloptwrap"))

mcm <- lmer(CpGoe ~ Motifs +  + (1|scaff),data = dc, control = lmerControl(optimizer ="nloptwrap"))
mcc <- lmer(CpGoe ~ GC_content +  + (1|scaff),data = dc, control = lmerControl(optimizer ="nloptwrap"))
mcT <- lmer(CpGoe ~ tpm +  + (1|scaff),data = dc, control = lmerControl(optimizer ="nloptwrap"))


GC <- ggplot(dc,aes(y=TE_ovlp,x=CpGoe, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "CpGoe of de novo transcripts",y= "Relative overlap with TEs")
GCx <- ggplot(dc,aes(y=TE_ovlp,x=CpGoe, color=Population, shape=Population,  linetype=Population, fill=Population)) + geom_point() +
	geom_smooth(method=lm) + facet_wrap(.~ Population,scales='free') +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "CpGoe of de novo transcripts",y= "Relative overlap with TEs")
	
GCa <- ggplot(dc,aes(y=TE_ovlp,x=CpGoe, color=scaff, shape=scaff,  linetype=scaff, fill=scaff)) + geom_point() +
	geom_smooth(method=lm)+ facet_wrap(.~ scaff,scales='free') +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "CpGoe of de novo transcripts",y= "Relative overlap with TEs")
	
GCb <- ggplot(dc,aes(y=TE_ovlp,x=CpGoe)) + geom_point() +
	geom_smooth(method=lm) +
	theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
	axis.title.x = element_text(size=24),axis.text.x = element_text(size=16),
	axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
	labs(x = "CpGoe of de novo transcripts",y= "Relative overlap with TEs")
	
GraphGG <- plot_grid(GGa, GGb, labels = c('A', 'B'), nrow=2, label_size = 24)
GraphG0 <- plot_grid(G0a, G0b, labels = c('A', 'B'), nrow=2, label_size = 24)
GraphG1 <- plot_grid(G1a, G1b, labels = c('A', 'B'), nrow=2, label_size = 24)

ggsave('CpGoe_TEovlp_Population.pdf',GC,width=800,height=400,units="mm",dpi=300)
ggsave('CpGoe_TEovlp_Population_detailed.pdf',GCx,width=800,height=400,units="mm",dpi=300)
ggsave('CpGoe_TEovlp_Chromosome.pdf',GCa,width=800,height=400,units="mm",dpi=300)
ggsave('CpGoe_TEovlp_General.pdf',GCb,width=800,height=400,units="mm",dpi=300)

