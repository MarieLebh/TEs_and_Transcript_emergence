#!/usr/bin/env Rscript

# R version 4.3.0
# Bioconducior 3.17

if (!require(dplyr)) {
    install.packages("dplyr")
    require(dplyr)
}
# DATA transformation

if (!require(boot)) {
    install.packages("boot")
    require(boot)
}
# for bootstrapping

if (!require(car)) {
    install.packages("car")
    require(car)
}
# for getting pvalues of lmer

if (!require(MASS)) {
    install.packages("MASS")
    require(MASS)
}

if (!require(tidyr)) {
    install.packages("tidyr")
    require(tidyr)
}
# Data management

if (!require(glmmTMB)) {
    install.packages("glmmTMB")
    require(glmmTMB)
}
# for glmm
if (!require(MuMIn)) {
    install.packages("MuMIn")
    require(MuMIn)
}
# for glmm model selection
 
if (!require(ggplot2)) {
    install.packages("ggplot2")
    require(ggplot2)
}
# Graphs

if (!require(cowplot)) {
    install.packages("cowplot")
    require(cowplot)
}
# make nice multigraphs

if (!require(emmeans)) {
    install.packages("emmeans")
    require(emmeans)
}
# to compute post-hoc comparisons of glmm within a categorical variable 

if (!require(lares)) {
    install.packages("lares")
    require(lares)
}
# to compute CI for graphs

######################################################################################################################
#######################      Now direct comparison  Transcript vs Homologs     #######################################
######################################################################################################################

gd <- read.csv('Filtered_de_novo_transcripts_homologs_BF.csv') ## Data modified to include all levels of TE classification

CountInOrtho <- gd %>% group_by(Number_Orthogroup) %>% summarize(count = n())
singletons <- as.list(CountInOrtho[CountInOrtho$count == 1,'Number_Orthogroup'])
############################ Data transformation ############################

gdf <- gd[gd$Transcription_status %in% c('Non_transcribed_Homolog','Transcript'),]
gdf$binT <- ifelse(gdf$Transcription_status == 'Transcript',1,0)
gdfg <- gdf[,c(2,3,5,6,9:12,14,15,17,18,20:37)]


## First only TEs from different gene regions (sequence, upstream, downstream)
subSeq <- gdfg[,c(1:6,21,10,19,20,28,15,16,17,18)]
subSeq$TE_ovlp <- subSeq[,6]
subSeq$TE_fam <- subSeq[,7]
subSeq$TE_num <- subSeq[,8]
subSeq$TE_cla <- subSeq[,9]
subSeq$TE_ord <- subSeq[,10]
subSeq$TE_aut <- subSeq[,11]
subSeq$TE_mot <- subSeq[,12]
subSeq$TE_moH <- subSeq[,13]
subSeq$TE_cor <- subSeq[,14]
subSeq$TE_coH <- subSeq[,15]
gSeq <- subSeq[,c(1:5,12:21)]
gSeq$reg <- rep('Sequence',length(gSeq[,1]))
subUps <- gdfg[,c(1:5,7,24,12,22,23,29,15:18)]
subUps$TE_ovlp <- subUps[,6]
subUps$TE_fam <- subUps[,7]
subUps$TE_num <- subUps[,8]
subUps$TE_cla <- subUps[,9]
subUps$TE_ord <- subUps[,10]
subUps$TE_aut <- subUps[,11]
subUps$TE_mot <- subUps[,12]
subUps$TE_moH <- subUps[,13]
subUps$TE_cor <- subUps[,14]
subUps$TE_coH <- subUps[,15]
gUps <- subUps[,c(1:5,12:21)]
gUps$reg <- rep('aUpstream',length(gSeq[,1]))
subDow <- gdfg[,c(1:5,8,27,14,25,26,30,15:18)]
subDow$TE_ovlp <- subDow[,6]
subDow$TE_fam <- subDow[,7]
subDow$TE_num <- subDow[,8]
subDow$TE_cla <- subDow[,9]
subDow$TE_ord <- subDow[,10]
subDow$TE_aut <- subDow[,11]
subDow$TE_mot <- subDow[,12]
subDow$TE_moH <- subDow[,13]
subDow$TE_cor <- subDow[,14]
subDow$TE_coH <- subDow[,15]
gDow <- subDow[,c(1:5,12:21)]
gDow$reg <- rep('Downstream',length(gSeq[,1]))
subTEo <- rbind(gSeq,gUps,gDow)
subTEo$binO <- ifelse(subTEo$TE_ovlp > 0, 1,0) ## data for checking TE only

dataTE <- subTEo[!subTEo$Number_Orthogroup %in% singletons$Number_Orthogroup,]

dataF <- read.csv('Transcript_Feature_Comparison.csv')
dataF$Population <- ifelse(dataF$population == 'FI','AK5',ifelse(dataF$population == 'DK','DK5',ifelse(dataF$population == 'TR','YE',ifelse(dataF$population == 'UA','UM',ifelse(dataF$population == 'ZI','Zamb',ifelse(dataF$population == 'ES','GI5','SW5'))))))
dataF$ID <- paste(dataF$Transcript_ID, dataF$Population, sep = '::')
TPMl <- dataF$TPM
names(TPMl) <- dataF$ID

dataTE$TPM <- ifelse(dataTE$Transcription_status == 'Transcript',TPMl[dataTE$Transcript_Id],0)
#dataTE$rTPM <- (dataTE$TPM * (nrow(dataTE) - 1) + 0.5) / nrow(dataTE)

dataM <- dataTE[!is.na(dataTE$Number_motifs),] ### remove all seq / entry with NA instead of motif number DATA for checking
dataMo <- dataM[!is.na(dataM$Number_core),] ### remove all seq / entry with NA instead of motif number DATA for checking
dataMot <- dataMo[dataMo$reg == 'aUpstream',] ### 

daTEMot <- dataMo[dataMo$reg == 'Sequence',]
daTEMotU <- dataMo[dataMo$reg == 'aUpstream',]

daTEMot$TE_fam2 <- ifelse(is.na(daTEMot$TE_fam),'None',daTEMot$TE_fam)
lstOrthoK <- ifelse(daTEMot$TE_cla == 'None' & daTEMot$Transcription_status == 'Transcript',daTEMot$Number_Orthogroup,NA)
orthoRem <- na.omit(lstOrthoK)
OnDaTE <- daTEMot[!daTEMot$Number_Orthogroup %in% orthoRem,]
OnTE <- daTEMot[!daTEMot$TE_cla == 'None',]
OnTE$TE_cla2 <- ifelse(OnTE$TE_cla == 'RNA','aRNA',OnTE$TE_cla)
NoTE <- daTEMot[daTEMot$TE_cla == 'None',]
NoTE2 <- daTEMot[daTEMot$Number_Orthogroup %in% orthoRem,]

TEUp <- cbind(daTEMotU$TE_ovlp,daTEMotU$TE_num,daTEMotU$binO)
colnames(TEUp) <- c('TE_ovlpU','TE_numU','binOu') 
daTEMSU <- cbind(daTEMot,TEUp)

####################################################################################
####################################################################################
########## models general TEs: Ovlp, Num, Bin Transcripts vs. Homologs #############
####################################################################################
####################################################################################
 
modTEou <- glmmTMB(as.factor(Transcription_status) ~ TE_cla + reg:TE_cla + TE_ovlp + TE_num + binO + TE_ovlp:reg + TE_num:reg + binO:reg + (1|Number_Orthogroup) + (1|Population), data=dataTE, family=binomial)

modTEoua <- glmmTMB((log10(TPM+1)+0.001) ~ TE_ovlp + TE_num + binO + TE_ovlp:reg + TE_num:reg + binO:reg + (1|Number_Orthogroup) + (1|Population),
	zi = ~ TE_ovlp + TE_num + binO + TE_ovlp:reg + TE_num:reg + binO:reg + (1|Number_Orthogroup) + (1|Population), data=dataTE, family=ziGamma)

options(na.action = "na.fail")
msTE1 <- dredge(modTEou, extra = alist(AIC,BIC,ICOMP,Cp))
model.sel(msTE1)
model.avg(msTE1, subset = delta < 4)
confset.95pTE <- get.models(msTE1, cumsum(weight) <= .95)
avgmod.95pTE <- model.avg(confset.95pTE)
summary(avgmod.95pTE)
confint(avgmod.95pTE)
## Best model:
summary(get.models(msTE1, 1)[[1]])
Anova(get.models(msTE1, 1)[[1]],type='III')
saveRDS(get.models(msTE1, 1)[[1]], file = "GLMM_binomial_TEgeneral_TvH.rda")


####### Graphs ########
GraON <- dataTE %>% group_by(Transcription_status) %>% summarise(mean_ovlp = mean(TE_ovlp), mean_num = mean(TE_num), ssdovlp = sd(TE_ovlp), ssdnum = sd(TE_num), count = n())  %>%  mutate(upper_ci_ovlp = ci_upper(mean_ovlp, ssdovlp, count, conf = 0.95), upper_ci_num = ci_upper(mean_num, ssdnum, count, conf = 0.95)) %>% mutate(ci_ovlp = upper_ci_ovlp - mean_ovlp,ci_num = upper_ci_ovlp - mean_num)

FigO <- ggplot(GraON,aes(x=Transcription_status,y=mean_ovlp)) + geom_bar(aes(fill =Transcription_status),stat='identity') + 
		geom_errorbar(aes(ymin=mean_ovlp-ci_ovlp, ymax=mean_ovlp+ci_ovlp,color= Transcription_status), width=.2) + 
		theme(legend.position='none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),axis.text.x = element_text(size=16,angle=90),
		axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(y= "Relative overlap with TEs")
FigNb <- ggplot(GraON,aes(x=Transcription_status,y=mean_num)) +  geom_bar(aes(fill =Transcription_status),stat='identity') + 
		geom_errorbar(aes(ymin=mean_num-ci_num, ymax=mean_num+ci_num,color= Transcription_status), width=.2)  + 
		theme(legend.position ='none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Transcription status", color = "Transcription status",y= "Number of TEs overlaping")

GraONr <- dataTE %>% group_by(Transcription_status,reg) %>% summarise(mean_ovlp = mean(TE_ovlp), mean_num = mean(TE_num), ssdovlp = sd(TE_ovlp), ssdnum = sd(TE_num), count = n())  %>%  mutate(upper_ci_ovlp = ci_upper(mean_ovlp, ssdovlp, count, conf = 0.95), upper_ci_num = ci_upper(mean_num, ssdnum, count, conf = 0.95)) %>% mutate(ci_ovlp = upper_ci_ovlp - mean_ovlp,ci_num = upper_ci_ovlp - mean_num)

figOr <- ggplot(GraONr,aes(x=Transcription_status,y=mean_ovlp)) + geom_point(aes(color = Transcription_status, shape =Transcription_status)) +
		facet_wrap(.~reg) + geom_errorbar(aes(ymin=mean_ovlp-ci_ovlp, ymax=mean_ovlp+ci_ovlp,color= Transcription_status), width=.2) + 
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),axis.text.x = element_text(size=16,angle=90),
		axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Transcription status", color = "Transcription status",y= "Relative overlap with TEs")

figNr <- ggplot(GraONr,aes(x=Transcription_status,y=mean_num)) + geom_point(aes(color = Transcription_status, shape =Transcription_status)) +  facet_wrap(.~reg) + geom_errorbar(aes(ymin=mean_num-ci_num, ymax=mean_num+ci_num,color= Transcription_status), width=.2)  + 
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Transcription status", color = "Transcription status",y= "Number of TEs overlaping")

GraBr <- dataTE %>% group_by(Transcription_status,reg) %>% summarise(sumB = sum(binO), count = n(), propB = sumB / count)
figBr <- ggplot(dataTE,aes(x=Transcription_status)) + geom_bar(aes(fill=as.factor(binO)),position='fill') + 
		facet_wrap(.~reg,nrow=1) + 
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Transcription status", color = "Transcription status",y= "Presence / Absence of TEs")

Graph_TEg <- plot_grid(FigNb,FigO,figOr,  labels = c('A', 'B','C'), nrow=1, label_size = 24)

ggsave('TEs_relovlp_num_TvH.pdf',Graph_TEg,width=800,height=400,units="mm",dpi=300)
ggsave('TEs_PresAbs_general_TvH.pdf',figBr,width=800,height=400,units="mm",dpi=300)

########   Models to account for  motifs & co   ################

options(na.action = "na.omit")

modFu1 <- glmmTMB(as.factor(Transcription_status) ~ Number_motifs +
	Number_motifs_high +
	Number_core +
	Number_core_high +
	Number_motifs:Number_core+
	Number_motifs:Number_core_high+
	Number_motifs_high:Number_core+
	Number_motifs_high:Number_core_high+
	(1|Number_Orthogroup)+
	(1|Population),
	data=dataMot,family=binomial)

modFu1a <- glmmTMB((log10(TPM+1)+0.001) ~ Number_motifs +
	Number_motifs_high +
	Number_core +
	Number_core_high +
	(1|Number_Orthogroup)+
	(1|Population),
	zi = ~ Number_motifs + Number_motifs_high + Number_core + Number_core_high + (1|Number_Orthogroup) + (1|Population),
	data=dataMot,family=ziGamma)

options(na.action = "na.fail")
msMot1 <- dredge(modFu1, extra = alist(AIC,BIC,ICOMP,Cp))
model.sel(msMot1)
model.avg(msMot1, subset = delta < 4)
confset.95pMot <- get.models(msMot1, cumsum(weight) <= .95)
avgmod.95pMot <- model.avg(confset.95pMot)
summary(avgmod.95pMot)
confint(avgmod.95pMot)
## Best model:
Anova(get.models(msMot1, 1)[[1]],type='III')
summary(get.models(msMot1, 1)[[1]])
saveRDS(get.models(msMot1, 1)[[1]], file = "GLMM_binomial_Motifs_TvH.rda")

#### Best model modFu1
saveRDS(modFu1, file = "GLMM_binomial_TEMotifs_TvH.rda")




## Graphs

figTEMo <- ggplot(Fgd) + geom_violin(aes(y = Number_motifs, x = cut_interval(TE_ovlp,4), color = Transcription_status)) + 
		theme(legend.position="none",
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Relative overlap with TEs (bins)", color = "Transcription status",y= "Number of Motifs")

figTECo <- ggplot(dataMot) + geom_violin(aes(y = Number_core, x = Transcription_status, color = Transcription_status)) + 
		theme(legend.position="none",
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Relative overlap with TEs (bins)", color = "Transcription status",y= "Number of cores")

figTEcMo <- ggplot(dataMot) + geom_boxplot(aes(y = Number_core, x = cut_interval(TE_num,4), color = Transcription_status)) + 
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Log2 of Number of TEs overlaping (bins)", color = "Transcription status",y= "Number of Motifs")
		
Graph_TEMo <- plot_grid(figTEMo,figTECo, figTEcMo,  labels = c('A', 'B', 'C'), nrow=1, label_size = 24)

ggsave('TEs_Motifs_TvH.pdf',Graph_TEMo,width=800,height=400,units="mm",dpi=300)


########   Models to account for TEs and motifs & co   ################

options(na.action = "na.omit")

modTEMot1 <- glmmTMB(as.factor(Transcription_status)~ TE_ovlp + 
	TE_num + 
	Number_motifs+
	Number_motifs_high+
	Number_core+
	Number_core_high+
	TE_ovlp:Number_motifs +
	TE_ovlp:Number_motifs_high+
	TE_ovlp:Number_core +
	TE_ovlp:Number_core_high +
	TE_num:Number_motifs +
	TE_num:Number_motifs_high +
	TE_num:Number_core +
	TE_num:Number_core_high +
	binO:Number_motifs +
	binO:Number_motifs_high +
	binO:Number_core +
	binO:Number_core_high +
	(1|Number_Orthogroup)+
	(1|Population),
	data=daTEMot,family=binomial)

modTEMot2 <- glmmTMB(as.factor(Transcription_status)~ TE_ovlp + 
	TE_num + 
	Number_core+
	binO:Number_motifs +
	(1|Number_Orthogroup)+
	(1|Population),
	data=daTEMot,family=binomial)

options(na.action = "na.fail")
msTEMot1 <- dredge(modTEMot1 , extra = alist(AIC,BIC,ICOMP,Cp))
model.sel(msTEMot1)
model.avg(msTEMot1, subset = delta < 4)
confset.95pTEMot <- get.models(msTEMot1, cumsum(weight) <= .95)
avgmod.95pTEMot <- model.avg(confset.95pTEMot)
summary(avgmod.95pTEMot)
confint(avgmod.95pTEMot)
## Best model:
summary(get.models(msTEMot1, 1)[[1]])

########## Graphs ###############

graphdataMT <- daTEMot %>% group_by(Transcription_status,binO) %>% summarise(mMot = mean(Number_motifs), mCoH = mean(Number_core_high), mCor = mean(Number_core), ssdMot = sd(Number_motifs), ssdCoH = sd(Number_core_high), ssdCor = sd(Number_core), count = n())  %>%  mutate(upper_ci_Mot = ci_upper(mMot, ssdMot, count, conf = 0.95), upper_ci_CoH = ci_upper(mCoH, ssdCoH, count, conf = 0.95), upper_ci_Cor = ci_upper(mCor, ssdCor, count, conf = 0.95),prop = count/sum(count)) %>% mutate(ci_Mot = upper_ci_Mot - mMot,ci_CoH = upper_ci_CoH - mCoH,ci_Cor = upper_ci_Cor - mCor)

figMT <- ggplot(graphdataMT,aes(x=binO,y=mMot)) + 
		geom_point(aes(color = Transcription_status, shape =Transcription_status),position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin=mMot-ci_Mot, ymax=mMot+ci_Mot,color= Transcription_status),position=position_dodge(width=0.5), width=.2)+
		theme(legend.position = 'none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Presence / absence TE insertion", color = "Transcription status",y= "Number of Motifs")

ggsave('TE_Motifs_PA_TvH.pdf',figMT,width=800,height=400,units="mm",dpi=300)

############################    Models Now accounting for type of TEs   ############################

modFCu <- glmmTMB(as.factor(Transcription_status)~ TE_cla + Number_core + TE_ovlp:TE_cla + TE_num:TE_cla + # Number_core:TE_aut +
	Number_core:TE_cla + # Number_core_high:TE_aut +
	(1|Number_Orthogroup)+
	(1|Population),
	data=OnDaTE,family=binomial)
modFCu <- glmmTMB(as.factor(Transcription_status)~ TE_cla + TE_cla:TE_num:Number_core  +# Number_core:TE_aut +
	Number_core:TE_cla + # Number_core_high:TE_aut +
	(1|Number_Orthogroup)+
	(1|Population),
	data=OnTE,family=binomial)

### best model modFCu only upstream regions !!!!!
saveRDS(modFCu, file = "GLMM_binomial_TEclassMotif_TvH.rda")

modFCuR <- glmmTMB(as.factor(Transcription_status)~ TE_cla +
	Number_motifs:TE_cla:reg +
	Number_core:TE_cla +
	Number_core_high + 
	(1|Number_Orthogroup)+
	(1|Population),
	data=nFgd,family=binomial)

saveRDS(modFCuR, file = "GLMM_binomial_TEclassMotif_UpstrSeq_TvH.rda")


########## Graphs ###############

graphdataC <- OnTE %>% group_by(Transcription_status,TE_cla) %>% summarise(mMot = mean(Number_motifs), mCoH = mean(Number_core_high), mCor = mean(Number_core), ssdMot = sd(Number_motifs), ssdCoH = sd(Number_core_high), ssdCor = sd(Number_core), count = n())  %>%  mutate(upper_ci_Mot = ci_upper(mMot, ssdMot, count, conf = 0.95), upper_ci_CoH = ci_upper(mCoH, ssdCoH, count, conf = 0.95), upper_ci_Cor = ci_upper(mCor, ssdCor, count, conf = 0.95),prop = count/sum(count)) %>% mutate(ci_Mot = upper_ci_Mot - mMot,ci_CoH = upper_ci_CoH - mCoH,ci_Cor = upper_ci_Cor - mCor)

graphdataCa <-daTEMot %>% group_by(TE_fam2,Transcription_status) %>% summarise(mMot = mean(Number_motifs), mCoH = mean(Number_core_high), mCor = mean(Number_core), ssdMot = sd(Number_motifs), ssdCoH = sd(Number_core_high), ssdCor = sd(Number_core), count = n())  %>%  mutate(upper_ci_Mot = ci_upper(mMot, ssdMot, count, conf = 0.95), upper_ci_CoH = ci_upper(mCoH, ssdCoH, count, conf = 0.95), upper_ci_Cor = ci_upper(mCor, ssdCor, count, conf = 0.95),prop = count/sum(count)) %>% mutate(ci_Mot = upper_ci_Mot - mMot,ci_CoH = upper_ci_CoH - mCoH,ci_Cor = upper_ci_Cor - mCor)

graphdataCR <- OnTE %>% group_by(Transcription_status,TE_cla,reg) %>% summarise(mMot = mean(Number_motifs), mCoH = mean(Number_core_high), mCor = mean(Number_core), ssdMot = sd(Number_motifs), ssdCoH = sd(Number_core_high), ssdCor = sd(Number_core), count = n())  %>%  mutate(upper_ci_Mot = ci_upper(mMot, ssdMot, count, conf = 0.95), upper_ci_CoH = ci_upper(mCoH, ssdCoH, count, conf = 0.95), upper_ci_Cor = ci_upper(mCor, ssdCor, count, conf = 0.95)) %>% mutate(ci_Mot = upper_ci_Mot - mMot,ci_CoH = upper_ci_CoH - mCoH,ci_Cor = upper_ci_Cor - mCor)

figprop1 <- ggplot(graphdataCa,aes(x=TE_fam2,y=prop)) + 
		geom_bar(aes(fill = Transcription_status),position='stack',stat='identity') +
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "TE class (DNA vs. RNA)", fill = "Transcription status",y= "Proportion of sequences with TE overlap")
figprop2 <- ggplot(graphdataC,aes(x=Transcription_status,y=prop)) + 
		geom_bar(aes(fill = TE_cla),position='stack',stat='identity') +
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_blank(),axis.text.y = element_text(size=16)) + 
		labs(x = "Transcription status", color = "TE class")
figprop2a <- ggplot(graphdataC,aes(x=TE_cla,y=prop)) + 
		geom_bar(aes(fill = Transcription_status),position='stack',stat='identity') +
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_blank(),axis.text.y = element_text(size=16)) + 
		labs(x = "TE class", color = "Transcription status")

Grap_TEprop <- plot_grid(figprop1, figprop2, labels = c('A', 'B'), nrow=1, label_size = 24)

figovlpC <- ggplot(graphdataCR,aes(x=TE_cla,y=mMot)) + 
		geom_point(aes(color = Transcription_status, shape =Transcription_status),position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin=mMot-ci_Mot, ymax=mMot+ci_Mot,color= Transcription_status),position=position_dodge(width=0.5), width=.2) + 
		facet_wrap(.~reg) + 
		theme(legend.position = 'none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "TE class (DNA vs. RNA)", color = "Transcription status",y= "Number of Motifs")

figovlpC3 <- ggplot(graphdataC,aes(x=TE_cla,y=mCor)) + 
		geom_point(aes(color = Transcription_status, shape =Transcription_status),position=position_dodge(width=0.5)) +
		geom_errorbar(aes(ymin=mCor-ci_Cor, ymax=mCor+ci_Cor,color= Transcription_status),position=position_dodge(width=0.5), width=.2) + 
		theme(legend.text = element_text(size = 20),legend.title = element_text(size = 24),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "TE class (DNA vs. RNA)", color = "Transcription status",y= "Number of cores") 

Grap_TEclMo <- plot_grid(figovlpC, figovlpC3, labels = c('A', 'B'), nrow=1, label_size = 24)

ggsave('TE_clas_Motifs_TvH.pdf',figovlpC3,width=800,height=400,units="mm",dpi=300)
ggsave('TE_clas_prop_TvH.pdf',figprop2,width=800,height=400,units="mm",dpi=300)

