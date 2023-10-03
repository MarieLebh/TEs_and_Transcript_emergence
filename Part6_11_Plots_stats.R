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

##########################################################
####################        TEs       ####################
##########################################################

dataTEa <- read.csv('Data_rel_TE_overlap_new.csv')
head(dataTEa)
tapply(dataTEa$rel_overlap,list(dataTEa$Type,dataTEa$Population),mean)
dataTEa$P_A <- ifelse(dataTEa$rel_overlap > 0,'Presence','Absence')
dataTEa$IdPop <- paste0(dataTEa$Id,':',dataTEa$Population)

lTEt <- unique(dataTEa[dataTEa$rel_overlap > 0.80 & dataTEa$Type %in% c("Transcript","Intergenic"),'IdPop'])

dTE <- dataTEa[!dataTEa$IdPop %in% lTEt ,] ### Need to decide how to filter data to remove TEs (rel_overlap > 0.8) everywhere? or just for transcript & intergenic as now? 

modTFzg_fd <- glmmTMB(rel_overlap ~ as.factor(Type) +
	(1|Population),
	zi = ~ as.factor(Type) + (1|Population),
	data=dataTEa,family=ziGamma)


modTFzg <- glmmTMB(rel_overlap ~ as.factor(Type) +
	(1|Population),
	zi = ~ as.factor(Type) + (1|Population),
	data=dTE,family=ziGamma)

emmTF = emmeans(modTFzg,specs = pairwise ~ Type,adjust='bonferroni',type='response')


## This modTFzg is the best model in terms of significance against other similar models and in terms of family distribution zero-inflated gamma
saveRDS(modTFzg, file = "GLMM_ziGamma_TE_TvA.rda")

modTNzg <- glmmTMB(rel_overlap ~ 1 +
	(1|Population),
	zi = ~ as.factor(Type),
	data=dTE,family=ziGamma)
	
	
	
########### Graph TEs #############

figTEovlp <- ggplot(dTE,aes(x = Type)) + 
		scale_y_continuous(name = "Proportion of TE overlaping (bars)", sec.axis = sec_axis(~.*1, name="Relative overleap with TEs (violins)")) + 
		#geom_boxplot(data=dTE[dTE$P_A=='Presence',],aes(x=rel_overlap, y = Type))
		geom_violin(data=dTE[dTE$P_A=='Presence',], aes(y=rel_overlap, x = Type)) +
  		geom_point(data=dTE[dTE$P_A=='Presence',], aes(y=rel_overlap, x = Type), position = "jitter",size=0.1) + 
  		geom_bar(aes(fill=P_A),position='fill',alpha=0.5) +
  		#geom_boxplot(width=0.1, color="grey", alpha=0.2) +
		coord_flip() + 
		theme(legend.text = element_text(size = 17),legend.title = element_text(size = 18),
		strip.text.x = element_text(size = 24),axis.text.x = element_text(size=14,angle=90),
		axis.title.x = element_text(size=24),axis.title.y = element_text(size=20),axis.text.y = element_text(size=16))

figTEovlp

ggsave('TE_overlap_TvA.pdf',figTEovlp,width=800,height=400,units="mm",dpi=300)

############ bootstrap GLMM #####################

ownTE <- function(formula, zif, data, indices) {
  subdat <- data %>% group_by(data$Type) %>% slice_sample(n=5000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, zi = zif, data=d,family=ziGamma)
  return(summary(fit)$coeff$cond[,4])
}

bootTE <- boot(data=dTE, statistic = ownTE, formula = rel_overlap ~ as.factor(Type) + (1|Population), zif =  ~ as.factor(Type) + (1|Population), R = 999)

resBTE <- summary(bootTE)
rownames(resBTE) <- names(bootTE$t0)
resBTE$CB <- rep("Ga",4)

ownTEzi <- function(formula, zif, data, indices) {
  subdat <- data %>% group_by(data$Type) %>% slice_sample(n=5000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, zi = zif, data=d,family=ziGamma)
  return(summary(fit)$coeff$zi[,4])
}

bootTEzi <- boot(data=dTE, statistic = ownTEzi, formula = rel_overlap ~ as.factor(Type) + (1|Population), zif =  ~ as.factor(Type) + (1|Population), R = 999)

resBTEzi <- summary(bootTEzi)
rownames(resBTEzi) <- names(bootTEzi$t0)
resBTEzi$CB <- rep("Zi",4)
resTEb <- rbind(resBTE,resBTEzi)
write.csv(resTEb,'results_bootpval_TE_TvA_modelZIGamma.csv')		

ownTEpc <- function(formula, zif, data, indices) {
  subdat <- data %>% group_by(data$Type) %>% slice_sample(n=5000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, zi = zif, data=d,family=ziGamma)
  emm <- emmeans(fit,specs = pairwise ~ Type,adjust='bonferroni',type='response')
  memm <- as.data.frame(emm$contrasts)
  lstEM <- memm$p.value
  names(lstEM) <- memm$contrast
  return(lstEM)
}

bootTEpc <- boot(data=dTE, statistic = ownTEpc, formula = rel_overlap ~ as.factor(Type) + (1|Population), zif =  ~ as.factor(Type) + (1|Population), R = 999)

resBTEpc <- summary(bootTEpc)
rownames(resBTEpc) <- names(bootTEpc$t0)
write.csv(resBTEpc,'results_bootpval_TE_TvA_modelZIGamma_posthoc_comp.csv')		


##########################################################
####################      Motifs      ####################
##########################################################

dataM <- read.csv('Data_motifs_new.csv')
head(dataM)

modMI <- glmmTMB(number_motifs_per_transcript ~ type_te + (1|population), data=dataM,family=poisson)

summary(modMI)

saveRDS(modMI, file = "GLMM_poisson_Motifs_TvA.rda")

emmMI = emmeans(modMI,specs = pairwise ~ type_te,adjust='bonferroni',type='response')

modMIh <- glmmTMB(number_motifs_high_per_transcript ~ type_te + (1|population), data=dataM,family=poisson)

emmMIh = emmeans(modMIh,specs = pairwise ~ type_te,adjust='bonferroni',type='response')

saveRDS(modMIh, file = "GLMM_poisson_Motifs_high_TvA.rda")

############## bootstrap #################

own <- function(formula, data, indices) {
  subdat <- data %>% group_by(data$type) %>% slice_sample(n=1000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, data=d,family=poisson)
  return(summary(fit)$coeff$cond[,4])
}

trial <- boot(data=dataM, statistic = own, formula = number_motifs_per_transcript ~ type_te + (1|population), R = 999)
trialh <- boot(data=dataM, statistic = own, formula = number_motifs_high_per_transcript ~ type_te + (1|population), R = 999)

#boot.ci(trial, type='bca', index = 20)

resMot <- summary(trial)
rownames(resMot) <- names(trial$t0)
write.csv(resMot,'results_bootpval_Motifs_TvA_modelpoisson.csv')		

resMoth <- summary(trialh)
rownames(resMoth) <- names(trialh$t0)
write.csv(resMoth,'results_bootpval_Motifs_high_TvA_modelpoisson.csv')		

ownpc <- function(formula, data, indices) {
  subdat <- data %>% group_by(data$type) %>% slice_sample(n=1000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, data=d,family=poisson)
  emm <- emmeans(fit,specs = pairwise ~ type_te,adjust='bonferroni',type='response')
  memm <- as.data.frame(emm$contrasts)
  lstEM <- memm$p.value
  names(lstEM) <- memm$contrast
  return(lstEM)
}

trialpc <- boot(data=dataM, statistic = ownpc, formula = number_motifs_per_transcript ~ type_te + (1|population), R = 999)
trialpch <- boot(data=dataM, statistic = ownpc, formula = number_motifs_high_per_transcript ~ type_te + (1|population), R = 999)

#boot.ci(trial, type='bca', index = 20)

resMotpc <- summary(trialpc)
rownames(resMotpc) <- names(trialpc$t0)
write.csv(resMotpc,'results_bootpval_Motifs_TvA_modelpoisson_posthoc_comp.csv')		

resMotpch <- summary(trialpch)
rownames(resMotpch) <- names(trialpch$t0)
write.csv(resMotpch,'results_bootpval_Motifs_high_TvA_modelpoisson_posthoc_comp.csv')		

##### Graph #####

GmotD <- dataM %>% group_by(type_te) %>% summarise(mMot = mean(number_motifs_per_transcript), mMoH = mean(number_motifs_high_per_transcript), ssdMot = sd(number_motifs_per_transcript), ssdMoH = sd(number_motifs_high_per_transcript), count = n())  %>%  mutate(upper_ci_Mot = ci_upper(mMot, ssdMot, count, conf = 0.95), upper_ci_MoH = ci_upper(mMoH, ssdMoH, count, conf = 0.95)) %>% mutate(ci_Mot = upper_ci_Mot - mMot,ci_MoH = upper_ci_MoH - mMoH)

FigMot <- ggplot(GmotD,aes(x=type_te,y=mMot)) + geom_point(aes(color = type_te)) + 
		geom_errorbar(aes(ymin=mMot-ci_Mot, ymax=mMot+ci_Mot,color= type_te), width=.2) + 
		theme(legend.position='none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=14,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Number of Motifs")
 
FigMoH <- ggplot(GmotD,aes(x=type_te,y=mMoH)) + geom_point(aes(color = type_te)) + 
		geom_errorbar(aes(ymin=mMoH-ci_MoH, ymax=mMoH+ci_MoH,color= type_te), width=.2) + 
		theme(legend.text = element_text(size = 17),legend.title = element_text(size = 18),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=14,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Number of Motifs high")
 
FigMo <- plot_grid(FigMot,FigMoH,labels=c("A","B"), nrow= 1, label_size = 24)
ggsave('Motifs_SeqTyp_TvA.pdf',FigMo,width=800,height=400,units="mm",dpi=300)
ggsave('Motifs_SeqTyp_TvA_a.pdf',FigMot,width=800,height=400,units="mm",dpi=300)
ggsave('Motifs_SeqTyp_TvA_b.pdf',FigMoH,width=800,height=400,units="mm",dpi=300)

##########################################################
####################      Cores      ####################
##########################################################

dataC <- read.csv('Data_core_new.csv')
head(dataC)

modC <- glmmTMB(number_core_per_transcript ~ type_te + (1|population), data=dataC,family=poisson)

summary(modC)

saveRDS(modC, file = "GLMM_poisson_cores_TvA.rda")

emmC = emmeans(modC,specs = pairwise ~ type_te,adjust='bonferroni',type='response')

modCh <- glmmTMB(number_core_high_per_transcript ~ type_te + (1|population), data=dataC,family=poisson)

emmCh = emmeans(modCh,specs = pairwise ~ type_te,adjust='bonferroni',type='response')

saveRDS(modCh, file = "GLMM_poisson_cores_high_TvA.rda")

############## bootstrap #################

own <- function(formula, data, indices) {
  subdat <- data %>% group_by(data$type) %>% slice_sample(n=1000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, data=d,family=poisson)
  return(summary(fit)$coeff$cond[,4])
}

trialc <- boot(data=dataC, statistic = own, formula = number_core_per_transcript ~ type_te + (1|population), R = 999)
trialch <- boot(data=dataC, statistic = own, formula = number_core_high_per_transcript ~ type_te + (1|population), R = 999)

#boot.ci(trial, type='bca', index = 20)

resCor <- summary(trialc)
rownames(resCor) <- names(trialc$t0)
write.csv(resCor,'results_bootpval_cores_TvA_modelpoisson.csv')		

resCorh <- summary(trialch)
rownames(resCorh) <- names(trialch$t0)
write.csv(resCorh,'results_bootpval_cores_high_TvA_modelpoisson.csv')		

ownpc <- function(formula, data, indices) {
  subdat <- data %>% group_by(data$type) %>% slice_sample(n=1000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, data=d,family=poisson)
  emm <- emmeans(fit,specs = pairwise ~ type_te,adjust='bonferroni',type='response')
  memm <- as.data.frame(emm$contrasts)
  lstEM <- memm$p.value
  names(lstEM) <- memm$contrast
  return(lstEM)
}

trialcpc <- boot(data=dataC, statistic = ownpc, formula = number_core_per_transcript ~ type_te + (1|population), R = 999)
trialcpch <- boot(data=dataC, statistic = ownpc, formula = number_core_high_per_transcript ~ type_te + (1|population), R = 999)

#boot.ci(trial, type='bca', index = 20)

resCorpc <- summary(trialcpc)
rownames(resCorpc) <- names(trialcpc$t0)
write.csv(resCorpc,'results_bootpval_cores_TvA_modelpoisson_posthoc_comp.csv')		

resCorpch <- summary(trialcpch)
rownames(resCorpch) <- names(trialcpch$t0)
write.csv(resCorpch,'results_bootpval_cores_high_TvA_modelpoisson_posthoc_comp.csv')		

##### Graph #####

GcorD <- dataC %>% group_by(type_te) %>% summarise(mCor = mean(number_core_per_transcript), mCoH = mean(number_core_high_per_transcript), ssdCor = sd(number_core_per_transcript), ssdCoH = sd(number_core_high_per_transcript), count = n())  %>%  mutate(upper_ci_Cor = ci_upper(mCor, ssdCor, count, conf = 0.95), upper_ci_CoH = ci_upper(mCoH, ssdCoH, count, conf = 0.95)) %>% mutate(ci_Cor = upper_ci_Cor - mCor,ci_CoH = upper_ci_CoH - mCoH)

FigCor <- ggplot(GcorD,aes(x=type_te,y=mCor)) + geom_point(aes(color = type_te)) + 
		geom_errorbar(aes(ymin=mCor-ci_Cor, ymax=mCor+ci_Cor,color= type_te), width=.2) + 
		theme(legend.position='none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=14,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Number of cores")
 
FigCoH <- ggplot(GcorD,aes(x=type_te,y=mCoH)) + geom_point(aes(color = type_te)) + 
		geom_errorbar(aes(ymin=mCoH-ci_CoH, ymax=mCoH+ci_CoH,color= type_te), width=.2) + 
		theme(legend.text = element_text(size = 17),legend.title = element_text(size = 18),
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=14,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Number of cores high")
 
FigCo <- plot_grid(FigCor,FigCoH,labels=c("A","B"), nrow= 1, label_size = 24)
ggsave('Cores_SeqTyp_TvA.pdf',FigCo,width=800,height=400,units="mm",dpi=300)
ggsave('Cores_SeqTyp_TvA_a.pdf',FigCor,width=800,height=400,units="mm",dpi=300)
ggsave('Cores_SeqTyp_TvA_b.pdf',FigCoH,width=800,height=400,units="mm",dpi=300)


###################################################################
####################      Genome Features      ####################
###################################################################

dataF <- read.csv('Transcript_Feature_Comparison.csv')
head(dataF)

dataF$ID <- paste(dataF$Transcript_ID, dataF$population, sep = '::')
TPMl <- dataF$TPM
names(TPMl) <- dataF$ID

modBinAalt <- glmmTMB(as.factor(type)~log10(TPM)+
	gc_content+
	spliced_length+
	exon_number+
	(1|chromosome)+
	(1|population),
	data=dataF,family=binomial)
## best model

saveRDS(modBinAalt, file = "GLMM_binomial_GenomFeat_TvA.rda")

################ bootstrap ###############

ownF <- function(formula, data, indices) {
  subdat <- data %>% group_by(data$type) %>% slice_sample(n=8000)
  d <- subdat[indices,] # allows boot to select sample
  fit <- glmmTMB(formula, data=d,family=binomial)
  return(summary(fit)$coeff$cond[,4])
}

trialF <- boot(data=dataF, statistic = ownF, formula = as.factor(type)~log10(TPM)+
	gc_content+
	spliced_length*
	exon_number+
	(1|chromosome)+
	(1|population), R = 999)

resGF <- summary(trialF)
rownames(resGF) <- names(trialF$t0)
write.csv(resGF,'results_bootpval_GenomFeat_TvA_modelbinomial.csv')		

############## Graphs ##################

fTPM <- ggplot(dataF) + geom_boxplot(aes(x= type,y = log10(TPM), color = type))  + 
		theme(legend.position = 'none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Log(TPM)")
fGCc <- ggplot(dataF) + geom_boxplot(aes(x= type,y = gc_content, color = type)) + 
		theme(legend.position = 'none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "GC content")
fsLEn <- ggplot(dataF) + geom_boxplot(aes(y = spliced_length, x = type, color = type)) + 
		theme(legend.position = 'none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Spliced length")

fsEXn <- ggplot(dataF) + geom_boxplot(aes(y = exon_number, x = type, color = type)) + 
		theme(legend.position = 'none',
		strip.text.x = element_text(size = 20),axis.title.x = element_text(size=24),
		axis.text.x = element_text(size=16,angle=90),axis.title.y = element_text(size=24),axis.text.y = element_text(size=16)) + 
		labs(x = "Type", color = "Type",y= "Exon number")

GraGenFt <- plot_grid(fTPM,fGCc,fsLEn,fsEXn,labels=c("A","B","C","D"), nrow= 2, label_size = 24)

ggsave("GenomFeat_SeqTyp_TvA.pdf",GraGenFt,width=800,height=400,units="mm",dpi=300)


