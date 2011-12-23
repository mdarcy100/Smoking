#################################################################################
##
## Load required libraries
##
#################################################################################
 ##########################################################################################
 #
 #     function to get full p-value can apply to rows of data
 #	    depends on the smokefullstrata and bmi_cat need to be the same length as x which is
 #      a row (gene)
 
 #
 ##########################################################################################
pextractor= function(x){
	#print("about to run lm")
#	x
	#dim(x)
	#smokefullstrata
	#dim(smokefullstrata)
	#dim(bmi_cat)
	#bmi_cat
	tmp<-lm(x~factor(age_cat)+factor(smokefullstrata)+factor(bmi_cat))
	#print("just ran lm")
   	sum<-summary(tmp)	
	coef<-as.matrix(sum$coefficients)
	smoke1p = coef[2,4] # this is the first factor of smoking
	print(smoke1p)
	smoke2p= coef[3,4] # this is the second factor of smoking
	print(smoke2p)
   	fstat<-sum$fstatistic
   	pvalue=pf(fstat[1],fstat[2],fstat[3],,lower.tail=FALSE)
	print(pvalue)
	#return(pvalue)
	ret<-c(pvalue,smoke1p,smoke2p)
	return(ret)
}
 ##########################################################################################
f.smoker<-factor(smoke.eset$smoking,levels=c(0,1,2))
tmp<-lm(exprs(smoke.eset)[1,]~f.smoker)
#(Intercept)    f.smoker1    f.smoker2  
#    -0.1367      -0.3759      -0.2096  
#-0.09155      0.21342      0.07528 
rownames(smoke.eset)=="A_23_P365412"
tmp<-lm(exprs(smoke.eset)[8891,]~factor(smoke.eset$smoking))
tmp<-lm(exprs(smoke.eset)[8891,]~f.smoker) 
smoke.limma.filtered
#$coefficients
#             (Intercept) factor(smoke.eset$smoking)1 factor(smoke.eset$smoking)2
#A_23_P100001 -0.13670130                  -0.3758908                 -0.20957648
#A_23_P100011 -0.09155455                   0.2134230                  0.07527677
#A_23_P100022  2.83404560                   0.4442044                 -0.60033749
#A_23_P100056  0.03856421                   0.5720700                 -0.28048088
#A_23_P100074  0.01301797                   0.2127978                 -0.40437908



#MD: needed to install biobase prior to using the library Biobase
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biobase")
#########

library(Biobase)
library(hgug4112a.db)
library(gdata)
library(qvalue)
library(survival)
library(gplots)
library(limma)
library(RColorBrewer)
library(made4)
###MD: needed to also install mclust, flexmix, modeltools,multcomp,mvtnorm for fpc 
library(fpc)
library(impute)

setwd("/Users/mdarcy100/Desktop/MTroester/")   
source("ToolsHelp/Functions.R")
source("AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("AgePaper/AgeManuscript_062011/pcaFuncs.R")

##########################################################################
# Analysis with all data (data with bmi, age, smoking status)
##########################################################################
#read in data
#all.rm.eset
load("SmokerAnalysisOct2011/all_rm_eset.RData")
dim(all.rm.eset)
#Features  Samples 
#   31642      133 

#filter out those without entrezids
## Basic prefiltering
## Reduce probes to those having an ENTREZ Id
entrezIds <- mget(featureNames(all.rm.eset), envir=hgug4112aENTREZID)
haveEntrezIds <- names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(all.rm.eset)) - length(haveEntrezIds)
smoke.eset <- all.rm.eset[haveEntrezIds,]

dim(smoke.eset)
#Features  Samples 
#   24066      133 

## Determine the interquartile range for each gene - will be used to filter results
## after applying the empirical Bayes step in Limma
## With limma

iqr.smoke.eset <- esApply(smoke.eset,1,IQR)

# create regression analysis model without bmi/with bmi/predicting 1,2,3 as fctors
# there are only 11 smokers with bmi and age information
# Regression on smoke.smoking (0= never, 1= past, 2 = current )
f.smoker<-factor(smoke.eset$smoking,levels=c(0,1,2))
design.smoker<-model.matrix(~0+f.smoker)
colnames(design.smoker) <- c("never","past","current")

#################################################################################
#       adjusted for bmi/age; not sure if they are truly confounding factors.  
#       (will also do a matched analysis)
####### confounding factors of bmi and age 

bmi_cat<-cut(smoke.eset$bmi,breaks=c(0,30,100),right=FALSE) #previously had been <25, 25-30, >30 (how many significant with15%)589 116
f.bmi<-factor(bmi_cat)
age_cat<-cut(smoke.eset$age,breaks=c(0,35,50,100),right=FALSE)
f.age<-factor(age_cat)
###################################################################################
design.smoke.full <- model.matrix(~ f.smoker+f.age+f.bmi)
colnames(design.smoke.full) <- c("Intercept","past_smoker","current_smoker","middle_age","old_age","bmi_high")
#smokefull.exprs<-exprs(smoke.eset)[,as.numeric(rownames(design.smoke.full))]
smokefull.eset<-smoke.eset[,as.numeric(rownames(design.smoke.full))]

smokefull.limma.fit1 <- lmFit(smokefull.eset,design.smoke.full)
smokefull.limma.ebayes <- eBayes(smokefull.limma.fit1)
#filter by variability (do top 25% -> not 50% - may change to top 25%)
smokefull.limma.filtered <- smokefull.limma.ebayes[iqr.smoke.eset>median(iqr.smoke.eset),]
dim(smoke.limma.filtered)

currentsmoke.limma.summary <- topTable(smokefull.limma.filtered,coef="current_smoker",adjust.method="fdr",num=Inf)
pastsmoke.limma.summary <- topTable(smokefull.limma.filtered,coef="past_smoker",adjust.method="fdr",num=Inf)
bmi_high.limma.sumamry<-topTable(smokefull.limma.filtered,coef="bmi_high",adjust.method="fdr",num=Inf)
middle_age.limma.summary<-topTable(smokefull.limma.filtered,coef="middle_age",adjust.method="fdr",num=Inf)
old_age.limma.sumamry<-topTable(smokefull.limma.filtered,coef="old_age",adjust.method="fdr",num=Inf)
smoke.limma.summary <- topTable(smokefull.limma.filtered,adjust.method="fdr",num=Inf)

## Distribution of p-values
pval.summary(currentsmoke.limma.summary$P.Value)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#  1335   2567   4276   5210   5967   6901   7533   8013   8435  12026 
pval.summary(currentsmoke.limma.summary$adj.P.Val)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#   208   1404   3216   4130   4999   5954   6606   7156   7627  12026 


#####only current smoking status seems to have an effect - very strange
pval.summary(bmi_high.limma.sumamry$adj.P.Val)

pval.summary(old_age.limma.sumamry$adj.P.Val)

pval.summary(middle_age.limma.summary$adj.P.Val)

pval.summary(pastsmoke.limma.summary$adj.P.Val)

currentsmokefulladjust.fold <- currentsmoke.limma.summary[abs(currentsmoke.limma.summary$logFC)>=1.5
                                    & currentsmoke.limma.summary$adj.P.Val  < 0.025,]

                                
dim(currentsmokefulladjust.fold)
# [1] 993   7                                  
currentsmokefulladjust.fold.probe <- currentsmokefulladjust.fold$ID
currentsmokefulladjust.fold.entrez <- unlist(mget(currentsmokefulladjust.fold$ID,env=hgug4112aENTREZID))
currentsmokefulladjust.fold.symbol <- unlist(mget(currentsmokefulladjust.fold$ID,env=hgug4112aSYMBOL))

sig_fold_identifier<-paste(currentsmokefulladjust.fold.entrez,currentsmokefulladjust.fold.symbol,sep="|")

exprs.currentsmokefulladjust.fold.clust <- exprs(smokefull.eset)[currentsmokefulladjust.fold$ID,]
dim(exprs.currentsmokefulladjust.fold.clust)
#[1] 993 77

#write out to a file
colnames(exprs.currentsmokefulladjust.fold.clust)<-paste(colnames(exprs.currentsmokefulladjust.fold.clust),"Smoke:",smokefull.eset$smoking, "Age:",smokefull.eset$age,"BMI:",smokefull.eset$bmi)
rownames(exprs.currentsmokefulladjust.fold.clust)<-sig_fold_identifier;
#clusterFoldFile<-cbind(sig_fold_identifier,as.matrix(exprs.currentsmoke.fold.clust))
write.table(exprs.currentsmokefulladjust.fold.clust,file=paste("SmokeSigFold_FullyAdjusted_", Sys.Date(),".txt",sep=""),sep="\t",col.names=NA)
#redo analysis adjusted for age groups (< 30, 30-39, >40) and bmi > 30, bmi < 30



#################################################################################
#################################################################################
#unadjusted smoking regression below
#################################################################################
#the code below gives similar results to the regression (plain linear regression)
design.smoke.smoking <- model.matrix(~ factor(smoke.eset$smoking))
colnames(design.smoke.smoking) <- c("Intercept","past_smoker","current_smoker")

smoke.limma.fit1 <- lmFit(smoke.eset,design.smoke.smoking)
smoke.limma.ebayes <- eBayes(smoke.limma.fit1)
#filter by variability (do top 25% -> not 50% - may change to top 25%)
smoke.limma.filtered <- smoke.limma.ebayes[iqr.smoke.eset>median(iqr.smoke.eset),]
dim(smoke.limma.filtered)
#12026     3

#pvalue of interest is located here 
smoke.limma.filtered[1,]$p.value[3]
smoke.limma.filtered[1,]$F
smoke.limma.filtered[1,]$F.p.value
smoke.limma.filtered[i,3]$p.value[3]

#look at the contrast matrix:  compare smokers to non-smokers, past smokers to non-smokers and current smokers to past smokers
#contrast.matrix<-makeContrasts(design.smoke.smoking$current_smoker~design.smoke.smoking$Intercept,past_smoker$design.smoke.smoking~design.smoke.smoking#$Intercept,design.smoke.smoking$current_smoker~design.smoke.smoking$past_smoker)
#smoke.limma.fit2<-contrasts.fit(smoke.limma.filtered,contrast.matrix)
#smoke.ebayes.fit2<-ebayes(smoke.limma.fit2)

#look at adjusted p-values, q-values for both past and current smokers
#currentsmoke.limma.summary <- topTable(smoke.limma.filtered,coef="factor.smoke.eset.smoking.2",adjust.method="fdr",num=Inf)
currentsmoke.limma.summary <- topTable(smoke.limma.filtered,coef="current_smoker",adjust.method="fdr",num=Inf)
pastsmoke.limma.summary <- topTable(smoke.limma.filtered,coef="past_smoker",adjust.method="fdr",num=Inf)
smoke.limma.summary <- topTable(smoke.limma.filtered,adjust.method="fdr",num=Inf)

#pastsmoke.limma.summary <- topTable(smoke.limma.filtered,coef="factor.smoke.eset.smoking.1",adjust.method="fdr",num=Inf)


## Distribution of p-values
pval.summary(currentsmoke.limma.summary$P.Value)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#   672   1566   3036   3884   4727   5717   6390   6950   7433  12026 
pval.summary(currentsmoke.limma.summary$adj.P.Val)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#    30    474   1773   2560   3315   4289   4985   5596   6080  12026 

#take all those where < 0.025

q.value.currentsmoker <- qvalue(p=currentsmoke.limma.summary$P.Value)
summary(q.value.currentsmoker)
#Cumulative number of significant calls:

#        <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#p-value    672   1566  3036   3884  4727 5717 12026
#q-value     96    830  2529   3528  4623 6005 12026


pval.summary(pastsmoke.limma.summary$P.Value)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     9    211   1800   3113   4358   5694   6531   7196   7733  12026 

pval.summary(pastsmoke.limma.summary$adj.P.Val)
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1  <0.15   <0.2  <0.25     <1 
#     0      0      0      0    649   3202   4620   5514   6218  12026 
q.value.pastsmoker <- qvalue(p=pastsmoke.limma.summary$P.Value)
summary(q.value.pastsmoker)
#        <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#p-value      9    211  1800   3113  4358 5694 12026
#q-value      0      0     0   1913  4369 6544 12026

##########################################################################################
# get significant genes and those with > fold change
##########################################################################################
currentsmoke.fold <- currentsmoke.limma.summary[abs(currentsmoke.limma.summary$logFC)>=1.0
                                    & currentsmoke.limma.summary$adj.P.Val  < 0.025,]
                                
dim(currentsmoke.fold)
# [1] 605   7                                   
currentsmoke.fold.probe <- currentsmoke.fold$ID
currentsmoke.fold.entrez <- unlist(mget(currentsmoke.fold$ID,env=hgug4112aENTREZID))
currentsmoke.fold.symbol <- unlist(mget(currentsmoke.fold$ID,env=hgug4112aSYMBOL))

sig_fold_identifier<-paste(currentsmoke.fold.entrez,currentsmoke.fold.symbol,sep="|")

exprs.currentsmoke.fold.clust <- exprs(smoke.eset)[currentsmoke.fold$ID,]
dim(exprs.currentsmoke.fold.clust)
#[1] 605 133

#write out to a file
colnames(exprs.currentsmoke.fold.clust)<-paste("Smoke:",smoke.eset$smoking)
rownames(exprs.currentsmoke.fold.clust)<-sig_fold_identifier;
clusterFoldFile<-cbind(sig_fold_identifier,as.matrix(exprs.currentsmoke.fold.clust))
write.table(exprs.currentsmoke.fold.clust,file=paste("ClusterableSmokeSigFold_", Sys.Date(),".txt",sep=""),sep="\t",col.names=NA)
#redo analysis adjusted for age groups (< 30, 30-39, >40) and bmi > 30, bmi < 30


summary(classifyTestsF(smoke.limma.filtered,p.value=0.0001))
summary(classifyTestsF(smoke.limma.filtered,p.value=0.00001))
summary(classifyTestsF(smoke.limma.filtered,p.value=0.000001))

results000001<-classifyTestsF(smoke.limma.filtered,p.value=0.000001)
vennDiagram(results000001,include="up")
vennDiagram(results000001,include="down")

#make pairwise comparisons between the two groups
results <- decideTests(smoke.limma.filtered,p.value=0.0001)
vennDiagram(results,"up")
vennDiagram(results,"down")

selected<- p.adjust(smoke.limma.filtered$F.p.value, method="fdr") < 0.001

up_results000001 <- results000001[results000001$Intercept == 0 & results000001$past_smoker == 0 & results000001$current_smoker == 1 &,]
                                    
sig.age.decade.fold.probe <- sig.age.decade.fold$ID
sig.age.decade.fold.entrez <- unlist(mget(sig.age.decade.fold$ID,env=hgug4112aENTREZID))
sig.age.decade.fold.symbol <- unlist(mget(sig.age.decade.fold$ID,env=hgug4112aSYMBOL))

up_results000001 <- 
results <- decideTests(smoke.limma.filtered,method="global",p=0.001)
summary(results)
heatDiagram(results,smoke.limma.filtered$coef,primary="Difference")

##########################################################################
# Analysis with matched data (matched on age, bmi)
##########################################################################