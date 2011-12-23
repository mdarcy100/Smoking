library(impute)
library(gdata); # to read excel data

setwd("/Users/mdarcy100/Desktop/MTroester/")   
source("ToolsHelp/Functions.R")
source("AgePaper/AgeManuscript_062011/age_paper_scripts2.R")
source("AgePaper/AgeManuscript_062011/pcaFuncs.R")

#read in file
fulldata<-read.delim("SmokerAnalysisOct2011/UMAOct2011.txt",sep="\t",header=T,fill=T)
dim(fulldata)

##########################################################################
# process and filter data
##########################################################################
#filter by missingness
#	1) rows where >= 90% of cells are missing data. this will filter out platform differences
# 	2) columns where >=30% of cells missing data  not as great as i would like but ok...)
#   3) rows where >= 20% of cells missing data.
filterdata<-basicfilter(fulldata,.9,.3,.2)
dim(filterdata)
#31643   161

#we need to remove the 'IR'  samples - they are explants and have been Ionized radiated 
#- hopefully there aren't any in this dataset - MD thinks she removed them all
indices<-getIR(filterdata)
length(indices)
filterdata<-filterdata[,-indices]
dim(filterdata)
#31643   160

#find smokers and do a matched analysis for greater power (do this after imputation)

#processing data
all.rm.exprs <- filterdata[2:nrow(filterdata),4:ncol(filterdata)]

probe.names <- stripAGIjunk(as.character(filterdata$CLID)[2:nrow(filterdata)])
rownames(all.rm.exprs) <- probe.names

all.rm.exprs <- data.matrix(all.rm.exprs)
dim(all.rm.exprs)
#31642   157

na.frac.arrays <- apply(all.rm.exprs,2,
                                function(x){sum(is.na(x))/length(x)})

#look for replicates
reps <- getRepsFast(all.rm.exprs)  #in age_paper_scripts2

## Check correlations
unique.reps <-  unique(reps$sampleId.reps)


cor.list.ALL.RMsamples <- vector("list",length=length(unique.reps))
names(cor.list.ALL.RMsamples) <- unique.reps
for(rep.samp in unique.reps){
    rep.samp.data <- vector("list",length=2)
    names(rep.samp.data) <- c("cor","na.frac")
    rep.samp.ind <- grep(paste("^",rep.samp,sep=""),colnames(all.rm.exprs))
    rep.samp.data[["cor"]] <- cor(all.rm.exprs[,rep.samp.ind],use="complete.obs")
    rep.samp.data[["na.frac"]] <- na.frac.arrays[rep.samp.ind]
    cor.list.ALL.RMsamples[[rep.samp]] <- rep.samp.data
}


all.rm.exprs.comp <- na.omit(all.rm.exprs)
PCAtext(all.rm.exprs.comp,colnames(all.rm.exprs.comp),
        cent=TRUE,scle=TRUE,
        type="s",col="red",adj=c(-.2,.5),radius=.5)
        
        
## 07065_rep1 is most likely an outlier - eliminate
## everybody else is ambiguous, keep them in and average

## This will eliminate only 07065_rep1
##pca.07065.rm.ind <- grep("07065_rep[0-9]",colnames(all.rm.exprs))
##all.rm.exprs <- all.rm.exprs[,-pca.07065.rm.ind]

### MD Fall 2011 thinks one of the smokers may be problematic (7187) (based on the PCA)
## Eliminate c("UMA07065_rep1","UMA07039_rep2_244K","UMA07028")
dim(all.rm.exprs)
elim.inds <- which(colnames(all.rm.exprs) %in% c("UMA07065_rep1","UMA07039_rep2_244K","UMA07028"))
all.rm.exprs <- all.rm.exprs[,-elim.inds]
dim(all.rm.exprs)


## Impute data, then average replicates
all.rm.exprs.imputed <- impute.knn(all.rm.exprs,k=5, rng.seed=98765)$data

## Check by PCA
reps.imputed <- getRepsFast(all.rm.exprs.imputed)

unique.reps.imputed <-  unique(reps.imputed$sampleId.reps)

PCAtext(all.rm.exprs.imputed[,reps.imputed$all.reps],colnames(all.rm.exprs.imputed[,reps.imputed$all.reps]),
        cent=TRUE,scle=TRUE,
        type="s",col="red",adj=c(-.2,.5))
        
   
      
iqr.val.rm.imputed <- apply(all.rm.exprs.imputed,1,IQR)

#check correlation of imputed replicates
cor.list.imputed <- vector("list",length=length(unique.reps.imputed))
names(cor.list.imputed) <- unique.reps.imputed
for(rep.samp in unique.reps.imputed){
    rep.samp.ind <- grep(paste("^",rep.samp,sep=""),colnames(all.rm.exprs.imputed))
    cor.list.imputed[[rep.samp]] <- cor(all.rm.exprs.imputed[,rep.samp.ind],use="complete.obs")
}


## Average the replicates - doing the same way JP did this since MDs was ugly
all.rm.avg <- all.rm.exprs.imputed[,-which(colnames(all.rm.exprs.imputed) %in% reps.imputed$all.reps)]

avg.data <- matrix(nrow=nrow(all.rm.avg),ncol=length(unique.reps.imputed))
rownames(avg.data) <- rownames(all.rm.avg)
colnames(avg.data) <- unique.reps.imputed

for(id in unique.reps.imputed){
    avg.samps <- all.rm.exprs.imputed[,grep(paste("^",id,sep=""),colnames(all.rm.exprs.imputed))]
    avg.data[,id] <- apply(avg.samps,1,mean,na.rm=TRUE)
}

all(rownames(avg.data) == rownames(all.rm.avg))

dim(all.rm.avg)
#31642   118

all.rm.avg <- cbind(all.rm.avg,avg.data)
## Change name of UMA07088_rep1 to UMA07088
colnames(all.rm.avg) <- gsub("^UMA07088_rep1$","UMA07088",colnames(all.rm.avg))

dim(all.rm.avg)

##########################################################################
# Build data.frame containing patient data
##########################################################################
patient.data <- read.xls("SmokerAnalysisOct2011/avondata082310_deidentified.xls",sheet="avon")

rownames(patient.data) <- paste("UMA0",patient.data$STUDYID,sep="")

## Determine patients for which there is demographic data
patient.id.match <- rownames(patient.data)[which(rownames(patient.data) %in% colnames(all.rm.avg))]

colnames(all.rm.avg)[!(colnames(all.rm.avg) %in% rownames(patient.data))]
## [1] "UMA07026" "UMA07035" "UMA07051" - there was array information but no demographics

all.rm.match <- all.rm.avg[,patient.id.match]
patient.match <- patient.data[patient.id.match,]

all(colnames(all.rm.match) == rownames(patient.match))

patient.match <- data.frame(patient.match)

all.rm.phenoData <- new("AnnotatedDataFrame", data = patient.match)

all.rm.eset <- new("ExpressionSet", exprs = all.rm.match,
                   phenoData = all.rm.phenoData, annotation = "hgug4112a")

dim(all.rm.eset)
#this should save to the smoker directory
save(all.rm.eset,file="SmokerAnalysisOct2011/all_rm_eset.RData")


