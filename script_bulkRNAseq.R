library("limma")
library("edgeR")
library("ggplot2")

########################################
##### Functions
########################################

###Get DE genes
getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)
  
  return(genes_de_sorted)
}

###Normalize per gene
normalizePerGene<-function(expMatrix){
  resultMatrix<-t(apply(expMatrix, 1, function(x)(x-min(x))/(max(x)-min(x))))
  
  return(resultMatrix)
}

################################################################################
######### LOAD DATA
################################################################################

getwd()

########################################
##### Load needed count data
########################################

##### Load raw counts exp1
countDataKia1<-read.table(file="/srv/data/RNAseqSamples/exp1/outputHTSeqCount/all_counts.txt",sep="\t", header=TRUE, 
                      stringsAsFactors=FALSE)

##### Load raw counts exp2
countDataKia2<-read.table(file="/srv/data/RNAseqSamples/exp2/outputHTSeqCount/all_counts.txt",sep="\t", header=TRUE, 
                          stringsAsFactors=FALSE)

### Only NMG and NCPM needed
colsNMG<-grep('NMG',colnames(countDataKia2))
colsNCPM<-grep('NCPM',colnames(countDataKia2))
countDataKia2<-countDataKia2[,c(colsNMG, colsNCPM)]

### Remove NMG_3
colsNMG3<-grep('NMG_S3',colnames(countDataKia2))
countDataKia2<-countDataKia2[,- colsNMG3]

##### Load raw counts of GSE117081
countDataCharlie<-read.table(file="/srv/data/RNAseqSamples/GSE117081/outputHTSeqCount/all_counts.txt",sep="\t", header=TRUE, 
                          stringsAsFactors=FALSE)

### Only AMF Cre- and KC Cre- needed (and only 4 replicates of each)
colsAMF<-grep('AMF_Creneg_[1-4]',colnames(countDataCharlie))
colsKC<-grep('KC_Creneg_[1-4]',colnames(countDataCharlie))
countDataCharlie<-countDataCharlie[,c(colsAMF, colsKC)]


########################################
##### Merge needed count data
########################################

dim(countDataKia1)
dim(countDataKia2)
dim(countDataCharlie)

cbind(head(rownames(countDataKia1)),head(rownames(countDataKia2)),head(rownames(countDataCharlie)))
cbind(tail(rownames(countDataKia1)),tail(rownames(countDataKia2)),tail(rownames(countDataCharlie)))

countData<-cbind(countDataKia1,countDataKia2,countDataCharlie)
dim(countData)
##24421    27

### Change colnames
newColnames<-gsub('_S1','_1',colnames(countData))
newColnames<-gsub('_S2','_2',newColnames)
newColnames<-gsub('_S3','_3',newColnames)
newColnames<-gsub('_S4','_4',newColnames)
newColnames<-gsub('_S5','_5',newColnames)

theIDs<-grep('^PM_',newColnames)
newColnames[theIDs]<-substr(newColnames[theIDs],1,nchar(newColnames[theIDs])-4)

cbind(colnames(countData),newColnames)
colnames(countData)<-newColnames


########################################
##### Load meta data
########################################

### Load meta data
colData<-read.table(file="documentation/metadata.txt",sep="\t", header=TRUE, stringsAsFactors=TRUE)
rownames(colData)<-colData$fileName
colData<-colData[,-1]
##remove the bad samples here as well
rowNMG3<-grep('NMG_3',rownames(colData))
rowAMF5<-grep('AMF_Creneg_5',rownames(colData))
rowsKC56<-grep('KC_Creneg_[5-6]',rownames(colData))
colData<-colData[- c(rowNMG3, rowAMF5, rowsKC56),]

### Reorder
countData<-countData[,rownames(colData)]
dim(countData)
# 24421    27

################################################################################
########## CREATE OBJECT
################################################################################

y <- DGEList(counts = countData)

################################################################################
########## FILTER DATA
################################################################################

##### Filter low count genes
## Do filtering
yNoFilter<-y
keep = rowSums(cpm(y)>1) >= 3
y = y[keep,]
dim(y)
##13 471
dim(yNoFilter)
##24 421

##### Reset lib sizes
y$samples$lib.size = colSums(y$counts)
y$samples

################################################################################
########## NORMALISATION
################################################################################

##### Scale normalisation
yNoNormalisation<-y
y <- calcNormFactors(y)

##### MDS-plot
theColors<-c(rep("#619cff",4),rep("#f8766d",4),rep("darkgreen",4),rep("#d39200",3),rep("#a97070",4),
             rep("deeppink",4),rep("darkorchid",4))

plotMDS(y, cex=1.4, col=theColors, pch=16)
plotMDS(y, cex=1.0, col=theColors)


################################################################################
########## LOG2 TRANSFORMATION
################################################################################

#### Create design matrix
TS <- paste(colData$cell, colData$condition, sep=".")
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

##### Do voom
v <- voom(y, design, plot = TRUE)

expTable<-v$E


##### Normalised counts
countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)


################################################################################
########## CHECK FILTERING AND NORMALISATION
################################################################################


#################### BARPLOT ####################

##### Barplot lib sizes raw counts
par(mar = c(9,3,3,1)) #more margin: bottom, left, top, right
bp<-barplot(yNoFilter$samples$lib.size*1e-6,axisnames=FALSE,main="Barplot lib sizes of raw counts",ylab="Library size (millions)")
axis(1, labels=rownames(yNoFilter$samples), at = bp, las=2, cex.axis=0.7)

##### Barplot lib sizes normalised counts
# countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
bp<-barplot(colSums(countData_norm)*1e-6,axisnames=FALSE,main="Barplot lib sizes of normalised counts",ylab="Library size (millions)")
axis(1, labels=colnames(countData_norm), at = bp, las=2, cex.axis=0.7)


#################### BOXPLOT ####################
col <- rainbow(nrow(colData))

y2<-y
y2$samples$norm.factors<-1
y2$counts[,1]<-ceiling(y2$counts[,1]*0.05)
y2$counts[,2]<-y2$counts[,2]*5

par(mfrow=c(1,2), mar = c(12,4,3,1)) #more margin: bottom, left, top, right
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Unnormalised data", ylab="log-cpm", col=col)

y2<-calcNormFactors(y2)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Normalised data", ylab="log-cpm", col=col)



#################### DENSITY PLOT ####################
col <- topo.colors(nrow(colData))

par(mfrow=c(1,2))
### Plot log2-CPM values of each sample before filtering
theCpmNoFilter<-cpm(yNoFilter, log=TRUE)

plot(density(theCpmNoFilter[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="raw data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpmNoFilter)[i]))
  
  den <- density(theCpmNoFilter[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
# legend("topright", rownames(y$samples), text.col=col, cex=0.6, bty="n")

### Plot log2-CPM values of each sample after filtering (and normalisation)
theCpm<-cpm(y, log=TRUE)

plot(density(theCpm[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="filtered data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpm)[i]))
  
  den <- density(theCpm[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
par(mfrow=c(1,1))



#################### HISTOGRAM OF EXPTABLE ####################

par(mfrow=c(1,2))
### Histogram
hist(expTable)

### Density plot
plot(density(expTable[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="Density of expTable", xlab="log2")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0("Add line for sample ",colnames(expTable)[i]))
  
  den <- density(expTable[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
# legend("topright", rownames(y$samples), text.col=col, cex=0.6, bty="n")
par(mfrow=c(1,1))

################################################################################
########## PCA
################################################################################
library("rgl")

### Calculate variance
variance<-apply(expTable, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.15*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-expTable[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-theColors

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab="pc1", ylab="pc2", zlab="pc3")
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-2), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1.3)



### Save 3D
dirToSave<-paste(getwd(), "/results/",sep="")
writeWebGL(dir = dirToSave, filename = file.path(dirToSave, "pca_noLabels.html"),
           template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
           prefix = "",
           snapshot = TRUE, font = "Arial")

### Take snapshot
rgl.snapshot("results/pca_noLabels.png")


################################################################################
########## CORRELATION HEATMAP SAMPLES
################################################################################
library("RColorBrewer")
library("gplots")

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

##heatmap 1: based on distance
distsRL <- dist(t(expTable),method="euclidean")
hc <- hclust(distsRL,method="ward.D")

heatmap.2(as.matrix(distsRL),
          Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=rev(hmcol),margin=c(13, 13), cexRow=0.9, cexCol=0.9)


##heatmap 2: based on correlation
cm=cor(expTable)

distsRL <- as.dist(1-cor(expTable))
hc <- hclust(distsRL,method="ward.D")


heatmap.2(cm, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=hmcol, margin=c(13, 13), cexRow=0.9, cexCol=0.9)



#just dendrogram
library('ape')
distsRL <- as.dist(1-cor(expTable))
hc <- hclust(distsRL,method="ward.D")

plot(as.phylo(hc), cex = 0.9, label.offset = 0.01)


##heatmap 2: based on correlation (with values)
library('dplyr')
cm=cor(expTable[,hc$order])
head(cm)
library(reshape2)
melted_cormat <- melt(cm)
head(melted_cormat)

### First make as.numeric of the factors => these will be the x- and y-coordinates
sample_plotdata <- melted_cormat %>% 
  mutate(x=as.numeric(Var1), y=as.numeric(Var2))%>% 
  mutate(Var1_cluster = gsub("(.*)_.*", "\\1", Var1), Var2_cluster = gsub("(.*)_.*", "\\1", Var2))
### Then group by Var1_cluster and Var2_cluster to do this: mean of x-coordinates, mean of y-coordinates and the mean of the value
cluster_plotdata <- sample_plotdata %>% 
  group_by(Var1_cluster, Var2_cluster) %>% 
  summarise(x=mean(x), y=mean(y), value=mean(value))

ggplot(sample_plotdata, aes(Var1, Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradientn(colors=hmcol) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text = element_text(color = "black", size = 12), 
        axis.title = element_blank(), legend.text=element_text(size=11), legend.title = element_blank()) +
  geom_text(aes(x, y, label=round(value, 2), color=value < 0.9), data=cluster_plotdata) +
  scale_color_manual(values=c(`TRUE`="grey10", `FALSE`="grey98"), guide=FALSE)

################################################################################
########## GET DE GENES
################################################################################

########################################
##### 1. MF.MhcHiDura vs others
########################################

### Create design matrix
TS <- c(rep("hiDura",4),rep("others",23))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=hiDura-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup1<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup1<-getDEgenes(allGenesGroup1,0.05,1)
dim(DEgenesGroup1)
##35

########################################
##### 2. MF.MhcIntDura vs others
########################################

### Create design matrix
TS <- c(rep("others",4), rep("intDura",4),rep("others",19))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=intDura-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup2<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup2<-getDEgenes(allGenesGroup2,0.05,1)
dim(DEgenesGroup2)
##166

########################################
##### 4. MF.brain vs others
########################################

### Create design matrix
TS <- c(rep("others",12), rep("brain",3),rep("others",12))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=brain-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup4<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup4<-getDEgenes(allGenesGroup4,0.05,1)
dim(DEgenesGroup4)
##1741

########################################
##### 5. MF.MhcHiCp vs others
########################################

### Create design matrix
TS <- c(rep("others",15), rep("HiCp",4),rep("others",8))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=HiCp-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup5<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup5<-getDEgenes(allGenesGroup5,0.05,1)
dim(DEgenesGroup5)
##596

########################################
##### 8. (MhcHiDura+MhcIntDura+MhcHiCp
#####       +brain) vs others
########################################

### Create design matrix
TS <- c(rep("wanted",8),rep("others",4),rep("wanted",7),rep("others",8))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=wanted-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup8<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup8<-getDEgenes(allGenesGroup8,0.05,1)
dim(DEgenesGroup8)
##4383

upGenesGroup8<-DEgenesGroup8[DEgenesGroup8$logFC > 1,]
tmp<-upGenesGroup8[upGenesGroup8$AveExpr > 0,]
upGenesGroup8<-tmp
dim(upGenesGroup8)
##1208

########################################
##### 9. (MhcHiDura+MhcIntDura+MhcHiCp)
#####       vs others
########################################

### Create design matrix
TS <- c(rep("wanted",8),rep("others",7),rep("wanted",4),rep("others",8))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=wanted-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup9<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup9<-getDEgenes(allGenesGroup9,0.05,1)
dim(DEgenesGroup9)
##3280

upGenesGroup9<-DEgenesGroup9[DEgenesGroup9$logFC > 1,]
tmp<-upGenesGroup9[upGenesGroup9$AveExpr > 0,]
upGenesGroup9<-tmp
dim(upGenesGroup9)
##557

########################################
##### 10. (MhcHiDura+MhcHiCp)
#####       vs others
########################################

### Create design matrix
TS <- c(rep("wanted",4),rep("others",11),rep("wanted",4),rep("others",8))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=wanted-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup10<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup10<-getDEgenes(allGenesGroup10,0.05,1)
dim(DEgenesGroup10)
##1592

upGenesGroup10<-DEgenesGroup10[DEgenesGroup10$logFC > 1,]
tmp<-upGenesGroup10[upGenesGroup10$AveExpr > 0,]
upGenesGroup10<-tmp
dim(upGenesGroup10)
##197

########################################
##### 20. (MhcHiDura+MhcIntDura)
#####       vs others
########################################

### Create design matrix
TS <- c(rep("wanted",8),rep("others",19))
cbind(colnames(expTable),TS)
TS <- factor(TS, levels=unique(TS))
designTmp <- model.matrix(~0+TS)
colnames(designTmp) <- levels(TS)
designTmp

#### Fit linear model on data
fit <- lmFit(v, designTmp)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=wanted-others,levels=designTmp)
cont.matrix

fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Adjust P-values via Benjamini-Hochberg
allGenesGroup20<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=1)
DEgenesGroup20<-getDEgenes(allGenesGroup20,0.05,1)
dim(DEgenesGroup20)
##2457

upGenesGroup20<-DEgenesGroup20[DEgenesGroup20$logFC > 1,]
tmp<-upGenesGroup20[upGenesGroup20$AveExpr > 0,]
upGenesGroup20<-tmp
dim(upGenesGroup20)
##320


################################################################################
########## FIND CORE GENES
################################################################################
colsMhcHiDura<-grep('MHCII_high',colnames(expTable))
colsMhcIntDura<-grep('MHCII_int',colnames(expTable))
colsPeri<-grep('^PM_',colnames(expTable))
colsBrain<-grep('NMG',colnames(expTable))
colsMhcHiCp<-grep('NCPM',colnames(expTable))
colsLung<-grep('AMF',colnames(expTable))
colsLiver<-grep('KC_',colnames(expTable))

########################################
##### 1. MhcHiDura
########################################
colsWanted<-colsMhcHiDura
colsOthers<-c(colsMhcIntDura,colsPeri,colsBrain,colsMhcHiCp,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup1)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##16

### Select markers with high expression level
mhcHiDuraGenes<-getBestMarkers(bestGenes, 100)
length(mhcHiDuraGenes)
##2


########################################
##### 2. MhcIntDura
########################################
colsWanted<-colsMhcIntDura
colsOthers<-c(colsMhcHiDura,colsPeri,colsBrain,colsMhcHiCp,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup2)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##53

### Select markers with high expression level
mhcIntDuraGenes<-getBestMarkers(bestGenes, 200)
length(mhcIntDuraGenes)
##12

### Manually add Gas6
mhcIntDuraGenes<-c(mhcIntDuraGenes, 'Gas6')
length(mhcIntDuraGenes)
##13

### Put genes in a certain order
tmp<-c("Lyve1","Colec12","Cd4","Tmem8","Ddx60","Tbc1d4","Plekhg5","Pla2g2d","Gpx3","Mmp9","Bmp2","Gas6","Cebpd")
setdiff(mhcIntDuraGenes, tmp)
mhcIntDuraGenes<-tmp
length(mhcIntDuraGenes)
##13


########################################
##### 3. Brain
########################################
colsWanted<-colsBrain
colsOthers<-c(colsMhcHiDura,colsMhcIntDura,colsPeri,colsMhcHiCp,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup4)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-0
blueBlocksAllowed<-0

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##549

### Select markers with high expression level
brainGenes<-getBestMarkers(bestGenes, 400)
length(brainGenes)
##76

brainGenesAll<-getBestMarkers(bestGenes, 200)
length(brainGenesAll)
##161

### Manually remove 'Mlec', 'Id2', 'Edem1', 'Ski', 'Bin2' and 'Ivns1abp'
### Manually add 'Fcrls'
brainGenes<-brainGenes[brainGenes!='Mlec']
brainGenes<-brainGenes[brainGenes!='Id2']
brainGenes<-brainGenes[brainGenes!='Edem1']
brainGenes<-brainGenes[brainGenes!='Ski']
brainGenes<-brainGenes[brainGenes!='Bin2']
brainGenes<-brainGenes[brainGenes!='Ivns1abp']
brainGenes<-c(brainGenes, 'Fcrls')
length(brainGenes)
##71


########################################
##### 4. MhcHiCp
########################################
colsWanted<-colsMhcHiCp
colsOthers<-c(colsMhcHiDura,colsMhcIntDura,colsPeri,colsBrain,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup5)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##311

### Select markers with high expression level
mhcHiCpGenes<-getBestMarkers(bestGenes, 200)
length(mhcHiCpGenes)
##12

### Manually remove 'Pmepa1' and 'Timp3'
mhcHiCpGenes<-mhcHiCpGenes[mhcHiCpGenes!='Pmepa1']
mhcHiCpGenes<-mhcHiCpGenes[mhcHiCpGenes!='Timp3']
length(mhcHiCpGenes)
##10

### Put genes in a certain order
tmp<-c("Kl","Cpd","Slc22a17","Enpp2","Ttr","Cyr61","Col14a1","Igfbp2","Clu","Cited2")
setdiff(mhcHiCpGenes, tmp)
setdiff(tmp,mhcHiCpGenes)
mhcHiCpGenes<-tmp
length(mhcHiCpGenes)
##10



########################################
##### 5. Combi1: (MhcHiDura+MhcIntDura+
#####         MhcHiCp+brain) vs others
########################################
colsWanted<-c(colsMhcHiDura,colsMhcIntDura,colsBrain,colsMhcHiCp)
colsOthers<-c(colsPeri,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup8)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##94

### Select markers with high expression level
combi1Genes<-getBestMarkers(bestGenes, 200)
length(combi1Genes)
##33

### Put genes in a certain order
tmp<-c("Serinc3","Marcks","Tgfbr2","Abca9","P2ry6","Sft2d2","Gna12","P2rx7","Ptpra","Evi2a","Tnfrsf11a","Spred1","Il21r",
       "Tspan3","Ptpro","Efhd2","Rab3il1","Pip4k2a","Galnt1","Sec14l1","Tlr7","Nisch","Pacsin2","Rassf2","Arhgef6","Tifa",
       "Dgkd","Ophn1","Btg2","Ier5","Mef2c","Lmo2","Nfatc2")
setdiff(combi1Genes, tmp)
combi1Genes<-tmp
length(combi1Genes)
##33


########################################
##### 6. Combi2: (MhcHiDura+MhcIntDura+
#####               MhcHiCp) vs others
########################################
colsWanted<-c(colsMhcHiDura,colsMhcIntDura,colsMhcHiCp)
colsOthers<-c(colsPeri,colsBrain,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup9)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##68

### Select markers with high expression level
combi2Genes<-getBestMarkers(bestGenes, 200)
length(combi2Genes)
##25

### Manually add 'Clec12a','Tgfbi','Pla2g7' and 'Aoah'
combi2Genes<-c(combi2Genes,c("Clec12a","Tgfbi","Pla2g7","Aoah"))
length(combi2Genes)
##29

### Put genes in a certain order
tmp<-c("C3ar1","Nrp1","Cd63","Ms4a7","Adrb2","Ms4a6b","Ms4a6c","Ifnar1","Ptger4","Lifr","Clec12a","Zfp36l2","Ehd4","Zfp36l1",
       "Myo5a","Swap70","Rasa4","B3galnt1","St8sia4","H1f0","Sesn1","Apoe","Cd14","Sdc4","Aoah","Pla2g7","Tgfbi","Ier2","Maf")
setdiff(combi2Genes, tmp)
combi2Genes<-tmp
length(combi2Genes)
##29


########################################
##### 7. Combi3: (MhcHiDura+MhcHiCp)
#####               vs others
########################################
colsWanted<-c(colsMhcHiDura,colsMhcHiCp)
colsOthers<-c(colsMhcIntDura,colsPeri,colsBrain,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup10)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##32

### Select markers with high expression level
combi3Genes<-getBestMarkers(bestGenes, 100)
length(combi3Genes)
##10

### Put genes in a certain order
tmp<-c("H2-Ab1","H2-Aa","H2-Eb1","Tmem176a","Gpr65","Cd74","H2-DMb1","Ciita","Cnn2","Fgl2")
setdiff(combi3Genes, tmp)
combi3Genes<-tmp
length(combi3Genes)
##10


########################################
##### 8. Combi4: (MhcHiDura+MhcIntDura)
#####               vs others
########################################
colsWanted<-c(colsMhcHiDura,colsMhcIntDura)
colsOthers<-c(colsMhcHiCp,colsPeri,colsBrain,colsLung, colsLiver)
DEgenes<-rownames(DEgenesGroup20)
meanMedianStep<-1.0
isRedCutOff<-0.5
redBlocksAllowed<-1
blueBlocksAllowed<-1

bestGenes<-getCoreGenes(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed)
length(bestGenes)
##68

### Select markers with high expression level
combi4Genes<-getBestMarkers(bestGenes, 200)
length(combi4Genes)
##25

### Put genes in a certain order
tmp<-c("Mrc1","D17H6S56E-5","Eps15","Wwp1","Clec10a","Prkacb","Lrp6","Cbr2","Ap1b1","Dnajb1","Snx6","Samd9l","Parp14",
       "Sbf2","Snx8","Map3k8","Dab2","AI607873","Sepp1","F13a1","Ccl7","Ccl8","Cp","Gbp7","Nfxl1")
setdiff(combi4Genes, tmp)
combi4Genes<-tmp
length(combi4Genes)
##25


########################################
##### Functions
########################################

##### Function getCoreGenes
getCoreGenes<-function(colsWanted, colsOthers, DEgenes, meanMedianStep, isRedCutOff, redBlocksAllowed, blueBlocksAllowed){
  ##### Select some genes with nice difference
  myTable<-as.data.frame(countData_norm[DEgenes,], stringsAsFactors=FALSE)
  
  myTable$mean1<-apply(myTable[,colsOthers],1,mean)
  myTable$mean2<-apply(myTable[,colsWanted],1,mean)
  
  myTable$median1<-apply(myTable[,colsOthers],1,median)
  myTable$median2<-apply(myTable[,colsWanted],1,median)
  
  myTable$diffMean<-myTable$mean2/myTable$mean1
  myTable$diffMedian<-myTable$median2/myTable$median1
  
  goodGenes1<-myTable[abs(myTable$diffMean)>meanMedianStep,]
  goodGenes2<-myTable[abs(myTable$diffMedian)>meanMedianStep,]
  
  goodGenes<-intersect(rownames(goodGenes1),rownames(goodGenes2))
  length(goodGenes)
  
  ##### Check the selected genes in the heatmap
  ### Normalize
  expProfiles<-countData_norm[goodGenes,]
  expProfilesNorm<-normalizePerGene(expProfiles)
  
  ### Check the genes
  bluePart<-expProfilesNorm[,colsOthers]
  redPart<-expProfilesNorm[,colsWanted]
  
  ### Find good 'blue' genes
  goodBlueGenes<-c()
  for(i in 1:nrow(bluePart)){
    testGene<-bluePart[i,]
    theSum<-sum(testGene[testGene>0]>isRedCutOff)
    if(theSum<=redBlocksAllowed){
      goodBlueGenes<-c(goodBlueGenes,rownames(bluePart)[i])
    }
  }
  
  ### Find good 'red' genes
  tmp2<-apply(redPart>0.5, 1, sum)
  cutoff<-length(colsWanted)-blueBlocksAllowed
  part2<-tmp2[tmp2>=cutoff]
  goodRedGenes<-names(part2)
  
  ### Final core genes
  bestGenes<-intersect(goodBlueGenes,goodRedGenes)
  
  return(bestGenes)
}


##### Function getBestMarkers
getBestMarkers<-function(bestGenes, cutoff){
  test<-as.data.frame(countData_norm[bestGenes,], stringsAsFactors = FALSE)
  test$maxValue<-apply(test,1,max)
  
  tmp<-test[test$maxValue >= cutoff,]
  return(rownames(tmp[order(tmp$maxValue, decreasing = T),]))
}



######################################################################
########## CREATE HEATMAP
######################################################################

########################################
##### Preparation
########################################

wantedSamples<-c(colsBrain,colsMhcIntDura,colsMhcHiDura,colsMhcHiCp,colsPeri,colsLung,colsLiver)
gapsCol<-c(3,7,11,15,19,23)
gapsRow<-NULL

########## 1. Heatmaps with all found marker genes ##########
toPlot<-normalizePerGene(countData_norm[brainGenes,wantedSamples])
##we want the genes ordered based on their expression in NCPM (MHCII hi Cp):
meanNCPM<-rowSums(toPlot[,grep('NCPM_',colnames(toPlot))])
toPlot<-toPlot[names(sort(meanNCPM)),]
cellSize<-8
fontSize<-8
## => heatmap_markers_brain.pdf


mergeGenes<-c(combi1Genes,combi2Genes,combi4Genes,mhcIntDuraGenes,mhcHiDuraGenes,combi3Genes,mhcHiCpGenes)
toPlot<-normalizePerGene(countData_norm[mergeGenes,wantedSamples])
cellSize<-5.6
fontSize<-5.6
gapsRow<-c(length(combi1Genes),
           length(combi1Genes)+length(combi2Genes),
           length(combi1Genes)+length(combi2Genes)+length(combi4Genes),
           length(combi1Genes)+length(combi2Genes)+length(combi4Genes)+length(mhcIntDuraGenes),
           length(combi1Genes)+length(combi2Genes)+length(combi4Genes)+length(mhcIntDuraGenes)+length(mhcHiDuraGenes),
           length(combi1Genes)+length(combi2Genes)+length(combi4Genes)+length(mhcIntDuraGenes)+length(mhcHiDuraGenes)+length(combi3Genes))
## => heatmap_combiOfHeatmaps.pdf



########################################
##### Create heatmap normV1
########################################
### Prepare heatmap
library('pheatmap')
library('grid')
myColorPalette2<-c("#125ba2", "#135da4", "#1661a7", "#1864aa", "#1967ad", 
                   "#1c6ab0", "#1f6db1", "#2372b4", "#2675b6", "#2a79b8", "#2e7ebb", "#3181bc", "#3685bf", "#3888c1", "#3d8dc3", 
                   "#4090c5", "#4794c7", "#4c98c9", "#539dcb", "#5ba1ce", "#60a5d0", "#67aad2", "#6cadd4", "#74b2d6", "#79b5d8", 
                   "#81badb", "#8bbfde", "#99c7e2", "#a8d0e6", "#b2d5e9", "#c1dded", "#cbe3f0", "#daebf4", "#e4f0f7", "#f3f8fb", 
                   "#fffefd", "#fff8f5", "#feefe9", "#fee9e1", "#fee0d4", "#fedacc", "#fdd1c0", "#fdcbb8", "#fdc2ac", "#fcb9a0", 
                   "#fcb398", "#fcab8f", "#fca68a", "#fc9e82", "#fc997c", "#fc9174", "#fc8c6e", "#fb8466", "#fb7d5d", "#fb7758", 
                   "#fb7050", "#f96a4d", "#f66348", "#f45e45", "#f15640", "#ee4e3b", "#ec4937", "#ea4133", "#e83c2f", "#e5342a", 
                   "#e32f27", "#dd2c25", "#d92924", "#d32622", "#cd2220")

### Create heatmap for normV1
pheatmap(as.matrix(toPlot),color=myColorPalette2,scale="none", cluster_rows=F, cluster_cols=F,
         border_color = "gray25", cellwidth = cellSize, cellheight = cellSize, gaps_col = gapsCol, gaps_row = gapsRow, fontsize=fontSize,
         show_rownames = T, show_colnames = T)




######################################################################
########## CREATE VOLCANOPLOT
######################################################################


##### 14. MF.MhcHiCp vs MF.MhcHiDura #####
toPlot<-allGenesGroup14
neededCols<-c(colsMhcHiCp,colsMhcHiDura)

# nonSig      sig sigExprs 
# 12176     1006      289 

##### 15. MF.MhcHiDura vs MF.MhcIntDura #####
toPlot<-allGenesGroup15
neededCols<-c(colsMhcHiDura,colsMhcIntDura)

# nonSig      sig sigExprs 
# 12930      372      169 


#################### Create Volcano plot ####################
##### Preparation #####
toPlot$threshold1<-toPlot$adj.P.Val < 0.01
toPlot$threshold2<-toPlot$logFC < -1 | abs(toPlot$logFC)  > 1
toPlot$sig<-toPlot$threshold1 & toPlot$threshold2
toPlot$sig[toPlot$sig == TRUE]<-"sig"
toPlot$sig[toPlot$sig == FALSE]<-"nonSig"
sigGenes<-rownames(toPlot[toPlot$sig=="sig",])
length(sigGenes)
table(toPlot$sig)

exprsGenes<-rownames(countData_norm)[rowSums(countData_norm[,neededCols] >= 30) >= 3]
sigExprsGenes<-intersect(exprsGenes, sigGenes)
# sigExprsGenesGroup12<-sigExprsGenes
# sigExprsGenesGroup13<-sigExprsGenes
# sigExprsGenesGroup21<-sigExprsGenes
toPlot[sigExprsGenes,which(colnames(toPlot)=="sig")]<-"sigExprs"
toPlot$sig<-as.factor(toPlot$sig)
table(toPlot$sig)

##### Create plot #####
p <- ggplot(data=toPlot, aes(x=logFC, y =-log10(adj.P.Val))) +
  geom_point(aes(colour=factor(sig), fill = factor(sig)), shape=21, size = 3.0, alpha=0.8) +
  scale_fill_manual(values=c("nonSig"="gray10", "sig"="#0832F0","sigExprs"="#EA1B03")) + 
  scale_colour_manual(values=c("black", "black", "black")) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  theme_bw() +
  theme(legend.position="right", panel.grid = element_blank())
print(p)


######################################################################
########## CREATE VENNDIAGRAM
######################################################################

##### Venn diagram #####
length(sigExprsGenesGroup12)
length(sigExprsGenesGroup13)
length(sigExprsGenesGroup21)

list1<-sigExprsGenesGroup12
list2<-sigExprsGenesGroup13
list3<-sigExprsGenesGroup21

write.table(sigExprsGenesGroup12,file="sigExprsGenesGroup12.txt")
write.table(sigExprsGenesGroup13,file="sigExprsGenesGroup13.txt")
write.table(sigExprsGenesGroup21,file="sigExprsGenesGroup21.txt")

### Put in "Venny tool" http://bioinfogp.cnb.csic.es/tools/venny/index.html

