# You might need to install some of these dependencies first:
library("R.utils")
library("utils")
library("graphics")
library("stats")
library("data.table")
library("mixtools")
library("GSEABase")
library("SummarizedExperiment")
library("methods")
library("Biobase")
library("zoo")
library("DT")
library("NMF")
library("plotly")
library("BiocStyle")
library("rmarkdown")
library("doMC")
library("doRNG")
library("doParallel")
library("foreach")
library("dynamicTreeCut")

# setRepositories(ind=1:2)
# install.packag

# BiocManager::install("GENIE3", version = "3.8")
# BiocManager::install("AUCell", version = "3.8")
# BiocManager::install("RcisTarget", version = "3.8")
# devtools::install_github("aertslab/SCENIC", dep = FALSE)

library("GENIE3")
library("AUCell")
library("RcisTarget")
library("RcisTarget.mm9.motifDatabases.20k")
library("SCENIC")

# install.packages('qsub')
#More info: https://cran.r-project.org/web/packages/qsub/readme/README.html
library("qsub")


options(stringsAsFactors=FALSE)

################################################################################
########## LOAD DATA
################################################################################

########################################
### Load sample and
### select clusters of interest

### Here on microglia clusters of
### Irf8 WT+KO whole brain sample
########################################

library('Seurat')
seuratObj<-readRDS(file="Robjects/seuratObj.rds")
TSNEPlot(seuratObj, do.label=T, label.size = 7, pt.size = 1)

## WT microglia = cl0+1
## KO microglia = cl2+3+4+9
cellsWT<-WhichCells(seuratObj, ident=c(0,1))
cellsKO_part1<-WhichCells(seuratObj, ident=c(2,3,9))
cellsKO_part2<-WhichCells(seuratObj, ident=c(4))

##Create new tSNE plot
seuratObjNew<-seuratObj
seuratObjNew<-SetIdent(seuratObjNew, cells.use = cellsWT, ident.use = 20)
seuratObjNew<-SetIdent(seuratObjNew, cells.use = cellsKO_part1, ident.use = 21)
seuratObjNew<-SetIdent(seuratObjNew, cells.use = cellsKO_part2, ident.use = 22)

TSNEPlot(seuratObjNew, do.label=T, label.size = 7, pt.size = 1)

clusterMatrix<-seuratObj@meta.data
clusterMatrix$resFinal<-clusterMatrix$res.0.8
clusterMatrix[clusterMatrix$resFinal %in% c(0,1), which(colnames(clusterMatrix)=="resFinal")]<-"20"
clusterMatrix[clusterMatrix$resFinal %in% c(2,3,9), which(colnames(clusterMatrix)=="resFinal")]<-"21"
clusterMatrix[clusterMatrix$resFinal %in% c(4), which(colnames(clusterMatrix)=="resFinal")]<-"22"

##### Get raw counts #####
dataAggr <- Read10X("filtered_gene_bc_matrices_mex/mm10/")
rawData<-as.matrix(dataAggr)
rm(dataAggr)

exprMat <- rawData[,c(cellsWT, cellsKO_part1, cellsKO_part2)]
dim(exprMat)
# 27998 9606
save(exprMat, file="data/exprMat.RData")

exprMat[1:5,1:5]

##### Get cluster info #####
clusterMatrixSlice<-clusterMatrix[clusterMatrix$resFinal %in% c(20,21,22),]
dim(clusterMatrixSlice)
# 9606    6

clusterMatrix<-clusterMatrixSlice
rm(clusterMatrixSlice)

################################################################################
########## 1. RUN SCENIC (part 1)
################################################################################

########################################
### Load meta data info for each cell
### Select colors per selected cluster
########################################

########## Cell info ##########
### Create a variable 'cellInfo' with cellbarcodes as rownames, 'orig.ident' as column1 and the 'level1class' (= the cluster number) as column2
#                           orig.ident level1class
# K10_AAACCTGAGATTACCC        K10           2
# K10_AAACCTGAGCAATATG        K10           2
# K10_AAACCTGAGTGCCAGA        K10           2
# K10_AAACCTGCACTCTGTC        K10           3
# K10_AAACCTGCAGCCAGAA        K10           2
# K10_AAACCTGCATCGATGT        K10           0

neededCols<-c(grep('orig.ident',colnames(clusterMatrix)), grep('resFinal',colnames(clusterMatrix)))
cellInfo<-clusterMatrix[,neededCols]
colnames(cellInfo)[2]<-"level1class"
dim(cellInfo)
# 9606    2
table(cellInfo$level1class)
# 20   21   22 
# 5044 3465 1097
save(cellInfo, file="data/cellInfo.RData")

##### Get colors #####
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(level1class=setNames(c("#D39200", "#00BA38","#7f00ff"), 
                                     c("20", "21","22")),
                orig.ident=setNames(c("#F8766D", "#00bfc4"), 
                                    c("JP34", "JP35")))
save(colVars, file="data/colVars.RData")

plot.new(); legend(0,1, fill=colVars$level1class, legend=names(colVars$level1class))

###Remove some variables
rm(seuratObj)
rm(seuratObjNew)
rm(rawData)
rm(cellsKO)
rm(cellsWT)

########## Load organisms TFs ##########
# Get genes in databases:
data(mm9_500bpUpstream_motifRanking) # or 10kbp, they should have the same genes
genesInDatabase <- mm9_500bpUpstream_motifRanking@rankings$rn

# Get TFS in databases:
data(mm9_direct_motifAnnotation)
allTFs <- mm9_direct_motifAnnotation$allTFs

########## Gene filter/selection ##########
nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     0.0     3.0   515.2   411.0  9606.0
summary(nCountsPerGene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0       0       3    1302     431 2739397
max(exprMat)
# 2320
sum(exprMat>0) / sum(exprMat==0)
# 0.05667556

###First filter: remove genes that are very lowly expressed and thus might be noise
minReads <- 3*.01*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)
# 7993

###Second filter: remove genes that are expressed in only a few cells
minSamples <- ncol(exprMat)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)
# 7993

###Select only the genes present in the RcisTarget databases
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)
# 7329

########## Do filtering ##########
exprMatrix_filtered <- exprMat[genesLeft_minCells_inDatabases, ]
dim(exprMatrix_filtered)
# 7329 9606
save(exprMatrix_filtered, file="int/1.1_exprMatrix_filtered.RData")

# Check whether any relevant gene / potential gene of interest is missing:
interestingGenes <- c("Neurod1", "Sox10", "Dlx1")
interestingGenes[which(!interestingGenes %in% rownames(exprMatrix_filtered))]

rm(exprMat)

########## Extract list of TFs from RcisTarget annotation database ##########
inputTFs <- allTFs[allTFs%in% rownames(exprMatrix_filtered)]
save(inputTFs, file="int/1.2_inputTFs.RData")

c(allTFs=length(allTFs), inputTFs=length(inputTFs))
# allTFs inputTFs 
# 1620      629

########## Run GENIE3 ##########
## Input is expresssion matrix and list of candidate regulators (= the list of transcription factors from RcisTarget databases)

load(file="int/1.1_exprMatrix_filtered.RData")
load(file="int/1.2_inputTFs.RData")

### Run via cluster
weightMatrix <- GENIE3(exprMatrix_filtered, regulators=inputTFs, nCores=24)

weightMatrix<-output$method_output
dim(weightMatrix)
# 629 7329
save(weightMatrix, file="int/1.3_GENIE3_weightMatrix.RData")


########## Get correlation of GENIE3 results ##########
load(file="int/1.1_exprMatrix_filtered.RData")
dim(exprMatrix_filtered)
# 7329 9606

##This takes a while!!
corrMat <- cor(t(exprMatrix_filtered), method="spearman")
dim(corrMat)
# 7329 7329
save(corrMat, file="int/1.4_corrMat.RData")


################################################################################
########## 1b. RUN SCENIC: Create co-expression modules (based on GENIE3 output)
################################################################################

# ########## Run SCENIC using the wrapper ##########
# load("data/esetMouseBrain.RData")
# exprMat <- exprs(esetMouseBrain)
# # Optional: add log for TF expression plot, it does not affect any other calculation
# exprMat <- log2(exprMat+1)
# dim(exprMat)
# # 19970  3005
# 
# load("data/colVars.RData")
# cellInfo <- pData(esetMouseBrain)[colnames(exprMat), names(colVars), drop=F]
# 
# 
# library(SCENIC)
# runSCENIC(exprMat=exprMat, org="mm9", cellInfo=cellInfo, colVars=colVars, nCores=4,  stepsToRun=c("1.2", "2", "3.1", "3.2", "4"))


########## Load output from GENIE3 ##########
library("GENIE3")

### Convert the weight matrix into links (list):
load("int/1.3_GENIE3_weightMatrix.RData")
linkList <- getLinkList(weightMatrix, threshold=0.001) # (slighly faster)
colnames(linkList) <- c("TF", "Target", "weight")

#Order by weight
linkList <- linkList[order(linkList[,"weight"], decreasing=TRUE),]
save(linkList, file="int/1.5_GENIE3_linkList.RData")

load("int/1.5_GENIE3_linkList.RData")
dim(linkList)
# 3045479       3

head(linkList)
#     TF Target    weight
# 1 Rps4x Rpl37a 0.1433998
# 2 Rps4x  Rpl32 0.1433057
# 3 Rps4x  Rps23 0.1421255
# 4 Rps4x  Rps29 0.1413821
# 5 Rps4x  Rps27 0.1407486
# 6 H2afz  Stmn1 0.1395439

########## Creating TF modules ##########
##### A. Keep only the TF-targets with weight > 0.001 #####
quantile(linkList$weight, probs=c(0.75, 0.90))
# 75%         90% 
#   0.002238495 0.003302075

plot(linkList$weight[1:1000000], type="l", ylim=c(0, max(linkList$weight)), main="Weight of the links",
     ylab="Weight", xlab="Links sorted decreasingly")
abline(h=0.001, col="blue") # Threshold

###How many percent of the weight matrix survives this treshold?
sum(linkList$weight>0.001)/nrow(linkList)*100
# 100

###Do the filtering
linkList_001 <- linkList[which(linkList[,"weight"]>0.001),]
# Number of links over the threshold: 
nrow(linkList_001) 
# 3045479 = all links

##### B. Create the geneSets for each TF #####
## For each TF, select:
##1.Targets with weight > 0.001
##2.Targets with weight > 0.005
##3.Top 50 targets (targets with highest weight)
##4.Targets for which the TF is within its top 5 regulators
##5.Targets for which the TF is within its top 10 regulators
##6.Targets for which the TF is within its top 50 regulators

tfModules <- list()
linkList_001$TF <- as.character(linkList_001$TF)
linkList_001$Target <- as.character(linkList_001$Target)
head(linkList_001)

#### Create TF-modules:
# 1: Weight > 0.001 (filtered in previous step) 
tfModules[["w001"]] <- split(linkList_001$Target, factor(linkList_001$TF))

# 2: Weight > 0.005
llminW <- linkList_001[which(linkList_001[,"weight"]>0.005),]
tfModules[["w005"]] <- split(llminW$Target, factor(llminW$TF))

# 3: Top 50 targets for each TF
# ("w001" should be ordered decreasingly by weight)
tfModules[["top50"]] <- lapply(tfModules[["w001"]], function(x) x[1:(min(length(x), 50))])

# 4-6: Top regulators per target 
# (linkList_001 should be ordered by weight!)
linkList_001_byTarget <- split(linkList_001, factor(linkList_001$Target))
save(linkList_001_byTarget, file="int/1.5_linkList_001_byTarget.RData")

nTopTfs <- c(5, 10, 50)
nTopTfs <- setNames(nTopTfs, paste("top", nTopTfs, "perTarget", sep=""))

library(reshape2); library(data.table)
topTFsperTarget <- lapply(linkList_001_byTarget, function(llt) {
  nTFs <- nTopTfs[which(nTopTfs <= nrow(llt))]
  melt(lapply(nTFs, function(x) llt[1:x,"TF"]))
})
topTFsperTarget <- topTFsperTarget[which(!sapply(sapply(topTFsperTarget, nrow), is.null))]
topTFsperTarget.asDf <-  data.frame(rbindlist(topTFsperTarget, idcol=TRUE))
head(topTFsperTarget.asDf)
colnames(topTFsperTarget.asDf) <- c("Target", "TF", "method")

# Merge the all the gene-sets:
tfModules.melted <- melt(tfModules)
colnames(tfModules.melted) <- c("Target", "TF", "method")
tfModules <- rbind(tfModules.melted, topTFsperTarget.asDf)

save(tfModules, file="int/1.6_tfModules.RData")

load("int/1.6_tfModules.RData")
# Basic counts:
rbind(nGeneSets=nrow(tfModules), 
      nTFs=length(unique(tfModules$TF)), 
      nTargets=length(unique(tfModules$Target)))
# nGeneSets 3659687
# nTFs          629
# nTargets     7329

########## Split into positive- and negative-correlated targets ##########
load("int/1.4_corrMat.RData")
### Keep only correlation between TFs and potential targets
tfs <- unique(tfModules$TF)
corrMat <- corrMat[tfs,]

### Split TF modules according to correlation
##if corr > 0.03 => get value of 1 (is activation)
##if corr < -0.03 => get value of -1 (is repression)
##otherwise => get value of 0 (not able to say if the link between TF and target is an activation or repression)
tfModules_byTF <- split(tfModules, factor(tfModules$TF))
tfModules_withCorr_byTF <- lapply(tfModules_byTF, function(tfGeneSets)
{
  tf <- unique(tfGeneSets$TF)
  targets <- tfGeneSets$Target
  cbind(tfGeneSets, corr=c(as.numeric(corrMat[tf,targets] > 0.03) - as.numeric(corrMat[tf,targets] < -0.03)))
})
tfModules_withCorr <- data.frame(rbindlist(tfModules_withCorr_byTF))
save(tfModules_withCorr, file="int/1.7_tfModules_withCorr.RData")

load("int/1.7_tfModules_withCorr.RData")
head(tfModules_withCorr)
#           Target            TF method corr
# 1           Cck 1810024B03Rik   w001    1
# 2          Iars 1810024B03Rik   w001    0
# 3           Axl 1810024B03Rik   w001    0
# 4         Ndor1 1810024B03Rik   w001    0
# 5 1700096K18Rik 1810024B03Rik   w001    0
# 6 1110051M20Rik 1810024B03Rik   w001    0
table(tfModules_withCorr$method)
# top10perTarget          top50 top50perTarget  top5perTarget           w001           w005 
# 73290          31450         366450          36645        3045479         106373
table(tfModules_withCorr$corr)
# -1       0       1 
# 21974 2488722 1148991

dim(tfModules_withCorr)
# 3659687       4


################################################################################
########## 2. RUN SCENIC: Get Regulons (direct TF targets via RcisTarget)
################################################################################

########## Motif databases ##########
library(RcisTarget.mm9.motifDatabases.20k)

# Motif rankings (genes x motifs)
### This takes a while!!!
data(mm9_500bpUpstream_motifRanking)
data(mm9_10kbpAroundTss_motifRanking)
motifRankings <- list()
motifRankings[["500bp"]] <- mm9_500bpUpstream_motifRanking
motifRankings[["10kbp"]] <- mm9_10kbpAroundTss_motifRanking
save(motifRankings, file="int/2.1.a_motifRankings.RData")

# Motif annotation (TFs)
data(mm9_direct_motifAnnotation)
direct_motifAnnotation <- mm9_direct_motifAnnotation
data(mm9_inferred_motifAnnotation) # optional
inferred_motifAnnotation <- mm9_inferred_motifAnnotation

# Get TFS in databases:
allTFs <- mm9_direct_motifAnnotation$allTFs
length(allTFs)
# 1620

########## Load gene-sets ##########
load("int/1.7_tfModules_withCorr.RData")

### Select the targets with positive correlation
tfModules_Selected <- tfModules_withCorr[which(tfModules_withCorr$corr==1),]

### Add a column with the geneSet name (TF_method)
tfModules_Selected <- cbind(tfModules_Selected, geneSetName=paste(tfModules_Selected$TF, tfModules_Selected$method, sep="_"))
head(tfModules_Selected)
#   Target            TF method corr        geneSetName
# 1        Cck 1810024B03Rik   w001    1 1810024B03Rik_w001
# 13     Fnbp4 1810024B03Rik   w001    1 1810024B03Rik_w001
# 16     Ercc3 1810024B03Rik   w001    1 1810024B03Rik_w001
# 25    Zfp358 1810024B03Rik   w001    1 1810024B03Rik_w001
# 31 Tnfrsf11a 1810024B03Rik   w001    1 1810024B03Rik_w001
# 41     Acot2 1810024B03Rik   w001    1 1810024B03Rik_w001

### Create list again (TF-modules, with several methods)
tfModules <- split(tfModules_Selected$Target, tfModules_Selected$geneSetName)
head(names(tfModules))
# [1] "1810024B03Rik_top50"          "1810024B03Rik_top50perTarget" "1810024B03Rik_w001"          
# [4] "1810024B03Rik_w005"           "2310011J03Rik_top50"          "2310011J03Rik_top50perTarget"

### Keep gene sets with at least 20 genes
tfModules <- tfModules[which(lengths(tfModules)>=20)]

### Add TF to the gene set (used in the following steps, careful if editing)
tfModules <- setNames(lapply(names(tfModules), function(gsn) {
  tf <- strsplit(gsn, "_")[[1]][1]
  unique(c(tf, tfModules[[gsn]]))
}), names(tfModules))
save(tfModules, file="int/2.1_tfModules_forMotifEnrichmet.RData")

### Explore
load("int/2.1_tfModules_forMotifEnrichmet.RData")
tfModulesSummary <- t(sapply(strsplit(names(tfModules), "_"), function(x) x[1:2]))
sort(table(tfModulesSummary[,2]))
# top5perTarget top10perTarget           w005          top50 top50perTarget           w001 
# 41             62             92            147            206            629


########## Run RcisTarget ##########
# load(file="int/2.1.a_motifRankings.RData")
load(file="int/2.1_tfModules_forMotifEnrichmet.RData")

##### A. Calculate motif enrichment for each TF-module #####
### 1.1 Get present motifs and calculate enrichment
##Run via cluster
motifs_AUC <- lapply(motifRankings, function(ranking) calcAUC(tfModules, ranking, aucMaxRank=0.01*nrow(ranking@rankings), nCores=12, verbose=FALSE))
save(motifs_AUC, file="int/2.2_motifs_AUC_500bp_10kbp.RData")


### 1.2 Convert to table, filter by NES & add the TFs to which the motif is annotated
# (For each database...)
### This takes a while!!
motifEnrichment <- lapply(motifs_AUC, function(aucOutput)
{
  # Extract the TF of the gene-set name (i.e. MITF_w001):
  tf <- sapply(setNames(strsplit(rownames(aucOutput), "_"), rownames(aucOutput)), function(x) x[[1]])
  
  # Calculate NES and add motif annotation (provide tf in 'highlightTFs'):
  addMotifAnnotation(aucOutput, highlightTFs=tf, nesThreshold=3.0, digits=3,
                     motifAnnot_direct=direct_motifAnnotation,
                     motifAnnot_inferred=inferred_motifAnnotation)
})

# Merge both tables, adding a column that contains the 'motifDb' 
motifEnrichment <- do.call(rbind, lapply(names(motifEnrichment), function(dbName){
  cbind(motifDb=dbName, motifEnrichment[[dbName]])
}))
save(motifEnrichment, file="int/2.3_motifEnrichment.RData")
cat("Number of motifs in the initial enrichment: ", nrow(motifEnrichment))
# 661269


### 1.3 Keep only the motifs annotated to the initial TF
motifEnrichment_selfMotifs <- motifEnrichment[which(motifEnrichment$TFinDB != ""),, drop=FALSE]
save(motifEnrichment_selfMotifs, file="int/2.4_motifEnrichment_selfMotifs.RData")
cat("Number of motifs annotated to the initial TF: ", nrow(motifEnrichment_selfMotifs))
# 6518

rm(motifEnrichment)


##### B. Prune targets = keep only the targets that clearly have the motif #####
motifEnrichment_selfMotifs_wGenes <- lapply(names(motifRankings), function(motifDbName){
  addSignificantGenes(resultsTable=motifEnrichment_selfMotifs[motifDb==motifDbName],
                      geneSets=tfModules,
                      rankings=motifRankings[[motifDbName]],
                      maxRank=5000, method="aprox", nCores=4)
})


### Save
library(data.table)
motifEnrichment_selfMotifs_wGenes <- rbindlist(motifEnrichment_selfMotifs_wGenes)
save(motifEnrichment_selfMotifs_wGenes, file="int/2.5_motifEnrichment_selfMotifs_wGenes.RData")

write.table(motifEnrichment_selfMotifs_wGenes, file="output/Step2_MotifEnrichment.tsv", sep="\t", quote=FALSE, row.names=FALSE)

###Explore
load("int/2.5_motifEnrichment_selfMotifs_wGenes.RData")
dim(motifEnrichment_selfMotifs_wGenes)
# 6518   12

motifEnrichment_selfMotifs_wGenes[order(NES,decreasing=TRUE)][1:5,-"enrichedGenes", with=F]


##########  Format regulons = create per TF only 1 list of targets ##########
library(data.table)
motifEnrichment.asIncidList <- apply(motifEnrichment_selfMotifs_wGenes, 1, function(oneMotifRow) {
  genes <- strsplit(oneMotifRow["enrichedGenes"], ";")[[1]]
  oneMotifRow <- data.frame(rbind(oneMotifRow), stringsAsFactors=FALSE)
  data.frame(oneMotifRow[rep(1, length(genes)),c("NES", "motif", "highlightedTFs", "TFinDB")], genes, stringsAsFactors = FALSE)
})
motifEnrichment.asIncidList <- rbindlist(motifEnrichment.asIncidList)
colnames(motifEnrichment.asIncidList) <- c("NES", "motif", "TF", "annot", "gene")
motifEnrichment.asIncidList <- data.frame(motifEnrichment.asIncidList, stringsAsFactors = FALSE)

# Get targets for each TF, but keep info about best motif/enrichment 
# (directly annotated motifs are considered better)
regulonTargetsInfo <- lapply(split(motifEnrichment.asIncidList, motifEnrichment.asIncidList$TF), function(tfTargets){
  # print(unique(tfTargets$TF))
  tfTable <- as.data.frame(do.call(rbind, lapply(split(tfTargets, tfTargets$gene), function(enrOneGene){
    directAnnot <- "**" %in% enrOneGene$annot
    enrOneGeneByAnnot <- enrOneGene
    if(directAnnot) enrOneGeneByAnnot <- enrOneGeneByAnnot[which(enrOneGene$annot == "**"),]
    bestMotif <- which.max(enrOneGeneByAnnot$NES)
    
    cbind(TF=unique(enrOneGene$TF), gene=unique(enrOneGene$gene), nMotifs=nrow(enrOneGene), 
          bestMotif=as.character(enrOneGeneByAnnot[bestMotif,"motif"]), NES=as.numeric(enrOneGeneByAnnot[bestMotif,"NES"]), 
          directAnnot=directAnnot)
  })), stringsAsFactors=FALSE)
  tfTable[order(tfTable$NES, decreasing = TRUE),]
})
regulonTargetsInfo <- rbindlist(regulonTargetsInfo)
colnames(regulonTargetsInfo) <- c("TF", "gene", "nMotifs", "bestMotif", "NES", "directAnnot")

# Optional: Add Genie3 score
load("int/1.5_GENIE3_linkList.RData")
linkList <- linkList[which(linkList$weight>=0.001),]
rownames(linkList) <- paste(linkList$TF, linkList$Target,sep="__")
regulonTargetsInfo <- cbind(regulonTargetsInfo, Genie3Weight=linkList[paste(regulonTargetsInfo$TF, regulonTargetsInfo$gene,sep="__"),"weight"])

save(regulonTargetsInfo, file="int/2.6_regulonTargetsInfo.RData")
write.table(regulonTargetsInfo, file="output/Step2_regulonTargetsInfo.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

### Split into regulons… (output: list TF –> targets)
regulonTargetsInfo_splitByAnnot <- split(regulonTargetsInfo, regulonTargetsInfo$directAnnot)
regulons <- sapply(split(regulonTargetsInfo_splitByAnnot[["TRUE"]], regulonTargetsInfo_splitByAnnot[["TRUE"]][,"TF"]), function(x) sort(as.character(unlist(x[,"gene"]))))
regulons_extended <- sapply(split(regulonTargetsInfo_splitByAnnot[["FALSE"]],regulonTargetsInfo_splitByAnnot[["FALSE"]][,"TF"]), function(x) unname(x[,"gene"]))
regulons_extended <- sapply(names(regulons_extended), function(tf) sort(unique(c(regulons[[tf]], regulons_extended[[tf]]))))
names(regulons_extended) <- paste(names(regulons_extended), "_extended", sep="")
regulons <- c(regulons, regulons_extended)
save(regulons, file="int/2.6_regulons_asGeneSet.RData")

#### Explore results
load("int/2.6_regulons_asGeneSet.RData")
# Number of regulons and summary of sizes:
length(regulons)
# 217

summary(lengths(regulons))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1      16      63     297     310    4172

class(regulons)
head(names(regulons))
head(regulons)

##########  Transform into an incidence matrix (i.e. network) ##########
##an incidence matrix = TFs as rows, genes as columns, and 0/1 as value indicating whether the TF regulates the gene

incidList <- melt(regulons)
incidMat <- table(incidList[,2], incidList[,1])
save(incidMat, file="int/2.6_regulons_asIncidMat.RData")
dim(incidMat)
# 217 6242


##########  Exploring regulons ##########
# How many TFs are self-regulating?
table(sapply(names(regulons), function(x) x %in% regulons[[x]]))
# FALSE  TRUE 
# 149    68 

# # Motifs associated to a TF (i.e. MITF regulons):
# selTF <- "Dlx5"
# subsetTable <- motifEnrichment_selfMotifs_wGenes[highlightedTFs %in% selTF][order(NES,decreasing=TRUE)][,-"enrichedGenes", with=F]
# subsetTable <- addLogo(subsetTable)
# library(DT)
# 
# datatable(subsetTable, escape=FALSE, filter="top", options=list(pageLength=5))

# # Gene enrichment plots
# geneSetName <- "Dlx5_top50"
# motifDbName <- "10kbp"
# selectedMotifs <- subsetTable[geneSet==geneSetName & motifDb==motifDbName, motif]
# selectedMotifs <- selectedMotifs[1:3]
# 
# pdf("int/2.8_RCC_selectedMotifs.pdf")
# par(mfrow=c(2,2))
# signifGenes_SelectedMotifs <- getSignificantGenes(tfModules[[geneSetName]], 
#                                                   motifRankings[[motifDbName]],
#                                                   signifRankingNames=selectedMotifs,
#                                                   plotCurve=TRUE, maxRank=5000, nCores=4, 
#                                                   genesFormat="geneList", method="aprox")
# 
# dev.off()
# 
# # Motif & number of genes:
# cbind(lengths(signifGenes_SelectedMotifs$enrichedGenes))


################################################################################
########## 3a. RUN SCENIC: Network activity in each cell
################################################################################


##########  Load regulons & expression matrix ##########
load(file="data/exprMat.RData")
dim(exprMat)
# 27998  9606

load("data/colVars.RData")
load("data/cellInfo.RData")
dim(cellInfo)
# 9606     2


### Keep only the regulons with >= 10 genes
load("int/2.6_regulons_asGeneSet.RData")
regulons <- regulons[order(lengths(regulons), decreasing=TRUE)]
regulons <- regulons[lengths(regulons)>=10]

### Add the TF itself to the regulon & rename the names of the regulons
regulons <- setNames(lapply(names(regulons), function(tf) sort(unique(c(gsub("_extended", "", tf), regulons[[tf]])))), names(regulons))
names(regulons) <- paste(names(regulons), " (",lengths(regulons), "g)", sep="")
save(regulons, file="int/3.0_regulons_forAUCell.RData")
length(regulons)
# 183

cbind(names(regulons)[1:10])

########## Calculate AUCell ##########

library(AUCell)
##### 1. Create rankings = per cell rank the genes based on their expression #####
aucellRankings <- AUCell.buildRankings(exprMat, nCores=4, plotStats=TRUE)
# Quantiles for the number of genes detected by cell: 
#   (Non-detected genes are shuffled at the end of the ranking. Keep in mind when choosing the threshold for calculating the AUC).
# min      1%      5%     10%     50%    100% 
# 554.00  848.05  998.00 1085.00 1470.00 3423.00 
aucellRankings
abline(v=aucellRankings@nGenesDetected["1%"], col="skyblue3", lwd=5, lty=3)
save(aucellRankings, file="int/3.1_aucellRankings.RData")

##### 2. Calculate AUC values #####
##the more the target genes of a regulon match with the highly expressed genes in a cell, the higher the AUC value

# regulonAUC <- AUCell.calcAUC(regulons, aucellRankings, aucMaxRank=aucellRankings@nGenesDetected["1%"], nCores=10)
regulonAUC <- AUCell.calcAUC(regulons, aucellRankings, nCores=6)
save(regulonAUC, file="int/3.2_regulonAUC.RData")

load("int/3.2_regulonAUC.RData")
### Order the modules by similarity, for easier exploration in the upcoming steps & save
variableRegulons <- names(which(apply(getAuc(regulonAUC), 1, sd) > 0))
reguDist <-as.dist(1-cor(t(getAuc(regulonAUC)[variableRegulons,]), method="spear"))
reguClust <- hclust(reguDist, method="ward.D2")
regulonClusters <- setNames(dynamicTreeCut::cutreeDynamic(reguClust, distM=as.matrix(reguDist), verbose = FALSE), reguClust$labels)
regulonOrder <- reguClust$labels[reguClust$order]
regulonOrder <- regulonOrder[order(regulonClusters[regulonOrder], decreasing = TRUE)]
regulonAUC@matrix <- regulonAUC@matrix[regulonOrder,]
save(regulonAUC, file="int/3.2_regulonAUC.RData")


########## Overview of cell states according to module activity (tSNE on AUC) ##########
# (It is recommended to try different perplexity values)
regulonAUC_subset <- subset(regulonAUC, onlyNonDirectExtended(rownames(regulonAUC)))

##### Calculate tSNE #####
# PCA-based t-SNE
set.seed(123)
tsneAUC <- Rtsne::Rtsne(t(getAuc(regulonAUC_subset)), initial_dims=10, perplexity=10)
rownames(tsneAUC$Y) <- colnames(regulonAUC_subset)
colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
save(tsneAUC, file="int/3.3_tsneRegulonAUC_PCA.RData")

# # Alternative: Distance-based t-SNE:
# corDist <- as.dist(1-cor(getAuc(regulonAUC_subset)))
# set.seed(123)
# tsneAUC <- Rtsne::Rtsne(corDist, is_distance=TRUE, perplexity=10)
# rownames(tsneAUC$Y) <- labels(corDist)
# colnames(tsneAUC$Y) <- c("tsne1", "tsne2")
# save(tsneAUC, file="int/3.3_tsneRegulonAUC_Dist.RData")

##### Plot tSNE #####
load("int/3.3_tsneRegulonAUC_PCA.RData")
tSNE <- tsneAUC$Y
par(mfrow=c(1,2))

# Number of genes detected:
nGenesPerCell <- apply(exprMat, 2, function(x) sum(x>0))
colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
cellColorNgenes <- setNames(adjustcolor(colorPal(10), alpha=.8)[as.numeric(cut(nGenesPerCell,breaks=10, right=F,include.lowest=T))], names(nGenesPerCell))

plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")

# Other known properties:
for(varName in names(colVars))
{
  cellColor <- setNames(colVars[[varName]][cellInfo[,varName]], rownames(cellInfo))
  plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
}


########## Plot AUC histograms ##########
## The distribuion of the AUC of a regulon across all the cells can provide important information about its activity.
## Regulons that are differentialy active across the cells will often show bimodal or skewed distributions.
### Takes a while!!!
library(SCENIC)
Cairo::CairoPDF("output/Step3_RegulonActivity_AUCtSNE.pdf", width=20, height=5)
par(mfrow=c(1,4))

# tSNE (colored by number of genes detected per cell)
plot(tSNE, col=cellColorNgenes[rownames(tSNE)], pch=16, main="nGenes", sub="t-SNE on regulon activity (AUC)")
plot(tSNE, col=cellColor[rownames(tSNE)], pch=16, main=varName, sub="t-SNE on regulon activity (AUC)")
plot.new(); plot.new()

# Plot module activity, thresholds & assignment:
cells_AUCellThresholds <- plot_aucTsne(tSNE=tSNE, exprMat=exprMat, regulonAUC=regulonAUC, alphaOff=0.1)
dev.off()
save(cells_AUCellThresholds, file="int/3.4_AUCellThresholds.RData")



########## Save thresholds as text ##########
load("int/3.4_AUCellThresholds.RData")

# Get cells assigned to each regulon
regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)

### Save threshold info as text (e.g. to edit/modify...)
trhAssignment <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$selected))
commentsThresholds <- sapply(cells_AUCellThresholds, function(x) unname(x$aucThr$comment))

table2edit <- cbind(regulon=names(trhAssignment), 
                    threshold=trhAssignment, 
                    nCellsAssigned=lengths(regulonsCells)[names(trhAssignment)],
                    AUCellComment=commentsThresholds, 
                    nGenes=gsub("[\\(g\\)]", "", regmatches(names(cells_AUCellThresholds), gregexpr("\\(.*?\\)", names(cells_AUCellThresholds)))),
                    clusteringOrder=1:length(trhAssignment), 
                    clusterGroup=regulonClusters[names(trhAssignment)], 
                    onlyNonDirectExtended=(names(trhAssignment) %in% onlyNonDirectExtended(names(trhAssignment))),
                    personalNotes="")
write.table(table2edit, file="int/3.5_1_AUCellThresholds.txt", row.names=F, quote=F, sep="\t")


################################################################################
########## 3b. RUN SCENIC: Binary network activity
################################################################################

########## Modify the regulon activity thresholds (OPTIONAL) ##########


########## Binary regulon activity matrix (Active regulons per cell) ##########
load("int/3.4_AUCellThresholds.RData")

##### Get cells assigned to each regulon #####
regulonsCells <- lapply(cells_AUCellThresholds, function(x) x$assignment)
length(regulonsCells)
# 183

##### Convert to matrix (regulons with zero assigned cells are lost) #####
regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[,1], regulonActivity[,2]))
class(binaryRegulonActivity) <- "matrix"
save(binaryRegulonActivity, file="int/3.6_BinaryRegulonActivity.RData")

dim(binaryRegulonActivity)
# 115 9606

binaryRegulonActivity[1:10,1:3]
#                         AAACCTGAGACCTTTG-1 AAACCTGAGACGCACA-1 AAACCTGAGAGCTGGT-2
# Atf4 (594g)                              0                  0                  0
# Atf4_extended (786g)                     0                  0                  0
# Atf7 (49g)                               0                  0                  0
# Batf (14g)                               0                  0                  0
# Batf_extended (21g)                      0                  0                  0
# Bhlhe40_extended (15g)                   0                  0                  1
# Bhlhe41_extended (206g)                  0                  0                  1
# Bmyc (535g)                              1                  1                  0
# Bmyc_extended (1002g)                    1                  1                  0
# Ccnt2 (13g)                              0                  0                  0



##### Remove duplicate regulons #####
## This matrix contains some duplicated regulons (e.g. for some TFs, there is a regulon based on direct annotation, and also
## the extended version) => save a filtered version, containing only “extended” regulons if there is not a regulon based on 
## direct annotation.
binaryRegulonActivity_nonDupl <- binaryRegulonActivity[which(rownames(binaryRegulonActivity) %in% onlyNonDirectExtended(rownames(binaryRegulonActivity))),]
save(binaryRegulonActivity_nonDupl, file="int/3.7_binaryRegulonActivity_nonDupl.RData")
dim(binaryRegulonActivity_nonDupl)
# 76 9606

binaryRegulonActivity_nonDupl[1:10,1:3]

##### Matrix overview #####
##How many regulons are assigned to at least 5 cells?
sum(rowSums(binaryRegulonActivity)>5)
# 98

cbind(nCellsOn=sort(rowSums(binaryRegulonActivity), decreasing=TRUE)[1:15])
#                         nCellsOn
# Relb (281g)                9604
# Relb_extended (357g)       9604
# Crem (13g)                 8470
# Runx1 (10g)                6868
# Zbtb7a_extended (14g)      6589
# Egr1_extended (4172g)      5496
# Elf1 (907g)                5448
# Nr3c1 (1138g)              5274
# Etv3 (310g)                5243
# Foxp1 (414g)               5228
# Foxp1_extended (551g)      5151
# Nr3c1_extended (1423g)     5110
# Etv3_extended (428g)       5089
# Elf1_extended (1329g)      5078
# Irf8 (158g)                5075

summary(rowSums(binaryRegulonActivity))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0   107.5   515.0  2118.9  4955.5  9604.0

##As boxplot
par(mfrow=c(1,2))
boxplot(rowSums(binaryRegulonActivity_nonDupl), main="nCells per regulon", 
        sub='number of cells \nthat have the regulon active',
        col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)
boxplot(colSums(binaryRegulonActivity_nonDupl), main="nRegulons per Cell", 
        sub='number of regulons \nactive per cell',
        col="darkolivegreen1", border="#001100", lwd=2, frame=FALSE)


########## Binary regulon activity heatmap ##########
load("data/exprMat.RData")
dim(exprMat)
# 27998 9606

load("data/colVars.RData")
load("data/cellInfo.RData")
dim(cellInfo)
# 9606     2

minCells <- ncol(exprMat) * .01


load("int/3.6_BinaryRegulonActivity.RData")
load("int/3.7_binaryRegulonActivity_nonDupl.RData")
regulonSelection <- list()

##### 1. All regulons #####
regulonSelection[["All regulons \n (including duplicated regulons)"]] <- rownames(binaryRegulonActivity)
length(rownames(binaryRegulonActivity))
# 115

##### 2. Active in > 1% cells #####
##regulon needs to be active in at least 1% of the cells (2974 cells)
regMinCells <- names(which(rowSums(binaryRegulonActivity_nonDupl) > minCells))
regulonSelection[["Regulons active in more than 1% of cells"]] <- regMinCells
length(regMinCells)
# 53

##### 3. Correlation across regulons (based on binary cell activity) #####
reguCor <- cor(t(binaryRegulonActivity_nonDupl[regMinCells,]))
diag(reguCor) <- 0
nrow(reguCor)
# 53

##### 4. Regulons that co-ocurr in similar cells. If a regulon is relevant by itself it will not be shown, also check the regulons ignored. #####
corrRegs <- names(which(rowSums(abs(reguCor) > 0.30) > 0))
regulonSelection[["Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)"]]  <- corrRegs
length(corrRegs)
# 33

missingRegs <- rownames(binaryRegulonActivity_nonDupl)[which(!rownames(binaryRegulonActivity_nonDupl) %in% corrRegs)]
regulonSelection[["Regulons no other regulons correlated\n with abs(cor)>0.30 \n or active in fewer than 1% of cells"]]  <- missingRegs

save(regulonSelection,file="int/3.8_regulonSelections.RData")

## Set regulon order (for plotting)
binaryRegulonOrder <- hclust(as.dist(1-reguCor[corrRegs,corrRegs]))
binaryRegulonOrder <- binaryRegulonOrder$labels[binaryRegulonOrder$order]
save(binaryRegulonOrder,file="int/3.9_binaryRegulonOrder.RData")


##### Plot heatmap #####
names(regulonSelection)
# [1] "All regulons \n (including duplicated regulons)"                                                       
# [2] "Regulons active in more than 1% of cells"                                                              
# [3] "Regulons with any other regulon correlated\n with abs(cor)>0.30 \n(and active in at least 1% of cells)"
# [4] "Regulons no other regulons correlated\n with abs(cor)>0.30 \n or active in fewer than 1% of cells"  

for(i in seq_len(length(regulonSelection)))
{
  print(paste0("Print heatmap ",i))
  selRegs <- names(regulonSelection)[i]
  if(length(regulonSelection[[selRegs]])>1)
  {
    binaryMat <- binaryRegulonActivity[regulonSelection[[selRegs]],,drop=FALSE]
    NMF::aheatmap(binaryMat, scale="none", revC=TRUE, main=selRegs,
                  annCol=cellInfo[colnames(binaryMat),, drop=FALSE],
                  annColor=colVars,
                  color = c("white", "black"),
                  filename=paste0("output/Step3.3_binaryRegulonActivity_Heatmap_",i,".pdf"))
  }
}


################################################################################
########## 4. RUN SCENIC: Binary network activity
################################################################################

########## Load cell info and binary activity matrix ##########
load("data/exprMat.RData")
dim(exprMat)
# 27998 9606

load("data/colVars.RData")
load("data/cellInfo.RData")
dim(cellInfo)
# 9606     2

load("int/3.7_binaryRegulonActivity_nonDupl.RData")
tBinaryAct <- t(binaryRegulonActivity_nonDupl)
tBinaryAct[1:5,1:5]
#                     Atf4 (594g) Atf7 (49g) Batf (14g) Bhlhe40_extended (15g)
# AAACCTGAGACCTTTG-1           0          0          0                      0
# AAACCTGAGACGCACA-1           0          0          0                      0
# AAACCTGAGAGCTGGT-2           0          0          0                      1
# AAACCTGAGCCAACAG-2           0          0          0                      1
# AAACCTGAGCTACCGC-2           0          0          0                      1

########## Calculate t-SNE on the binary regulon activity ##########
library(Rtsne)

# PCA based t-SNE
set.seed(123)
tBinaryAct_jitter <- jitter(tBinaryAct, factor=1)
tsneBinaryActivity_PCA <- Rtsne(tBinaryAct_jitter, initial_dims=5, perplexity=30)
rownames(tsneBinaryActivity_PCA$Y) <- rownames(tBinaryAct_jitter)
colnames(tsneBinaryActivity_PCA$Y) <- c("tsne2", "tsne1")
tsneBinaryActivity_PCA$Y <- tsneBinaryActivity_PCA$Y[,c("tsne1", "tsne2")]
save(tsneBinaryActivity_PCA, file="int/4.1_tsneBinaryActivity_5PC.RData")

# PCA based t-SNE
set.seed(123)
tBinaryAct_jitter <- jitter(tBinaryAct, factor=1)
tsneBinaryActivity_PCA <- Rtsne(tBinaryAct_jitter, initial_dims=50, perplexity=30)
rownames(tsneBinaryActivity_PCA$Y) <- rownames(tBinaryAct_jitter)
colnames(tsneBinaryActivity_PCA$Y) <- c("tsne2", "tsne1")
tsneBinaryActivity_PCA$Y <- tsneBinaryActivity_PCA$Y[,c("tsne1", "tsne2")]
save(tsneBinaryActivity_PCA, file="int/4.1_tsneBinaryActivity_50PC.RData")

# Distance-based t-SNE
corDist <- as.dist(1-cor(t(tBinaryAct)))
set.seed(123)
tsneBinaryActivity_Dist <- Rtsne(corDist, is_distance=TRUE, perplexity=30)
rownames(tsneBinaryActivity_Dist$Y) <- labels(corDist)
colnames(tsneBinaryActivity_Dist$Y) <- c("tsne1", "tsne2")
save(tsneBinaryActivity_Dist, file="int/4.1_tsneBinaryActivity_Dist.RData")


########## Plot t-SNEs coloured by known cell properties ##########
tSNEs_binary <- list()
load("int/4.1_tsneBinaryActivity_Dist.RData")
tSNEs_binary[["Dist"]] <- tsneBinaryActivity_Dist$Y
load("int/4.1_tsneBinaryActivity_5PC.RData")
tSNEs_binary[["5PC"]] <- tsneBinaryActivity_PCA$Y
load("int/4.1_tsneBinaryActivity_50PC.RData")
tSNEs_binary[["50PC"]] <- tsneBinaryActivity_PCA$Y

for(tsneName in names(tSNEs_binary))
{
  tSNE_binary <- tSNEs_binary[[tsneName]]
  
  # Density
  library(KernSmooth)
  library(RColorBrewer)
  dens2d <- bkde2D(tSNE_binary, 1)$fhat
  
  Cairo::CairoPDF(paste0("output/Step4.1_tsneModuleActivity_",tsneName,".pdf"), width=15, height=5)
  par(mfrow=c(1,3))
  # nGenes
  plot(tSNE_binary, col=cellColorNgenes[rownames(tSNE_binary)], pch=20, cex = .5)
  # density
  image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
  contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)
  
  # Known phenotype:
  if(!is.null(cellInfo))
  {
    nVars <- ncol(cellInfo)
    for(varName in colnames(cellInfo))
    {
      cellColor <- setNames(colVars[[varName]][as.character(cellInfo[,varName])], 
                            rownames(cellInfo))
      plot(tSNE_binary, col=cellColor[rownames(tSNE_binary)], pch=20, cex = .5, 
           main=varName, sub="t-SNE on Binary regulon activity",
           xlab="", ylab="",axes = FALSE)
    }
  }
  # legend(10, 25, names(colVars[[varName]]), fill=colVars[[varName]], cex=.7, bty="n")
  for(i in seq_len(3 - ((nVars+2) %% 3))) # fill remaining slots in page
  {
    plot.new()
  }
  dev.off()
}



########## Plot t-SNEs coloured by regulon activity ##########
# Choose a t-SNE
load("int/4.1_tsneBinaryActivity_Dist.RData")
tSNE_binary <- tsneBinaryActivity_Dist$Y
tSNEname <- "tsneBinaryActivity_Dist"

# Load...
load("int/3.6_BinaryRegulonActivity.RData")
load("int/3.9_binaryRegulonOrder.RData")
load("int/3.2_regulonAUC.RData")
load("int/3.4_AUCellThresholds.RData")

regOrder<- binaryRegulonOrder[which(binaryRegulonOrder %in% rownames(tBinaryAct))]
Cairo::CairoPDF(paste0("output/Step4.2_",tSNEname,"_BinaryRegulons.pdf"), width=20, height=15)
par(mfrow=c(4,6))
cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat,
                                    regulonAUC=t(tBinaryAct)[binaryRegulonOrder,], cex=1.5, plots="binary", thresholds=0)
dev.off()

regOrder<- binaryRegulonOrder[which(binaryRegulonOrder %in% rownames(regulonAUC))]
Cairo::CairoPDF(paste0("output/Step4.2_",tSNEname,"_AUCRegulons.pdf"), width=20, height=15)
par(mfrow=c(4,6))
cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, 
                                    regulonAUC=regulonAUC[regOrder,], cex=1.5, plots="AUC", thresholds=0)
dev.off()

Cairo::CairoPDF(paste0("output/Step4.2_",tSNEname,"_allPlots.pdf"), width=20, height=5)
par(mfrow=c(1,4))
cells_trhAssignment <- plot_aucTsne(tSNE=tSNE_binary, exprMat=exprMat, 
                                    regulonAUC=regulonAUC[regOrder,],
                                    alphaOff=0.1, thresholds=cells_AUCellThresholds[regOrder])
dev.off()




