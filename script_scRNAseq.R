library('Seurat')
library('dplyr')
library("scater")
library("scran")
library('gridExtra')
library('rlist')

################################################################################
########## LOAD DATA
###############################################################################
getwd()

dataSample <- Read10X("filtered_gene_bc_matrices/mm10/")
rawData<-as.matrix(dataSample)
rm(dataSample)

diagnostics<-list()

########################################
##### Characteristics of data
########################################
rawData[1:5,1:5]

dim(rawData)

zeroInflated<-sum(rawData==0)/(nrow(rawData)*ncol(rawData))*100
zeroInflated

nrCells<-apply(rawData,1,function (x){sum(x>0)})
length(nrCells[nrCells<3])

##### Add to diagnostics #####
diagnostics[['nrGenes']]<-nrow(rawData)
diagnostics[['nrCells']]<-ncol(rawData)
diagnostics[['zeroInflation']]<-zeroInflated

diagnostics[['minGenesPerCell']]<-min(nrCells)
diagnostics[['maxGenesPerCell']]<-max(nrCells)
diagnostics[['meanGenesPerCell']]<-mean(nrCells)
diagnostics[['medianGenesPerCell']]<-median(nrCells)


################################################################################
########## QC: CELLS
################################################################################

########################################
########## Calculate metrics
########################################

##### Create object #####
library("scater")
sce<-SingleCellExperiment(list(counts=rawData))
dim(sce)

##### Get spike inns #####
is.spike <- grepl("^ERCC", rownames(sce))
sum(is.spike)

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
dim(colData(sce))
colnames(colData(sce))

##### Create metaData matrix #####
metaData<-data.frame("staticNr"=colnames(rawData),"orig.ident"=label_1,
                     "nGene"=sce$total_features_by_counts,"nUMI"=sce$total_counts,"percent.mito"=sce$pct_counts_Mt, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

metaData[grep('-2',rownames(metaData)), which(colnames(metaData)=="orig.ident")]<-label_2

table(metaData$orig.ident)

########################################
########## Get outliers
########################################

##### CHANGE PARAMETERS HERE! #####
### Change nmad
###################################

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##nGene
feature.drop.low <- isOutlier(sce$total_features_by_counts, nmads=4, type="lower", log=TRUE)
sum(feature.drop.low)

feature.drop.high <- isOutlier(sce$total_features_by_counts, nmads=4, type="higher", log=TRUE)
sum(feature.drop.high)

feature.drop<-as.logical(feature.drop.low + feature.drop.high)
sum(feature.drop)

##nUMI
libsize.drop.low <- isOutlier(sce$total_counts, nmads=4, type="lower", log=TRUE)
sum(libsize.drop.low)

libsize.drop.high <- isOutlier(sce$total_counts, nmads=4, type="higher", log=TRUE)
sum(libsize.drop.high)

libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
sum(libsize.drop)

##% mitochondrial genes
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=4, type="higher")
sum(mito.drop)


##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop

##### Add to diagnostics #####
diagnostics[['feature.drop']]<-sum(feature.drop)
diagnostics[['libsize.drop']]<-sum(libsize.drop)
diagnostics[['mito.drop']]<-sum(mito.drop)

########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###F8766D=red
###00BFC4=cyan
###7CAE00=green
###C77CFF=purple

toPlot<-metaData
##nGene
png(file="results/plots/1a_nGene.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
dev.off()

##nUMI
png(file="results/plots/1b_nUMI.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
dev.off()

##percent.mito
png(file="results/plots/1c_percMito.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
dev.off()

########################################
########## Create violinPlots
########################################

### Before filtering
toPlot<-metaData
drawVlnPlot(toPlot, "results/plots/2a_beforeFiltering.png",c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, "results/plots/2a_beforeFiltering_splitted.png",c('nGene.drop','nUMI.drop','mito.drop'))


### After filtering
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, "results/plots/2b_afterFiltering.png",c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, "results/plots/2b_afterFiltering_splitted.png",c('nGene.drop','nUMI.drop','mito.drop'))


########################################
########## Remove outliers
########################################

sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
dim(sce)

### Number of cells removed
nrow(metaData)-ncol(sce)

##### Add to diagnostics #####
diagnostics[['firstRemove']]<-nrow(metaData)-ncol(sce)



########################################
########## Create PCA
########################################
library('mvoutlier')
selected_variables <- c("pct_counts_in_top_100_features", 
                        "total_features_by_counts", "pct_counts_feature_control", 
                        "total_features_by_counts_feature_control", "log10_total_counts_endogenous", 
                        "log10_total_counts_feature_control")
setdiff(selected_variables, colnames(colData(sce)))

varsToUse<-selected_variables
setdiff(varsToUse, colnames(colData(sce)))

##### Detect bad cells #####
sceNew<-runPCA(sce,use_coldata=T, detect_outliers=T, selected_variables=varsToUse)
table(sceNew$outlier)
outs<-colnames(sceNew)[sceNew$outlier]
### Add to metaData
metaData$pca.drop<-metaData$final.drop
metaData[outs,which(colnames(metaData)=="pca.drop")]<-TRUE

##### Color bad cells on PCA plot #####
colorDF<-as.data.frame(cbind(colnames(sceNew),"1"), stringsAsFactors=F)
rownames(colorDF)<-colorDF[,1]
colorDF[outs,2]<-"2"
colorDF[,2]<-as.factor(colorDF[,2])
tmp<-colorDF[,2,drop=F]

png(file="results/plots/3a_pca.png",  width = 850, height = 642)
plotReducedDim(sceNew, use_dimred = "PCA_coldata", colour_by='outlier',shape_by='outlier') + labs(title="PCA with outliers colored")
dev.off()


#### Remove the bad cells based on the PCA plot ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
sum(pca.drop)

##### Create violinplots ####
##Before
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, "results/plots/3b_beforePcaFiltering_2.png",c(rep('pca.drop',3)))

##After
toPlot<-metaData[! metaData$pca.drop,]
drawVlnPlot(toPlot, "results/plots/3c_afterPcaFiltering.png",c('nGene.drop','nUMI.drop','mito.drop'))


##### Remove outlier cells ####
sce <- sce[,!(pca.drop)]
dim(sce)

rm(sceNew)

##### Add to diagnostics #####
diagnostics[['pcaRemove']]<-sum(pca.drop)

################################################################################
########## QC: GENES
################################################################################

ave.counts <- calcAverage(sce)
png(file="results/plots/4_QCgenes.png", width = 850, height = 642)
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=-3,col="blue",lwd=2,lty=2)
dev.off()


##### CHANGE PARAMETERS HERE! #####
### Change cutoff
###################################

cutoff<-10^-3
demo.keep <- ave.counts >= cutoff
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
sce<-filtered.sce

rm(filtered.sce)
rm(ave.counts)
rm(demo.keep)

saveRDS(sce, file="Robjects/sceFiltered.rds")


################################################################################
########## NORMALIZATION
################################################################################
library('limSolve')
dim(sce)

##### Method 2 #####
# While the deconvolution approach is robust to the high frequency of zeroes in scRNA-seq data, it will eventually fail if too many counts are zero. 
# This manifests as negative size factors, which are obviously nonsensical. To avoid this, the computeSumFactors function will automatically remove 
# low-abundance genes prior to the calculation of size factors. Genes with an average count below a specified threshold (min.mean) are ignored. 
# For read count data, the default value of 1 is usually satisfactory. For UMI data, counts are lower so a threshold of 0.1 is recommended.
sce <- computeSumFactors(sce, min.mean=0.1)
summary(sizeFactors(sce))

png(file="results/plots/5_normalization.png", width = 850, height = 642)
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", ylab="Library size (millions)", xlab="Size factor")
dev.off()


sce <- normalize(sce)
mat<-exprs(sce)

saveRDS(sce, file="Robjects/sce.rds")
saveRDS(mat, file="Robjects/mat.rds")

################################################################################
########## HVG DETECTION
################################################################################

var.fit<-trendVar(sce, parametric=TRUE, use.spikes=FALSE)
var.out<-decomposeVar(sce, var.fit)
dim(var.out)

hvg.out<-var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.01),]
hvg.out<-hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)


########################################
########## Get IEG genes
########################################

DEgenes_exp1<-readRDS(file="JP30-JP31/Robjects/listDEgenes.rds")
DEgenes_exp2<-readRDS(file="JP32-JP33/Robjects/listDEgenes.rds")

tmp<-DEgenes_exp1$JP31_vs_JP30
DEgenes_exp1<-as.character(tmp[tmp$avg_logFC > 0,which(colnames(tmp)=="gene")])

tmp<-DEgenes_exp2$JP33_vs_JP32
DEgenes_exp2<-as.character(tmp[tmp$avg_logFC > 0,which(colnames(tmp)=="gene")])

dissociationGenes<-intersect(DEgenes_exp1,DEgenes_exp2)
length(dissociationGenes)
## 224

##### Remove IEG genes #####
hvgGenes<-setdiff(rownames(hvg.out),dissociationGenes)
length(hvgGenes)

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

dim(mat)

##### Create object #####
seuratObj <- CreateSeuratObject(raw.data = mat, min.cells = 3, min.genes = 200, project = "seuratObj")
dim(seuratObj@raw.data)

##### Add to diagnostics #####
diagnostics[['nrGenes_afterFiltering']]<-nrow(seuratObj@raw.data)
diagnostics[['nrCells_afterFiltering']]<-ncol(seuratObj@raw.data)

##### Fill var.genes (HVG genes) #####
seuratObj@var.genes<-rownames(hvg.out)
length(seuratObj@var.genes)
##or if correcting for IEG
seuratObj@var.genes<-hvgGenes
length(seuratObj@var.genes)

##### Fill norm.data #####
seuratObj@data<-mat


################################################################################
########## FILTER DATA
################################################################################

mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")

head(seuratObj@meta.data)

png(file="results/plots/6a_afterNorm.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
metaDataTable[grep('-1',rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-label_1
metaDataTable[grep('-2',rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-label_2
metaDataTable$orig.ident<-factor(metaDataTable$orig.ident)
seuratObj@meta.data<-metaDataTable


################################################################################
########## NORMALIZE
################################################################################


##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj@data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, "results/plots/6b_afterNorm_splitted.png")

################################################################################
########## SCALE DATA
################################################################################
seuratObj <- ScaleData(object = seuratObj)

##### Check per group #####
head(metaDataTable)
metaDataTable$nUMI<-colSums(as.matrix(seuratObj@scale.data))
metaDataTable$nGene<-apply(as.matrix(seuratObj@scale.data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, "results/plots/7_afterScale_splitted.png")


################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 10, pcs.compute = 40)

########################################
########## PCA PLOT
########################################
pdf(file="results/plots/8a_PCA.pdf", width = 10)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 3)
dev.off()

########################################
########## HEATMAP OF PCs
########################################
pdf(file="results/plots/9a_selectPC.pdf")
PCHeatmap(seuratObj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 37:40, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

# PCHeatmap(seuratObj, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

# seuratObj <- JackStraw(seuratObj, num.replicate = 100, num.pc = 40)
# pdf(file="results/plots/9c_jackStrawPlot.pdf")
# JackStrawPlot(seuratObj, PCs = 1:40)
# dev.off()

png(file="results/plots/9b_selectPC.png", width = 850, height = 642)
PCElbowPlot(seuratObj, num.pc = 40)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################

##### CHANGE PARAMETERS HERE! #####
### Change listPCs
###################################

##Final PC
listPCs<-c(17)

for(PCs in listPCs){
  print(paste0("WORKING ON PCS ", PCs))
  seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:PCs, resolution = 0.8, print.output = 0, save.SNN = F)
  
  fileName<-paste0("results/plots/10_tSNE_1-",PCs,"_res08.png")
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:PCs, do.fast = TRUE)
  p1<-TSNEPlot(seuratObj, do.label=T, label.size = 7, pt.size = 1)
  
  ##### Split plot
  clusterMatrix<-seuratObj@meta.data
  tsneTable<-as.data.frame(seuratObj@dr$tsne@cell.embeddings, stringsAsFactors = F)
  p2 <- ggplot()+
    geom_point(aes(x=tSNE_1,y=tSNE_2, colour=clusterMatrix$orig.ident), data=tsneTable, size=1) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.text=element_text(size=12))
  print(p2)
  ggsave(grid.arrange(p1, p2,ncol=2),file=fileName, dpi=200, width = 25, height = 10)
}

##### Save object
saveRDS(seuratObj, file="Robjects/seuratObj.rds")


