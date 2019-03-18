library('Seurat')
library('dplyr')
library('Matrix')
library("scater")
library("scran")
library('gridExtra')

################################################################################
########## LOAD DATA
################################################################################
getwd()

dataK11 <- Read10X("filtered_gene_bc_matrices/mm10/")
rawData<-as.matrix(dataK11)
rm(dataK11)

########################################
##### Characteristics of data
########################################
rawData[1:5,1:5]

dim(rawData)

sum(rawData==0)

nrCells<-apply(rawData,1,function (x){sum(x>0)})
min(nrCells)
max(nrCells)
length(nrCells[nrCells<3])
##12775
##27998-12775 = 15223 genes to keep
length(nrCells[nrCells==0])
##10816


nrGenes<-apply(rawData,2,function (x){sum(x>0)})
min(nrGenes)
max(nrGenes)
##min: 85
##max: 6632
length(nrGenes[nrGenes<200])
##7
## 4055-7= 4048 cells to keep

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
# 27998  4055

##### Get spike inns #####
is.spike <- grepl("^ERCC", rownames(sce))
sum(is.spike)
##0

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)
##13

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
dim(colData(sce))
# 4055   28
colnames(colData(sce))

##### Create metaData matrix #####
metaData<-data.frame("staticNr"=colnames(rawData),"nGene"=sce$total_features,"nUMI"=sce$total_counts,"percent.mito"=sce$pct_counts_Mt, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

########################################
########## Get outliers
########################################

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop <- isOutlier(sce$total_features, nmads=4, type="lower", log=TRUE)
sum(feature.drop)
#15

##same as UMI in Seurat pipeline
libsize.drop <- isOutlier(sce$total_counts, nmads=4, type="lower", log=TRUE)
sum(libsize.drop)
#0

mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=4, type="higher")
sum(mito.drop)
#114

##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop

########################################
########## Remove outliers
########################################

sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
dim(sce)
# 27998  3940
###4055-3940 = 115 cells removed


########################################
########## Create violinPlots
########################################

toPlot<-metaData
toPlot<-metaData[! metaData$final.drop,]
toPlot<-metaData[! metaData$pca.drop,]

p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nGene.drop)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.drop)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.drop)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))

grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="results/plots/1a_beforeFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="results/plots/2_afterFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="results/plots/3b_beforePcaFiltering.png")
ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file="results/plots/4_afterPcaFiltering.png")


########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))
###F8766D=red
###00BFC4=cyan
###7CAE00=green
###C77CFF=purple

##nGene
png(file="results/plots/1b_nGene.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
dev.off()

##nUMI
png(file="results/plots/1v_nUMI.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
dev.off()

##percent.mito
png(file="results/plots/1d_percMito.png", width=850)
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
dev.off()


########################################
########## Create PCA
########################################
colnames(colData(sce))
##on position 4, 22 and 18
colnames(colData(sce))[colnames(colData(sce))=="log10_total_counts"]<-"log10_total_counts_feature_controls"
colnames(colData(sce))[colnames(colData(sce))=="pct_counts_feature_control"]<-"pct_counts_feature_controls"
colnames(colData(sce))[colnames(colData(sce))=="total_features_feature_control"]<-"total_features_feature_controls"

png(file="results/plots/3a_pca.png",  width = 850, height = 642)
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(sce, pca_data_input="pdata", detect_outliers=TRUE) + fontsize
dev.off()

##### Detect bad cells #####
sceNew<-runPCA(sce, pca_data_input = "pdata", detect_outliers=TRUE)
table(sceNew$outlier)
# FALSE  TRUE 
# 3896    44 
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

scater::plotPCA(sceNew, pca_data_input="pdata", return_SCESet = FALSE, colour_by=tmp, shape_by=tmp, theme_size=14)

##### Create violin plots again ####
##change color of geom_jitter() into 'pca.drop'

##### Remove the bad cells based on the PCA plot ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
sum(pca.drop)
##44

sce <- sce[,!(pca.drop)]
dim(sce)
## 27998  3896


rm(sceNew)

################################################################################
########## QC: GENES
################################################################################

ave.counts <- calcAverage(sce)
png(file="results/plots/5_QCgenes.png", width = 850, height = 642)
hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=-2.7,col="blue",lwd=2,lty=2)
dev.off()

cutoff<-10^-2.7
demo.keep <- ave.counts >= cutoff
filtered.sce <- sce[demo.keep,]
summary(demo.keep)
# Mode   FALSE    TRUE 
# logical   14816   13182 
sce<-filtered.sce

rm(filtered.sce)
rm(ave.counts)
rm(demo.keep)

########## Keep cells with at least 200 genes
to.keep<-which(sce@colData$total_counts>=200)
length(to.keep)
# 3896 = all the cells


################################################################################
########## NORMALIZATION
################################################################################
library('limSolve')
dim(sce)
# 13182  3896

##### Method 2 #####
# While the deconvolution approach is robust to the high frequency of zeroes in scRNA-seq data, it will eventually fail if too many counts are zero. 
# This manifests as negative size factors, which are obviously nonsensical. To avoid this, the computeSumFactors function will automatically remove 
# low-abundance genes prior to the calculation of size factors. Genes with an average count below a specified threshold (min.mean) are ignored. 
# For read count data, the default value of 1 is usually satisfactory. For UMI data, counts are lower so a threshold of 0.1 is recommended.
sce <- computeSumFactors(sce, min.mean=0.1)
summary(sizeFactors(sce))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2902  0.6639  0.8465  1.0000  1.1756  3.4444 


png(file="results/plots/6_normalization.png", width = 850, height = 642)
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy", ylab="Library size (millions)", xlab="Size factor")
dev.off()


sce <- normalize(sce)
mat<-exprs(sce)

saveRDS(sce, file="Robjects/sce.rds")
saveRDS(mat, file="Robjects/mat.rds")

################################################################################
########## HVG DETECTION
################################################################################

var.fit<-trendVar(sce, parametric=TRUE, use.spikes=FALSE, span=0.2)
var.out<-decomposeVar(sce, var.fit)
dim(var.out)
# 13182     6

hvg.out<-var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.01),]
hvg.out<-hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)
##1149


hist(hvg.out$bio[hvg.out$bio < 0.60], breaks = 30)
hist(hvg.out$bio[hvg.out$bio < 0.20], breaks = 30)


################################################################################
########## CREATE SEURAT OBJECT
################################################################################

dim(mat)
# 13182  3896

##### Create object #####
seuratObj <- CreateSeuratObject(raw.data = mat, min.cells = 3, min.genes = 0, project = "seuratObj")
dim(seuratObj@raw.data)
# 13177  3896
### 5 extra genes removed

##### Fill var.genes (HVG genes) #####
seuratObj@var.genes<-rownames(hvg.out)
length(seuratObj@var.genes)
##1149

##### Fill scaled.data (normalized data) #####
seuratObj@scale.data<-mat
seuratObj@data<-mat

################################################################################
########## FILTER DATA
################################################################################

mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")

head(seuratObj@meta.data)

png(file="results/plots/7_afterNormalization.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

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

PCHeatmap(seuratObj, pc.use = 2, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

# seuratObj <- JackStraw(seuratObj, num.replicate = 100, do.print = FALSE, num.pc = 30)
# JackStrawPlot(seuratObj, PCs = 1:30)

png(file="results/plots/9b_selectPC.png", width = 850, height = 642)
PCElbowPlot(seuratObj, num.pc = 30)
dev.off()


################################################################################
########## CLUSTER THE CELLS
################################################################################

seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:20, resolution = 0.8, print.output = 0, save.SNN = F)

##### Create tSNE plot
seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:20, do.fast = TRUE)
png(file="results/plots/10_tSNE_1-20_res08.png", width = 850, height = 642)
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)
dev.off()

##### Save object
saveRDS(seuratObj, file="Robjects/seuratObj.rds")


################################################################################
########## DETECT THE CELLS AFFECTED BY DISSOCIATION
################################################################################
dissocationGenes<-c("Actg1","Ankrd1","Arid5a","Atf3","Atf4","Bag3","Bhlhe40","Brd2","Btg1","Btg2","Ccnl1","Ccrn4l","Cebpb","Cebpd",
                    "Cebpg","Csrnp1","Cxcl1","Cyr61","Dcn","Ddx3x","Ddx5","Des","Dnaja1","Dnajb1","Dnajb4","Dusp1","Dusp8","Egr1",
                    "Egr2","Eif1","Eif5","Erf","Errfi1","Fam132b","Fos","Fosb","Fosl2","Gadd45a","Gcc1","Gem","H3f3b","Hipk3","Hsp90aa1",
                    "Hsp90ab1","Hspa1a","Hspa1b","Hspa5","Hspa8","Hspb1","Hsph1","Id3","Idi1","Ier2","Ier3","Ifrd1","Il6","Irf1","Irf8",
                    "Itpkc","Jun","Junb","Jund","Klf2","Klf4","Klf6","Klf9","Litaf","Lmna","Maff","Mafk","Mcl1","Midn","Mir22hg","Mt1",
                    "Mt2","Myadm","Myc","Myd88","Nckap5l","Ncoa7","Nfkbia","Nfkbiz","Nop58","Nppc","Nr4a1","Odc1","Osgin1","Oxnad1",
                    "Pcf11","Pde4b","Per1","Phlda1","Pnp","Pnrc1","Ppp1cc","Ppp1r15a","Pxdc1","Rap1b","Rassf1","Rhob","Rhoh","Ripk1",
                    "Sat1","Sbno2","Sdc4","Serpine1","Skil","Slc10a6","Slc38a2","Slc41a1","Socs3","Sqstm1","Srf","Srsf5","Srsf7","Stat3",
                    "Tagln2","Tiparp","Tnfaip3","Tnfaip6","Tpm3","Tppp3","Tra2a","Tra2b","Trib1","Tubb4b","Tubb6","Ubc","Usp2","Wac",
                    "Zc3h12a","Zfand5","Zfp36","Zfp36l1","Zfp36l2","Zyx","Gadd45g","Hspe1","Ier5","Kcne4")
setdiff(dissocationGenes, rownames(seuratObj@scale.data))

normData<-seuratObj@scale.data
dissocationGenes<-intersect(dissocationGenes, rownames(normData))
setdiff(dissocationGenes, rownames(normData))

########## Start the code from the paper ##########
Selection<-normData[dissocationGenes,]
Selection[is.na(Selection)]<-0

#Generate  a  new dataframe  (DataPercentages)  containing  for  each  cell:  t-SNE  dimensions  ("V1"  and "V2"), 
#sum of reads from all dissociation-affected genes ("Sums") and percentage of transcriptome of that cell that maps to 
#dissociation affected reads ("Percentage"; this equals "Sums"-column divided by total read count of that cell):
tSNECoordinates <-tsneTable
head(cbind(rownames(tSNECoordinates),colnames(normData)))
tail(cbind(rownames(tSNECoordinates),colnames(normData)))

Sums <-colSums(Selection)
TotalSums <-colSums(normData)
DataPercentages <-merge(tSNECoordinates, Sums, by="row.names", all=T)
row.names(DataPercentages) <-DataPercentages$Row.names
DataPercentages <-DataPercentages[,-1]
colnames(DataPercentages)[3] <-"Sums"
DataPercentages$Percentage <-DataPercentages$Sums*100/TotalSums

#Make a histogram showing the distribution of the metric "Percentage of dissociation-affected reads per cell":
png(file="results/plots/13a_histogramDissociationAffectedCells.png")
hist(DataPercentages$Percentage, breaks = 100, col = "lightgrey", main = "Expression level dissociation-affected genes", 
     xlab = "Sum expression level of dissociation-affected genes", ylab = "Number of cells")
abline(v=3.9,col="blue",lwd=2,lty=2)
dev.off()

#Based  on  this  histogram,  choose  an  appropriate  cutoff  value.  All  cells  with  a  percentage  equal  to  or above 
#this value are affected by the dissociation procedure, will be labelled as "Dissociation-affected" in the "DataPercentages" dataframe:
SetCutoff <-3.9
DataPercentages$Dissociation_affected <-ifelse(DataPercentages$Percentage >= SetCutoff,1,0)

#Calculate  what percentage of the cells  will be annotated as dissociation-affected cells  with this cutoff value:
PercentageDissociationAffected   <-round((nrow(DataPercentages[DataPercentages$Percentage   >= SetCutoff, ]))/(nrow(DataPercentages))*100, digits = 2)
print(c("Percentage of cells annotated as dissociation-affected cells is", PercentageDissociationAffected))
##62.78

cellsAffected<-rownames(DataPercentages[DataPercentages$Dissociation_affected == "1",])
length(cellsAffected)
# 2446
saveRDS(cellsAffected, file="Robjects/cellsAffected.rds")


########## Color the affected cells on the tSNE ##########
###F8766D=red
###00BFC4=cyan
tmp<-clusterMatrix
tmp$affected<-"FALSE"
tmp[cellsAffected,which(colnames(tmp)=="affected")]<-"TRUE"

p <- ggplot()+
  geom_point(aes(x=tSNE_1,y=tSNE_2, colour=tmp$affected), data=tsneTable, size=2, shape=20) +
  scale_color_manual(values=c("FALSE"="#00BFC4", "TRUE"="#F8766D")) + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
print(p)
ggsave(p, file="results/plots/13b_colorDissociationAffectedCells.png")


################################################################################
########## 2. CREATE SEURAT OBJECT
################################################################################

dim(mat)
# 13182  3896

rm(seuratObj)

##### Create object #####
seuratObj <- CreateSeuratObject(raw.data = mat, min.cells = 3, min.genes = 0, project = "seuratObj")
dim(seuratObj@raw.data)
# 13177  3896
### 5 extra genes removed

##### Fill var.genes (HVG genes) #####
hvgGenes<-setdiff(rownames(hvg.out),dissocationGenes)

seuratObj@var.genes<-hvgGenes
length(seuratObj@var.genes)
##1056

##### Fill scaled.data (normalized data) #####
seuratObj@scale.data<-mat

################################################################################
########## FILTER DATA
################################################################################

mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuratObj@data), value = TRUE)
percent.mito <- Matrix::colSums(seuratObj@raw.data[mito.genes, ])/Matrix::colSums(seuratObj@raw.data)
seuratObj <- AddMetaData(object = seuratObj, metadata = percent.mito, col.name = "percent.mito")

head(seuratObj@meta.data)

png(file="results/plots/clean_tSNE/7_afterNormalization.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features.plot = c("nGene", "nUMI","percent.mito"), nCol = 3)
dev.off()

################################################################################
########## REGRESS OUT
################################################################################

##### Preparation #####
cellsAffected<-readRDS(file="Robjects/cellsAffected.rds")
length(cellsAffected)
# 2446
tmp<-seuratObj@meta.data
tmp$affected<-FALSE
tmp[cellsAffected, which(colnames(tmp)=="affected")]<-TRUE

toAdd<-tmp$affected
names(toAdd)<-rownames(tmp)

seuratObj <- AddMetaData(object = seuratObj, metadata = toAdd, col.name = "affected")
head(seuratObj@meta.data)
table(seuratObj@meta.data$affected)
# FALSE  TRUE 
# 1450  2446

##### Regress out #####
seuratObj <- ScaleData(object = seuratObj, vars.to.regress = c("affected"), do.scale = FALSE, do.center = FALSE)

################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, pc.genes = seuratObj@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 10, pcs.compute = 40)

PrintPCA(object = seuratObj, pcs.print = 1:7, genes.print = 10, use.full = FALSE)

##### Check presence of dissociaton genes
genesPCA<-seuratObj@dr$pca@gene.loadings
topGenes<-c()
for(i in 1:20){
  testGenes<-genesPCA[,i]
  testGenesOrdered<-testGenes[order(testGenes, decreasing=T)]
  topGenes<-c(topGenes, head(testGenesOrdered,25), tail(testGenesOrdered,25))
}
intersect(names(topGenes), dissocationGenes)

########################################
########## PCA PLOT
########################################
pdf(file="results/plots/clean_tSNE/8a_PCA.pdf", width = 10)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 2)
PCAPlot(object = seuratObj, dim.1 = 1, dim.2 = 3)
dev.off()

########################################
########## HEATMAP OF PCs
########################################
pdf(file="results/plots/clean_tSNE/9a_selectPC.pdf")
PCHeatmap(seuratObj, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 13:24, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 25:36, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
PCHeatmap(seuratObj, pc.use = 37:40, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE) 
dev.off()

PCHeatmap(seuratObj, pc.use = 5, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

# seuratObj <- JackStraw(seuratObj, num.replicate = 100, do.print = FALSE, num.pc = 30)
# JackStrawPlot(seuratObj, PCs = 1:30)

png(file="results/plots/clean_tSNE/9b_selectPC.png", width = 850, height = 642)
PCElbowPlot(seuratObj, num.pc = 40)
dev.off()

################################################################################
########## CLUSTER THE CELLS
################################################################################

seuratObj <- FindClusters(object = seuratObj, reduction.type = "pca", dims.use = 1:20, resolution = 0.8, print.output = 0, save.SNN = F)

##### Create tSNE plot
seuratObj <- RunTSNE(object = seuratObj, dims.use = 1:20, do.fast = TRUE)
png(file="results/plots/clean_tSNE/10_tSNE_1-20_res08.png", width = 850, height = 642)
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)
dev.off()


##### Save object
saveRDS(seuratObj, file="Robjects/seuratObj_cleanRegressOut.rds")


