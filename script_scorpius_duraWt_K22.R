library('SCORPIUS')
library('Seurat')

################################################################################
########## LOAD DATA
################################################################################
getwd()

##### Load normalised counts #####
load(file="sample_K22_Dura/Robjects/seuratObj.Robj")
TSNEPlot(seuratObj, do.label=T, label.size = 8, pt.size = 2)

clusterMatrix<-seuratObj@meta.data
head(clusterMatrix)
wantedCells<-rownames(clusterMatrix[clusterMatrix$res.0.8 %in% c(2,11,1,0),])
length(wantedCells)
# 2433

normData<-seuratObj@data[,wantedCells]
dim(normData)
# 14331  2433

rm(seuratObj)
rm(seuratObjDclean)

##### Load meta data #####
sampleInfo<-as.data.frame(cbind(colnames(normData),"cl2"), stringsAsFactors = F)
colnames(sampleInfo)<-c("sampleName","group")
sampleInfo[sampleInfo$sampleName %in% rownames(clusterMatrix[clusterMatrix$res.0.8 == '11',]),2]<-"cl11"
sampleInfo[sampleInfo$sampleName %in% rownames(clusterMatrix[clusterMatrix$res.0.8 == '1',]),2]<-"cl1"
sampleInfo[sampleInfo$sampleName %in% rownames(clusterMatrix[clusterMatrix$res.0.8 == '0',]),2]<-"cl0"
sampleInfo$group<-as.factor(sampleInfo$group)
dim(sampleInfo)
# 2433    2
table(sampleInfo$group)
# cl0  cl1 cl11  cl2 
# 1001  754  184  494 
table(clusterMatrix$res.0.8)

##### Preparation for Scorpius #####
expression<-t(normData)
dim(expression)
# 2433 14331
group_name <- sampleInfo$group
length(group_name)
# 2433

################################################################################
########## REDUCE DIMENSIONALITY
################################################################################

dist <- correlation_distance(expression)
dim(dist)
# 1432 1432

plot(density(dist))

space <- reduce_dimensionality(dist)
dim(space)
# 1432    3

draw_trajectory_plot(space)


## Save manually as '1_dimensionalityReduction.png'
draw_trajectory_plot(space, progression_group = group_name)

################################################################################
########## OUTLIER FILTERING
################################################################################

library("ggplot2")
pdf(file="sample_K22_Dura/results/scorpius/1_dimensionalityReduction.pdf")
draw_trajectory_plot(space[, c(1, 3)]) + labs(y = "Component 3")
draw_trajectory_plot(space[, c(1, 3)], progression_group = group_name) + labs(y = "Component 3")

draw_trajectory_plot(space[, c(2, 3)]) + labs(x = "Component 2", y = "Component 3")
draw_trajectory_plot(space[, c(2, 3)], progression_group = group_name) + labs(x = "Component 2", y = "Component 3")
dev.off()

### Filter out outliers
filt <- outlier_filter(dist)
expression <- expression[filt, ]
dim(expression)
# 1424 14331
# 8 cells removed
saveRDS(filt, file="sample_K22_Dura/Robjects/scorpius/filt_cl2-11-1.rds")

removedCells<-setdiff(colnames(normData),rownames(expression))

group_name <- group_name[filt]
length(group_name)
# 2430

### Show the removed cells
tmpSpace<-as.data.frame(space[,c(1,2)], stringsAsFactors = F)
tmpSpace$isGood<-"yes"
tmpSpace[removedCells,3]<-"no"

p <- ggplot()+
  geom_point(aes(x=Comp1,y=Comp2, colour=tmpSpace$isGood), data=tmpSpace, size=3) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
print(p)
ggsave(p, file="sample_K22_Dura/results/scorpius/2_removedCell.png")

################################################################################
########## REDO REDUCE DIMENSIONALITY
################################################################################

dist <- dist[filt, filt]
dim(dist)
# 2430 2430
space <- reduce_dimensionality(dist)
dim(space)
# 2430    3

saveRDS(dist, file="sample_K22_Dura/Robjects/scorpius/dist_cl2-11-1.rds")
saveRDS(space, file="sample_K22_Dura/Robjects/scorpius/space_cl2-11-1.rds")

pdf(file="sample_K22_Dura/results/scorpius/3_dimensionalityReductionAfterFiltering.pdf")
draw_trajectory_plot(space[, c(1, 2)])
draw_trajectory_plot(space, progression_group = group_name)

draw_trajectory_plot(space[, c(1, 3)]) + labs(y = "Component 3")
draw_trajectory_plot(space[, c(1, 3)], progression_group = group_name) + labs(y = "Component 3")

draw_trajectory_plot(space[, c(2, 3)]) + labs(x = "Component 2", y = "Component 3")
draw_trajectory_plot(space[, c(2, 3)], progression_group = group_name) + labs(x = "Component 2", y = "Component 3")
dev.off()



##Use the colors from the heatmap
draw_trajectory_plot(space[, c(1, 3)], progression_group = group_name) + labs(y = "Component 3") +
  scale_color_manual(values=c("cl0"="#e41a1c", "cl1"="#377eb8","cl11"="#4daf4a","cl2"="#984ea3"))

################################################################################
########## INFERRING TRAJECTORY
################################################################################

traj <- infer_trajectory(space)

png(file="sample_K22_Dura/results/scorpius/4_findTrajectory.png", width=850)
draw_trajectory_plot(space, progression_group = group_name, path = traj$path)
dev.off()


### Based on other PCs of PCA plot
traj <- infer_trajectory(space[,c(1,3)])
draw_trajectory_plot(space[,c(1,3)], progression_group = group_name, path = traj$path)


##Use the colors from the heatmap + put line thicker
draw_trajectory_plot(space[,c(1,3)], progression_group = group_name) +
  scale_color_manual(values=c("cl0"="#e41a1c", "cl1"="#377eb8","cl11"="#4daf4a","cl2"="#984ea3")) +
  geom_path(aes(Comp1, Comp2), data.frame(traj$path), size=1.5)


################################################################################
########## FIND MARKER GENES
################################################################################

gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:100,]
expr_sel <- expression[,gene_sel$gene]

traj <- infer_trajectory(expr_sel)

# Save manually as pdf '5a_heatmapMarkers.pdf'
draw_trajectory_heatmap(expr_sel, traj$time, group_name, show_labels_row = T)


# Save manually as pdf '5b_heatmapMarkersWithModules.pdf'
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = F)
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = T)






