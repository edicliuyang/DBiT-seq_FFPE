library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)

dir <- "C:/Users/EDIC-THINKPAD/Desktop/FFPEs/FFPE-2"  
setwd(dir)

#change filename1 to name of txt file you want to load
data1 <- read.table("Filtered_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
data2 <- t(data1)
sample1.name <- "ffpe2"
matrix1.data <- Matrix(as.matrix(data2), sparse = TRUE)
ffpe2          <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)

ffpe2 <- PercentageFeatureSet(ffpe2, pattern = "^MT-", col.name = "percent.mt")
ffpe2 <- SCTransform(ffpe2, vars.to.regress = "percent.mt", verbose = FALSE)
ffpe2 <- RunPCA(ffpe2, verbose = FALSE)
ffpe2 <- RunUMAP(ffpe2, dims = 1:10, verbose = FALSE)

ffpe2 <- FindNeighbors(ffpe2, dims = 1:10, verbose = FALSE)
ffpe2 <- FindClusters(ffpe2, resolution = 0.8, verbose = FALSE)
DimPlot(ffpe2, label = TRUE) + NoLegend()

ffpe2.markers <- FindAllMarkers(ffpe2, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.01) # the min.pct and logfc.threshold need to be tuned based on data.

ffpe2 <-ScaleData(object=ffpe2, features = rownames(ffpe2))

top10 <- ffpe2.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ffpe2, features = top10$gene) + scale_fill_gradientn(colors = c("red", "black","green")) #plot the top 10 DEG for each cluster

write.table(top10, "top10.txt",sep="\t")

go <- ffpe2.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_logFC)
write.table(go, "go.txt",sep="\t") # write top 1000 DEGs for each cluster into a table, for GO analysis


# generate the logFC plot for each cluster
x=ffpe2.markers$cluster
aa = ffpe2.markers$p_val_adj >= 0.01
ffpe2.markers$pp = as.numeric(aa)
ffpe2.markers$pp <- as.factor(ffpe2.markers$pp)
ggplot(ffpe2.markers, aes(x=cluster, y=avg_logFC, color = pp)) + 
  geom_violin(color = NA, fill = NA) + geom_jitter(shape=16, size = 1, position=position_jitter(0.2)) + 
  ylab("average logFC") + xlab("Cluster") +
  geom_hline(yintercept = 0, size = 1) +
  scale_color_manual(values=c("red", "black")) +
  geom_text(aes(label=ifelse(!(avg_logFC < 0.25 & avg_logFC >-0.25) & pp == 0,as.character(gene),'')),position=position_jitter(width=0.5,height=0.3),colour = "black") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text=element_text(colour= "black",size=15),
        axis.title=element_text(colour= "black",size=15,face="bold"),
        legend.text=element_text(colour= "black", size=15),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line.y = element_line(color="black", size = 1),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank())


#visualize marker genes for each cluster
VlnPlot(ffpe2.markers)

deg <- read.table("go.tsv", header = TRUE, sep = "\t")

top5 <- ffpe2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = p_val)
VlnPlot(ffpe2, features = c("Igf2", "Prrx1"))

DotPlot(ffpe2, features = unique(top10$gene),cols = c("blue", "red"), dot.scale = 8,)+ RotatedAxis()



