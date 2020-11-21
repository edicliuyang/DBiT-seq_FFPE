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
library(org.Hs.eg.db)
library(clusterProfiler)
library(Hmisc)


mylist.names <- c("X0", "X1", "X2","X3","X4","X5","X6","X7","X8")
mylist <- vector("list", length(mylist.names))
names(mylist) <- mylist.names

for (i in 0:8) {
  x <- deg$gene[deg$cluster == i]
  gene.df <- clusterProfiler::bitr(toupper(x), fromType = "SYMBOL",
                                   toType = c("ENSEMBL", "ENTREZID"),
                                   OrgDb = org.Hs.eg.db)
  mylist[[i+1]] = gene.df$ENTREZID
}

mylist

ck <- clusterProfiler::compareCluster(geneCluster = mylist, fun = "enrichGO",OrgDb='org.Hs.eg.db')
head(as.data.frame(ck))
clusterProfiler::dotplot(ck)
ck@compareClusterResult$Description = NULL



cc <- gofilter(ck, level = 1)

sc@compareClusterResult$Description

sc <- simplify(ck)

clusterProfiler::dotplot(cc)
ck@compareClusterResult$Description


gene.df

mapIds(org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL')
keys()
