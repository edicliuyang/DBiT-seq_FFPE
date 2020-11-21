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


dir <- "C:/Users/EDIC-THINKPAD/Desktop/"  
setwd("C:/Users/EDIC-THINKPAD/Desktop/")

data1 <- read.table("go1.tsv", header = FALSE, sep = "\t")
levels(data1$V1)
data1$V6 <- as.factor(data1$V6)
temp <- data1$V1
levels(temp)
data1$V1 <- factor(data1$V1,levels = c("forebrain development (GO:0030900)",
                                       "telencephalon development (GO:0021537)", 
                                       "central nervous system development (GO:0007417)",
                                       "eye morphogenesis (GO:0048592)",
                                       "limb morphogenesis (GO:0035108)" ,
                                       "gland development (GO:0048732)",
                                       "muscle structure development (GO:0061061)",
                                       "nervous system development (GO:0007399)" ,
                                       "artery development (GO:0060840)" ,
                                       "mesenchyme development (GO:0060485)" ,
                                       "epithelium development (GO:0060429)" ,
                                       "skeletal system development (GO:0001501)",
                                       "generation of neurons (GO:0048699)",
                                       "neurogenesis (GO:0022008)",
                                       "regulation of mesenchymal cell proliferation (GO:0010464)",
                                       "neuron development (GO:0048666)",
                                       "blood vessel development (GO:0001568)",
                                       "vasculature development (GO:0001944)",
                                       "tube development (GO:0035295)",
                                       "circulatory system development (GO:0072359)" ,
                                       "ganglion development (GO:0061548)",
                                       "olfactory lobe development (GO:0021988)" ))

S1 <- ggplot(data1, aes(x = V6, y = V1 , color = V5, size = V3)) + 
  geom_point() +xlab("Cluster") + 
  ylab(NULL) + 
  labs(colour = "P.adjust", size = "GeneRatio") +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=12,colour ="black"),
        axis.title=element_text(size=14,face="bold",colour ="black"),
        legend.text=element_text(size=12),
        legend.title = element_text(colour="black", size=14, face="bold"),
        axis.line = element_line(colour = "black", size =1),
        panel.grid.major = element_line(colour = "lightgray", size =1.5), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

S1 = S1+scale_color_gradient(low = "red2",  high = "mediumblue", space = "Lab", limit = c(0.000000000000000000000000000000000000000000000004, 0.03))
S1+scale_size(range = c(2, 8))
S1
