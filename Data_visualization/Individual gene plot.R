library(ggplot2)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(dplyr)

dir <- "C:/Users/EDIC-THINKPAD/Desktop/replot/fig2/gene"  
setwd(dir)


##read expression matrix
data_filtered <- read.table(file = 'cL50.tsv', sep = '\t', header = TRUE, row.names = 1)
data2 <- t(data_filtered)
sample1.name <- "mouse"
matrix1.data <- Matrix(as.matrix(data2))
pbmc          <- CreateSeuratObject(matrix1.data, min.cells = 10, project = sample1.name)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

genedata <- pbmc[["RNA"]]@data
genedata <- t(genedata)
gene <-as.data.frame(as.matrix(genedata))
gene$X = row.names(gene)
gene$X
test <- gene %>% separate(X, c("A", "B"),  sep = "x")



#UMI heatmap
pdf(file = paste("Epcam.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=test$Epcam)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),
                        oob = scales::squish) +
  ggtitle("Epcam") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 40, face = "bold"),
        axis.text=element_text(colour="black",size=30),
        axis.title=element_text(colour="black",size=30,face="bold"),
        legend.text=element_text(colour="black",size=30),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(size =1, fill = NA))
dev.off()



