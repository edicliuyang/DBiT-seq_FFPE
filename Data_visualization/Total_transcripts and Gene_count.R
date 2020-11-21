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

# set the expression matrix containing folder as the working directory
dir <- "/FFPEs/FFPE-2"  
setwd(dir)

##read in the coordinates of points lying on top of the tissue.position.txt is generated from matlab script "Pixel_identification.m". 
location <- read.table("position.txt", sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])
x = x[-1]

##read expression matrix and generate the Filtered_matrix.tsv, which contains only the useful pixels
my_data <- read.table(file = 'FFPE-2_exp_matrix.tsv', sep = '\t', header = TRUE, stringsAsFactors=FALSE)
data_filtered <- my_data[my_data$X %in% x,]
write.table(data_filtered, file = 'Filtered_matrix.tsv', sep = '\t',col.names=TRUE, row.names = FALSE,quote = FALSE)

##calculate the total UMI count and Gene count
count <- rowSums(data_filtered[,2:ncol(data_filtered)])
data_filtered_binary <- data_filtered[,2:ncol(data_filtered)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)

##UMI Count 
region <- 2500  #change the x axis maxium, need to adjust based on different sample
test <- data_filtered %>% separate(X, c("A", "B"),  sep = "x")
df <- data.frame(number=1, c=count)
pdf(file = paste("UMI.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=region/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="UMI",limits = c(0,region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
  #xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

##Gene Count
df <- data.frame(number=1, c=gene_count)
region = 1500 #change the x axis maxium, need to adjust based on different sample
pdf(file = paste("Gene.pdf",sep =""), width=8.6, height=8.6)
ggplot(df,aes(x=c),color='blue', xlab="Gene") + 
  geom_histogram(aes(y=..density..),binwidth=region/20, color="black", fill="white",size=1)+ 
  geom_density(alpha=.2, fill="#FF6666",size=1,color ="red") +
  scale_x_continuous(name="Gene",limits = c(0,region)) + 
  scale_y_continuous(name="Density", expand = c(0, 0)) + 
#xlim(0,4000) +
  #expand_limits(x = 0, y = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(colour="black",size=20),
        axis.title=element_text(colour="black",size=25,face="bold"),
        legend.text=element_text(colour="black",size=20),
        legend.title = element_text(colour="black", size=20, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
dev.off()

#imported_raster=OpenImageR::readImage("ventricle.jpg")     #if you want the microscope image under the heatmap, then uncomment this line.
g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

#UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number.
pdf(file = paste("UMI_heatmap.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=count)) +
  #scale_color_gradientn(colours = c("black", "green")) + 
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,1000),
                        oob = scales::squish) +
  ggtitle("UMI") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

#Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number.
pdf(file = paste("Gene_heatmap.pdf",sep =""), width=8.6, height=8.6)
ggplot(test, aes(x = as.numeric(A), y = as.numeric(B), color=gene_count)) +
  scale_color_gradientn(colours = c("blue","green", "red"),limits=c(0,1000),
                        oob = scales::squish) +
  ggtitle("Gene") +
  #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  #if you want the microscope image under the heatmap, then uncomment this line.
  guides(colour = guide_colourbar(barwidth = 1, barheight = 30)) +
  geom_point(shape = 15, size = 3)+
  expand_limits(x = 0, y = 0) +
  scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
  scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
  coord_equal(xlim=c(0,51),ylim=c(51,1)) +
  theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        #legend.title = element_text(colour="black", size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()



