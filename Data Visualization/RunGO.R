require(clusterProfiler)
require(ReactomePA)
library(mygene)
library(tidyverse)
genetable <-  data.frame(cluster = go$cluster, gene = go$gene)
temp = genetable[genetable$cluster == 0,] 

mylist.names <- c("X0","X1", "X2", "X3","X4","X5","X6","X7","X8")
mylist <- vector("list", length(mylist.names))
names(mylist) <- mylist.names

for (i in c(0:8)) {
  temp = genetable[genetable$cluster == i,] 
  xli = temp$gene
  result <- queryMany(xli, scopes="symbol", fields=c("uniprot", "ensembl.gene", "reporter"), species="mouse")
  mylist[[i+1]] = result$`_id`
}
mylist <- lapply(mylist, function(x) x[!is.na(x)])
mylist <- lapply(mylist, function(x) str_remove(x,"ENSMUSG000000"))
mylist

res <- compareCluster(mylist, fun="enrichPathway")
dotplot(res)
