### 4. Constructing homeostasis and transplantation cells differentiation tsne map

rm(list=ls())
pwd <- getwd()
library(Seurat)
### load expression and meta table of cells under homeostasis and transplantation 
load(paste0(pwd,"/input/01.QC.Comb.cell.count.tpm.meta.Rdata"))
### load PCA top genes and clusters of homeostasis cells
load(paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))


comb.pc <- prcomp(t(comb.data.tpm[pc.genes$Gene,rownames(comb.data.meta)]))
comb.dist <- dist(comb.pc$x[,1:10])
set.seed(123)
comb.tsne <- Rtsne::Rtsne(comb.dist, is_distance=TRUE, perplexity=30, check_duplicates=F,
                            verbose = TRUE,max_iter=1500, theta=0)
plot(comb.tsne$Y, pch=20)
cell.names <- attr(comb.dist, "Labels")


# save data
save(cell.names, comb.tsne, file = paste0(pwd,"/input/04.Comb_tsne.homeostasis_transplantation.RData"))