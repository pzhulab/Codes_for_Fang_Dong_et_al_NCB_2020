### 5. Classification of transplantation cells to homeostasis cell types

rm(list=ls())
pwd <- getwd()
### note we used more than 1,031 transplantation cells as training data
### set to perform this step, the result have been slightly influenced by 
### the number of training cells, so we can see a few cells changed their 
### cell type. However, it didn't affect the conclusions. 


### load expression and meta table of homeostasis and transplantation 
load(paste0(pwd,"/input/01.QC.Comb.cell.count.tpm.meta.Rdata"))
### load PCA top genes and clusters of homeostasis cells
load(paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))
### load dimensionality reduction results produced by the forward step 
load(paste0(pwd,"/input/04.Comb_tsne.homeostasis_transplantation.RData"))


homeostasis.train <- comb.data.tpm[pc.genes$Gene, row.names(comb.data.meta)[comb.data.meta$Ho_Tx == "Ho"]]
transplantation.test <- comb.data.tpm[pc.genes$Gene, row.names(comb.data.meta)[comb.data.meta$Ho_Tx == "Tx"]]
dim(homeostasis.train) # 1044 1270
dim(transplantation.test) # 1044 1031


Cell_type.v2 <- paste0("t", pc.cells$Celltype, pc.cells$Cluster)
Cell_type.v2.order <- unlist(sapply(unique(pc.cells$Celltype), function(x){
  paste0("t", x, names(table(pc.cells$Cluster[pc.cells$Celltype == x ])))
}), use.names = F)
Cell_type.v2 <- Cell_type.v2[match(colnames(homeostasis.train), pc.cells$Cell)]


# KNN classification
library(class)
# nearest number of cells
Knum <- 10
# by sub cell type
transplantation.classification <- knn(train = t(homeostasis.train), test = t(transplantation.test), 
                                cl = Cell_type.v2, prob = T, k = Knum)


# prepare data frame for transplantation cells
tsne.df <- data.frame( x= comb.tsne$Y[,1], y=comb.tsne$Y[,2])
tsne.df <- cbind(tsne.df, comb.data.meta)
tsne.df$predict_cell_type <- c(Cell_type.v2, as.character(transplantation.classification))
tsne.df$predict_cell_type.family <- tsne.df$predict_cell_type
tsne.df$predict_cell_type.family <- factor(tsne.df$predict_cell_type.family, levels=Cell_type.v2.order, ordered=T)
time.order <- c("Ho", "1Day", "3Day", "5Day", "7Day")
tsne.df$Time <- factor(tsne.df$Time, levels = time.order, ordered = T)


## save data 
save(tsne.df,file = paste0(pwd,"/input/05.Cluster.tsne.homeostasis_transplantation.RData"))


# renumber tMPP subset
# cell subset of tMPPs were renumbered that cell cluster 1, 2, 3, 4 and 5 
# represent tMPP5, tMPP4, tMPP1, tMPP2 and tMPP3 respectively.
New.Celltype <- sapply(tsne.df$predict_cell_type, function(x){
  if( x == "tMPP1" ){
    return("tMPP5")
  } else if (x == "tMPP2"){
    return("tMPP4")
  } else if (x == "tMPP3"){
    return("tMPP1")
  } else if (x == "tMPP4"){
    return("tMPP2")
  } else if (x == "tMPP5"){
	return("tMPP3")
  } else {
    return(x)
  }
})

tsne.df$New.Celltype <- New.Celltype

## save data 
save(tsne.df,file = paste0(pwd,"/input/05.Cluster.renumber.tMPPs.tsne.homeostasis_transplantation.RData"))