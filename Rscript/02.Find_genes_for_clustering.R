### 2. find genes contributing most in distinguish homeostasis cells

rm(list=ls())
pwd <- getwd()
library(pheatmap)
### load homeostasis cell expression and meta table
load(paste0(pwd,"/input/01.Homeostasis.Cells.UMI_TPM_metadata.RData"))
  
  
pc.genes <- data.frame("Gene"=NULL, "Cluster"=NULL, "Celltype"=NULL, stringsAsFactors = F)
pc.cells <- data.frame("Cell"=NULL, "Cluster"=NULL, "Celltype"=NULL, stringsAsFactors = F)
group <- c("HSC","MPP","CP","ME","GM","Lym")
for(i in group){
     if(i=="HSC"){
	 cell.type <- "HSC"
	 cell.subtype <- c("LT_HSC", "Fraction I", "Fraction III", "ESLAM", "ESLAMSK") 
	 cut.row <- 4
	 cut.col <- 3}
	 if(i=="MPP"){
     cell.type <- "MPP"
	 cell.subtype <- c("ST_HSC", "LMPP", "MPP1", "MPP2", "MPP3", "MPP4", "Fraction II", "HPC2", "HPC3") 
	 cut.row <- 5
	 cut.col <- 5}
	 if(i=="CP"){
	 cell.type <- "CP"
	 cell.subtype <- c("CMP", "GMP", "MEP", "CLP")
	 cut.row <- 3
	 cut.col <- 3}
	 if(i=="ME"){
	 cell.type <- "ME"
	 cell.subtype <- c("MK", "EryA", "EryB")
	 cut.row <- 3
	 cut.col <- 3}
	 if(i=="GM"){
	 cell.type <- "GM"
	 cell.subtype <- c("Macrophage", "Monocyte", "Granulocyte")
	 cut.row <- 3
	 cut.col <- 3}
	 if(i=="Lym"){
	 cell.type <- "Lym"
	 cell.subtype <- c("NK cell", "B cell", "CD4T", "CD8T")
	 cut.row <- 3
	 cut.col <- 4}
	 input.cells <- tpm.Ho[, meta.Ho$phenotype %in% cell.subtype]
     input.meta <- meta.Ho[ colnames(input.cells), ]
	 input.cells.pca <- prcomp(t(input.cells))
     pc.topGenes <- unique(as.character(apply(input.cells.pca$rotation[,1:3], 2, function(x){ 
     index <- sort(abs(x), decreasing=T, index.return =T)
	 row.names(input.cells.pca$rotation)[index$ix[1:100]]
	 })))
	 save(pc.topGenes, file = paste0(pwd,"/input/02.pc.topGenes.", cell.type, ".RData"))
	 input.mtx <- input.cells[pc.topGenes,]
     input.dist.row <- as.dist(1-cor(t(input.mtx), method = "pearson"))
     input.clust.row <- hclust(input.dist.row, method = "ward.D2")
     input.dist.col <- as.dist(1-cor(input.mtx, method = "pearson"))
     input.clust.col <- hclust(input.dist.col, method = "ward.D2")
	 col.panel <- colorRampPalette(colors = c("black", "gold"))
	 input.cutree.row <- cutree(input.clust.row, k = cut.row)
	 input.cutree.col <- cutree(input.clust.col, k = cut.col)
	 input.meta$clusters <- input.cutree.col
	 input.meta <- input.meta[order(input.meta$clusters),]
	 input.mtx <- tpm.Ho[pc.topGenes, row.names(input.meta)]
	 pdf(paste0(pwd,"/output/", cell.type, "_groups_comparison.heatmap.pdf", sep = ""))
	 pheatmap(input.mtx, color=col.panel(100), cluster_rows = input.clust.row,
         cluster_cols = input.clust.col, legend = T, show_colnames = F, show_rownames = F,
         annotation_col = data.frame("Cell type" = input.meta$phenotype, row.names= row.names(input.meta))
         )
	 dev.off()
     this.genes <- pc.topGenes
     this.genes.cluster <- as.numeric(input.cutree.row)
     this.cells <- names(input.cutree.col)
     this.cells.cluster <- as.numeric(input.cutree.col)
     this.celltype <- i
     overlap.idx <- which(this.genes %in% pc.genes$Gene)
     if(length(overlap.idx) > 0){
         this.genes <- this.genes[-overlap.idx]
         this.genes.cluster <- this.genes.cluster[-overlap.idx]
         }
     this.genes.data <- data.frame("Gene"=this.genes, "Cluster"=this.genes.cluster, "Celltype"=rep(this.celltype, length(this.genes)), stringsAsFactors = F)
     pc.genes <- rbind(pc.genes, this.genes.data)
     this.cells.data <- data.frame("Cell"=this.cells, "Cluster"=this.cells.cluster, "Celltype"=rep(this.celltype, length(this.cells)), stringsAsFactors = F)
     pc.cells <- rbind(pc.cells, this.cells.data)
     }
dim(pc.genes) # 1044    3
dim(pc.cells) # 1270    3


## save data 
save(pc.genes,pc.cells,file = paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))