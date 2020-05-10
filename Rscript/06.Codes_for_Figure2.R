### Generating figures 2
### note: cell subset of tMPPs were renumber that cell cluster 1, 2, 3, 4 and 5 represent tMPP5, 
### tMPP4, tMPP1, tMPP2 and tMPP3 respectively.

### Figure 2a
rm(list=ls())
pwd <- getwd()
library(pheatmap)
### load MPP PCA top genes
load(paste0(pwd,"/input/02.pc.topGenes.MPP.RData"))
### load expression and meta data of homeostasis and transplantation cells
load(paste0(pwd,"/input/01.QC.Comb.cell.count.tpm.meta.Rdata"))
### load homeostasis cell cluster 
load(paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))


cell.type <- "MPP"
cell.subtype <- c("ST_HSC", "LMPP", "MPP1", "MPP2", "MPP3", "MPP4", "Fraction II", "HPC2", "HPC3")
# number of gene sets distinguishing cell type subsets 
max(max(pc.genes[pc.genes$Celltype == cell.type,]$Cluster)) # 5
cut.row <- 5
# number of subsets within HSC cell type 
max(max(pc.cells[pc.cells$Celltype == cell.type,]$Cluster)) # 5
cut.col <- 5


input.meta <- comb.data.meta[pc.cells[pc.cells$Celltype == cell.type,]$Cell,]
input.mtx <- comb.data.tpm[pc.topGenes, rownames(input.meta)]
input.dist.row <- as.dist(1-cor(t(input.mtx), method = "pearson"))
input.clust.row <- hclust(input.dist.row, method = "ward.D2")
input.dist.col <- as.dist(1-cor(input.mtx, method = "pearson"))
input.clust.col <- hclust(input.dist.col, method = "ward.D2")
col.panel <- colorRampPalette(colors = c("black", "gold"))
input.cutree.row <- cutree(input.clust.row, k = cut.row)
input.cutree.col <- cutree(input.clust.col, k = cut.col)
input.meta$clusters <- input.cutree.col
input.meta <- input.meta[order(input.meta$clusters),]


input.heatmap <- pheatmap(input.mtx, color=col.panel(100), cluster_rows = T,
                          cluster_cols = T, legend = T, show_colnames = F, show_rownames = F,
                          clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
                          clustering_method = "ward.D2", 
                          annotation_col = data.frame("Cell clusters" = factor(input.meta$clusters), "Cell type" = input.meta$phenotype, row.names= row.names(input.meta)),
                          annotation_row = data.frame("Gene clusters" = factor(input.cutree.row), row.names = pc.topGenes),
                          cutree_rows = max(input.cutree.row), cutree_cols = max(input.cutree.col))
dev.off()


# visualised expression distribution
library(data.table)
library(ggplot2)
ggplot(reshape2::melt(input.mtx),aes(x=value))+geom_line(stat="density",colour="black")
ggsave(file=paste0(pwd,"/output/Expression level distribution.",cell.type,".pdf"),width=4,height=3,device = "pdf")
dev.off()


# set genes expression level higher than 8 equally to 8, and visualise clustering
input.mtx[ input.mtx > 8] <- 8
pdf(paste0(pwd,"/output/Figure.1b.", cell.type, ".clustering.heatmap.pdf"), onefile = F)
pheatmap(input.mtx, color=col.panel(100), cluster_rows = input.heatmap$tree_row,
         cluster_cols = input.heatmap$tree_col, legend = T, show_colnames = F, show_rownames = F,
         annotation_col = data.frame("Cell clusters" = factor(input.meta$clusters), "Cell type" = input.meta$phenotype, row.names= row.names(input.meta)),
         annotation_row = data.frame("Gene clusters" = factor(input.cutree.row), row.names = pc.topGenes),
         cutree_rows = max(input.cutree.row), cutree_cols = max(input.cutree.col),
         main = paste(cell.type, "clustering"))
         dev.off()


### figure 2b
rm(list=ls())
pwd <- getwd()
library(plyr)
library(ggplot2)
### load meta table contain transcriptional cell type have renumbered tMPPs subset
load(paste0(pwd,"/input/05.Cluster.renumber.tMPPs.tsne.homeostasis_transplantation.RData"))


cell.type <- "MPP"
cell.order <- c("ST_HSC", "LMPP", "MPP1", "MPP2", "MPP3", "MPP4", "Fraction II", "HPC2", "HPC3")
cell.order.predict <- c("tMPP1","tMPP2","tMPP3","tMPP4","tMPP5")
tMPPs <- tsne.df[tsne.df$New.Celltype %in% cell.order.predict & tsne.df$Ho_Tx == "Ho",]


count.data <- ddply(tMPPs, .(New.Celltype, phenotype), summarize, counts= length(phenotype))
percent.data <- ddply(count.data, .(New.Celltype), summarise, phenotype=phenotype, counts = counts, percent = counts/sum(counts))
percent.data$phenotype <- factor(percent.data$phenotype,levels=cell.order,ordered=T)
percent.data$New.Celltype <- factor(percent.data$New.Celltype,levels=cell.order.predict,ordered=T)
percent.data <- percent.data[ order(percent.data$New.Celltype, percent.data$phenotype), ]


ggplot(data=percent.data, aes(fill=phenotype)) + 
  geom_bar(aes(x=New.Celltype, y=percent), stat = "identity") + 
  theme(legend.position="top") + 
  scale_fill_manual(name="celltype",values=c("ST_HSC"="#F5AB53","LMPP"="#97C627","MPP1"="#ED7A92","MPP2"="#D59EC5",
    "MPP3"="#BED970","MPP4"="#81CCDA","Fraction II"="#008B8C","HPC2"="#C8C8E4","HPC3"="#38AE37",breaks=cell.order)) + 
	guides(fill=guide_legend(keywidth=0.2,keyheight=0.1,default.unit="inch")) +
    xlab("celltype") + ylab("Percentages") 
ggsave(filename = paste0(pwd,"/output/Figure.2b.",cell.type,".phenotype.percentage.pdf"),width=5,height=5,device = "pdf")


### figure 2c
count.data <- ddply(tMPPs, .(phenotype, New.Celltype), summarize, counts= length(phenotype))
percent.data <- ddply(count.data, .(phenotype), summarise, New.Celltype=New.Celltype, counts=counts, percent=counts/sum(counts))
percent.data$phenotype <- factor(percent.data$phenotype,levels=cell.order,ordered=T)
percent.data$New.Celltype <- factor(percent.data$New.Celltype,levels=cell.order.predict,ordered=T)
percent.data <- percent.data[ order(percent.data$New.Celltype, percent.data$phenotype), ]


ggplot(data=percent.data, aes(fill=New.Celltype)) + 
  geom_bar(aes(x=phenotype, y=percent), stat = "identity") + 
  theme(legend.position="right") + 
  scale_fill_manual(name="celltype",values=c("tMPP1"="#ED7975","tMPP2"="#F28F19",
    "tMPP3"="#008B8C","tMPP4"="#CCDF7C","tMPP5"="#AE86BA",breaks=cell.order.predict)) + 
	guides(fill=guide_legend(keywidth=0.2,keyheight=0.2,default.unit="inch")) +
    xlab("phenotype") + ylab("Percentages") 
ggsave(filename = paste0(pwd,"/output/Figure.2c.",cell.type,".predict.celltype.percentage.pdf"),width=7,height=4,device = "pdf")


### figure 2d
rm(list=ls())
pwd <- getwd()
library(monocle) 
### load expression data of homeostasis cells
load(paste0(pwd,"/input/01.Homeostasis.Cells.UMI_TPM_metadata.RData"))
### load PCA top genes and homeostasis cell clusters
load(paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))
### load meta have renumber tMPP subsets
load(paste0(pwd,"/input/05.Cluster.renumber.tMPPs.tsne.homeostasis_transplantation.RData"))


homeostasis.exp <- tpm.Ho[pc.genes$Gene,]
homeostasis.meta <- tsne.df[tsne.df$Ho_Tx == "Ho",]


Cell_type.v2 <- homeostasis.meta$New.Celltype
Cell_type.v2.order <- unlist(sapply(unique(pc.cells$Celltype), function(x){
  paste0("t", x, names(table(pc.cells$Cluster[pc.cells$Celltype == x ])))
}), use.names = F)
Cell_type.v2 <- Cell_type.v2[match(colnames(homeostasis.exp), pc.cells$Cell)]
homeostasis.meta$New.Celltype <- factor(homeostasis.meta$New.Celltype,levels = Cell_type.v2.order, ordered = T)
dim(homeostasis.exp) # 1044 1270
dim(homeostasis.meta) # 1270   16


# Store Data in a CellDataSet Object
pd <- new("AnnotatedDataFrame", data = homeostasis.meta)
fd <- new("AnnotatedDataFrame", data = data.frame("gene_short_name"=row.names(homeostasis.exp), row.names = row.names(homeostasis.exp)))
cds <- newCellDataSet(as(data.matrix(homeostasis.exp), "sparseMatrix"), 
                      phenoData = pd, featureData = fd, 
                      expressionFamily=gaussianff())
dim(cds) # 1044     1270


# Trajectory step 1: choose genes that define a cell's progress
cds_traj <- cds
cds_traj <- setOrderingFilter(cds_traj, pc.genes$Gene)
# Trajectory step 2: reduce data dimensionality
cds_traj.rd <- reduceDimension(cds_traj, max_components = 4, norm_method = "none",
                               reduction_method = 'DDRTree')
# Trajectory step 3: order cells along the trajectory
library(RColorBrewer)
cds_traj.path <- orderCells(cds_traj.rd)
plot_cell_trajectory(cds_traj.path, color_by = "New.Celltype", show_branch_points = F, theta = 45) +
  facet_wrap(~New.Celltype, nrow=4) + 
  scale_fill_manual(name="Cell type", values = colorRampPalette(brewer.pal(n=9, name = "Set1"))(length(Cell_type.v2.order))) + 
  theme(legend.position="none") +
  theme(text = element_text(size=20))


homeostasis.meta$x <- cds_traj.path@reducedDimS[1,]
homeostasis.meta$y <- cds_traj.path@reducedDimS[2,]


# add cell type label
{
  cellTypeLabel <- function(cellType, tsne, meta){
    as.data.frame(t(sapply(cellType, function(x){      
      x1 <- median(tsne$x[meta == x])
      y1 <- median(tsne$y[meta == x])
      return(c(x1, y1))
    })))
  }
  cell_type.label <- cellTypeLabel(Cell_type.v2.order, homeostasis.meta, as.vector(homeostasis.meta$New.Celltype))
  cell_type.label$cellType <- Cell_type.v2.order
  names(cell_type.label) <- c("x", "y", "Celltype")
}


# backbone 
reduced_dim_coords <- reducedDimK(cds_traj.path)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(1,2), ]))
colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)
ica_space_df$sample_state <- row.names(ica_space_df)
dp_mst <- minSpanningTree(cds_traj.path)
if (is.null(dp_mst)) {
  stop("You must first call orderCells() before using this function")
}
library(igraph)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c("source", "target")
edge_df <- merge(ica_space_df, edge_list, by.x = "sample_name", 
                 by.y = "source", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "source_prin_graph_dim_2"))
edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                           "prin_graph_dim_1", "prin_graph_dim_2")], by.x = "target", 
                 by.y = "sample_name", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "target_prin_graph_dim_2"))
S_matrix <- reducedDimS(cds_traj.path)
data_df <- data.frame(t(S_matrix[c(1, 2), ]))
sample_state <- pData(cds_traj.path)$State
data_df <- cbind(data_df, sample_state)
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- row.names(data_df)
lib_info_with_pseudo <- pData(cds_traj.path)
data_df <- merge(data_df, lib_info_with_pseudo, by.x = "sample_name", 
                 by.y = "row.names")
return_rotation_mat <- function(theta) {
  theta <- theta/180 * pi
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
         nrow = 2)
}
theta <- 0
tmp <- return_rotation_mat(theta) %*% t(as.matrix(data_df[, 
                                                          c(2, 3)]))
data_df$data_dim_1 <- tmp[1, ]
data_df$data_dim_2 <- tmp[2, ]
tmp <- return_rotation_mat(theta = theta) %*% t(as.matrix(edge_df[, 
                                                                  c("source_prin_graph_dim_1", "source_prin_graph_dim_2")]))
edge_df$source_prin_graph_dim_1 <- tmp[1, ]
edge_df$source_prin_graph_dim_2 <- tmp[2, ]
tmp <- return_rotation_mat(theta) %*% t(as.matrix(edge_df[, 
                                                          c("target_prin_graph_dim_1", "target_prin_graph_dim_2")]))
edge_df$target_prin_graph_dim_1 <- tmp[1, ]
edge_df$target_prin_graph_dim_2 <- tmp[2, ]
bg.plot <- data.frame(x = homeostasis.meta$x, y = homeostasis.meta$y, celltype = homeostasis.meta$New.Celltype)
library(ggplot2)
library(RColorBrewer)
col.values = colorRampPalette(brewer.pal(n=9, name = "Set1"))(length(Cell_type.v2.order))
name <- sapply(Cell_type.v2.order, function(x){
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
names(col.values) <- name

# save color values
save(col.values, file = paste0(pwd,"/input/Colors.21Cells.Rdata"))


ggplot(data=homeostasis.meta, aes(x = x, y = y)) + 
     geom_point(data=bg.plot, aes(x = x, y = y,fill = celltype), shape = 21, size=1, colour="grey", alpha=I(0.3)) +
     geom_segment(aes_string(x = "source_prin_graph_dim_1", 
         y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
         yend = "target_prin_graph_dim_2"), size = 0.75, 
         linetype = "solid", na.rm = TRUE, data = edge_df) +
     geom_point(aes(fill=factor(New.Celltype)), shape=21, colour="black", size=3) +
     scale_fill_manual(name="Cell type", values = col.values) + 
     scale_alpha(guide=F)+ theme(text = element_text(size=20)) + 
     facet_wrap(~New.Celltype, nrow=4) + 
     theme(legend.position="right", legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
     theme(panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black")) +
     xlab("Dimension 1") + ylab("Dimension 2")
     ggsave(filename = paste0(pwd,"/output/monocle_single-cell_tracjectory.facetByCelltype2.pdf"),device = "pdf")


ggplot(data=homeostasis.meta, aes(x = x, y = y)) + 
     geom_segment(aes_string(x = "source_prin_graph_dim_1", 
         y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
         yend = "target_prin_graph_dim_2"), size = 0.75, 
          linetype = "solid", na.rm = TRUE, data = edge_df) +
     geom_point(aes(fill = factor(New.Celltype)), shape=21, colour="grey", size=3) +
     scale_fill_manual(name = "Cell type", values = col.values) + 
     geom_label(data = cell_type.label, aes(x = x, y =  y,label = Celltype, alpha=0.1), size=4) + 
     scale_alpha(guide=F)+ theme(text = element_text(size=20)) + 
     theme(legend.position="right", legend.key.height = grid::unit(0.35, "in")) + theme(legend.key = element_blank()) + 
     theme(panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black")) +
     xlab("Dimension 1") + ylab("Dimension 2")
      ggsave(filename = paste0(pwd,"/output/Figure.2d.monocle_single-cell_tracjectory.facetByCelltype.pdf"),width=8,height=5,device = "pdf")


# save data
save(homeostasis.meta, edge_df, Cell_type.v2.order, cell_type.label, bg.plot, cds_traj.path, file = paste0(pwd,"/input/Figure.2d.homeostasis.monocle.RData"))