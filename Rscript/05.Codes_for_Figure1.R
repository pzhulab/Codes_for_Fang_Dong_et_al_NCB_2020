### Generating figures 1b-d and 1f

### figure 1b
rm(list=ls())
pwd <- getwd()
library(pheatmap)
library(grid)
### load expression and meta data of homeostasis and transplantation cell
load(paste0(pwd,"/input/01.QC.Comb.cell.count.tpm.meta.Rdata"))
### load HSC PCA top genes and homeostasis cell clusters
load(paste0(pwd,"/input/02.pc.topGenes.HSC.RData"))
load(paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))


cell.type <- "HSC"
cell.subtype <- c("LT_HSC", "Fraction I", "Fraction III", "ESLAM", "ESLAMSK") 
# number of gene sets distinguishing cell type subsets 
max(max(pc.genes[pc.genes$Celltype == cell.type,]$Cluster)) # 4
cut.row <- 4
# number of subsets within HSC cell type 
max(max(pc.cells[pc.cells$Celltype == cell.type,]$Cluster)) # 3
cut.col <- 3


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

### figure 1c
### note:figure 1c showed slightly different from published one for uwot
### packages had updated and same parameters we had used were missing.
library(Seurat)
library(uwot)
input.umi <- comb.data.counts[,rownames(input.meta)]
raw.obj <- CreateSeuratObject(
	raw.data = input.umi, min.cells = 3, min.genes = 0, 
	meta.data = input.meta, project = "HSC")
raw.obj <- NormalizeData(object = raw.obj, normalization.method = "LogNormalize", scale.factor = 10000)
raw.obj <- ScaleData(raw.obj, display.progress = F)
raw.obj <- FindVariableGenes(object = raw.obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.25, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = raw.obj@var.genes)
high.var.gene = raw.obj@var.genes[1:1101]
raw.obj <- RunPCA(object = raw.obj,pc.genes =high.var.gene, do.print = TRUE, pcs.print = 1:5, genes.print = 5)


# extract PCA for uwot
data.plot <- as.matrix(raw.obj@dr$pca@cell.embeddings)
set.seed(123)
umap <- umap(data.plot[,1:11],n_neighbors = 15, n_components = 3, metric = "euclidean",
     n_epochs = NULL, scale = FALSE, init = "random", min_dist=0.5, spread=0.7,
	 set_op_mix_ratio = 1,local_connectivity = 1, bandwidth = 1, repulsion_strength = 1,
     target_weight = 0.5)
umap <- as.data.frame(umap)
colnames(umap) <- paste0("UMAP",c(1:3))
rownames(umap) <- rownames(raw.obj@meta.data)
umap <- cbind(umap[,c(1:2)],raw.obj@meta.data)
umap$celltype <- paste0("tHSC",umap$cluster)


theme_imfo <- theme(legend.position="right",legend.spacing.y = unit(0.1, 'cm'),legend.text = element_text(colour="black", size=6), 
	axis.text=element_text(size=8,color="black"),axis.text.x=element_text(size=8),axis.title.x = element_text(size = 8),
	axis.text.y = element_text(size = 8,color="black"),axis.title.y = element_text(size = 8),panel.grid.minor = element_blank())
ggplot(data=umap, aes(x=UMAP1, y=UMAP2,colour=celltype)) + geom_point(aes(colour=celltype),alpha=1,size=1) + 
	 labs(fill="Celltype")+ guides(fill=guide_legend(keywidth=0.1,keyheight=0.1,default.unit="inch")) + 
     xlab("UMAP1") + ylab("UMAP2") + theme_imfo
     ggsave(file=paste(pwd,"/output/Figure.1c.HSC.UMAP.pdf", sep=""),width=4,height=3,device = "pdf")


### figure 1d 
library(ggpubr)
library(data.table)
genes <- c("Egr1","Nr4a1","Sh3gl1","S100a9","Cd79a","Blnk")
input.mtx <- comb.data.tpm[genes, row.names(input.meta)]
input.meta$celltype <- paste0("tHSC",input.meta$clusters)
cell.info <- as.data.frame(input.meta[,"celltype"],row.names=rownames(input.meta),stringsAsFactors=F)
Mat <- cbind(cell.info,t(input.mtx))
mMat <- reshape2::melt(Mat)
colnames(mMat) <- c("celltype","gene","Expr")


ggplot(mMat,aes(x=factor(celltype,levels=c("tHSC1","tHSC2","tHSC3")),y=Expr)) + 
         geom_boxplot(aes(fill = celltype)) +
         xlab("celltype")+ ylab("log2(TPM+1)") + 
		 ggtitle("HSC specific genes") + coord_cartesian(ylim = c(0, 10)) + 
         theme(panel.background = element_blank(),panel.grid=element_line(color="black",size=0.5))	+	 
		 facet_wrap(~gene,ncol=2)
         ggsave(file=paste0(pwd,"/output/Figure.1d.HSC_specificgene_genes.pdf"),width=6,height=6,device = "pdf") 


### figure 1e 
library(corrplot)
library(ggplot2)
library(plyr)
### load data
load(paste0(pwd,"/input/05.Cluster.renumber.tMPPs.tsne.homeostasis_transplantation.RData"))


# calculate each phenotype cell percentages in predicted cell types
tsne.df.hsc <- tsne.df[tsne.df$predict_cell_type %in% c("tHSC1","tHSC2","tHSC3") & tsne.df$Time=="Ho",]
count.data <- ddply(tsne.df.hsc, .(predict_cell_type, phenotype), summarize, counts= length(phenotype))
percent.data <- ddply(count.data, .(predict_cell_type), summarise, phenotype=phenotype, counts=counts, percent=counts/sum(counts))
cell.order <- c("LT_HSC", "Fraction I", "Fraction III", "ESLAM", "ESLAMSK") 
cell.order.predict <- c("tHSC1","tHSC2","tHSC3")
percent.data$phenotype <- factor(percent.data$phenotype,levels=cell.order,ordered=T)
percent.data$predict_cell_type <- factor(percent.data$predict_cell_type,levels=cell.order.predict,ordered=T)
percent.data <- percent.data[ order(percent.data$predict_cell_type, percent.data$phenotype), ]


ggplot(data=percent.data, aes(fill=phenotype)) + 
     geom_bar(aes(x=predict_cell_type, y=percent), stat = "identity") + 
     theme(legend.position="right") + 
     scale_fill_manual(name="celltype",values=c("LT_HSC"="#AFD0E8","Fraction I"="#E08BB8","Fraction III"="#5183C4","ESLAM"="#33AE3E",
     "ESLAMSK"="#ED7A92",breaks=cell.order)) + 
     guides(fill=guide_legend(keywidth=0.2,keyheight=0.1,default.unit="inch")) +
     xlab("celltype") + ylab("Percentages") 
     ggsave(filename = paste0(pwd,"/output/Figure.1e.",cell.type,".phenotype.percentages.pdf"),width=6,height=6,device = "pdf")
     dev.off()

	 
### figure 1f
# calculate each predicted cell type percentages in phenotype cells
count.data <- ddply(tsne.df.hsc, .(phenotype, predict_cell_type), summarize, counts= length(predict_cell_type))
percent.data <- ddply(count.data, .(phenotype), summarise, predict_cell_type=predict_cell_type, counts=counts, percent=counts/sum(counts))
percent.data$phenotype <- factor(percent.data$phenotype,levels=cell.order,ordered=T)
percent.data$predict_cell_type <- factor(percent.data$predict_cell_type,levels=cell.order.predict,ordered=T)
percent.data <- percent.data[ order(percent.data$predict_cell_type, percent.data$phenotype), ]


Cellfreq <- data.frame("phenotype" = NULL, "predict_cell_type" = NULL, "counts"= NULL,"percent"= NULL,stringsAsFactors = F)
for(i in cell.order.predict){
	 thisfreq <- data.frame("phenotype" = cell.order,"predict_cell_type" = i,"counts" = 0, "percent" = 0)
	 filter <- percent.data[percent.data$predict_cell_type==i,]
	 thisfreq <- rbind(filter,thisfreq)
	 thisfreq <- thisfreq[!duplicated(thisfreq$phenotype),] 
	 thisfreq$phenotype <- factor(thisfreq$phenotype,levels=cell.order,ordered=T)
	 thisfreq <- thisfreq[order(thisfreq$phenotype),]
	 Cellfreq <- rbind(Cellfreq,thisfreq)
	 }


## convert percentage value into a numeric matrix
Mat <- matrix(NA,nrow=length(cell.order),ncol=length(cell.order.predict))
rownames(Mat) <- cell.order
colnames(Mat) <- cell.order.predict
for(i in 1:length(cell.order.predict)){
     j <- cell.order.predict[i]
     Mat[,j] <- Cellfreq[Cellfreq$predict_cell_type==j,"percent"]
	 }
dat <- apply(Mat,2,as.numeric)
rownames(dat) <- rownames(Mat)
col <- colorRampPalette(c( 'white','orange','red3'))


pdf(paste0(pwd,"/output/Figure.1f.HSC_percent.circle.pdf"),3,4)
     corrplot(dat, is.corr = FALSE,tl.cex = 0.8,addCoef.col = "black",
     addrect = 2,
     rect.col = "white",
     method="circle",
     cl.lim=c(0,1),
	 cl.length=3,
     addCoefasPercent = FALSE,
     col=col(100),
     cl.pos="b",
     cl.cex = 0.8,
     cl.ratio = 0.15,
	 number.cex=0.8,
     number.font=1,
     tl.col="black")
     dev.off()