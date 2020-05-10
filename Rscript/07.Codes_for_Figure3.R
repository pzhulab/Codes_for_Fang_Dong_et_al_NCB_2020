### Generating figure 3c-g

### figure 3b
rm(list=ls())
pwd <- getwd()
### load data
load(paste0(pwd,"/input/05.Cluster.renumber.tMPPs.tsne.homeostasis_transplantation.RData"))
load(paste0(pwd,"/input/02.Homeostasis.Cells_cluster_pc.topGenes.RData"))
load(paste0(pwd,"/input/Colors.21Cells.Rdata"))


# cell type coordinate
# add cell type label
{
Cell_type.v2.order <- unlist(sapply(unique(pc.cells$Celltype), function(x){
     paste("t", x, names(table(pc.cells$Cluster[pc.cells$Celltype == x ])), sep = "")
}), use.names = F)
cellTypeLabel <- function(cellType, tsne, meta){
      as.data.frame(t(sapply(cellType, function(x){      
      x1 <- median(tsne$x[meta == x])
      y1 <- median(tsne$y[meta == x])
      return(c(x1, y1))
     })))
 }
cell_type.label <- cellTypeLabel(Cell_type.v2.order,tsne.df,tsne.df$New.Celltype)
cell_type.label$cellType <- Cell_type.v2.order
names(cell_type.label) <- c("x", "y", "Celltype")
}
tsne.df$New.Celltype <- factor(tsne.df$New.Celltype,levels=Cell_type.v2.order)


### homeostasis and transplantation cells
library(ggplot2)
ggplot(data = tsne.df, aes(x = x, y = y)) + 
     geom_point(aes(fill=factor(New.Celltype)), shape=21, colour="grey", size=3) +
     theme_classic() +
     scale_fill_manual(name="Cell type", values = col.values) + 
     geom_label(data=cell_type.label, aes(x=x,y=y,label=Celltype, alpha=0.1), size=4) + 
     scale_alpha(guide=F)+ theme(text = element_text(size=20)) + 
     xlab("tSNE1") + ylab("tSNE2")
     ggsave(filename = paste0(pwd,"/output/Figure.3c-d.total.cell.Predict_cell_types.tsne.pdf"), width=10)


### figure 3c for homeostasis cells
ggplot(data = tsne.df, aes(x = x,y = y)) + 
     geom_point(fill="grey", shape = 21, colour="white", size=4) +
     geom_point(data=subset(tsne.df, Time %in% "Ho"), aes(fill=factor(New.Celltype)), shape=21, colour="white", size=4) + 
     theme_classic() +
     scale_fill_manual(name="Cell type", values = col.values) + 
     scale_alpha(guide=F)+ theme(text = element_text(size=20)) + 
     xlab("tSNE1") + ylab("tSNE2") +
     ggtitle("Homeostasis")
     ggsave(filename = paste0(pwd,"/output/Figure.3c.Predict_cell_types.tsne.Homeostasis.pdf"), width=10)


### figure 3d for transplantation
ggplot(data = tsne.df, aes(x = x,y = y)) + 
     geom_point(fill="grey", shape = 21, colour="white", size=4) +
     geom_point(data=subset(tsne.df, Time %in% c("1Day", "3Day", "5Day", "7Day")), aes(fill=factor(New.Celltype)), shape=21, colour="white", size=4) + 
     theme_classic() +
     scale_fill_manual(name="Cell type", values = col.values) + 
     scale_alpha(guide=F)+ theme(text = element_text(size=20)) + 
     xlab("tSNE1") + ylab("tSNE2") +
     ggtitle("Transplantation")
     ggsave(filename = paste0(pwd,"/output/Figure.3d.Predict_cell_types.tsne.pdf"), width=10)


### figure 3e for different time point
time.order <- levels(tsne.df$Time)
for( i in time.order){
     if( i == "Ho"){
	     dat <- subset(tsne.df, Time == i & New.Celltype %in% c("tHSC1","tHSC2","tHSC3") & phenotype %in% c("ESLAMSK"))
	     } 
	 else{
         dat <- subset(tsne.df, Time == i)
			 }
     ggplot(data = tsne.df, aes(x = x,y = y)) + 
         geom_point(fill = "grey", shape = 21, size = 2, colour = "white", alpha = I(0.3)) +
         geom_point(data = dat, aes(fill = factor(New.Celltype)), shape = 21, colour = "white",size = 7) + 
         theme_classic() +
         scale_fill_manual(name = "Cell type", values = col.values) + 
         scale_alpha(guide = F)+ theme(text = element_text(size=20)) + 
         xlab("tSNE1") + ylab("tSNE2") +
         ggtitle(i)
         ggsave(filename = paste0(pwd,"/output/Figure.3e.Predict_cell_types.tsne.", i, ".pdf"), width=8)
     }


### figure 3f for dynamic frequencies of each cell type at day 1, 3, 5 and 7 after transplantation
library(plyr)
meta <- subset(tsne.df, Ho_Tx == "Tx" | phenotype %in% c("ESLAMSK"))
count.data <- ddply(meta,.(Time, New.Celltype), summarize, counts= length(New.Celltype))
percent.data <- ddply(count.data, .(Time), summarise, Celltype=New.Celltype, counts=counts, percent=counts/sum(counts))

cell.order <- Cell_type.v2.order[Cell_type.v2.order %in% unique(percent.data$Celltype)]
Cellfreq <- data.frame("Time" = NULL, "Celltype" = NULL, "counts" = NULL, "percent"=NULL,stringsAsFactors = F)
for(i in time.order){
	 thisfreq <- data.frame("Time"=i,"Celltype"=cell.order,"counts"=0,"percent"=0)
	 filter <- percent.data[percent.data$Time==i,]
	 thisfreq <- rbind(filter,thisfreq)
	 thisfreq <- thisfreq[!duplicated(thisfreq$Celltype),] 
	 thisfreq$Celltype <- factor(thisfreq$Celltype,levels=cell.order,ordered=T)
	 thisfreq <- thisfreq[order(thisfreq$Celltype),]
	 Cellfreq <- rbind(Cellfreq,thisfreq)
	 }


# plot 
library(ggplot2)
col <- col.values[names(col.values) %in% cell.order]
Cellfreq$Celltype <- factor(Cellfreq$Celltype,levels=rev(cell.order))
ggplot(Cellfreq, aes_string(x="Time", y="percent", group="Celltype", fill="Celltype", col=NULL)) +
     geom_area() + 
	 scale_fill_manual(values=col)+
	 theme_bw()+ theme(axis.text.x=element_text(angle=0, hjust=1),legend.position="right") +
	 xlab("Proportion") + ylab(" ") +
	 ggsave(filename = paste0(pwd,"/output/Figure.3f.CellFreq.Muller.plot.pdf"), width=10)


### figure 3g for dynamic frequencies of different subgroups in homeostasis
### cells at different time points
Lineage <- sapply(percent.data$Celltype, function(x){
  if( x %in% c("tHSC1", "tHSC2", "tHSC3")){
    return("tHSCs")
  } else if (x %in% c("tMPP1", "tMPP2", "tMPP3", "tMPP4", "tMPP5")){
    return("tMPPs")
  } else if (x %in% c("tCP1", "tCP2","tCP3")){
    return("tCPs")
  } else {
    return("tLineages")
  }
})


percent.data$Lineage <- Lineage
percent.data$Lineage <- factor(percent.data$Lineage,levels=c("tLineages","tCPs","tMPPs","tHSCs"))
ggplot(data= percent.data, aes(fill = Lineage)) + 
     geom_bar(aes(x = Time, y = percent), stat = "identity") + 
     theme(panel.background = element_blank(),axis.line= element_line(color='black'),legend.position="right")+
     scale_fill_manual(name = " ",values = c("tHSCs"="#EA5951","tMPPs" = "#30B283","tCPs" = "#78CCF2","tLineages" = "#F5AB53")) + 
     guides(fill = guide_legend(keywidth=0.2,keyheight = 0.1,default.unit = "inch")) +
	 xlab(" ") + ylab("Frequency") 
     ggsave(filename = paste0(pwd,"/output/Figure.3g.Lineage.percent.pdf"),width=5,height=5,device = "pdf")
