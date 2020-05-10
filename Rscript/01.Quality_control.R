### 1. Quality contol on cells and genes

rm(list=ls())
pwd <- getwd()
### the number of transplantation cells we had used for this step 
### was more than 2879. This is the reason the cell number remain 
### on this step is not equally to 2301, the number of cells retained
### for further analysis.
### 1. Quality control on cells and genes 
### for homeostasis and transplantation cells respectively
### load raw cell and meta data
load(paste0(pwd,"/input/01.RawCell_UMI_TPM_Meta.Rdata")) 


## create output directory
dir.create(paste0(pwd,"/output"))


### homeostasis data
metadata.Ho <- raw_cells[raw_cells$Ho_Tx == "Ho",]
cells.idx.Ho <- match(rownames(metadata.Ho), colnames(raw.counts))
all.counts.Ho <- as.data.frame(raw.counts[, cells.idx.Ho])
dim(all.counts.Ho) # 23154  1569
 
 
### Transpantation data
metadata.Tx <- raw_cells[raw_cells$Ho_Tx == "Tx",]
cells.idx.Tx <- match(rownames(metadata.Tx), colnames(raw.counts))
all.counts.Tx <- as.data.frame(raw.counts[, cells.idx.Tx])
dim(all.counts.Tx) # 23154  1310


### Quality contol 
# Report total counts and expressed genes in each cell
SCQCstat <- function(data=data, is.expr=1){
  total_counts <- colSums(data, na.rm = T)
  total_genes <- apply(data, 2, function(x){
    sum(x>=is.expr)
  })
  return(list("total_counts"=total_counts, "total_genes"=unlist(total_genes)))
}


QC <- function(x,label){
     if(x=="Ho"){
	 counts.input <- all.counts.Ho
     meta.input <- metadata.Ho
	 } else {
	 counts.input <- all.counts.Tx
     meta.input <- metadata.Tx
	 }
     # quality control for cells and genes
     counts.qc <- SCQCstat(counts.input)
     size.drop <- scater::isOutlier(counts.qc$total_counts, nmads=3, type="both", log=F, batch = meta.input$phenotype)
     gene.drop <- scater::isOutlier(counts.qc$total_genes, nmads=3, type="both", log=F, batch = meta.input$phenotype)
     counts.filter <- counts.input[, !(size.drop | gene.drop)]
     counts.filter.qc <- SCQCstat(counts.filter)
     # filtering with gene number and umi counts
     gene.count <- 1000
     umi.count <- 20000
     counts.filter <- counts.filter[, counts.filter.qc$total_counts >= umi.count & counts.filter.qc$total_genes >= gene.count]
     dim(counts.filter) # Ho:23154 1270; Tx:23154  1059
     counts.filter.qc <- SCQCstat(counts.filter)
     # filtering low abundance genes
     ave.cut <- 0.2
     ave.counts <- rowMeans(counts.filter)
     counts.filter <- counts.filter[ log10(ave.counts) >= log10(ave.cut),]
	 if(x=="Ho"){
     filter.counts.Ho <- counts.filter
     filter.tpm.Ho <- raw.tpm[row.names(counts.filter), colnames(counts.filter)]
     meta.Ho <- meta.input[colnames(counts.filter),]
	 counts.Ho <- counts.input[,colnames(counts.filter)]
	 tpm.Ho <- raw.tpm[,colnames(counts.filter)]
	 save(filter.tpm.Ho, filter.tpm.Ho, counts.Ho, tpm.Ho, meta.Ho, file = paste0(pwd,"/input/01.Homeostasis.Cells.UMI_TPM_metadata.RData"))
	 return(meta.Ho)
	 } else {
     filter.counts.Tx <- counts.filter
     filter.tpm.Tx <- raw.tpm[row.names(counts.filter), colnames(counts.filter)]
     meta.Tx <- meta.input[colnames(counts.filter),]	 
	 counts.Ho <- counts.input[,colnames(counts.filter)]
	 tpm.Ho <- raw.tpm[,colnames(counts.filter)]
	 save(filter.counts.Tx, filter.tpm.Tx, counts.Ho, tpm.Ho, meta.Tx,file = paste0(pwd,"/input/01.Transplantation.Cells.UMI_TPM_metadata.RData"))  
	 return(meta.Tx)
	 }
}


meta.Ho <- QC("Ho")
meta.Tx <- QC("Tx")
## homeostasis cells remain 
dim(meta.Ho) # 1270   10
## transplantation cells remain 
dim(meta.Tx) # 1059   10 
