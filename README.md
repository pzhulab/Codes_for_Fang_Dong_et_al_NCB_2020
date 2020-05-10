# Codes_for_Fang_Dong_et_al_NCB_2020

This repository provides the codes used to perform data anlsysis in our recently published article: 

F. Dong et al., Differentiation of transplanted haematopoietic stem cells tracked by single-cell transcriptomic analysis. Nat Cell Biol. 125, 1–10 (2020).

If you use the scripts, please cite the paper!

Note: Due to the size of the expression data, we are not allowed to upload the data files to github. Therefore this repository only show you the codes we used to analyze the sequencing data in the paper. We are attempting to upload the data files by dividing them into small files, then you will completely reproduce the reuslts in our paper. Thanks for your interesting in our work!

script "Rscript/01.Quality_control.R" is used for quality control on cells and genes.

script "Rscript/02.Find_genes_for_clustering.R" is used for finding genes contributing most in distinguish homeostasis cell clusters.

script "Rscript/03.Construct_tSNE_map.R" is used for constructing homeostasis and transplantation cells differentiation tsne map.

script "Rscript/04.Classify_transplantation_cells.R" is used for classifying of transplantation cells to homeostasis cell clusters.

script "Rscript/05.Codes_for_Figure1.R" is used for generating Figures 1b-d and 1f.

script "Rscript/06.Codes_for_Figure2.R" is used for generating Figure 2.

script "Rscript/07.Codes_for_Figure3.R" is used for generating Figure 3c-g.

The file containing raw counts, normalized expression of genes on all cell before and after quality control is stored in directory named input.
