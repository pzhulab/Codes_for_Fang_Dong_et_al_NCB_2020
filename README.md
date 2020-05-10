# Codes_for_Fang_Dong_et_al_NCB_2020

This repository provides the codes used to perform data anlsysis in our recently published article: 

F. Dong et al., Differentiation of transplanted haematopoietic stem cells tracked by single-cell transcriptomic analysis. Nat Cell Biol. 125, 1–10 (2020).

If you use the scripts, please cite the paper!

Note: Due to the size of the expression data, we are not allowed to upload the data files to github. Therefore this repository only shows you the codes we used to analyze the sequencing data in the paper. We are attempting to upload the data files by dividing them into small files, then you will completely reproduce the reuslts in our paper. Thanks for your interest in our work!

Details of scripts:

- Quality control: `01.Quality_control.R`

- Gene selection for cell clustering: `02.Find_genes_for_clustering.R`

- Dimension reduction and tSNE visualization of single cell data: `03.Construct_tSNE_map.R`

- Classification of cells after transplantation based on transcriptome reference constructed in homeostasis: `4.Classify_transplantation_cells.R`

- Figure 1b-d,f: `05.Codes_for_Figure1.R`

- Figure 2: `06.Codes_for_Figure2.R`

- Figure 3c-g: `07.Codes_for_Figure3.R`

Contact: Caiying Zhu (zhucaiying[At]ihcams.ac.cn)
