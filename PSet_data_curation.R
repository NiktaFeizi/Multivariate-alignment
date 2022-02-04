options(encoding = "UTF-8")
work_dir = "~/Data"
setwd(work_dir)
library(qs)
library(dplyr)
library(PharmacoGx)

# Curate expression data from PSets
data_extractor <- function(path){
  
  pset <- readRDS(path)
  
  
  # Both RNA-seq and microarray data:
  marray = NULL
  if("rna" %in% mDataNames(pset)){
    marray <- summarizeMolecularProfiles(pset, mDataType="rna") %>% assay() %>% t() 
    marray <- marray[rowSums(is.na(marray)) != ncol(marray), ] # Removing rows in which the gene expression is NA for all genes
    gene_ann_arr <- featureInfo(pset, mDataType = "rna") %>% data.frame() 
    
    if (!identical(grep("AFFX", colnames(marray)), integer(0))) { 
      marray <- marray[ , -grep("AFFX", colnames(marray))]
      gene_ann_arr <- gene_ann_arr[-grep("AFFX", rownames(gene_ann_arr)) , ]
      
    } 
    
  }
  
  # Removing duplications and filtering tissue types -> duplicates have "PAR_Y" at the end of the column name 
  exp_mat <- summarizeMolecularProfiles(pset, mDataType = "Kallisto_0.46.1.rnaseq") %>% assay() %>% t() 
  exp_mat <- exp_mat[ , -grep("PAR_Y", colnames(exp_mat))]
  exp_mat <- exp_mat[rowSums(is.na(exp_mat)) != ncol(exp_mat), ] # Removing rows in which the gene expression is NA for all genes
  
  gene_ann <- as.data.frame(featureInfo(pset, mDataType = "Kallisto_0.46.1.rnaseq")) 
  gene_ann <- gene_ann[colnames(exp_mat),] # Removing "PAR_Y" genes annotations
  
  # Stripping off the version number from ENSEMBLE IDs
  colnames(exp_mat) <- sub("\\..*", "", colnames(exp_mat))
  rownames(gene_ann) <- sub("\\..*", "", rownames(gene_ann))
  
  if(!is.null(marray)) {
    
    # Keeping only the tissue type with at least 20 cell lines
    tumor_types <- cellInfo(pset) %>%  filter(cellid %in% rownames(marray)) %>% group_by(tissueid) %>% dplyr::summarise(n = n()) %>% filter(n>=20) #tissue type with at least 20 cell lines
    cell_ann_marray <- cellInfo(pset) %>% filter(tissueid %in% tumor_types$tissueid) 
    marray <- marray[intersect(rownames(cell_ann_marray) , rownames(marray)), ] 
    cell_ann_marray <- cell_ann_marray[intersect(rownames(cell_ann_marray), rownames(marray)), ]
    
    tumor_types <- cellInfo(pset) %>%  filter(cellid %in% rownames(exp_mat)) %>% group_by(tissueid) %>% dplyr::summarise(n = n()) %>% filter(n>=20) #tissue type with at least 20 cell lines
    cell_ann_seq <- cellInfo(pset) %>% filter(tissueid %in% tumor_types$tissueid) 
    exp_mat <- exp_mat[intersect(rownames(cell_ann_seq ) , rownames(exp_mat)), ] 
    cell_ann_seq <- cell_ann_seq[intersect(rownames(cell_ann_seq), rownames(exp_mat)), ]
    
    
    return(list(marray_exp = marray, marray_gene_ann = gene_ann_arr, cell_ann_marray = cell_ann_marray, 
                rna_seq_exp = exp_mat, rna_seq_gene_ann = gene_ann, cell_ann_seq = cell_ann_seq))
  }
  
  else { 
    tumor_types <- cellInfo(pset) %>%  filter(cellid %in% rownames(exp_mat)) %>% group_by(tissueid) %>% dplyr::summarise(n = n()) %>% filter(n>=20) #tissue type with at least 20 cell lines
    cell_ann_seq <- cellInfo(pset) %>% filter(tissueid %in% tumor_types$tissueid) 
    exp_mat <- exp_mat[intersect(rownames(cell_ann_seq ) , rownames(exp_mat)), ] 
    cell_ann_seq <- cell_ann_seq[intersect(rownames(cell_ann_seq), rownames(exp_mat)), ]
    
    return(list(rna_seq_exp = exp_mat, rna_seq_gene_ann = gene_ann, cell_ann_seq = cell_ann_seq))
    
  }
}


# All PSets downloaded from ORCESTRA 

pset_names <- c("GDSCv2_u133a.rds", "GDSCv1_u133a.rds", "GDSCv2.rds", "GDSCv1.rds", "CCLE.rds","gCSI2.rds", "UHNBreast.rds") 
tissues <- c()

for (name in pset_names){
  
  res <- data_extractor( path = paste("~/PSets-ORCESTRA", name, sep = "/"))
  print(name)
  tissues <- c(tissues, unique(res[["cell_ann_seq"]]$tissueid) , unique(res[["cell_ann_marray"]]$tissueid))
  
  if(length(res) == 6){
    saveRDS(res[["marray_exp"]] , paste(work_dir, sub("\\..*", "_marray_exp.rds", name), sep="/"))
    saveRDS(res[["marray_gene_ann"]] , paste(work_dir, sub("\\..*", "_marray_gene_ann.rds", name), sep="/"))
    saveRDS(res[["cell_ann_marray"]] , paste(work_dir, sub("\\..*", "_cell_ann_marray.rds", name), sep="/"))
    saveRDS(res[["rna_seq_exp"]] , paste(work_dir, sub("\\..*", "_rna_seq_exp.rds", name), sep="/"))
    saveRDS(res[["rna_seq_gene_ann"]] , paste(work_dir, sub("\\..*", "_rna_seq_gene_ann.rds", name), sep="/"))
    saveRDS(res[["cell_ann_seq"]] , paste(work_dir, sub("\\..*", "_cell_ann_seq.rds", name), sep="/"))
    
    print(name)
    print(dim (res[[1]]))
    print(dim (res[[2]]))
    print(dim (res[[3]]))
    print(dim (res[[4]]))
    print(dim (res[[5]]))
    print(dim (res[[6]]))
    
  } else {
    
    saveRDS(res[["rna_seq_exp"]] , paste(work_dir, sub("\\..*", "_rna_seq_exp.rds", name), sep="/"))
    saveRDS(res[["rna_seq_gene_ann"]] , paste(work_dir, sub("\\..*", "_rna_seq_gene_ann.rds", name), sep="/"))
    saveRDS(res[["cell_ann_seq"]] , paste(work_dir, sub("\\..*", "_cell_ann_seq.rds", name), sep="/"))
    
    print(name)
    print(dim (res[[1]]))
    print(dim (res[[2]]))
    print(dim (res[[3]]))
  }
  
}

unique(tissues)


########## Removing breast from CCLE:

#### RNA-seq ####
ccle_ann <- readRDS("CCLE_cell_ann_seq.rds") %>% dplyr::filter(!tissueid == "Breast")
ccle_mat <- readRDS("CCLE_rna_seq_exp.rds")

ccle_mat <- ccle_mat[intersect(rownames(ccle_ann), rownames(ccle_mat)), ]
ccle_ann <- ccle_ann[intersect(rownames(ccle_ann), rownames(ccle_mat)), ]


saveRDS(ccle_ann, "CCLE_cell_ann_seq.rds")
saveRDS(ccle_mat, "CCLE_rna_seq_exp.rds")

#### Microarray ####
ccle_ann <- readRDS("CCLE_cell_ann_marray.rds") %>% dplyr::filter(!tissueid == "Breast")
ccle_mat <- readRDS("CCLE_marray_exp.rds")

ccle_mat <- ccle_mat[intersect(rownames(ccle_ann), rownames(ccle_mat)), ]
ccle_ann <- ccle_ann[intersect(rownames(ccle_ann), rownames(ccle_mat)), ]


saveRDS(ccle_ann, "CCLE_cell_ann_marray.rds")
saveRDS(ccle_mat, "CCLE_marray_exp.rds")


