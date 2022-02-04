options(encoding = "UTF-8")
work_dir = "~/Data"
setwd(work_dir)
# library(data.table)
# library(magrittr)

# =========== load_data function =========== 
load_data <- function(data_dir , tumor_file, cell_line_file, 
                      tcga_cell_ann, hgnc_file, tumor_types , pset_cell_ann , type="microarray") {
  
  hgnc.complete.set <- data.table::fread(hgnc_file) %>% as.data.frame() # Gene annotation (ensemble ids)
  
  TCGA_mat <- qread(file.path(data_dir, tumor_file))
  pset_mat <-  readRDS(file.path(data_dir, cell_line_file)) 
  
  # Removing genes with any NAs:
  if(any(is.na(TCGA_mat))) {
    TCGA_mat <- TCGA_mat[, -which(apply(TCGA_mat , 2, function(x) any(is.na(x))))]
  }
  
  if(any(is.na(pset_mat))) {
    pset_mat <- pset_mat[, -which(apply(pset_mat , 2, function(x) any(is.na(x))))]
  }
  
  TCGA_ann <- qread(file.path(data_dir, tcga_cell_ann)) %>% as.data.frame()
  pset_ann <- readRDS(file.path(data_dir, pset_cell_ann)) %>% as.data.frame()
  
  ann <- data.frame(sampleID = c(TCGA_ann$sampleID, pset_ann$cellid),
                    lineage = c(TCGA_ann$tissue_name, pset_ann$tissueid),
                    subtype = NA,
                    `Primary/Metastasis` =NA,
                    type = c(rep('tumor', nrow(TCGA_ann)), rep('CL', nrow(pset_ann))))
  
  
  TCGA_ann <- dplyr::filter(ann, type=='tumor')
  pset_ann <- dplyr::filter(ann, type=='CL')
  
  func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id # hgnc.complete.set downloaded through "TCGA_data_curation.R"
  genes_used <- intersect(colnames(TCGA_mat), colnames(pset_mat))
  genes_used <- intersect(genes_used, func_genes)
  
  TCGA_mat <- TCGA_mat[,genes_used]
  pset_mat <- pset_mat[,genes_used]
  
  # Keeping TCGA tumors intersectig with tumor types in the pset:
  # TCGA_ann  <- dplyr::filter(TCGA_ann , lineage %in% tumor_types)
  TCGA_ann  <- dplyr::filter(TCGA_ann , lineage %in% intersect(pset_ann$lineage , TCGA_ann$lineage))
  TCGA_mat  <- TCGA_mat[intersect(rownames(TCGA_mat), TCGA_ann $sampleID) , ] 
  
  if (type == "rnaseq"){ # Log2 is required only for RNAseq
    TCGA_mat  <- log2(TCGA_mat+0.001) }
  
  TCGA_ann  <- filter(TCGA_ann, sampleID %in% intersect(rownames(TCGA_mat), TCGA_ann $sampleID))
  
  # Only keeping intersecting tissue types with TCGA data
  pset_ann  <- dplyr::filter(pset_ann , lineage %in% intersect(pset_ann$lineage , TCGA_ann$lineage))
  pset_mat  <- pset_mat[intersect(rownames(pset_mat), pset_ann $sampleID) , ] 
  #pset_ann  <- filter(pset_ann, sampleID %in% intersect(rownames(pset_mat), pset_ann $sampleID))
  
  return(list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, pset_mat = pset_mat, pset_ann = pset_ann))
}

# ================================== RNA-seq ==================================
# NOTE: ALL GDSC PSETs HAVE AN IDENTICAL RNA-SEQ CONTENT. ONLY MICROARRAY IS DEIFERENT BETWEEN THE TWO ARRAY VRERSIONS; V1 AND V2 REFFER TO 
# DIFFERENT SENSITIVITY PROFILES ONLY.

CCLE_rnaseq <- load_data (data_dir = "~/Data", 
                          tumor_file = "TCGA_RNA_seq.qs", 
                          cell_line_file = "CCLE_rna_seq_exp.rds", 
                          tcga_cell_ann = "TCGA_all_sample_annot.qs", 
                          hgnc_file = "~/Data/hgnc_complete_set.txt",
                          pset_cell_ann = "CCLE_cell_ann_seq.rds",
                          type = "rnaseq")# Breast is removed 


GDSC_rnaseq <- load_data (data_dir = "~/Data", 
                          tumor_file = "TCGA_RNA_seq.qs", 
                          cell_line_file = "GDSCv1_rna_seq_exp.rds", 
                          tcga_cell_ann = "TCGA_all_sample_annot.qs", 
                          hgnc_file = "~/Data/hgnc_complete_set.txt",
                          pset_cell_ann = "GDSCv1_cell_ann_seq.rds",
                          type = "rnaseq")

gCSI_rnaseq <- load_data (data_dir = "~/Data", 
                          tumor_file = "TCGA_RNA_seq.qs", 
                          cell_line_file = "gCSI2_rna_seq_exp.rds", 
                          tcga_cell_ann = "TCGA_all_sample_annot.qs", 
                          hgnc_file = "~/Data/hgnc_complete_set.txt",
                          pset_cell_ann = "gCSI2_cell_ann_seq.rds",
                          type = "rnaseq") 

UHNBreast_rnaseq <- load_data (data_dir = "~/Data", 
                               tumor_file = "TCGA_RNA_seq.qs", 
                               cell_line_file = "UHNBreast_rna_seq_exp.rds", 
                               tcga_cell_ann = "TCGA_all_sample_annot.qs", 
                               hgnc_file = "~/Data/hgnc_complete_set.txt",
                               pset_cell_ann = "UHNBreast_cell_ann_seq.rds",
                               type = "rnaseq")



# ================================== Save RNA-seq ==================================
qsave(CCLE_rnaseq, "~/Data/CCLE_tcga_rnaseq.qs")
qsave(GDSC_rnaseq, "~/Data/GDSC_tcga_rnaseq.qs")
qsave(gCSI_rnaseq, "~/Data/gCSI_tcga_rnaseq.qs")
qsave(UHNBreast_rnaseq, "~/Data/UHNBreast_tcga_rnaseq.qs")

# ================================== Microarray ==================================
# identical(GDSCv1@molecularProfiles[["rna"]], GDSCv2@molecularProfiles[["rna"]]) #TRUE

CCLE_marray <- load_data(data_dir = "~/Data",
                         tumor_file = "TCGA_marray.qs",
                         cell_line_file = "CCLE_marray_exp.rds", 
                         tcga_cell_ann = "TCGA_all_sample_annot.qs",
                         hgnc_file = "~/Data/hgnc_complete_set.txt",
                         pset_cell_ann = "CCLE_cell_ann_marray.rds") 

GDSC_U219_marray <- load_data(data_dir = "~/Data",
                              tumor_file = "TCGA_marray.qs",
                              cell_line_file = "GDSCv1_marray_exp.rds", 
                              tcga_cell_ann = "TCGA_all_sample_annot.qs",
                              hgnc_file = "~/Data/hgnc_complete_set.txt",
                              pset_cell_ann = "GDSCv1_cell_ann_marray.rds") 

GDSC_U133a_marray <- load_data(data_dir = "~/Data",
                               tumor_file = "TCGA_marray.qs",
                               cell_line_file = "GDSCv1_u133a_marray_exp.rds", 
                               tcga_cell_ann = "TCGA_all_sample_annot.qs",
                               hgnc_file = "~/Data/hgnc_complete_set.txt",
                               pset_cell_ann = "GDSCv1_u133a_cell_ann_marray.rds") 

# ================================== Save Microarray ==================================
qsave(CCLE_marray, "~/Data/CCLE_tcga_marray.qs")
qsave(GDSC_U219_marray, "~/Data/GDSC_U219_tcga_marray.qs")
qsave(GDSC_U133a_marray, "~/Data/GDSC_U133a_tcga_marray.qs")

# ================================== Merge CCLE with UHNBreast for RNAseq ==================================
CCLE <- qread( "~/Data/CCLE_tcga_rnaseq.qs")
UHN <- qread( "~/Data/UHNBreast_tcga_rnaseq.qs")

TCGA_mat = rbind(CCLE[["TCGA_mat"]] , UHN[["TCGA_mat"]]) 
TCGA_ann = rbind(CCLE[["TCGA_ann"]] , UHN[["TCGA_ann"]]) 
pset_mat = rbind(CCLE[["pset_mat"]] , UHN[["pset_mat"]]) 
pset_ann = rbind(CCLE[["pset_ann"]] , UHN[["pset_ann"]]) 

CCLE_UHN <- list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, pset_mat = pset_mat, pset_ann = pset_ann)
qsave(CCLE_UHN, "~/Data/CCLE_UHN_tcga_rnaseq.qs")

# ================================== ccle- Keep mutual tumour types between ccle and gCSI ==================================
CCLE <- qread( "~/Data/CCLE_tcga_rnaseq.qs")
gCSI <- qread( "~/Data/gCSI_tcga_rnaseq.qs")

TCGA_ann <- filter(CCLE[["TCGA_ann"]], lineage %in% gCSI[["TCGA_ann"]]$lineage)
TCGA_mat <- CCLE[["TCGA_mat"]][TCGA_ann$sampleID, ]

pset_ann <- filter(CCLE[["pset_ann"]], lineage %in% gCSI[["pset_ann"]]$lineage)
pset_mat <- CCLE[["pset_mat"]][pset_ann$sampleID, ]

CCLE_mutual_gCSI <- list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, pset_mat = pset_mat, pset_ann = pset_ann)
qsave(CCLE_mutual_gCSI, "~/Data/CCLE_mutual_gCSI_tcga_rnaseq.qs")







