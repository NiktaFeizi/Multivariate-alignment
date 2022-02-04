#library(here)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(qs)

# ====== Global_params ======
global <- list(
  n_genes = 'all', # set to 'all' to use all protein coding genes found in both datasets
  umap_n_neighbors = 10, # num nearest neighbors used to create UMAP plot
  umap_min_dist = 0.5, # min distance used to create UMAP plot
  mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line data
  mnn_k_tumor = 50, # number of nearest neighbors of cell lines in the tumor data
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  # remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data
  distance_metric = 'euclidean', # distance metric used for the UMAP projection
  mod_clust_res = 5, # resolution parameter used for clustering the data
  mnn_ndist = 3, # ndist parameter used for MNN
  n_PC_dims = 70, # number of PCs to use for dimensionality reduction
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 10 # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
)

# ========================================== figure 1a ==========================================
plot_uncorrected_data <- function(org_dat, before_plot) {
  
  comb_ann <- cbind.data.frame(`sampleID` = c(rownames(org_dat$TCGA_mat), rownames(org_dat$pset_mat)),
                               `type` = c(rep('tumor', nrow(org_dat$TCGA_mat)), rep("CL", nrow(org_dat$pset_mat))))
  
  # using Seurat object to run cPCA and UMAP
  original_combined_obj <-  Seurat::CreateSeuratObject(t(rbind(org_dat$TCGA_mat,
                                                               org_dat$pset_mat)),
                                                       min.cells = 0,
                                                       min.features = 0,
                                                       meta.data = comb_ann %>%
                                                         magrittr::set_rownames(comb_ann$sampleID))
  
  original_combined_obj <- Seurat::ScaleData(original_combined_obj, 
                                             features = rownames(Seurat::GetAssayData(original_combined_obj)), 
                                             do.scale = F)
  
  original_combined_obj %<>% Seurat::RunPCA(assay='RNA',
                                            features = rownames(Seurat::GetAssayData(original_combined_obj)),
                                            npcs = global$n_PC_dims,
                                            verbose = F)
  
  original_combined_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                             reduction = 'pca',
                                             n.neighbors = global$umap_n_neighbors,
                                             min.dist = global$umap_min_dist,
                                             metric = global$distance_metric,
                                             verbose=F)
  
  uncorrected_alignment <- Seurat::Embeddings(original_combined_obj, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID')
  
  ggplot2::ggplot(uncorrected_alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, type=='tumor'), alpha=0.6, size=0.5, pch=21, color='white', aes(fill=type)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, type=='CL'), alpha=0.6, size=0.6, pch=3, aes(color=type), stroke=0.5) +
    ggplot2::scale_color_manual(values=c(CL="#F8766D")) +
    ggplot2::scale_fill_manual(values=c(tumor="#00BFC4")) +
    ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='bottom',
                   text = ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size=6),
                   axis.title = ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0), 
                   legend.box.margin=ggplot2::margin(-10,-30,-10,-30),
                   axis.line = ggplot2::element_line(size = .3))
  
  ggsave(paste(paste("~/Results/", before_plot, sep="") ,".png" , sep=""), width = 7, height = 5)
  
  
}   

# ========================================== figure 2  ==========================================
Celligner_alignment_plot <- function(aligned_data, after_plot) {
  
  aligned_data %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                    reduction = 'pca',
                                    n.neighbors = global$umap_n_neighbors,
                                    min.dist = global$umap_min_dist,
                                    metric = global$distance_metric,
                                    verbose=F)
  
  comb_ann <- aligned_data[[]] # This contains all the metadata
  
  lineage_map <- data.frame("lineage"= c("BLCA","BRCA", "COAD", "DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP",
                                         "LAML","LGG","LIHC","LUAD","LUSC","OV","PAAD","PCPG","PRAD","READ","SARC",
                                         "SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS"),
                            
                            "tissue_name"= c("Bladder/Urinary Tract", "Breast", "Bowel","Lymphoid","Esophagus/Stomach",
                                             "CNS/Brain","Head and Neck","Kidney","Kidney","Kidney","Myeloid","CNS/Brain",
                                             "Liver","Lung","Lung","Ovary/Fallopian Tube","Pancreas","Peripheral Nervous System",
                                             "Prostate","Bowel","Bone","Skin","Esophagus/Stomach","Testicular","Thyroid","Thymoma",
                                             "Uterus","Uterus" ) )
  
  alignment <- Seurat::Embeddings(aligned_data, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID') %>%
    left_join(lineage_map, by = 'lineage')
  
  alignment$tissue_name <- ifelse(is.na(alignment$tissue_name), alignment$lineage, alignment$tissue_name)                              
  
  
  ggplot2::ggplot(alignment,
                  ggplot2::aes(UMAP_1, UMAP_2, fill=tissue_name, size=type, color = type)) +
    ggplot2::geom_point(pch=21, alpha=0.7)  +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=1, `tumor`=0.75)) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'bottom', 
                   text=ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0)) +
    #ggplot2::guides(fill=FALSE, color=FALSE) +
    #ggplot2::scale_fill_manual(values= seq(1: length(unique(alignment$tissue_name)))) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  ggsave(paste(paste("~/Results/", after_plot, sep="") ,".png" , sep=""), width = 8, height = 5)

}

# ============================= Both figures =============================

figs <- function( org_dat, before_plot, aligned_data, after_plot){
  
  plot_uncorrected_data ( org_dat, before_plot)
  Celligner_alignment_plot(aligned_data, after_plot) 
}

# ============================= ccle+uhn_rnaseq =============================
# original data from Task1_data_curation
# aligned data from Task1_data_Method

figs ( org_dat = qread("CCLE_UHN_tcga_rnaseq.qs")  , 
       before_plot = "ccle_before_alignment", 
       aligned_data = qread("ccle_uhn_rnaseq_celligner.qs") , 
       after_plot = "ccle_after_alignment")

# ============================= gCSI_rnaseq =============================

figs ( org_dat = qread("gCSI_tcga_rnaseq.qs")  , 
       before_plot = "gCSI_before_alignment", 
       aligned_data = qread("gCSI_rnaseq_celligner.qs") , 
       after_plot = "gCSI_after_alignment")

# ============================= GDSC_U133a _ microarray =============================

figs ( org_dat = qread("GDSC_U133a_tcga_marray.qs")  , 
       before_plot = "GDSC_U133a_before_alignment", 
       aligned_data = qread("GDSC_U133_marray_celligner_pca2.qs") , 
       after_plot = "GDSC_U133a_after_alignment_2pca")

# ============================= GDSC_U219 _ microarray =============================

figs ( org_dat = qread("GDSC_U219_tcga_marray.qs")  , 
       before_plot = "GDSC_U219_before_alignment", 
       aligned_data = qread("GDSC_U219_marray_celligner.qs") , 
       after_plot = "GDSC_U219_after_alignment")




