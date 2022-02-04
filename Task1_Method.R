# Task 1 
# https://github.com/broadinstitute/Celligner_ms/blob/master/src/Celligner_methods.R

options(encoding = "UTF-8" , timeout = 1000)
work_dir = "~/Data"
setwd(work_dir)
library(data.table)
library(magrittr)

# ========================================== github codes ==========================================
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

# ====== calculate gene stats function ======
calc_gene_stats <- function(dat, data_dir, hgnc_file) { 
  common_genes <- intersect(colnames(dat$TCGA_mat), colnames(dat$pset_mat))
  
  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  hgnc.complete.set <- hgnc.complete.set %>% 
    dplyr::select(Gene = ensembl_gene_id, Symbol = symbol) %>%
    dplyr::filter(Gene %in% common_genes)  # Nikta edited
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$Gene),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$Gene
  hgnc.complete.set <- hgnc.complete.set[common_genes,]  
  
  gene_stats <- data.frame(
    Tumor_SD = apply(dat$TCGA_mat, 2, sd, na.rm=T),
    pset_SD = apply(dat$pset_mat, 2, sd, na.rm=T),
    Tumor_mean = colMeans(dat$TCGA_mat, na.rm=T),
    pset_mean = colMeans(dat$pset_mat, na.rm=T),
    Gene = common_genes,
    stringsAsFactors = F) %>% 
    dplyr::mutate(max_SD = pmax(Tumor_SD, pset_SD, na.rm=T)) #add avg and max SD per gene
  
  gene_stats <- dplyr::left_join(hgnc.complete.set, gene_stats, by = "Gene") # Nikta edited
  
  return(gene_stats)
  
}

# ====== create Seurat object of expression data function ======
create_Seurat_object <- function(exp_mat, ann, type = NULL) { ######## I GET TWO WARNINGS HERE!!!!!!!!!!!!!!!!!
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat),
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann %>%
                                          magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  
  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = rownames(Seurat::GetAssayData(seu_obj)),
                              npcs = global$n_PC_dims, verbose = F)
  
  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                               reduction = 'pca',
                               n.neighbors = global$umap_n_neighbors,
                               min.dist =  global$umap_min_dist,
                               metric = global$distance_metric, verbose=F)
  
  return(seu_obj)
}


# ====== Seurat clustering function ======
cluster_data <- function(seu_obj) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = 'pca',
                                   dims = 1:global$n_PC_dims,
                                   k.param = 20, 
                                   force.recalc = TRUE,
                                   verbose = FALSE)
  
  seu_obj %<>% Seurat::FindClusters(reduction = 'pca', 
                                    resolution = global$mod_clust_res)
  
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  
  return(seu_obj)
  
}

# ====== Differentially expressed genes function ======
# Estimate linear-model stats for a matrix of data with respect to a group of phenotype variables
# using limma with empirical Bayes moderated F-stats for p-values
run_lm_stats_limma_group <- function (mat, phenos, covars = NULL, weights = NULL, target_type = "Gene", 
                                      limma_trend = FALSE) 
{
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  udata <- rownames(mat) %>% intersect(rownames(phenos))
  if (!is.null(covars)) {
    udata %<>% intersect(rownames(covars))
  }
  form <- as.formula(paste("~", paste0(colnames(phenos), collapse = " + ")))
  design <- model.matrix(form, data = phenos[udata, , drop = F])
  if (!is.null(covars)) {
    covars <- data.frame(covars)
    form <- as.formula(paste("~", paste0(colnames(covars), 
                                         collapse = " + ")))
    Cdesign <- model.matrix(form, data = covars[udata, , 
                                                drop = F])
    Cdesign <- Cdesign[, setdiff(colnames(Cdesign), "(Intercept)"), 
                       drop = FALSE]
    stopifnot(length(intersect(colnames(Cdesign), colnames(design))) == 
                0)
    design %<>% cbind(Cdesign)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  design <- design[, colSums(design) > 2, drop = FALSE]
  targ_coefs <- setdiff(colnames(design), "(Intercept)")
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- which(colnames(design) %in% targ_coefs)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf, 
                             sort.by = "F", genelist = colnames(mat))
  results %<>% tibble::rownames_to_column(var = target_type)
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(AveExpr = "Avg", 
                                                             F = "F_stat", P.Value = "p.value", adj.P.Val = "q.value"))) %>% 
    na.omit() %>% dplyr::select(-ProbeID)
  return(results)
}

# ====== find_differentially_expressed_genes ====== 

find_differentially_expressed_genes <- function(seu_obj) {
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(
      t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
      seu_obj@meta.data %>% dplyr::select(seurat_clusters),
      limma_trend = TRUE) %>%
      dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
                                       seu_obj@meta.data$cluster,
                                       limma_trend = TRUE) %>%
      dplyr::mutate(gene_stat = abs(t_stat)) %>%
      dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }
  
  return(cur_DE_genes)
  
}

# ====== cPCA function  ======
run_cPCA_analysis <- function(TCGA_dat, pset_dat, tumor_cluster_df, CL_cluster_df, pc_dims=NULL) {
  tumor_clust_avgs <- get_cluster_averages(TCGA_dat, tumor_cluster_df)
  CL_clust_avgs <- get_cluster_averages(pset_dat, CL_cluster_df)
  
  TCGA_subtype_ms <- TCGA_dat - tumor_clust_avgs[tumor_cluster_df$seurat_clusters,]
  pset_subtype_ms <- pset_dat - CL_clust_avgs[CL_cluster_df$seurat_clusters,]
  
  TCGA_cov <- cov(TCGA_subtype_ms)
  pset_cov <- cov(pset_subtype_ms)
  
  if(!is.null(pc_dims)) {
    cov_diff_eig <- irlba::prcomp_irlba(TCGA_cov - pset_cov, n = pc_dims)
  } else {
    cov_diff_eig <- eigen(TCGA_cov - pset_cov)
  }
  return(cov_diff_eig)
}

# ====== calculate the average expression per cluster  ======
get_cluster_averages <- function(mat, cluster_df) {
  n_clusts <- nlevels(cluster_df$seurat_clusters)
  clust_avgs <- matrix(NA, nrow = n_clusts, ncol = ncol(mat)) %>% 
    magrittr::set_colnames(colnames(mat)) %>% 
    magrittr::set_rownames(levels(cluster_df$seurat_clusters))
  for (ii in levels(cluster_df$seurat_clusters)) {
    clust_avgs[ii,] <- colMeans(mat[cluster_df$seurat_clusters == ii,], na.rm=T)
  }
  return(clust_avgs)
}

# ====== run_cPCA ======
run_cPCA <- function(TCGA_obj, pset_obj, pc_dims = NULL) {
  cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data')), 
                                    t(Seurat::GetAssayData(pset_obj, assay='RNA', slot='scale.data')), 
                                    TCGA_obj@meta.data, pset_obj@meta.data, pc_dims=pc_dims)
  return(cov_diff_eig) 
}

# ====== run mutual nearest neighbors batch correction  ======
# Allows for separate k values per dataset, and simplifies some of the IO and doesn't use PCA reduction
modified_mnnCorrect <- function(ref_mat, targ_mat, k1 = 20, k2 = 20, 
                                ndist = 3, subset_genes = NULL) {
  if (is.null(subset_genes)) {
    subset_genes <- colnames(ref_mat) 
  }  
  
  sets <- batchelor::findMutualNN(ref_mat[, subset_genes], 
                                  targ_mat[, subset_genes], 
                                  k1 = k2, k2 = k1, 
                                  BPPARAM = BiocParallel::SerialParam())
  mnn_pairs <- as.data.frame(sets) %>% 
    dplyr::mutate(ref_ID = rownames(ref_mat)[first],
                  targ_ID = rownames(targ_mat)[second],
                  pair = seq(nrow(.))) %>% 
    dplyr::select(-first, -second)
  
  # Estimate the overall batch vector.
  ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  overall.batch <- colMeans(ave.out$averaged)
  
  #remove variation along the overall batch vector
  ref_mat <- .center_along_batch_vector(ref_mat, overall.batch)
  targ_mat <- .center_along_batch_vector(targ_mat, overall.batch)
  
  # Recompute correction vectors and apply them.
  re.ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  targ_mat <- .tricube_weighted_correction(targ_mat, re.ave.out$averaged, re.ave.out$second, k=k2, ndist=ndist, subset_genes, BPPARAM=BiocParallel::SerialParam())
  
  final <- list(corrected = targ_mat, 
                pairs = mnn_pairs)
  return(final)
}

# ====== .average_correction  ======
# Copied from dev version of scran (2018-10-28) with slight modifications as noted
#https://github.com/MarioniLab/scran
.average_correction <- function(refdata, mnn1, curdata, mnn2)
  # Computes correction vectors for each MNN pair, and then
  # averages them for each MNN-involved cell in the second batch.
{
  corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
  corvec <- rowsum(corvec, mnn2)
  npairs <- table(mnn2)
  stopifnot(identical(names(npairs), rownames(corvec)))
  corvec <- unname(corvec)/as.vector(npairs)
  list(averaged=corvec, second=as.integer(names(npairs)))
}

# ====== .center_along_batch_vector  ======
.center_along_batch_vector <- function(mat, batch.vec) 
  # Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
  # This removes any variation along the overall batch vector within each matrix.
{
  batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
  batch.loc <- as.vector(mat %*% batch.vec)
  central.loc <- mean(batch.loc)
  mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
  return(mat)
}

# ====== .tricube_weighted_correction  ======
#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, subset_genes, BNPARAM=NULL, BPPARAM=BiocParallel::SerialParam())
  # Computing tricube-weighted correction vectors for individual cells,
  # using the nearest neighbouring cells _involved in MNN pairs_.
  # Modified to use FNN rather than queryKNN for nearest neighbor finding
{
  cur.uniq <- curdata[in.mnn,,drop=FALSE]
  safe.k <- min(k, nrow(cur.uniq))
  # closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
  closest <- FNN::get.knnx(cur.uniq[, subset_genes], query=curdata[, subset_genes], k=safe.k)
  # weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
  weighted.correction <- .compute_tricube_average(correction, closest$nn.index, closest$nn.dist, ndist=ndist)
  curdata + weighted.correction
}

# ====== .compute_tricube_average ======
.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3) 
  # Centralized function to compute tricube averages.
  # Bandwidth is set at 'ndist' times the median distance, if not specified.
{
  if (is.null(bandwidth)) {
    middle <- ceiling(ncol(indices)/2L)
    mid.dist <- distances[,middle]
    bandwidth <- mid.dist * ndist
  }
  bandwidth <- pmax(1e-8, bandwidth)
  
  rel.dist <- distances/bandwidth
  rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
  tricube <- (1 - rel.dist^3)^3
  weight <- tricube/rowSums(tricube)
  
  output <- 0
  for (kdx in seq_len(ncol(indices))) {
    output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
  }
  
  if (is.null(dim(output))) {
    matrix(0, nrow(vals), ncol(vals))
  } else {
    output
  }
}


# ====== run_MNN ======
run_MNN <- function(pset_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                    subset_genes) {
  mnn_res <- modified_mnnCorrect(pset_cor, TCGA_cor, k1 = k1, k2 = k2, ndist = ndist, 
                                 subset_genes = subset_genes)
  
  return(mnn_res)
}


# ====== cor between tumors and cell lines function  ======
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- dplyr::filter(Celligner_info, type=='tumor')$sampleID
  cl_samples <- dplyr::filter(Celligner_info, type=='CL')$sampleID
  tumor_CL_cor <- cor(t(Celligner_aligned_data[tumor_samples,]), t(Celligner_aligned_data[cl_samples,]),
                      use='pairwise')
  
  
  return(tumor_CL_cor)
}


# ====== run all Celligner methods  ======

run_Celligner <- function(data_dir, dat, remove_cPCA_dims = c(1,2,3,4), plotname) {
  
  dat <- qread(dat)
  gene_stats <- calc_gene_stats(dat, data_dir, hgnc_file = 'hgnc_complete_set.txt') # UPDATED hgnc FILE -> EDIT BY NIKTA
  
  comb_ann <- rbind(
    dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary.Metastasis`) %>%
      dplyr::mutate(type = 'tumor'),
    dat$pset_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary.Metastasis`) %>%
      dplyr::mutate(type = 'CL')
  )
  
  TCGA_obj <- create_Seurat_object(exp_mat=dat$TCGA_mat, ann=dat$TCGA_ann, type='tumor')
  pset_obj <- create_Seurat_object(dat$pset_mat, dat$pset_ann, type='CL')
  
  TCGA_obj <- cluster_data(seu_obj = TCGA_obj)
  pset_obj <- cluster_data(pset_obj)
  
  tumor_DE_genes <- find_differentially_expressed_genes(seu_obj = TCGA_obj) # Warning message: Zero sample variances detected, have been offset away from zero 
  CL_DE_genes <- find_differentially_expressed_genes(pset_obj) # Warning message: Zero sample variances detected, have been offset away from zero 
  
  DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(
      tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
    dplyr::left_join(gene_stats, by = 'Gene')
  
  # take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < global$top_DE_genes_per) %>%
    .[['Gene']]
  
  
  cov_diff_eig <- run_cPCA(TCGA_obj, pset_obj, global$fast_cPCA)
  
  jpeg(paste(paste("~/nikta/Roche/Task1/Results/", plotname, sep="") ,".png" , sep=""))
  png <- barplot(cov_diff_eig$sdev, main = plotname)
  dev.off()
  
  if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, remove_cPCA_dims, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, remove_cPCA_dims, drop = FALSE]
  }
  
  rownames(cur_vecs) <- colnames(dat$TCGA_mat)
  TCGA_cor <- resid(lm(t(dat$TCGA_mat[complete.cases(dat$TCGA_mat),]) ~ 0 + cur_vecs)) %>% t() #nikta edited
  pset_cor <- resid(lm(t(dat$pset_mat) ~ 0 + cur_vecs)) %>% t()
  
  mnn_res <- run_MNN(pset_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist,
                     subset_genes = DE_gene_set)
  
  combined_mat <- rbind(mnn_res$corrected, pset_cor)
  
  comb_obj <- create_Seurat_object(combined_mat, comb_ann)
  comb_obj <- cluster_data(seu_obj = comb_obj)
  
  return(comb_obj) 
}

# ========================================== RNA-seq ==========================================
#===== CCLE =====
CCLE_UHN_rnaseq <- run_Celligner(data_dir = work_dir , dat = "CCLE_UHN_tcga_rnaseq.qs",
                                 remove_cPCA_dims = 1 , plotname = "CCLE_rnaseq")

qsave(CCLE_UHN_rnaseq, "ccle_uhn_rnaseq_celligner.qs")

#===== gCSI =====
gCSI_rnaseq <- run_Celligner(data_dir = work_dir , dat = "gCSI_tcga_rnaseq.qs", remove_cPCA_dims = 1 ,
                             plotname = "gCSI_rnaseq")

qsave(gCSI_rnaseq, "gCSI_rnaseq_celligner.qs")


# 1000 significant genes to be used for kmeans => run cell-ligner function up to "DE_gene_set" and save the genes
saveRDS(DE_gene_set[c(1:1000)], "1000_DE_gcsi_genes.rds")


# ========================================== Microarray ==========================================
#===== GDSC_U133 ==========
GDSC_U133_marray <- run_Celligner(data_dir = work_dir , dat = "GDSC_U133a_tcga_marray.qs",
                                  remove_cPCA_dims = c(1,2) , plotname = "GDSC_U133_marray")

qsave(GDSC_U133_marray , "GDSC_U133_marray_celligner.qs")

#===== GDSC2 =====
GDSC_U219_marray <- run_Celligner(data_dir = work_dir , dat = "GDSC_U219_tcga_marray.qs",
                                  remove_cPCA_dims = 1 , plotname = "GDSC_U219_marray")

qsave(GDSC_U219_marray , "GDSC_U219_marray_celligner.qs")



