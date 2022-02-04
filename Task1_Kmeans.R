# Running kmeans for different distances

## packages
library(factoextra)
library(fpc)
library(NbClust)


ploter= function(mat, ann, data_name = "gSCI_bfore",  distance = "euclidean", genes = "common_genes"){
  
  km.res <- eclust(mat, "kmeans", k = length(unique(ann[["lineage"]])), nstart = 25, graph = FALSE, hc_metric = distance) 
  
  cluster = merge(data.frame(km.res$cluster),  ann %>% `row.names<-`(., NULL)  %>% tibble::column_to_rownames('sampleID'), by = "row.names", all.x=T)
  cluster$km.res.cluster <- factor(cluster$km.res.cluster)
  
  plot_title <- paste(paste(data_name, distance, sep="_"), genes, sep= "_")
  ggplot(data = cluster, aes(y = km.res.cluster)) +
    geom_bar(aes(fill = lineage), width = 0.5) +
    ggtitle(plot_title) + xlab("Sample") + ylab("Cluster")
  theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x=element_text(size=12,  face="bold"),
          axis.title=element_text(size=12,face="bold"),
          axis.text.y=element_text(size=10, face = "bold"),
          strip.text = element_text(size=10, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position="bottom",
          legend.text = element_text(size = 8, face="bold"),
          legend.title = element_blank()) 
  ggsave(paste(paste("~/Results/Kmeans/", plot_title, sep="") ,"kmeans.png" , sep="_"), width = 7, height = 5)
  
  
  fviz_silhouette(km.res, main = plot_title) 
  ggsave(paste(paste("~/Results/Kmeans/", plot_title, sep="") ,"silhouette.png" , sep="_"), width = 7, height = 5)
  
}



data_loader <- function(before_align, after_align, data_name ,  distance , genes, sig_genes = NA) {
  
  before_data = qread(before_align) 
  after_data = qread(after_align) 
  
  
  
  #================ Before alignment plots
  
  if(! is.na(sig_genes)) {
    
    before_data$TCGA_mat = before_data$TCGA_mat [,sig_genes]
    before_data$pset_mat = before_data$pset_mat [,sig_genes]
    
  }
  
  pset_before <- ploter (mat = before_data$pset_mat , 
                         ann = before_data$pset_ann ,
                         data_name = paste(data_name, "before", sep="_") ,
                         distance= distance ,
                         genes = genes)
  
  
  TCGA_before <- ploter (mat = before_data$TCGA_mat , 
                         ann = before_data$TCGA_ann ,
                         data_name = paste(data_name, "TCGA_before", sep="_")  ,
                         distance= distance ,
                         genes = genes)
  
  
  all_before <- ploter ( mat = rbind(before_data$TCGA_mat , before_data$pset_mat)  , 
                         ann = rbind(before_data$TCGA_ann , before_data$pset_ann) ,
                         data_name = paste(data_name, "all_before", sep="_")  ,
                         distance= distance ,
                         genes = genes)
  
  gc()
  #================ After alignment plots
  
  all_mat <- after_data@assays[["RNA"]]@counts %>% as.data.frame() 
  all_ann <- after_data[[c( "sampleID", "lineage", "subtype")]]  %>% as.data.frame() 
  
  
  if(! is.na(sig_genes)) { all_mat = all_mat [sig_genes, ] }
  
  pset_after <- ploter ( mat = t(all_mat[, !grepl("TCGA", colnames(all_mat))]) , 
                         ann = all_ann[!grepl("TCGA", all_ann$sampleID),] ,
                         data_name = paste(data_name, "after", sep="_")  ,
                         distance= distance ,
                         genes = genes)
  
  TCGA_after <- ploter ( mat = t(all_mat[, grepl("TCGA", colnames(all_mat))]) , 
                         ann = all_ann[grepl("TCGA", all_ann$sampleID),] ,
                         data_name = paste(data_name, "TCGA_after", sep="_") ,
                         distance= distance ,
                         genes = genes)
  
  all_after <- ploter ( mat = t(all_mat) , 
                        ann = all_ann,
                        data_name = paste(data_name, "all_after", sep="_")  ,
                        distance= distance ,
                        genes = genes)
}


# Before alignment data from Task1_data_curation
# After alignment data from Task1_data_Method

data_loader(before_align = "gCSI_tcga_rnaseq.qs", 
            after_align = "gCSI_rnaseq_celligner.qs",
            data_name = "gCSI",  distance = "euclidean", genes = "common_genes") 


data_loader(before_align = "gCSI_tcga_rnaseq.qs", 
            after_align = "gCSI_rnaseq_celligner.qs",
            data_name = "gCSI",  distance = "pearson", genes = "common_genes") 


# Run for a subset of genes (genes with highest variability) # created in Task1_Method
DE_gene_set <- readRDS("1000_DE_gcsi_genes.rds")

data_loader(before_align = "gCSI_tcga_rnaseq.qs", 
            after_align = "gCSI_rnaseq_celligner.qs",
            data_name = "gCSI",  distance = "euclidean", genes = "1000_genes", sig_genes = DE_gene_set) 




