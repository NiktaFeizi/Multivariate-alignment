options(encoding = "UTF-8")
library(tidyr)
library(PharmacoGx)

# Download from zenodo (Uploaded to zenodo by Minoru)
diseaseCodes = c("ACC" , "BLCA",  "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH" , "KIRC", "KIRP", "LAML",
"LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "STAD", "SARC", "SKCM", "TGCT", "THCA", "THYM","UCEC","UCS", "UVM")

for (dis in diseaseCodes){
     print(dis)
     web_path <- paste("https://zenodo.org/record/5731017/files/TCGA_", paste(dis, ".rds?download=1", sep="") , sep="")
     save_path <- paste("~/nikta/Roche/Task1/Data/TCGA_Minoru", paste(dis, "_tcga.rds", sep="") , sep="/")
     download.file(web_path, save_path)
  }

# Download gene annotation (ensemble ids) 
url <-"http://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"
download.file(url, "hgnc_complete_set.txt")
hgnc.complete.set <- data.table::fread("hgnc_complete_set.txt") %>% as.data.frame() 

# Extracting expression profiles from TCGA data
data_extractor_tcga <- function(rds_file , hgnc = hgnc.complete.set , data_dir = "~/nikta/Roche/Task1/Data/TCGA_Minoru"){
  
  dataset <- readRDS(paste(data_dir, rds_file, sep="/"))
  
  all_types <- list(rnaseq = paste(sub("\\_.*", "",rds_file) , "_RNASeq2Gene-20160128", sep=""),
                    marray = paste(sub("\\_.*", "",rds_file) , "_mRNAArray-20160128", sep=""))
  
  res <- list()
  
  for (type in all_types){
    
    if (!is.null(dataset@ExperimentList@listData[[type]])){
      
      tcga_mat <- dataset@ExperimentList@listData[[type]] %>% 
        assay() %>% as.matrix() %>% 
        WGCNA::transposeBigData() 
      
      common_genes <- intersect(colnames(tcga_mat), hgnc[,"symbol"])
      tcga_mat <- tcga_mat[,common_genes]
      gene_ann <- dplyr::filter(hgnc, symbol %in% common_genes) 
      gene_ann <- gene_ann[!duplicated(gene_ann$symbol),]
      rownames(gene_ann) <- gene_ann$symbol
      gene_ann <- gene_ann[common_genes,]
      colnames(tcga_mat) <- ifelse(gene_ann$ensembl_gene_id == "", gene_ann$symbol, gene_ann$ensembl_gene_id) # where nsembl_gene_id is not available; gene symbol is used
      
      
      res[[type]] <- list(expression_mat = tcga_mat , gene_annotation = gene_ann)
    }              
  }
  return(res)
}

###################### All TCGA data in one matrix ###################### 
# Each project has potentially different platforms for microarrays (and some no microarrays at all).

# Create a seed matrix to append the data from other tissues:
temp <- data_extractor_tcga (rds_file = "READ_tcga.rds") # this one has both rna-seq and marray

TCGA_RNA_seq <- matrix(NA, ncol= ncol(temp[[1]][["expression_mat"]])) 
colnames(TCGA_RNA_seq) <- colnames(temp[[1]][["expression_mat"]])

TCGA_marray <- matrix(NA, ncol= ncol(temp[[2]][["expression_mat"]])) 
colnames(TCGA_marray) <- colnames(temp[[2]][["expression_mat"]])

rm(temp)

# Sample annotation: 
all_sample_annot <- data.frame(sampleID =NA , lineage=NA )


for (name in list.files(path = "~/nikta/Roche/Task1/Data/TCGA_Minoru")){
  
  res <- data_extractor_tcga (rds_file = name)
  
  TCGA_RNA_seq <- rbind(TCGA_RNA_seq, res[[1]][["expression_mat"]])
  
  sample_annot <- data.frame(sampleID = rownames(res[[1]][["expression_mat"]]) , lineage = sub("\\_.*", "", name))
  all_sample_annot <- rbind(all_sample_annot, sample_annot)
  
  if (length(res) == 2){
    
    TCGA_marray <- rbind(TCGA_marray, res[[2]][["expression_mat"]]) 
    sample_annot <- data.frame(sampleID = rownames(res[[2]][["expression_mat"]]) , lineage = sub("\\_.*", "", name)) # added 12/21/2021
    all_sample_annot <- rbind(all_sample_annot, sample_annot) # added 12/21/2021
    
  }
  
  print(name)
  
}

# Removing the first row which is all NA:
TCGA_RNA_seq <- TCGA_RNA_seq[-1, ]
TCGA_marray <- TCGA_marray[-1, ]
all_sample_annot <- all_sample_annot[-1, ]
all_sample_annot <- unique(all_sample_annot)


lineage_map <- data.frame("lineage"= c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH",
                                       "KIRC","KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD",
                                       "READ","SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"),
                          
                          "tissue_name"= c("Adrenal Gland", "Bladder/Urinary Tract", "Breast", "Cervix", "Biliary Tract", "Bowel",
                                           "Lymphoid", "Esophagus/Stomach", "CNS/Brain", "Head and Neck", "Kidney", "Kidney", "Kidney",
                                           "Myeloid", "CNS/Brain", "Liver", "Lung", "Lung", NA, "Ovary/Fallopian Tube", "Pancreas",
                                           "Peripheral Nervous System", "Prostate", "Bowel", "Bone", "Skin", "Esophagus/Stomach", "Testis",
                                           "Thyroid", "Thymus", "Uterus", "Uterus","Eye") )



all_sample_annot <- merge(lineage_map , all_sample_annot , by = "lineage", all=T)

qsave( TCGA_RNA_seq , "TCGA_RNA_seq.qs")
qsave( TCGA_marray, "TCGA_marray.qs")
qsave( all_sample_annot, "TCGA_all_sample_annot.qs")















