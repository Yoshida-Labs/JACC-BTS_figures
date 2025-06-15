                                
#                                                                      
# 1.  User parameters --------------------------------------------------
DATA_DIR  <- "~/analysis_flexRNA/outs/per_sample_outs"  # raw 10X dirs
OUT_DIR   <- "~/AZD4547_scRNA_results"                                       # all output
SAMPLE_TBL <- data.frame(
  sample_id  = c("Ctrl","NS","AZD"),
  condition  = c("Control","Pump+NS","Pump+AZDtreated"),
  path       = file.path(DATA_DIR,c("FAC989","FAC990","FAC991"),
                         "count","sample_filtered_feature_bc_matrix"),
  stringsAsFactors = FALSE
)

# 2.  Libraries & helpers ---------------------------------------------
library(Seurat)      ; packageVersion("Seurat")
library(tidyverse)   ; library(patchwork)
set.seed(1234)

dir.create(OUT_DIR,showWarnings = FALSE,recursive = TRUE)
fig_dir <- file.path(OUT_DIR,"figures") ; dir.create(fig_dir,showWarnings=FALSE)

qc_plot <- function(obj,prefix){
  png(file.path(fig_dir,paste0(prefix,"_QC_violin.png")),800,500)
  print(VlnPlot(obj,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3))
  dev.off()
  png(file.path(fig_dir,paste0(prefix,"_QC_scatter.png")),800,500)
  p1<-FeatureScatter(obj,feature1="nCount_RNA",feature2="percent.mt")
  p2<-FeatureScatter(obj,feature1="nCount_RNA",feature2="nFeature_RNA")
  print(p1+p2)
  dev.off()
}

umap_plot <- function(obj,prefix){
  png(file.path(fig_dir,paste0(prefix,"_UMAP.png")),800,500)
  print(DimPlot(obj,reduction="umap"))
  dev.off()
}

# 3.  Per‑sample preprocessing ----------------------------------------
#     (Creates <sample_id>_SCT.rds in OUT_DIR)
process_sample <- function(row){
  message("\n▶ Processing ",row$sample_id)
  counts <- Read10X(row$path)
  seu <- CreateSeuratObject(counts=counts,project=row$sample_id,
                           min.cells = 3,min.features = 200)
  seu$condition <- row$condition
  seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern="^mt-")
  qc_plot(seu,row$sample_id)_bfr
  seu <- subset(seu,subset = nFeature_RNA>200 & nFeature_RNA<5000 & percent.mt<30)
  qc_plot(seu,row$sample_id)_aft
  seu <- SCTransform(seu,vst.flavor="v2",verbose=FALSE) %>%
          RunPCA(npcs=30,verbose=FALSE) %>%
          RunUMAP(reduction="pca",dims=1:30,verbose=FALSE)
  umap_plot(seu,row$sample_id)
  saveRDS(seu,file.path(OUT_DIR,paste0(row$sample_id,"_SCT.rds")))
  invisible(NULL)
}
apply(SAMPLE_TBL,1,process_sample)


