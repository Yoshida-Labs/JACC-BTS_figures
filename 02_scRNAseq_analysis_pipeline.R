#===============================================================================
# Single‐cell RNA‐seq Analysis Pipeline for Publication
# QC → Doublet removal → SCT integration → Clustering → Figures 7/8 & S13–S17
#===============================================================================

#--- 1. Load required packages -------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(DoubletFinder)
library(lisi)
library(CellChat)
library(clusterProfiler)
library(org.Mm.eg.db)
options(future.globals.maxSize = 50 * 1024^3)

#--- 2. Paths and file names ---------------------------------------------------
base.dir <- "your_path_to_directory"
out.dir  <- file.path(base.dir, "results")
dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)

#--- 3. Load individual SCT‐normalized Seurat objects ---------------------------
heart_c   <- readRDS(file.path(base.dir, "heart_scRNA_SCT_c.rds"))
heart_ns  <- readRDS(file.path(base.dir, "heart_scRNA_SCT_ns.rds"))
heart_azd <- readRDS(file.path(base.dir, "heart_scRNA_SCT_azd.rds"))

# assign metadata
heart_c$condition   <- "Control"
heart_ns$condition  <- "Pump+NS"
heart_azd$condition <- "Pump+AZDtreated"
heart_c$sample_id   <- "C"
heart_ns$sample_id  <- "NS"
heart_azd$sample_id <- "AZD"

#--- 4. Merge and SCTransform integration --------------------------------------
merged <- merge(heart_c, y=list(heart_ns, heart_azd), project="Heart")
merged <- SCTransform(merged, verbose=FALSE)

# initial PCA/UMAP for QC visualization
merged <- RunPCA(merged, verbose=FALSE)
merged <- RunUMAP(merged, dims=1:30, verbose=FALSE)

#--- 5. Doublet & low‐quality filtering ----------------------------------------
# mitochondrial %
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern="^MT-")

# compute LISI for batch mixing
lisi_scores <- compute_lisi(embeddings=merged@reductions$umap@cell.embeddings,
                            meta_data=data.frame(condition=merged$condition),
                            label_colnames="condition")
merged$LISI <- lisi_scores$condition
merged$gene_lisi_ratio <- merged$nFeature_RNA / merged$LISI

# DoubletFinder
annotations <- merged@meta.data$SCT_snn_res.0.6
hom.prop <- modelHomotypic(annotations)
nExp     <- round(0.075 * ncol(merged))
nExp.adj <- round(nExp * (1 - hom.prop))
merged   <- doubletFinder(merged, PCs=1:30, pN=0.10, pK=0.06,
                          nExp=nExp.adj, sct=TRUE)

# filter: singlets + QC thresholds
merged <- subset(merged,
                 subset = DF.classifications_0.1_0.06_.* == "Singlet" &
                          nFeature_RNA > 200 & nFeature_RNA < 5000 &
                          percent.mt < 30)

#--- 6. Re‐run PCA/UMAP/clustering on filtered data -----------------------------
merged <- RunPCA(merged, verbose=FALSE)
merged <- RunUMAP(merged, dims=1:30, verbose=FALSE)
merged <- FindNeighbors(merged, dims=1:30, verbose=FALSE)
merged <- FindClusters(merged, resolution=0.6, verbose=FALSE)

# save cleaned object
saveRDS(merged, file.path(out.dir, "integrated_data_filtered.rds"))

#--- 7. Generate Figure 7 ------------------------------------------------------
# 7A: UMAP colored by celltype & condition
merged_heart <- readRDS("integrated_data_with_QC_res0.6.rds")
new_cluster_labels <- c(
  "EC1","Macro1","FB1","Peri","CM1","CM-EC","FB2","Macro2",
  "lowQC-CM1","EC2","lowQC-CM2","EC3","CM2","EC-FB","EC4",
  "LC1","EC5","SM","LC2","EC6","LC3","CM3","LC4","CM-FB",
  "Macro3","EC7","CM-Peri"
)
names(new_cluster_labels) <- as.character(0:26)
merged_heart <- RenameIdents(merged_heart, new_cluster_labels)

# Subset out the two low-QC clusters
merged_no_lowQC <- subset(
  merged_heart,
  idents = c("lowQC-CM1", "lowQC-CM2"),
  invert = TRUE
)

# Plot UMAP on the filtered object
umap_no_lowQC <- DimPlot(
  merged_no_lowQC,
  reduction = "umap",
  label = TRUE,
    repel = TRUE,
  label.size = 12
) + theme(legend.position = "right")

# Save to PNG
png("Fig7A_UMAP.png", width = 1000, height = 800)
print(umap_no_lowQC)
dev.off()



# 2. Figure 7B–D: Fibroblast subset reclustering & plots
load("integrated_data_pump+AZD_FBsubset.Rdata")  # gives subset_fibroblast.data

# 7B: UMAP of fibroblast subclusters
umap_fb <- DimPlot(subset_fibroblast.data, reduction="umap", label=TRUE, label.size=10) +
  theme(
    plot.title=element_text(size=24),
    axis.title=element_text(size=24),
    axis.text=element_text(size=24),
    legend.title=element_text(size=24),
    legend.text=element_text(size=24)
  )
png("Fig7B_UMAP_FB.png", width=500, height=500)
print(umap_fb)
dev.off()

# 7C: 100% stacked bar of subcluster distribution by condition
df_fb <- subset_fibroblast.data@meta.data %>%
  count(condition, seurat_clusters) %>%
  group_by(condition) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  mutate(condition = factor(condition, levels=c("Control","Pump+NS","Pump+AZDtreated")))
bar_fb <- ggplot(df_fb, aes(condition, freq, fill=factor(seurat_clusters))) +
  geom_bar(stat="identity", position="fill") +
  scale_y_continuous(labels=percent) +
  labs(x="Condition", y="Percentage", fill="Clusters") +
  theme_minimal() +
  theme(
    axis.title=element_text(size=18),
    axis.text=element_text(size=18),
    legend.title=element_text(size=18),
    legend.text=element_text(size=18)
  )
png("Fig7C_Bar_FB.png", width=500, height=600)
print(bar_fb)
dev.off()

# 7D: Violin + boxplots for fibrosis genes
Idents(subset_fibroblast.data) <- subset_fibroblast.data$condition
genes_fb <- c("Col1a1","Col1a2","Mmp2")
for(g in genes_fb){
  dat <- FetchData(subset_fibroblast.data, vars=c("condition", g))
  p_num <- ggviolin(dat, x="condition", y=g, color="condition", add="boxplot",
                    add.params=list(fill="white",width=0.05)) +
    stat_summary(fun.y=mean, geom="point", color="red") +
    ggtitle(paste0(g,": numeric p-values")) +
    theme(axis.title=element_text(size=18), axis.text=element_text(size=16),
          strip.text=element_text(size=16))
  p_star <- ggviolin(dat, x="condition", y=g, color="condition", add="boxplot",
                     add.params=list(fill="white",width=0.05)) +
    stat_compare_means(comparisons=list(c("Control","Pump+NS"),c("Pump+NS","Pump+AZDtreated"),c("Control","Pump+AZDtreated")),
                       method="t.test", label="p.signif", size=6) +
    ggtitle(paste0(g,": significance")) +
    theme(axis.title=element_text(size=18), axis.text=element_text(size=16),
          strip.text=element_text(size=16))
  ggsave(paste0("Fig7D_",g,"_numeric.png"), plot=p_num, width=6, height=10, units="in", dpi=300)
  ggsave(paste0("Fig7D_",g,"_asterisk.png"), plot=p_star, width=6, height=10, units="in", dpi=300)
}


# 7E–G: Cardiomyocyte subset
cm <- subset(merged, idents=grep("^CM", levels(Idents(merged)), value=TRUE))
cm <- RunUMAP(cm, dims=1:30)
cm <- FindClusters(cm, resolution=0.7)
pE <- DimPlot(cm, reduction="umap", label=TRUE)
dfE <- cm@meta.data %>% count(condition, seurat_clusters) %>%
       group_by(condition) %>% mutate(frac = n/sum(n))
pF <- ggplot(dfE, aes(condition, frac, fill=seurat_clusters)) +
      geom_bar(stat="identity", position="fill") + ylab("Proportion")

# violin of Mmp2 & Fgfr1 across major cell types
topcell <- subset(merged, idents=c("FB","EC","Macro","CM"))
pG1 <- VlnPlot(topcell, features="Mmp2", group.by="condition", split.by="ident", pt.size=0)
pG2 <- VlnPlot(topcell, features="Fgfr1", group.by="condition", split.by="ident", pt.size=0)
ggsave("Fig7E-G_CMs.pdf", pE + pF + (pG1|pG2), path=out.dir, width=15, height=5)

#--- 8. Generate Figure 8 ------------------------------------------------------
# prepare CellChat objects for NS vs AZD
cc_ns <- readRDS("cellchat_heart_scRNA_SCT_Pump+NS.rds")
cc_azd<- readRDS("cellchat_heart_scRNA_SCT_Pump+AZDtreated.rds")
cc_list <- list(NS=cc_ns, AZD=cc_azd)
cellchat <- mergeCellChat(cc_list, add.names=names(cc_list))

# 8A: network diagram
p1 <- compareInteractions(cellchat, show.legend=FALSE, group=c(1,2)) +
  ggtitle("Number of Interactions") +
  theme(plot.title=element_text(size=20,hjust=0.5), axis.title=element_text(size=18), axis.text=element_text(size=16))

p2 <- compareInteractions(cellchat, measure="weight", show.legend=FALSE, group=c(1,2)) +
  ggtitle("Interaction Strength") +
  theme(plot.title=element_text(size=20,hjust=0.5), axis.title=element_text(size=18), axis.text=element_text(size=16))
png("Fig8A_Interactions.png", width=800, height=400)
print(p1 + p2)
dev.off()

# 8B: top10 FGF→CM and NPR1→FB
pFGF <- rankNet(cellchat, mode="comparison", measure="weight",
                sources.use="Ventricular CM", targets.use=NULL)
pNPR <- rankNet(cellchat, mode="comparison", measure="weight",
                sources.use="Fibroblast", targets.use=NULL)
ggsave("Fig8B_FGF_NPR1.pdf", pFGF+pNPR, path=out.dir, width=10, height=5)

# 8C–D: chord diagrams for FGF & NPR1
for(sig in c("FGF","NPR1")){
  wmax <- getMaxWeight(cellchat, slot.name="netP", attribute=sig)
  png(file.path(out.dir, paste0("Fig8C-D_",sig,".png")), width=800, height=400)
  par(mfrow=c(1,2), xpd=TRUE)
  for(i in 1:2){
    netVisual_aggregate(object.list(cellchat)[[i]], signaling=sig,
                        layout="circle", edge.weight.max=wmax,
                        signaling.name=paste(sig,names(object.list(cellchat))[i]))
  }
  dev.off()
}

# 8E: violin of Nppa,Nppb,Npr1 by celltype/condition
genes8E <- c("Nppa","Nppb","Npr1")
p8E <- VlnPlot(merged, features=genes8E, group.by="celltype", split.by="condition",
               pt.size=0) & theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("Fig8E_Nppa_Nppb_Npr1.pdf", p8E, path=out.dir, width=8, height=6)

#--- 9. Supplementary S13–S17 -------------------------------------------------
# S13: FeaturePlots in FB & CM subsets
genes_S13_FB <- c("Dcn","Col1a1","Postn","Mmp2")
genes_S13_CM <- c("Myh6","Tnnt2","Mybpc3","Ryr2")
png(file.path(out.dir,"FigS13_FB.png"), width=900, height=900)
print(FeaturePlot(fib, genes_S13_FB, ncol=2))
dev.off()
png(file.path(out.dir,"FigS13_CM.png"), width=900, height=900)
print(FeaturePlot(cm, genes_S13_CM, ncol=2))
dev.off()

# S14 & S15: GO BP enrichment bubble plots for each cluster
go_bubble_plot <- function(markers, title){
  ego <- compareCluster(geneCluster=markers, fun="enrichGO",
                        OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH")
  dotplot(ego, showCategory=10) + ggtitle(title)
}
# prepare list of upregulated genes per subcluster
markers_FB <- lapply(levels(fib), function(cl){
  Idents(fib)<-cl
  rownames(FindMarkers(fib, ident.1=cl, logfc.threshold=0, only.pos=TRUE))
})
names(markers_FB) <- levels(fib)
pS14 <- go_bubble_plot(markers_FB, "FB subclusters GO BP")
ggsave("FigS14_GO_FB.pdf", pS14, path=out.dir, width=12, height=8)

markers_CM <- lapply(levels(cm), function(cl){
  Idents(cm)<-cl
  rownames(FindMarkers(cm, ident.1=cl, logfc.threshold=0, only.pos=TRUE))
})
names(markers_CM) <- levels(cm)
pS15 <- go_bubble_plot(markers_CM, "CM subclusters GO BP")
ggsave("FigS15_GO_CM.pdf", pS15, path=out.dir, width=12, height=10)

# S16: Comparative bubble plot of information flow
png(file.path(out.dir,"FigS16_bubble.png"), width=1200, height=800)
netVisual_bubble(cellchat, sources.use=NULL, targets.use=NULL,
                 comparison=c(1,2), angle.x=45)
dev.off()

# S17: violin of 9 fibrosis/signaling genes across 7 cell types
genes_S17 <- c("Mmp2","Spry1","Sorbs3","Smad3","Rhoa",
               "Hrh2","Vim","Fgfr1","Ctnnb1")
pS17 <- VlnPlot(merged, features=genes_S17, group.by="celltype",
                split.by="condition", pt.size=0) &
        theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave("FigS17_violins.pdf", pS17, path=out.dir, width=12, height=8)

