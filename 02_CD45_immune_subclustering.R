# ==============================================================================
# CD45-positive immune-cell subclustering
# ==============================================================================

setwd("./2sample")
library(Seurat)

load("scdata_anno_res1.5.RData")
table(scdata$celltype)

immune <- subset(scdata,celltype %in% c("Macrophage","CD4 Tcell","CD8 Tcell","Other Tcell","Granulocytes","NK","DC","B cell"))

CD45pos <- subset(immune, subset = Ptprc > 0)
DimPlot(CD45pos, reduction = "umap",group.by = "celltype",  pt.size = .3)

CD45pos <- NormalizeData(CD45pos)
CD45pos <- FindVariableFeatures(CD45pos)
CD45pos <- ScaleData(CD45pos)
CD45pos <- RunPCA(CD45pos)
CD45pos <- FindNeighbors(CD45pos, reduction = "harmony", dims = 1:30)
CD45pos <- FindClusters(CD45pos, resolution = 0.6)
CD45pos <- RunUMAP(CD45pos, reduction = "harmony", dims = 1:30)
CD45pos

DimPlot(object = CD45pos, reduction = "umap",label = T,pt.size = 0.01) #+ NoLegend()
DimPlot(CD45pos,reduction = "umap", group.by = "celltype",label = T,pt.size = 0.001,label.size=3)
DimPlot(CD45pos,reduction = "umap", group.by = "orig.ident",label = T,pt.size = 0.001)

saveRDS(CD45pos,"CD45pos.RDS")

# FeaturePlot(CD45pos,"Cd8a")
# FeaturePlot(CD45pos,"Cd4")

CD45pos <- readRDS("CD45pos.RDS")
DimPlot(CD45pos,label = T)


table(CD45pos$RNA_snn_res.0.6)
CD45pos$immune_cluster <- CD45pos$RNA_snn_res.0.6
Idents(CD45pos) <- CD45pos$immune_cluster
DimPlot(CD45pos,label = T)
DimPlot(CD45pos,label = T,group.by = "celltype",pt.size = 0.01)


# Immune-cluster metadata and abundance summaries
celltype_counts <- table(CD45pos$orig.ident, CD45pos$immune_cluster)
celltype_counts
celltype_counts <- table(CD45pos$group, CD45pos$immune_cluster)
celltype_counts
table(CD45pos$group)


features1 <- c("Cd3d","Cd4", "Cd8a",# t cells
               "Cd79a","Ms4a1","Mzb1",# b_cells
               "Ncr1","Eomes", #NK
               "S100a8", "S100a9", # Granulocytes
               "Apoe","Cd68", "Csf1r",  # tam-like
               "Clec10a", "Clec9a"  # dc
)
CD45pos$celltype <- factor(CD45pos$celltype,levels = c("DC","Macrophage", "Granulocytes","NK", "B cell","Other Tcell","CD8 Tcell","CD4 Tcell"))

Idents(CD45pos) <- CD45pos$celltype
DotPlot(CD45pos, features = features1, dot.scale = 5) + RotatedAxis() + theme(axis.text.x = element_text(size = 12,angle = 90))+ theme(axis.text.y = element_text(size = 12))


