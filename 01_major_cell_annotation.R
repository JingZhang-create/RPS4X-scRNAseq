# ==============================================================================
# Major_cell_annotation
# ==============================================================================

setwd("./2sample")
library(Seurat)
library(dplyr)
library(ggplot2)

# Major cell-type annotation
load("scdata_preanno.RData")

features1 <- c("Krt8","Krt18", "Krt19",  # Epithelial cells
               "Pecam1", "Cldn5", "Vwf", # Endothelial cells
               "Col1a1","Col1a2", "Dcn", # Fibroblasts
               "Cd79a","Ms4a1","Mzb1",   # B cells
               "Cd3d","Cd4", "Cd8a",     # T cells
               "Ncr1","Eomes","Il2rb",   # NK cells
               "S100a8", "S100a9",       # Granulocytes
               "Apoe","Cd68", "Csf1r", "Cd86","Mrc1","Arg1", # Macrophage
               "Clec10a", "Clec9a"  # Dendritic cells
)


library("ggplot2")
DotPlot(scdata, features = features1, dot.scale = 3) + 
  RotatedAxis() + 
  theme(axis.text.x = element_text(size = 8))+ theme(axis.text.y = element_text(size = 8))

new.cluster.ids <- c(
  "Epithelial cell",   # Cluster 0
  "Epithelial cell",   # Cluster 1
  "Epithelial cell",   # Cluster 2
  "Epithelial cell",   # Cluster 3
  "Granulocytes",      # Cluster 4
  "Macrophage",        # Cluster 5
  "Epithelial cell",   # Cluster 6
  "Epithelial cell",   # Cluster 7
  "Epithelial cell",   # Cluster 8
  "CD4 Tcell",         # Cluster 9
  "Macrophage",        # Cluster 10
  "Epithelial cell",   # Cluster 11
  "Epithelial cell",   # Cluster 12
  "CD8 Tcell",         # Cluster 13
  "Granulocytes",      # Cluster 14
  "Other Tcell",       # Cluster 15
  "CD8 Tcell",         # Cluster 16
  "CD8 Tcell",         # Cluster 17
  "CAF",               # Cluster 18
  "Epithelial cell",   # Cluster 19
  "DC",                # Cluster 20
  "NK",                # Cluster 21
  "B cell",            # Cluster 22
  "Endothelial cell"   # Cluster 23
)

names(new.cluster.ids) <- levels(scdata)
scdata<- RenameIdents(scdata, new.cluster.ids)
DimPlot(scdata, reduction = "umap", label = F, label.size=3,pt.size = .05)+ theme(
  legend.text = element_text(size = 10))


# Cell-type abundance
scdata$celltype <- scdata@active.ident
scdata$celltype <- factor(scdata$celltype,levels = c("Epithelial cell","Endothelial cell","CAF","CD4 Tcell","CD8 Tcell","Other Tcell", "B cell","Granulocytes","Macrophage","NK","DC"))
scdata@active.ident <- scdata$celltype

celltype_counts <- table(scdata$group, scdata$celltype)
celltype_counts
celltype_proportions <- prop.table(celltype_counts, margin = 1)
celltype_proportions <- as.data.frame(celltype_proportions)
colnames(celltype_proportions) <- c("group", "celltype", "proportion")
str(celltype_proportions)
celltype_proportions$group <- factor(celltype_proportions$group,levels = c("NC","KD"))
ggplot(celltype_proportions, aes(x = group, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Proportion") +theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size = 10,colour = 'black'))

save(scdata,file="scdata_anno_res1.5.RData")

scdata$celltype <- scdata@active.ident
scdata$celltype <- as.character(Idents(scdata))
scdata$celltype[
  scdata$celltype %in% c("Other Tcell", "CD8 Tcell", "CD4 Tcell")
] <- "T cell"
table(scdata$celltype)
scdata$celltype <- factor(scdata$celltype,levels = c("DC","Macrophage", "Granulocytes","NK", "B cell","T cell",
                                                     "CAF","Endothelial cell",  "Epithelial cell"))

features1 <- c("Krt8","Krt18","Krt14",  # epithelial cells
               "Pecam1", "Cldn5", "Vwf", # endothelial cells
               "Col1a1", "Col1a2","Dcn", # fibroblasts
               "Cd3d","Cd4", "Cd8a",# t cells
               "Cd79a","Ms4a1","Mzb1",# b_cells
               "Ncr1","Eomes","Il2rb", #NK
               "S100a8", "S100a9", # Granulocytes
               "Apoe","Cd68", "Csf1r",  # tam-like
               "Clec10a", "Clec9a"  # dc
)


scdata@active.ident <- scdata$celltype
table(scdata$celltype)
DotPlot(scdata, features = features1, dot.scale = 5) + RotatedAxis() + theme(axis.text.x = element_text(size = 12,angle = 90))+ theme(axis.text.y = element_text(size = 12))

