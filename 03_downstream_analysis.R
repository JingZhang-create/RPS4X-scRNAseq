# ==============================================================================
# Downstream_analysis
# ==============================================================================

library(Seurat)
library(ggplot2)
library(pheatmap)

#Differential expression in epithelial cells 
load("scdata_anno_res1.5.RData")
epi <- subset(scdata,celltype == "Epithelial cell")
DotPlot(epi,c("Rps4x","Stat1","Irf1","H2-K1","H2-D1","Tap1","Tap2","Psmb9","Tapbp","B2m"),group.by = "group") + RotatedAxis()



#Average-expression heatmap for immune-cell clusters
CD45pos <- readRDS("CD45pos.RDS")
features1 <- c("Cd3e","Cd4","Cd8a","Pdcd1", "Havcr2","Lag3","Tigit","Ctla4","Foxp3","Tox",
               "Cd40lg","Cd69","Icos","Prf1","Gzmb","Nkg7","Gzmk","Ifng","Mki67","Top2a","Il7r","Ccr7","Sell","Il2ra",
               "Apoe","Lyz2","C1qb","Nos2","Cd80","Cd86","Mrc1","Arg1","Cd163","Fcgr3","Csf3r","S100a8",
               "S100a9","Mmp9","Cxcr4","Prok2","Itgam","Ly6c1","Itga2","Eomes","Il2rb","Ncr1","Cst3","Clec9a","Clec10a","Cd19",
               "Ms4a1","Fcer2a","Slamf6","Tbx21","Tcf7")

avg_exp <- AverageExpression(
  CD45pos,
  assays = "RNA", 
  features = features1,
  group.by = "RNA_snn_res.0.6" 
)$RNA 


colnames(avg_exp)
avg_exp <- avg_exp[,c(4,5,9,1,12,7,3,6,2,8,10,13,11)]
avg_exp <- as.data.frame(avg_exp)
avg_exp_ordered <- avg_exp[intersect(features1, rownames(avg_exp)), ]
colnames(avg_exp_ordered) <- c("3","4","8","0","11","6","2","5","1","7","9","12","10")


pdf("plot.pdf",width = 4,height = 10)
pheatmap(
  avg_exp_ordered,
  scale = "row",        
  cluster_rows = F,   
  cluster_cols = F,   
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  angle_col = 45,       
  fontsize_row = 7,    
  fontsize_col = 10     
)
dev.off()

