library(Seurat)
library(dplyr)
library(patchwork)
library(ggsci)
library(ggplot2)
library(ggrepel)


dir.create('./04_differential_expression_analysis')


merged_object_newer = readRDS('./01_rds_file/integrated_seurat.rds')


merged_object_newer@meta.data$group[merged_object_newer@meta.data$Sample == 'OZ_3_B'] = 'S'
merged_object_newer@meta.data$group[merged_object_newer@meta.data$Sample == 'OZ-1-NFP-A'] = 'NS'


merged_object_newer$celltype.group = paste0(merged_object_newer$seurat_clusters, '_', merged_object_newer$group)
level_1 = paste0(rep(levels(merged_object_newer$seurat_clusters), each = 2), c('_S', '_NS'))
merged_object_newer$celltype.group = factor(merged_object_newer$celltype.group, levels = level_1)


merged_object_newer$celltype = Idents(merged_object_newer)
Idents(merged_object_newer) = "celltype.group"


cell_DEGs = FindMarkers(object = merged_object_newer,
                      ident.1 = '14_S',
                      ident.2 = '14_NS',
                      test.use = "wilcox",
                      logfc.threshold = 0,
                      min.pct = 0.01,
                      verbose = FALSE)



save_csv = paste0("./04_differential_expression_analysis/", "cluster_14_DEGs", ".csv")
write.csv(cell_DEGs, file = save_csv)



markers = cell_DEGs %>% mutate(Difference = pct.1 - pct.2)

markers[['genes']] = ifelse(markers[['p_val_adj']] <0.05 & abs(markers[['avg_log2FC']]) >0.25, rownames(markers), NA)

markers[['significance']] = ifelse(markers[['p_val_adj']] <0.05 & abs(markers[['avg_log2FC']]) >0.25, ifelse(markers[['avg_log2FC']] >0.25, 'up', 'down'), 'ns')


P = ggplot()+
  geom_point(data = markers,
             aes(x=avg_log2FC,
                 y=-log10(p_val_adj),
                 # color=significance,
                 color=-log10(p_val_adj),
                 ))+ 
  scale_color_gradientn(values = seq(0, 1, 0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # scale_color_manual(values = c(down = "#1C71B6", ns = "grey60", up = "#EE2A29"))+
  geom_text_repel(data = markers,
                  aes(x=avg_log2FC, y=-log10(p_val_adj),label = genes),
                  size = 1.5,
                  force = 0.5)+
  geom_vline(xintercept = -0.25, linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "#999999")+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  theme_test()


save_pdf = paste0("./04_differential_expression_analysis/", "Volcanic_plot", ".pdf")
ggsave(P, file = save_pdf, width = 4, height = 5)


rm(list=ls())

dev.off()