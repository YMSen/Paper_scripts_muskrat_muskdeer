
library(Seurat)
library(dplyr)
library(patchwork)
library(ggsci)
library(ggplot2)


dir.create('./03_gene_expression_analysis')


merged_object_newer = readRDS('./01_rds_file/integrated_seurat.rds')


merged_object_newer@meta.data$group[merged_object_newer@meta.data$Sample == 'OZ_3_B'] = 'S'
merged_object_newer@meta.data$group[merged_object_newer@meta.data$Sample == 'OZ-1-NFP-A'] = 'NS'


merged_object_newer$celltype.group = paste0(merged_object_newer$seurat_clusters, '_', merged_object_newer$group)
level_1 = paste0(rep(levels(merged_object_newer$seurat_clusters), each = 2), c('_S', '_NS'))
merged_object_newer$celltype.group = factor(merged_object_newer$celltype.group, levels = level_1)


merged_object_newer$celltype = Idents(merged_object_newer)
Idents(merged_object_newer) = "celltype.group"

genes = c('Ozibchr11G0599','Ozibchr2G1080')

P = DotPlot(object = merged_object_newer,
        assay = "RNA",
        features = genes,
        group.by = "celltype.group"
        )+
  theme_classic()+
  coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))+
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = rev(c("#D93F49","#E28187","#EBBFC2","#D5E1E3","#AFC9CF","#8FB4BE")))+
  labs(x = NULL, y = NULL)+
  guides(size = guide_legend(order=3))

save_pdf = paste0('./03_gene_expression_analysis/', 'gene_expression', '.pdf')
ggsave(P, file = save_pdf, width = 10, height = 2)



rm(list=ls())

quit()