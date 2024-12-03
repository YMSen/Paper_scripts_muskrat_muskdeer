
library(Seurat)
library(dplyr)
library(patchwork)
library(ggsci)
library(ggplot2)


dir.create('./02_sample_t-SNE')


merged_object_newer = readRDS('./01_rds_file/integrated_seurat.rds')


table(merged_object_newer@meta.data$Sample)



p1 = DimPlot(merged_object_newer,
             reduction = "tsne",
             group.by = "Sample",
             label = FALSE,
             repel = TRUE,
             raster=FALSE,
             pt.size = 0.0001)+
  scale_color_manual(values = c("OZ_3_B" = "#FF7B6E", "OZ-1-NFP-A" = "#4F82FF"))+
  theme_test()+
  xlab('tSNE_1')+
  ylab('tSNE_2')+
  labs('')


out_name_1 = paste0('./02_sample_t-SNE/', '/diff_sample_merge_tSNE.pdf')
ggsave(out_name_1, plot=p1, width = 4.5, height =3)




p2 = DimPlot(merged_object_newer,
             reduction = "tsne",
             group.by = "Sample",
             split.by = "Sample",
             label = FALSE,
             repel = TRUE,
             raster=FALSE,
             pt.size = 0.0001)+
  scale_color_manual(values = c("OZ_3_B" = "#FF7B6E", "OZ-1-NFP-A" = "#4F82FF"))+
  theme_test()+
  xlab('tSNE_1')+
  ylab('tSNE_2')+
  labs('')


out_name_2 = paste0('./02_sample_t-SNE/', '/diff_sample_tSNE.pdf')
ggsave(out_name_2, plot=p2, width = 6.5, height =3)



p3 = DimPlot(merged_object_newer,
             reduction = "tsne",
             split.by = "Sample",
             label = TRUE,
             repel = TRUE,
             raster=FALSE,
             pt.size = 0.0001)+
  scale_color_manual(values = c('#5E8CBF', '#385487', '#E0C73C', '#E25C60', '#EFC0B9', '#EB999B',
                                '#49A27B', '#64BFCA', '#D15AA1', '#91306C', '#B1A859', '#6DA2D6',
                                '#462C5C', '#F0DE83', '#E1544A', '#E37C7E', '#EFB1AC', '#D0D3DE',
                                '#D08DBD', '#576A8B', '#3AB3C4', '#75C6BB', '#F3E4C8'))+
  theme_test()+
  xlab('tSNE_1')+
  ylab('tSNE_2')+
  labs('')


out_name_3 = paste0('./02_sample_t-SNE/', 'diff_sample_show_clusters_tSNE.pdf')
ggsave(out_name_3, plot=p3, width = 7, height =3)



rm(list=ls())

quit()
