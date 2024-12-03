library(dplyr)
library(pheatmap)
library(ggplot2)

file_name_1 = paste0("./00_files/List_of_tissue-specific_genes.txt")
table_1 = read.table(file_name_1, header = T)

file_name_2 = paste0("./00_files/All_samples_tpm_1.0_expression_filter.txt")
table_2 = read.table(file_name_2, header = T, row.names = 1)


table_3 = table_2[table_1[["id"]], ]

table_3 = table_3[ , c("OZ_1_brain",'OZ_2_brain',"OZ_3_brain","OZ_NFP_1_brain",'OZ_NFP_2_brain',"OZ_NFP_3_brain","OZ_NFP_4_brain",'OZ_NFP_5_brain',"OZ_NFP_6_brain",
                       "OZ_NFP_1_eyeball","OZ_NFP_2_eyeball","OZ_NFP_3_eyeball","OZ_NFP_4_eyeball","OZ_NFP_5_eyeball","OZ_NFP_6_eyeball",
                       "OZ_NFP_1_fat","OZ_NFP_2_fat","OZ_NFP_3_fat","OZ_NFP_4_fat","OZ_NFP_5_fat","OZ_NFP_6_fat",
                       "OZ_NFP_1_gland","OZ_NFP_2_gland","OZ_NFP_3_gland","OZ_NFP_4_gland","OZ_NFP_5_gland","OZ_NFP_6_gland",
                       "OZ_NFP_1_heart","OZ_NFP_2_heart",'OZ_NFP_3_heart',"OZ_NFP_4_heart","OZ_NFP_5_heart","OZ_NFP_6_heart",
                       "OZ_NFP_1_kidney","OZ_NFP_2_kidney","OZ_NFP_3_kidney","OZ_NFP_4_kidney","OZ_NFP_5_kidney","OZ_NFP_6_kidney",
                       "OZ_NFP_1_liver","OZ_NFP_2_liver","OZ_NFP_3_liver","OZ_NFP_4_liver","OZ_NFP_5_liver","OZ_NFP_6_liver",
                       "OZ_NFP_1_lung","OZ_NFP_2_lung","OZ_NFP_3_lung","OZ_NFP_4_lung","OZ_NFP_5_lung","OZ_NFP_6_lung",
                       "OZ_NFP_1_muscle","OZ_NFP_2_muscle","OZ_NFP_3_muscle","OZ_NFP_4_muscle",'OZ_NFP_5_muscle',"OZ_NFP_6_muscle",
                       "OZ_NFP_1_ovary","OZ_NFP_2_ovary","OZ_NFP_3_ovary","OZ_NFP_4_ovary","OZ_NFP_5_ovary","OZ_NFP_6_ovary",
                       "OZ_NFP_1_spleen","OZ_NFP_2_spleen","OZ_NFP_3_spleen","OZ_NFP_4_spleen",'OZ_NFP_5_spleen','OZ_NFP_6_spleen',
                       "OZ_1_testis","OZ_2_testis","OZ_3_testis","OZ_NFP_1_testis","OZ_NFP_2_testis",'OZ_NFP_3_testis',"OZ_NFP_4_testis","OZ_NFP_5_testis","OZ_NFP_6_testis",
                       "OZ_NFP_1_uterus","OZ_NFP_2_uterus","OZ_NFP_3_uterus","OZ_NFP_4_uterus",'OZ_NFP_5_uterus',"OZ_NFP_6_uterus"
)]

info = data.frame(sample = colnames(table_3), condition = as.character(lapply(strsplit(colnames(table_3), "_"), tail, 1)))

table_3 = t(table_3)

table_zscore = as.data.frame(table_3)
for (gene_id in colnames(table_zscore)) {
  table_zscore[[gene_id]] = scale(as.numeric(table_zscore[[gene_id]]), center=T, scale=T)  # 数据标准化
}


data = table_zscore

annotation_col = data.frame(Group = factor(info[["condition"]]))
annotation_row = data.frame(Group = factor(table_1[["Tissue"]]))

rownames(annotation_col) = info[["sample"]]
rownames(annotation_row) = colnames(data)

ann_colors = list(Group = c(brain = '#f1c933', eyeball = '#a9eee6', fat = '#626f92', gland = "#f5c7f7", heart = "#a393eb",
                            kidney = "#ffaa64", liver = '#7098da', lung = '#3fbac2', muscle = '#de95ba', ovary = "#cd595a",
                            spleen = "#ff9900", testis = "#6a7efc", uterus = "#7cbd1e"))

pdf("Tissue-Specific_Gene_Z-Score_Heatmap.pdf", height = 16.5, width = 5.5)

pheatmap(t(data),
         #clustering_distance_rows=sampleDists_2,
         #clustering_distance_cols=sampleDists_2,
         breaks=seq(-1,1,0.05),
         color = colorRampPalette(c("#788cb6", "white", "#c68143"))(40),
         cluster_cols = F,
         cluster_rows = F,
         border_color= NA,
         show_rownames = F,
         show_colnames = F,
         main = "",
         display_numbers=F,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors[1],
         angle_col = 90,
         fontsize_row = 12,
         fontsize_col = 12
)+theme_test()

dev.off()

rm(list=ls())