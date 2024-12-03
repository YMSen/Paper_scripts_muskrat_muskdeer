library(ggplot2)

data_exp = read.table(file = "./00_files/All_samples_tpm_1.0_expression_filter.txt", header = T, check.names = F, row.names = 1)
data_exp = data_exp[ , c("OZ_1_brain",'OZ_2_brain',"OZ_3_brain","OZ_NFP_1_brain",'OZ_NFP_2_brain',"OZ_NFP_3_brain","OZ_NFP_4_brain",'OZ_NFP_5_brain',"OZ_NFP_6_brain",
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

info = data.frame(sample = colnames(data_exp), condition = as.character(lapply(strsplit(colnames(data_exp), "_"), tail, 1)))


data_average_exp = as.data.frame(matrix(ncol = 13, nrow = nrow(data_exp)))
colnames(data_average_exp) = unique(info[["condition"]])
rownames(data_average_exp) = rownames(data_exp)

for (tissue in unique(info[["condition"]])) {
  data_one_tissue = data_exp[ , grep(tissue, colnames(data_exp))]
  data_average_exp[ , tissue] = rowMeans(data_one_tissue)
}

data = as.data.frame(matrix(ncol = 5, nrow = 1))
colnames(data) = c("Order","Proportion_cumulative_expression","Sample","Tissue","Expression")

info_2 = data.frame(sample = colnames(data_average_exp), condition = as.character(lapply(strsplit(colnames(data_average_exp), "_"), tail, 1)))

for (order in c(1:nrow(info_2))) {
  data_1 = data_average_exp[ , info_2[order,"sample"]]
  data_1 = data_1[order(data_1, decreasing = T)]
  data_one_tissue = as.data.frame(matrix(ncol = 5, nrow = length(data_1)))
  colnames(data_one_tissue) = c("Order","Proportion_cumulative_expression","Sample","Tissue","Expression")
  data_one_tissue[["Expression"]] = data_1
  data_one_tissue[["Order"]] = c(1:nrow(data_one_tissue))
  data_one_tissue[["Sample"]] = info_2[order,"sample"]
  data_one_tissue[["Tissue"]] = info_2[order,"condition"]
  for (number in c(1:nrow(data_one_tissue))) {
    if (number == 1){
      data_one_tissue[["Proportion_cumulative_expression"]][number] = data_one_tissue[["Expression"]][number]/sum(data_one_tissue[["Expression"]])*100
    }
    if (number != 1){
      data_one_tissue[["Proportion_cumulative_expression"]][number] = data_one_tissue[["Expression"]][number]/sum(data_one_tissue[["Expression"]])*100 + data_one_tissue[["Proportion_cumulative_expression"]][number-1]
    }
  }
  data = rbind(data, data_one_tissue)
  data = na.omit(data)
}

write.table(data, file = "Summary_of_cumulative_percentage_of_average_gene_expression_across_all_tissues.txt", col.names = T, row.names = F, sep = "\t", quote = F)


data_plot = read.table(file = "Summary_of_cumulative_percentage_of_average_gene_expression_across_all_tissues.txt", header = T, check.names = F)
P = ggplot(data_plot, aes(Order, Proportion_cumulative_expression))+ 
  geom_line(aes(colour = Tissue), linewidth = 0.5)+
  scale_color_manual(values=c(brain = '#f1c933', eyeball = '#a9eee6', fat = '#626f92', gland = "#f5c7f7", heart = "#a393eb",
                              kidney = "#ffaa64", liver = '#7098da', lung = '#3fbac2', muscle = '#de95ba', ovary = "#cd595a",
                              spleen = "#ff9900", testis = "#6a7efc", uterus = "#7cbd1e"))+
  theme_test()+
  #geom_hline(yintercept = 50, lty=2, col="grey", lwd=0.5)+ 
  #geom_vline(xintercept = 100, lty=2, col="grey", lwd=0.5)+
  geom_vline(xintercept = 1000, lty=2, col="grey", lwd=0.25)+
  annotate(geom = "text", label = "N=1000", x = 2000, y = 90, size = 4, colour = "black")+
  scale_x_continuous(expand = expansion(add = c(1, 500)))+
  scale_y_continuous(expand = expansion(add = c(0, 10)))+
  #scale_x_continuous(limits = c(1, 15000))+
  #scale_y_continuous(limits = c(0, 100))+
  ylab("Proportion of all transcriptional outputs (%)")+
  xlab("Rank order")+
  labs(title = "")

ggsave(P, filename = "Transcription complexity analysis of all tissue.pdf", width = 6, height = 3.5)

rm(list = ls())
dev.off()