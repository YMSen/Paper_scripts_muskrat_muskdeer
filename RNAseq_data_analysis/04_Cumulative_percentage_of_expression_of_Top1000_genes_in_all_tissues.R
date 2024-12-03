library(ggplot2)
library(reshape2)
library(RColorBrewer)

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


data = as.data.frame(matrix(ncol = 3, nrow = nrow(info)))
colnames(data) = c("Sample","Tissue","Cumulative_proportion")
data[["Sample"]] = info[["sample"]]
data[["Tissue"]] = info[["condition"]]


for (order in c(1:nrow(info))) {
  data_1 = data_exp[ , info[order,"sample"]]
  data_1 = data_1[order(data_1, decreasing = T)]
  data[["Cumulative_proportion"]][order] = sum(data_1[1:1000])/sum(data_1)*100
}


data_plot_1 = data

mean <- aggregate(data_plot_1[["Cumulative_proportion"]], by=list(data_plot_1[["Tissue"]]), FUN=mean)
colnames(mean) = c("Tissue","mean")
sd <- aggregate(data_plot_1[["Cumulative_proportion"]], by=list(data_plot_1[["Tissue"]]), FUN=sd)
colnames(sd) = c("Tissue","sd")

data_plot_2 = cbind(mean,sd)[ ,c(1,2,4)]

data_plot_2 = data_plot_2[order(data_plot_2[["mean"]], decreasing = T), ]

data_plot_2[["Tissue"]] = factor(data_plot_2[["Tissue"]], levels = data_plot_2[["Tissue"]], ordered = T)
data_plot_1[["Tissue"]] = factor(data_plot_1[["Tissue"]], levels = data_plot_2[["Tissue"]], ordered = T)


P = ggplot()+
  geom_bar(data = data_plot_2, aes(x=Tissue, y=mean, colour=Tissue), stat = "identity", fill="white", width=0.4, position=position_dodge(.6))+
  geom_errorbar(data = data_plot_2, aes(x=Tissue, ymin=mean-sd, ymax=mean+sd, colour=Tissue), position=position_dodge(.6), width=.17)+
  geom_jitter(data = data_plot_1, aes(x=Tissue, y=Cumulative_proportion, fill=Tissue), size = 1.2, shape = 21,
              stroke = 0.01, show.legend = FALSE, 
              position = position_jitter(width = 0.2, height = 0.0, seed = 123))+
  scale_color_manual(values=c(brain = '#f1c933', eyeball = '#a9eee6', fat = '#626f92', gland = "#f5c7f7", heart = "#a393eb",
                              kidney = "#ffaa64", liver = '#7098da', lung = '#3fbac2', muscle = '#de95ba', ovary = "#cd595a",
                              spleen = "#ff9900", testis = "#6a7efc", uterus = "#7cbd1e"))+
  scale_fill_manual(values=c(brain = '#f1c933', eyeball = '#a9eee6', fat = '#626f92', gland = "#f5c7f7", heart = "#a393eb",
                             kidney = "#ffaa64", liver = '#7098da', lung = '#3fbac2', muscle = '#de95ba', ovary = "#cd595a",
                             spleen = "#ff9900", testis = "#6a7efc", uterus = "#7cbd1e"))+
  theme_test()+
  scale_y_continuous(expand = expansion(add = c(0, 25)))+
  ylab("Cumulative proportion (%)")+
  xlab("")+
  labs(title = "")+
  theme(axis.text.x=element_text(angle=35, hjust = 1, colour="black", family="Times", size=10),
        panel.grid.minor = element_blank())

ggsave(P, filename = "Cumulative proportion of Top1000 gene expression in all tissue.pdf", width = 6, height = 3.5)


rm(list = ls())
dev.off()
