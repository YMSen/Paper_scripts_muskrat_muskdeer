library(ggsignif)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(reshape2
library(RColorBrewer)
library(nycflights13)
library(dplyr)
library(ggforce)
library(patchwork)

data = read.table(file = "./00_files/All_samples_TPM_0.1_expression_filter_miRNA.txt", header = T, check.names = F, row.names = 1)

data_gland_gene = data[c("miRNA_00054","miRNA_00290","miRNA_03298"), ]

data_gland_gene = as.data.frame(t(data_gland_gene))

data_gland_gene[["tissue"]] = as.character(lapply(strsplit(rownames(data_gland_gene), "_"), tail, 1))

data_gland_gene[["sample"]] = rownames(data_gland_gene)

data_plot = gather(data_gland_gene, gene, data, "miRNA_00054":"miRNA_03298")


P = ggplot(data = data_plot, mapping = aes(y = tissue,
                                           x = log2(data+1), 
                                           fill = tissue, 
                                           label = tissue))+
  geom_violin(scale = 'width', adjust = 1, trim = TRUE)+
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = '', times = length(x) - 2), x[length(x) - 1], ''))+
  facet_grid(cols = vars(gene), scales = 'free')+
  cowplot::theme_cowplot(font_family = '')+
  scale_fill_manual(values=c('#f1c933', '#a9eee6', '#626f92', "#f5c7f7", "#a393eb",
                             "#ffaa64", '#7098da', '#3fbac2', '#de95ba', "#cd595a",
                             "#ff9900", "#6a7efc", "#7cbd1e"))+
  xlab('Expression Level')+
  ylab('Tissue')+
  theme(
    panel.spacing = unit(x = 0, units = 'lines'),
    axis.line = element_blank(),
    panel.background = element_rect(fill = NA, color = 'black'),
    strip.background = element_blank(),
    strip.text = element_text(color = 'black', size = 10, angle = 285),
    axis.text.x = element_text(color = 'black', size = 7),
    axis.text.y = element_text(color = 'black', size = 8),
    axis.title.x = element_text(color = 'black', size = 11),
    axis.title.y = element_text(color = 'black', size = 11),
    axis.ticks.x = element_line(color = 'black', lineend = 'round'),
    axis.ticks.y = element_blank(),
    axis.ticks.length = unit(x = 0.1, units = 'cm'))

ggsave(P, filename = "Gland_specific_miRNA_violin_heatmap.pdf", width = 4, height = 6)


rm(list = ls())
dev.off()