
library(ggsignif)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(nycflights13)
library(dplyr)
library(ggforce)
library(patchwork)

data = read.table(file = "./00_files/All_samples_TPM_1.0_expression_filter_lnc.txt", header = T, check.names = F, row.names = 1)

data_gland_gene = data[c("G26977","G33936","G2915","G20176","G26992","G64566","G58402","G29813","G31515","G68842"), ]

data_gland_gene = as.data.frame(t(data_gland_gene))

data_gland_gene[["tissue"]] = as.character(lapply(strsplit(rownames(data_gland_gene), "_"), tail, 1))

data_gland_gene[["sample"]] = rownames(data_gland_gene)

data_plot = gather(data_gland_gene, gene, data, "G26977":"G68842")


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

ggsave(P, filename = "Gland_specific_lncRNA_violin_heatmap.pdf", width = 8, height = 6)


rm(list = ls())
dev.off()