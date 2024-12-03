library(DESeq2)
library(ggplot2)
library(ggrepel)

countsData = read.table("./00_files/All_samples_counts_1.0_expression_filter.txt", header=T, row.names=1)


info = data.frame(sample = colnames(countsData), condition = as.character(lapply(strsplit(colnames(countsData), "_"), tail, 1)))
write.table(info, file = "info.txt", col.names = T, row.names = F, quote = F, sep = '\t')
condition = as.factor(info[["condition"]])

infoData = data.frame(row.names = info[["sample"]], condition)

dds = DESeqDataSetFromMatrix(countsData, infoData, design= ~ condition)
dds = DESeq(dds)
sample_combine = combn(levels(condition), 2, simplify=F)

table_diff_gene_number = as.data.frame(matrix(ncol =13, nrow = 13))
colnames(table_diff_gene_number) = levels(condition)
rownames(table_diff_gene_number) = levels(condition)

for (i in sample_combine) {
  message("start_", i[1], "_vs_", i[2])
  result01 = results(dds, contrast = c("condition", i[1], i[2]))
  result02 = result01[order(result01$pvalue),]
  filename1 = paste("diff_", i[1], "_vs_", i[2], ".csv", sep = "")
  write.csv(result02, file = filename1)
  
  diff_gene_deseq2 = subset(result02, padj < 0.01 & (log2FoldChange >= 2 | log2FoldChange <= -2))
  
  filename2 = paste("signification_padj0.01_logFold2_diff_", i[1], "_vs_", i[2], ".csv", sep = "")
  write.csv(diff_gene_deseq2, file = filename2)

  

  result02$threshold = factor(ifelse(result02$padj < 0.01 & abs(result02$log2FoldChange) >= 2, ifelse(result02$log2FoldChange >= 2, 'Up', 'Down'), 'Non-Sig'), levels=c('Up', 'Down', 'Non-Sig'))
  filename3 = paste("signification_padj0.01_logFold2_diff_", i[1], "_vs_", i[2], "_marked.csv", sep = "")
  write.csv(result02, file = filename3)
  
  Diff_gene_number = nrow(subset(result02, threshold == "Up" | threshold == "Down"))
  table_diff_gene_number[i[1], i[2]] = Diff_gene_number
  
  diff_gene_VolcanoPlot = ggplot(as.data.frame(result02), aes(x = log2FoldChange, y = -log10(padj), color = threshold, alpha(0.5), ))+
    geom_point(size = 1)+
    scale_color_manual(values = c("#DC143C", "#00008B", "#808080"))+
    theme_test()+
    theme(legend.title = element_blank())+
    ylab('-log10(p-adj)')+
    xlab('log2(FoldChange)')+
    geom_vline(xintercept=c(-2, 2), lty=3, col="black", lwd=0.5) +
    geom_hline(yintercept = -log10(0.01), lty=3, col="black", lwd=0.5)+
    theme_test()

  VolcanoPlot_name = paste("diff_", i[1], "_vs_", i[2], "_0.01_2_VolcanoPlot", ".pdf", sep = "")
  ggsave(diff_gene_VolcanoPlot, filename = VolcanoPlot_name, width = 5, height = 5)
  message("end_", i[1], "_vs_", i[2])
}

filename4 = paste("./signification_padj0.01_logFold2_diff_gene", "_number.csv", sep = "")
write.csv(table_diff_gene_number, file = filename4, col.names = T, row.names = T)


rm(list=ls())
dev.off()