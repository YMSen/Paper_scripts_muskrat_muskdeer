library(ggplot2)
library(ggforce)
library(Rtsne)

data = read.table(file = "./00_files/All_samples_tpm_1.0_expression_filter.txt", header = T, check.names = F, row.names = 1)
info = data.frame(sample = colnames(data), condition = as.character(lapply(strsplit(colnames(data), "_"), tail, 1)))

data = log2(data+1)

genes = names(head(sort(apply(data,1,function(x) 100*sd(x)/mean(x)),decreasing=T), 20024))

data = as.data.frame(t(data[genes,]))
data$class = info$condition
data = as.data.frame(data)
data_unique = unique(data)
data_matrix = as.matrix(data_unique[,1:20024])

set.seed(125)
tsne_out = Rtsne(data_matrix, pca=T, dims=2, perplexity=23, theta=0.05)

tsne_res = as.data.frame(tsne_out$Y)
colnames(tsne_res) = c("tSNE1", "tSNE2")
tsne_res$sample = info$sample
tsne_res$condition = info$condition


P = ggplot(tsne_res, aes(tSNE1, tSNE2))+ 
  geom_point(aes(shape = condition, colour = condition), size = 1.5)+
  geom_mark_ellipse(aes(color = condition), expand = unit(0, "mm"))+
  scale_fill_manual(values=c('#f1c933', '#a9eee6', '#626f92', "#f5c7f7", "#a393eb",
                             "#ffaa64", '#7098da', '#3fbac2', '#de95ba', "#cd595a",
                             "#ff9900", "#6a7efc", "#7cbd1e"))+
  scale_color_manual(values=c('#f1c933', '#a9eee6', '#626f92', "#f5c7f7", "#a393eb",
                              "#ffaa64", '#7098da', '#3fbac2', '#de95ba', "#cd595a",
                              "#ff9900", "#6a7efc", "#7cbd1e"))+
  scale_shape_manual(values = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))+
  theme_test()+
  geom_hline(yintercept = 0, lty=2, col="grey", lwd=0.3)+ 
  geom_vline(xintercept = 0, lty=2, col="grey", lwd=0.3)+
  xlab("t-SNE 1")+
  ylab("t-SNE 2")+
  labs(title = "")

ggsave(P, filename = "t-SNE plot.pdf", width = 6.8, height = 6)


rm(list = ls())
dev.off()