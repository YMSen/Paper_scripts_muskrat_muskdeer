library('pheatmap')
average_matrix = read.csv('./05_top50_cluster_heatmap/whole.top50.vst.genes.matrix.csv', row.names = 1)

dat = t(average_matrix)
dat = scale(dat,
            center=TRUE,
            scale=TRUE)
dat = t(dat)
dat[dat > 1] = 1
dat[dat < -1] = -1

min = min(dat)
max = max(dat)

out_name_8 = paste0("./05_top50_cluster_heatmap/", "top50_cluster_heatmap.pdf")

pheatmap(dat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cellwidth = 12,
         cellheight = 12,
         show_rownames = TRUE,
         show_colnames = TRUE,
         scale = "none",
         border_color = 'NA',
         breaks=seq(min, max, (max - min)/50),
         color = colorRampPalette(c("#1C71B6", "white", "#EE2A29"))(50),
         fontsize = 9,
         fontsize_row = 9,
         fontsize_col = 9,
         angle_col = 45,
         main = "",
         filename = out_name_8)