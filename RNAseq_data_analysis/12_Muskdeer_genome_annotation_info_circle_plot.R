library(RCircos)
library(RIdeogram)


karyotype <- read.table("./00_files/muskdeer/muskdeer_karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)

muskdeer_gff_file <- read.table("./00_files/muskdeer/Musk_deer.gff", sep = "\t", header = T, stringsAsFactors = F)


gene_density <- GFFex(input = "./00_files/muskdeer/muskdeer_karyotype.txt", feature = "gene", window = 1000000)

gene_density[["Value"]] = log2(gene_density[["Value"]] + 1)


miRNA_density <- GFFex(input = "./00_files/muskdeer/miRNA_musk_deer.gff3", karyotype = "./00_files/muskdeer/muskdeer_karyotype.txt", feature = "miRNA", window = 1000000)

miRNA_density[["Value"]] = log2(miRNA_density[["Value"]] + 1)


snRNA_density <- GFFex(input = "./00_files/muskdeer/snRNA_musk_deer.gff3", karyotype = "./00_files/muskdeer/muskdeer_karyotype.txt", feature = "snRNA", window = 1000000)

snRNA_density[["Value"]] = log2(snRNA_density[["Value"]] + 1)


rRNA_density <- GFFex(input = "./00_files/muskdeer/rRNA_musk_deer.gff3", karyotype = "./00_files/muskdeer/muskdeer_karyotype.txt", feature = "rRNA", window = 1000000)

rRNA_density[["Value"]] = log2(rRNA_density[["Value"]] + 1)


tRNA_density <- GFFex(input = "./00_files/muskdeer/tRNA_musk_deer.gff3", karyotype = "./00_files/muskdeer/muskdeer_karyotype.txt", feature = "tRNA", window = 1000000)

tRNA_density[["Value"]] = log2(tRNA_density[["Value"]] + 1)


cyto.info = read.table("./00_files/muskdeer/muskdeer_karyotype.txt", header = T)

cyto.info[["Band"]] <- NA

cyto.info[["Stain"]] <- cyto.info[["Chr"]]

chr.exclude <- NULL

tracks.inside <- 5

tracks.outside <- 0

RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)


rcircos.params = RCircos.Get.Plot.Parameters()

rcircos.params[["heatmap.color"]] = "BlueWhiteRed"

rcircos.params[["scatter.color"]] = "green"

rcircos.params[["track.background"]] = "white"

rcircos.params[["grid.line.color"]] = "white"

rcircos.params[["line.color"]] = "skyblue"

rcircos.params[["hist.color"]] = "brown"

RCircos.Reset.Plot.Parameters(rcircos.params)


pdf(file = "RCircos_Muskdeer_Genome_annotation_distribution.pdf", height=8, width=8)

RCircos.Set.Plot.Area()

RCircos.Chromosome.Ideogram.Plot()


data.col = 5

track.num = 1

side = "in"

gene_density[["Gene_Name"]] = NA

gene_density = gene_density[ ,c("Chr","Start","End","Gene_Name","Value")]

RCircos.Heatmap.Plot(gene_density, data.col, track.num, side)


rcircos.params = RCircos.Get.Plot.Parameters()

rcircos.params[["heatmap.color"]] = "BlackOnly"

RCircos.Reset.Plot.Parameters(rcircos.params)


data.col = 5

track.num = 2

side = "in"

tRNA_density[["tRNA_Name"]] = NA

tRNA_density = tRNA_density[ ,c("Chr","Start","End","tRNA_Name","Value")]

RCircos.Heatmap.Plot(tRNA_density, data.col, track.num, side)


data.col = 5

track.num = 3

side = "in"

by.fold = 1

miRNA_density[["miRNA_Name"]] = NA

miRNA_density = miRNA_density[ ,c("Chr","Start","End","miRNA_Name","Value")]

RCircos.Scatter.Plot(miRNA_density, data.col, track.num, side, by.fold)


data.col = 5

track.num = 4

side = "in"

snRNA_density[["snRNA_Name"]] = NA

snRNA_density = snRNA_density[ ,c("Chr","Start","End","snRNA_Name","Value")]

RCircos.Line.Plot(snRNA_density, data.col, track.num, side)


data.col = 5

track.num = 5

side = "in"

rRNA_density[["rRNA_Name"]] = NA

rRNA_density = rRNA_density[ ,c("Chr","Start","End","rRNA_Name","Value")]

RCircos.Histogram.Plot(rRNA_density, data.col, track.num, side)


rm(list = ls())
dev.off()