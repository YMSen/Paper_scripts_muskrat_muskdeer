library(ggplot2)

data = read.table(file = "./00_files/Summary_of_the_lengths_of_the_three_versions_musk_deer.txt", header = T, check.names = F)

data_plot = as.data.frame(matrix(ncol = 3, nrow = 1))
colnames(data_plot) = c("Length","Version","Position_percentage")

for (version in c("GCA_006459085.1","GCA_022376915.1","This_study")) {
  data_one = data[data[["Version"]]==version, ]
  data_one_version = as.data.frame(matrix(ncol = 3, nrow = nrow(data_one)*2))
  colnames(data_one_version) = c("Length","Version","Position_percentage")
  data_one_version[["Length"]] = rep(data_one[["Length"]], each = 2)
  data_one_version[["Version"]] = version
  for (num in 1:nrow(data_one)) {
    if(num == 1) {
      data_one_version[["Position_percentage"]][2*num-1] = 0
      data_one_version[["Position_percentage"]][2*num] = data_one[["Position_percentage"]][num]
    }
    if(num != 1) {
      data_one_version[["Position_percentage"]][2*num-1] = data_one[["Position_percentage"]][num-1]
      data_one_version[["Position_percentage"]][2*num] = data_one[["Position_percentage"]][num]
    }
  }
  data_plot = rbind(data_plot, data_one_version)
  data_plot = na.omit(data_plot)
}



P = ggplot(data_plot, aes(Position_percentage, Length/1000000))+ 
  geom_line(aes(colour = Version), linewidth = 0.5)+
  scale_color_manual(values=c("GCA_006459085.1" = "#56a0d3",
                              "GCA_022376915.1" = "#8ec06c",
                              "This_study" = '#ecb731'))+
  theme_test()+
  #geom_hline(yintercept = 50, lty=2, col="grey", lwd=0.5)+ 
  geom_vline(xintercept = 50, lty=2, col="grey", lwd=0.5)+
  geom_vline(xintercept = 90, lty=2, col="grey", lwd=0.5)+
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 202))+
  xlab("Tally of scaffold (%)")+
  ylab("Scaffold length (Mb)")+
  labs(title = "")

ggsave(P, filename = "Comparison of different versions of genome scaffold N50.pdf", width = 5, height = 3.5)


rm(list = ls())
dev.off()
