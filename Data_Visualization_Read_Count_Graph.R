######### R script to generate read count vs. position graph for skim-seq population ######

f4 <- fread("DNA191014P04_B11.txt", header = T, check.names = T, data.table = F)
colnames(df4)[1:5] = c('rawread', 'chr', 'pos', 'nread', 'sample' ) ## rawread = read count; chr = chromosome; pos = position; nread = normalized read
head(df4)
p <- ggplot(data = subset(df4, df4$nread<230), aes(pos,nread)) + geom_point(size = 0.4, colour="blue") + facet_wrap( ~ chr, ncol=1 ) +   ## ,scales = "free_x"
  scale_x_continuous(labels=function(x)x/1000000, breaks = seq(0,900000000, by = 100000000),expand = c(0.01, 0)) + xlab("Genomic Position (Mb)") + ylab("Normalized Read Count")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.001, linetype = 'solid',colour = "gray88"),   ##
 axis.line = element_line(colour = "black")
  )
ggsave("DNA191014P04_B11.col.new.github2.pdf", p, width = 6, height = 5, units = "in", limitsize = F) 
