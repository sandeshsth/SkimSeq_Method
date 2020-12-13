df4 <- fread("DNA191014P01_G11.txt", header = T, check.names = T, data.table = F) # DNA191014P01_G11.txt is input file 
colnames(df4)[1:5] = c('raw.read_count', 'chr', 'pos', 'nread', 'sample' ) ## chr= chromosome; pos= position; nread= normalized read b
head(df4)
p <- ggplot(data = subset(df4, df4$nread<230), aes(pos,nread)) + geom_point(size = 0.8, colour="black") + facet_wrap( ~ chr, ncol=2 ) +   ## ,scales = "free_x"
  scale_x_continuous(labels=function(x)x/1000000, breaks = seq(0,900000000, by = 100000000),expand = c(0.01, 0)) + xlab("Genomic Position (Mb)") + ylab("Normalized Read Count")+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.01, linetype = 'solid',colour = "gray88"),   ##
        axis.line = element_line(colour = "black")
  )
ggsave("DNA191014P01_G11.txt.two.row.pdf", p, width = 9, height = 3, units = "in", limitsize = F) 
