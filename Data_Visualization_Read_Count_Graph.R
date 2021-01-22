######### R script to generate read count vs. position graph for skim-seq population ######

df4 <- fread("DNA191014P04_B11.TA5798.txt", header = T, check.names = T, data.table = F) ## wheat-barley recombinant- homozygous translocation
colnames(df4)[1:5] = c('rawread', 'chr', 'pos', 'nread', 'sample' )
#head(df4)
p <- ggplot(data = subset(df4, df4$nread<200), aes(pos,nread)) + geom_point(size = 0.8, colour="blue") + facet_wrap( ~ chr, ncol=1, strip.position = "right" ) +   ## ,scales = "free_x"
  scale_x_continuous(labels=function(x)x/1000000, breaks = seq(0,900000000, by = 100000000),expand = c(0.01, 0)) + xlab("Genomic Position (Mb)") + ylab("Normalized Read Count")+
  theme(panel.background = element_rect(fill="white"),
        panel.spacing = unit(1.5, "lines"),
        strip.text.y.right = element_text(angle = 0, face="bold", size="14"),
        panel.grid.major = element_line(size = 0.0001, linetype = 'solid',colour = "gray95"), 
    panel.border=element_rect(colour="black",size=0.1, fill=NA), 
 axis.text = element_text( size = 10, face="bold" ),
 axis.title = element_text( size = 14, face = "bold" ),
 strip.text = element_text(size = 20)

  )
ggsave("DNA191014P04_B11.TA5798.pdf", p, width = 8, height = 5.5, units = "in", limitsize = F) 
