######### R script to generate read count vs. position graph for skim-seq population ######

df4 <- fread("sample_DNA200317P02_G01.in.manuscript.txt", header = T, check.names = T, data.table = F)
colnames(df4)[1:5] = c('rawread', 'chr', 'pos', 'nread', 'sample' )
head(df4)
p <- ggplot(data = subset(df4, df4$nread<2), aes(pos,nread)) + geom_point(size = 1.8, colour="blue") + facet_wrap( ~ chr, ncol=1, strip.position = "right" ) +   ## ,scales = "free_x"
  scale_x_continuous(labels=function(x)x/1000000, breaks = seq(0,900000000, by = 100000000),expand = c(0.01, 0)) + xlab("\nGenomic Position (Mb)") + ylab("Normalized coverage (x)\n")+
  theme(
    
    panel.grid.major.x = element_line(size = 0.08, linetype = 'solid',colour = "gray68"),
    panel.grid.major.y = element_line(size = 0.08, linetype = 'solid',colour = "gray68"), 
    panel.background = element_rect(fill="white"),
        panel.spacing = unit(1.5, "lines"),
        #   strip.text=element_text(vjust=0),
        
        strip.text.y.right = element_text(angle = 0, face="bold", size="14"),
        panel.border=element_rect(colour="black",size=0.4, fill=NA), ##
        #axis.line = element_line(colour = "black")
        axis.text = element_text( size = 10, face="bold" ),
        #axis.text.x = element_text( size = 12 ),
        axis.title = element_text( size = 14, face = "bold" ),
        # legend.position="none",
        # The new stuff
        strip.text = element_text(size = 20)
        
  )
ggsave("sample_DNA200317P02_G01.Figure4.pdf", p, width = 13.5, height = 5, units = "in", limitsize = F) 
