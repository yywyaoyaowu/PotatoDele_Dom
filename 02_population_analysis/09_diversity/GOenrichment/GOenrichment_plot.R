library(ggplot2)
Args <- commandArgs(TRUE)

#input
input_file <- Args[1]
prefix <- strsplit(input_file, '.', fixed = T)[[1]][1]

Enrich=read.table(input_file,header=T,sep='\t')
Enrich$log.P.value= -log10(as.numeric(Enrich$classicFisher))
dat=Enrich[Enrich$Significant>=5 & Enrich$FC>=2 & Enrich$classicFisher <= 0.01,]
write.table(dat, file=paste0(prefix, "_Sig5_FC2"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

dat$Term=factor(dat$Term,levels=c(dat$Term))
p=ggplot(dat, aes(Term,FC )) +
  geom_point(aes(colour = log.P.value, size =Significant )) +
  scale_colour_gradient(low = "blue", high = "red")+
  coord_flip() + theme_bw()
p.axis=p + theme(axis.text.x = element_text(size = 22, color = "black", vjust =0.5, hjust =0.5, angle =0))+ theme(axis.text.y = element_text(size = 22))
##change xlab and ylab
p.lab=p.axis + xlab("GO term") + ylab("Fold change")+ theme(axis.title.x = element_text(size = 25, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))+ theme(axis.title.y = element_text(size = 25, color = "black"))
p.lenged=p.lab+theme(legend.text = element_text(size = 20))+theme(legend.title = element_text(size = 20))

ggsave(p.lenged,filename=paste0(prefix, "_Sig5_FC2.pdf"),height =8,width = 14)




#p <- ggplot(dat, aes(Term, FC)) +
#  geom_point(aes(colour = pvalue, size = FC)) +
#  scale_color_gradientn(colours =rainbow(10), breaks = c(0, 0.01,0.02,0.03,0.04,0.05), limits = c(0,0.05)) +
#  scale_size_continuous(limits =c(0,700) ,breaks =c(200,400,600)) +
#  #scale_color_continuous(,                         clours - rainbow(6)) +
#  #scale_color_distiller(palette = "Spectral") +
#  #scale_color_gradientn(values = c(0,0.01,0.05,1),
#  #                      colours = c('red','orange','yellow','blue')) +
#  xlab('Go Term')+
#  coord_flip() + theme_bw() +facet_wrap(~face)
#
#pdf('keggsmall.pdf', height = 4, width = 7)
#print(p)
#dev.off()


