library(CMplot)

Args <- commandArgs(T)
input <- Args[1]

data<-read.table(input,header=F)

CMplot(data, plot.type="d", bin.size=1e6, col=c("#FFA12C","#D3632C", "#AD1414"), file="pdf", dpi=300) 

