Args <- commandArgs(T)
recombRate <- Args[1]
recombRate_sort <- Args[2]

recombRate <- read.table(recombRate, header=T,sep="\t")

recombRate.sorted <- recombRate[order(recombRate$recombRate, decreasing = TRUE), ]
write.table(recombRate.sorted, recombRate_sort, quote = F, sep = "\t", row.names = F)


