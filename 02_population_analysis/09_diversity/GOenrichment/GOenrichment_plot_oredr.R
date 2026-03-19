library(ggplot2)
Args <- commandArgs(TRUE)

# input
input_file <- Args[1]
prefix <- strsplit(input_file, '.', fixed = T)[[1]][1]

Enrich <- read.table(input_file, header = T, sep = '\t')
Enrich$log.P.value <- -log10(as.numeric(Enrich$classicFisher))
dat <- Enrich[Enrich$Significant >= 5 & Enrich$FC >= 2 & Enrich$classicFisher <= 0.01, ]
write.table(dat, file = paste0(prefix, "_Sig5_FC2"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

if (length(Args) >= 2 && Args[2] == "FC") {
    dat <- dat[order(-dat$FC), ]
    dat$Term <- factor(dat$Term, levels = dat$Term)
} else if (length(Args) >= 2 && Args[2] == "Pvalue") {
    dat <- dat[order(dat$classicFisher), ]
    dat$Term <- factor(dat$Term, levels = dat$Term)
} else if (length(Args) >= 2 && Args[2] == "Significant") {
    dat <- dat[order(-dat$Significant), ]
    dat$Term <- factor(dat$Term, levels = dat$Term)
} else if (length(Args) >= 2 && Args[2] == "file" && length(Args) >= 3) {
    order_file <- Args[3]
    term_order <- read.csv(order_file, header = FALSE, stringsAsFactors = FALSE)$V1
    term_order <- rev(term_order[term_order %in% dat$Term])
    dat$Term <- factor(dat$Term, levels = term_order)
    dat <- dat[!is.na(dat$Term), ]
} else {
    dat$Term <- reorder(dat$Term, -dat$FC)
}

p <- ggplot(dat, aes(Term, FC)) +
  geom_point(aes(colour = log.P.value, size = Significant)) +
  scale_colour_gradient(low = "blue", high = "red") +
  coord_flip() + theme_bw()

p.axis <- p + theme(axis.text.x = element_text(size = 22, color = "black", vjust = 0.5, hjust = 0.5, angle = 0)) + 
  theme(axis.text.y = element_text(size = 22))

p.lab <- p.axis + xlab("GO term") + ylab("Fold change") + 
  theme(axis.title.x = element_text(size = 25, color = "black", vjust = 0.5, hjust = 0.5, angle = 0)) + 
  theme(axis.title.y = element_text(size = 25, color = "black"))

p.lenged <- p.lab + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20))

ggsave(p.lenged, filename = paste0(prefix, "_Sig5_FC2.pdf"), height = 8, width = 14)

