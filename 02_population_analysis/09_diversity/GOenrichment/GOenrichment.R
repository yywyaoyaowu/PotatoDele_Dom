library("topGO")
Args <- commandArgs(TRUE)

# Input
geneID2GO_file <- Args[1]
Geneset_file <- Args[2]

# Output
GeneGoInfor_file <- paste0(strsplit(Geneset_file, ".", fixed = T)[[1]][1], "_GoInfo")

geneID2GO <- readMappings(geneID2GO_file)
Geneset <- read.table(Geneset_file, header = F)
geneNames <- names(geneID2GO)
myInterestingGenes <- as.character(Geneset[, 1])
geneList_nscore <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList_nscore) <- geneNames

GOdata_nscore_BP <- new("topGOdata",
                        ontology = "BP",
                        allGenes = geneList_nscore,
                        annot = annFUN.gene2GO, gene2GO = geneID2GO,
                        nodeSize = 5)

result <- runTest(
  GOdata_nscore_BP,
  algorithm = "classic",
  statistic = "fisher"
)

gtFis <- GenTable(GOdata_nscore_BP, classicFisher = result, orderBy = "classic", ranksOf = "classicFisher", topNodes = 100)

go_terms <- Term(gtFis$GO.ID)
gtFis$Term <- go_terms[gtFis$GO.ID]

fdr <- p.adjust(p = gtFis[, "classicFisher"], method = "fdr")

all <- paste(as.character(gtFis[, 'Annotated']), '/', as.character(length(geneScore(GOdata_nscore_BP, use.names = FALSE))), 
             ' ', 
             '(',
             as.character(round(gtFis[, 'Annotated'] * 100 / length(geneScore(GOdata_nscore_BP, use.names = FALSE)), 2)),
             '%)', sep = '')
dis <- paste(as.character(gtFis[, 'Significant']), '/', as.character(numSigGenes(GOdata_nscore_BP)), 
             ' ', 
             '(',
             as.character(round(gtFis[, 'Significant'] * 100 / numSigGenes(GOdata_nscore_BP), 2)),
             '%)', sep = '')

r <- cbind(gtFis[, 1], gtFis[, 2], dis, all, gtFis[, 6], fdr, deparse.level = 0)

write.table(r, file = GeneGoInfor_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

gtFis$classicFisher <- as.numeric(gtFis$classicFisher)
gtFis$classicFisher[is.na(gtFis$classicFisher)] <- 1e-30

gtFis$P.adjust = fdr
gtFis$AllGoGenesNum = length(geneScore(GOdata_nscore_BP))
gtFis$InputGenesNum = numSigGenes(GOdata_nscore_BP)
gtFis$BackGroudRatio = gtFis$Annotated / gtFis$AllGoGenesNum
gtFis$InputRatio = gtFis$Significant / gtFis$InputGenesNum
gtFis$FC = gtFis$InputRatio / gtFis$BackGroudRatio

write.table(gtFis, file = paste0(GeneGoInfor_file, "_GoEnrichment.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
