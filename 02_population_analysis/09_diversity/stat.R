genome.pi <- read.table("Lan_Can.windowed.pi", header=T,sep="\t")

genome.pi$Can_PI.Lan_PI <- genome.pi$Can_PI / genome.pi$Lan_PI
write.table(genome.pi, "Lan_Can.windowed.pi.dom", quote = F, sep = "\t", row.names = F)
genome.pi.sorted <- genome.pi[order(genome.pi$Can_PI.Lan_PI, decreasing = TRUE), ]
write.table(genome.pi.sorted, "Lan_Can.windowed.pi.dom.sorted", quote = F, sep = "\t", row.names = F)


