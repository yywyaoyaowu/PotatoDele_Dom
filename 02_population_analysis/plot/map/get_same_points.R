dat <- read.table('map_info.txt', header = T, sep = "\t")
dat2 <- dat
dat2[dat2$Group == "New Candolleanum", ]$Group = "New_Candolleanum"
dat2[dat2$Group == "Candolleanum", ]$Group = "Previous_Candolleanum"
dat2[dat2$Group == "phureja", ]$Group <- 'Phureja'
dat2[dat2$Group == "ajanhuiri", ]$Group <- 'Ajanhuiri'
dat2[dat2$Group == "stenotomum", ]$Group <- 'Stenotomum'
dat2[dat2$Group == "goniocalyx", ]$Group <- 'Goniocalyx'

dat_diff <- dat2[dat2$Country != dat$State, ]
dat_same <- dat2[dat2$Country == dat$State, ]
write.table(dat_same[, c(1:4)], "dat_same.txt", quote = F, row.names = F, sep='\t')
write.table(dat_diff, "dat_diff.txt", quote = F, row.names = F, sep='\t')


