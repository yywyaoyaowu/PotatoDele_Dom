library("ggplot2")
library("ggmap")
library("sp")
library("maptools")
library("maps")
library("ggpubr")
setwd("/Volumes/外置硬盘/Users/starry_sky/文件/马铃薯重测序-有害突变/dom_del_final/材料地图")

dat_same <- read.table('dat_same_points.txt', header = F, sep = "\t")
dat_diff <- read.table('dat_diff.txt', header = T, sep = "\t")
colnames(dat_same) <- colnames(dat_diff)

dat4 <- rbind(dat_diff, dat_same)
dat3 <- rbind(dat_diff, dat_same)
dat3[dat3$Group == "Stenotomum", ]$Group <- 'Landrace'
dat3[dat3$Group == "Goniocalyx", ]$Group <- 'Landrace'
dat3[dat3$Group == "Phureja", ]$Group <- 'Landrace'
dat3[dat3$Group == "Ajanhuiri", ]$Group <- 'Landrace'

dat_diff <- dat3[dat3$Country != dat3$State, ]
dat_same <- dat3[dat3$Country == dat3$State, ]

dat_diff_New_Can <- dat_diff[dat_diff$Group == "New_Candolleanum", ]
dat_diff_Previous <- dat_diff[dat_diff$Group != "New_Candolleanum", ]
dat_same_New_Can <- dat_same[dat_same$Group == "New_Candolleanum", ]
dat_same_Previous <- dat_same[dat_same$Group != "New_Candolleanum", ]

mp<-NULL
mapworld<-borders("world",colour = "white",fill="#DCDCDC",alpha=0.9)
mp <- ggplot() + mapworld + ylim(-60,15) + xlim(-85,-30)

pal <- c("palegreen4", "#F0A419","#F3C97F")
dat3$Group <- factor(dat3$Group,levels=c("Landrace","New_Candolleanum","Previous_Candolleanum"))

p <- mp + geom_point(data=dat_same_Previous, aes(x=Longitude, y=Latitude, colour = Group), size=1.4, alpha=0.9, shape = 16)
p <- p + geom_point(data=dat_same_New_Can, aes(x=Longitude, y=Latitude, colour = Group), size=1.4, alpha=0.9, shape = 16)

p <- p + geom_point(data=dat_diff_Previous, aes(x=Longitude, y=Latitude, colour = Group), size=1.4, alpha=0.9, shape = 16,
                     position = position_jitter(width = 0.8, height = 0.8, seed = 0))
p + geom_point(data=dat_diff_New_Can, aes(x=Longitude, y=Latitude, colour = Group), size=1.4, alpha=0.9, shape = 16,
               position = position_jitter(width = 0.8, height = 0.8, seed = 0)) +
  scale_color_manual(values=c("New_Candolleanum" = "#F0A419",
                              "Previous_Candolleanum" = "#F3C97F",
                              "Landrace" = "palegreen4")) + theme_classic() +
  theme(legend.key.size = unit(0.6,'cm'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = c(0.8,0.9)) +
  guides(color = guide_legend(title = "Group", override.aes = list(size = 5)))
ggsave("map_Landrace1_alpha0.9.pdf",height = 16,width=10.65,units="cm")


##### Landrace split 3
dat4_diff <- dat4[dat4$Country != dat4$State, ]
dat4_same <- dat4[dat4$Country == dat4$State, ]

mp<-NULL
mapworld<-borders("world",colour = "white",fill="#DCDCDC",alpha=0.9)
mp <- ggplot() + mapworld + ylim(-60,15) + xlim(-85,-30)

dat4$Group <- factor(dat4$Group,levels=c("Ajanhuiri", "Stenotomum", "Goniocalyx", "Phureja", "New_Candolleanum", "Previous_Candolleanum"))

p <- mp + geom_point(data=dat4_same, aes(x=Longitude, y=Latitude, colour = Group), size=1, alpha=1)
p + geom_point(data=dat4_diff, aes(x=Longitude, y=Latitude, colour = Group), size=1, alpha=1,
               position = position_jitter(width = 0.8, height = 0.8, seed = 2)) +
  scale_color_manual(values=c("New_Candolleanum" = "#F0A419",
                              "Previous_Candolleanum" = "#F3C97F",
                              "Stenotomum" = "#8CC63E",
                              "Goniocalyx" = "#8CA0CB",
                              "Phureja" = "#56B4E9",
                              "Ajanhuiri" = "#009E73")) + theme_classic() +
  theme(legend.key.size = unit(0.6,'cm'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.position = c(0.8,0.9)) +
  guides(color = guide_legend(title = "Group", override.aes = list(size = 5)))
ggsave("map_Landrace3.pdf",height = 18,width=18,units="cm")




###### pie
data <- as.data.frame(table(dat3$Group))
colnames(data) <- c('class', 'count')
data$percent <- data$count / sum(data$count)
colnames(data) <- c('class','count','percent')
data$percent<-paste0(round(data$percent,2)*100,"%")
data$class1<-paste0(data$class,"(",data$percent,")")

ggpie(data, "count", fill = 'class', color="white",
         label = 'percent', lab.font = c(5, "black"), lab.pos = "in")+
  scale_fill_manual(values = c("New_Candolleanum" = "#F0A419",
                               "Previous_Candolleanum" = "#F3C97F",
                               "Landrace" = "palegreen4")) +
  theme(legend.position = "none", legend.text = element_text(size=10))+
  labs(fill="Group", family='Arial')
ggsave("pie_Landrace1.pdf",height = 6,width=6,units="cm")


##### pie Landrace split 3
data <- as.data.frame(table(dat4$Group))
colnames(data) <- c('class', 'count')
data$percent <- data$count / sum(data$count)
colnames(data) <- c('class','count','percent')
data$percent<-paste0(round(data$percent,2)*100,"%")
data$class1<-paste0(data$class,"(",data$percent,")")

data$class <- factor(data$class, levels = c("Ajanhuiri", "Stenotomum", "Goniocalyx", "Phureja", "New_Candolleanum", "Previous_Candolleanum"))

ggpie(data, "count", fill = 'class', color="white",
      label = 'percent', lab.font = c(5, "black"), lab.pos = "in")+
  scale_fill_manual(values = c("New_Candolleanum" = "#F0A419",
                               "Previous_Candolleanum" = "#F3C97F",
                               "Stenotomum" = "#8CC63E",
                               "Goniocalyx" = "#8CA0CB",
                               "Phureja" = "#56B4E9",
                               "Ajanhuiri" = "#009E73")) +
  theme(legend.position = "none")+
  labs(fill="Group", family='Arial')
ggsave("pie_Landrace3.pdf", height = 6, width=6, units="cm")
 






