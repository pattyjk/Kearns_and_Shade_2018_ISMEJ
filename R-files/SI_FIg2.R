####################################################################
##############################SI FIGURE 2###########################
####################################################################

#analysis performed in R V 3.4
#Correlation between PICRUSt and rrnDB in Centralia
picrust_rrndb <- read.delim("C:/Users/patty/Documents/picrust_rrndb.txt")
picrust_rrndb<-picrust_rrndb[,-4]
library(reshape2)
picrust_rrndb_m<-melt(picrust_rrndb)
library(plyr)
cent_cor<-ddply(picrust_rrndb_m, c("Type", "variable"), summarise, N= length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
library(ggplot2)
cor_pos<-position_dodge(width=0.3)
ggplot(cent_cor, aes(Type, mean, colour=variable))+
  geom_point(aes(cex=1.2), position=cor_pos)+
  scale_size(guide='none')+
  theme_bw()+
  ylab("Weighted Mean Copy Number")+
  xlab("")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2), position=cor_pos)
  
  

cor.test(picrust_rrndb$Weighted_Copy_PC, picrust_rrndb$Weighted_Copy_rrn, method='spearman')

Spearman's rank correlation rho

data:  picrust_rrndb$Weighted_Copy_PC and picrust_rrndb$Weighted_Copy_rrn
S = 7798.7, p-value = 3.152e-09
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.8627374
