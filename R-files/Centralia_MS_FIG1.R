####################################################################
##############################FIGURE 1##############################
####################################################################

#analysis was performed in R version 3.4.1

#PiCRUSt 16S operon data
#read in PiCRUSt generated operon data
library(readr)
cent_copy <- read_delim("~/cent_copy.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#summarize data
library(plyr)
cent_copy_sum<-ddply(cent_copy, c("Type"), summarize, N=length(Weighted_Copy), mean=mean(Weighted_Copy), sd=sd(Weighted_Copy), se=sd/sqrt(N))


centralia_pi_copy<-ggplot(cent_copy_sum, aes(Type, mean, shape=Type))+
  geom_point(aes(cex=1.5),position=dodge2)+
  scale_y_continuous(trans='log10')+
  xlab("")+
  ylab("Weighted Mean Ribosomal Copy number")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2, cex=0.2), position=dodge2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_shape_manual(values=c(1, 19, 15))

#####PiCRUSt dormancy data
#read in data
library(readr)
cent_dorm <- read_delim("~/cent_dorm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#make it ggplot compatible
library(reshape2)
cent_dorm_m<-melt(cent_dorm)

#summarize the data with plyr
library(plyr)
cent_dorm_sum<-ddply(cent_dorm_m, c("Type", "variable"), summarize, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

#define a position gap for ease of viewing
dodge2 <- position_dodge(width=1)
#plot it
cent_pi_dorm<-ggplot(cent_dorm_sum, aes(variable, mean, shape=Type))+
  geom_point(aes(cex=1.5),position=dodge2)+
  scale_y_continuous(trans='log10')+
  xlab("")+
  ylab("Normalized Gene Abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2, cex=0.2), position=dodge2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_shape_manual(values=c(1, 19, 15))


######Centralia metagenome dormancy genes
#read in data
library(readr)
cent_dorm_normalized <- read_delim("~/cent_dorm_normalized.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#melt it
library(reshape2)
norm_m<-melt(cent_dorm_normalized)

#summarize data
library(plyr)
cent_dormM_sum<-ddply(norm_m, c("Status", "variable"), summarize, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))

dodge2 <- position_dodge(width=1)
cent_meta_dorm<-  ggplot(cent_dormM_sum, aes(variable, mean, shape=Status))+
  geom_point(aes(cex=1.5),position=dodge2)+
  scale_y_continuous(trans='log10')+
  xlab("")+
  ylab("Normalized Gene Abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2, cex=0.2), position=dodge2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_shape_manual(values=c(1, 19, 15))

###Centralia operon count per rplB
#read in data
library(readr)
cent_tRNA <- read_delim("~/cent_tRNA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#summarize data
library(plyr)
cent_trna_sum<-ddply(cent_tRNA, c("Classification"), summarize, N=length(tRNA_rel), mean=mean(tRNA_rel), sd=sd(tRNA_rel), se=sd/sqrt(N))

#plot it
cent_trna<-ggplot(cent_trna_sum, aes(Classification, mean, shape=Classification))+
  geom_point(aes(cex=1.5),position=dodge2)+
  xlab("")+
  ylab("Normalized tRNA Abundance")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.2, cex=0.2), position=dodge2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  scale_shape_manual(values=c(1, 19, 15))

library(gridExtra)
grid.arrange(centralia_pi_copy, cent_pi_dorm, cent_trna, cent_meta_dorm, ncol=2)
#export as PDF 15 by 10


####################################################################
############################Statistics##############################
####################################################################


###Centralia PiCRUSt operon

bartlett.test(Weighted_Copy ~ Type, data=cent_copy)

#Bartlett's K-squared = 20.655, df = 2, p-value = 3.272e-05
#kruskal time

cent_copy$Type<-as.factor(cent_copy$Type)
kruskal.test(Weighted_Copy ~ Type, data=cent_copy)
##Kruskal-Wallis chi-squared = 12.071, df = 2, p-value = 0.002392

#significance, pair-wise comparisons.
library(dunn.test)
dunn.test(cent_copy$Weighted_Copy, cent_copy$Type, method="bh")

#Col Mean-|
#Row Mean |   FireAffe   Recovere
#---------+----------------------
# Recovere |  -1.791798
#|     0.0366
#|
# Referenc |  -3.361201  -2.150941
#|    0.0012*    0.0236*

#All comparisons significant

####Centralia PiCRUSt dormancy
#read in parsed out data to make things easier
centP_dorm_resus <- read_delim("~/centP_dorm_resus.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
centP_dorm_spor <- read_delim("~/centP_dorm_spor.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
centP_dorm_toxin <- read_delim("~/centP_dorm_toxin.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


bartlett.test(Toxin ~ Type, data=centP_dorm_toxin)
#Bartlett's K-squared = 3.625, df = 2, p-value = 0.1632
bartlett.test(resusitation ~ Type, data=centP_dorm_resus)
#Bartlett's K-squared = 5.1753, df = 2, p-value = 0.0752
bartlett.test(sporulation ~ Type, data=centP_dorm_spor)
#Bartlett's K-squared = 41.407, df = 2, p-value = 1.02e-09


toxin_aov<-aov(Toxin ~ Type, data=centP_dorm_toxin)
summary(toxin_aov)
#Df    Sum Sq   Mean Sq F value Pr(>F)
#Type         2 8.785e+09 4.393e+09    0.27  0.767
#Residuals   17 2.768e+11 1.628e+10
#no significance, won't need a Tukey HSD

resus_aov<-aov(resusitation ~ Type, data=centP_dorm_resus)
summary(resus_aov)
#             Df  Sum Sq Mean Sq F value   Pr(>F)    
#Type         2 5777864 2888932   10.92 0.000892 ***
#Residuals   17 4498025  264590 

#pairwise comparison
TukeyHSD(resus_aov)

#diff        lwr       upr     p adj
#Bfire affected-Areference -1425.875 -2233.9467 -617.8033 0.0008266
#Crecovered-Areference     -1223.000 -2031.0717 -414.9283 0.0032492
#Crecovered-Bfire affected   202.875  -456.9128  862.6628 0.714734

#all are difference from fire affected, ref and recovered not significantly different.




centP_dorm_spor$Type<-as.factor(centP_dorm_spor$Type)
kruskal.test(sporulation ~ Type, data=centP_dorm_spor)
#Kruskal-Wallis chi-squared = 12.268, df = 2, p-value = 0.002168
library(dunn.test)
dunn.test(centP_dorm_spor$sporulation , centP_dorm_spor$Type, method="bh")
#Col Mean-|
# Row Mean |   Areferen   Bfire af
#---------+----------------------
# Bfire af |   3.105295
#|    0.0029*
#  |
#  Crecover |   0.862581  -2.746751
#|     0.1942    0.0045*

#all significantly different

#global significance test
cent_dorm_m
bartlett.test(value ~ variable + value, data=cent_dorm_m)

kruskal.test(value ~ variable + value, data=cent_dorm_m)

#####Centralia metagenome dormancy genes
#read in data split by type of dormancy strategy
library(readr)
cent_dorm_resus <- read_delim("~/cent_dorm_resus.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cent_dorm_spor <- read_delim("~/cent_dorm_spor.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
cent_dorm_toxins <- read_delim("~/cent_dorm_toxins.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#remove reference
cent_dorm_resus <- cent_dorm_resus[-12,]
cent_dorm_spor <- cent_dorm_spor[-12,]
cent_dorm_toxins <-cent_dorm_toxins[-12,] 

#test for homogeneity
bartlett.test(Toxin ~ Status, data=cent_dorm_toxins)
#Bartlett test of homogeneity of variances

#data:  Toxin by Status
#Bartlett's K-squared = 4.1343, df = 1, p-value = 0.04202
#not significant, can use regular ANOVA

meta_tox_aov<-aov(Toxin ~ Status, data=cent_dorm_toxins)
summary(meta_tox_aov)


bartlett.test(Resusitation ~ Status, data=cent_dorm_resus)
#Bartlett test of homogeneity of variances

#data:  Resusitation by Status
#Bartlett's K-squared = 4.8899, df = 1, p-value = 0.07701
#not significantly different

meta_resus_aov<-aov(Resusitation ~ Status, data=cent_dorm_resus)
summary(meta_resus_aov)
#Df   Sum Sq   Mean Sq F value Pr(>F)
#Status       1 0.001027 0.0010269   1.282  0.287
#Residuals    9 0.007210 0.0008012 

bartlett.test(Sporulation ~ Status, data=cent_dorm_spor)
#not significant
#Bartlett test of homogeneity of variances

#data:  Sporulation by Status
#Bartlett's K-squared = 0.34452, df = 1, p-value = 0.5572

meta_dorm_spor<-aov(Sporulation ~ Status, data=cent_dorm_spor)
summary(meta_dorm_spor)

#Df    Sum Sq   Mean Sq F value   Pr(>F)    
#Status       1 6.158e-05 6.158e-05   27.18 0.000554 ***
#  Residuals    9 2.039e-05 2.270e-06 

######tRNA operon
#remove reference since only 1 rep
cent_tRNA2<-cent_tRNA[-12,]
bartlett.test(tRNA_rel ~ Classification, data=cent_tRNA2)

#Bartlett test of homogeneity of variances

#data:  tRNA_rel by Classification
#Bartlett's K-squared = 0.37801, df = 1, p-value = 0.5387

#data is parametric, use an regular ANOVA
tRNA_aov<-aov(tRNA_rel ~ Classification, data=cent_tRNA2)
summary(tRNA_aov)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#Classification  1 127.19  127.19   12.96 0.00575 **
# Residuals       9  88.35    9.82  


