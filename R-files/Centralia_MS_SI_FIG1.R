####################################################################
##############################SI FIGURE 1###########################
####################################################################
library(readr)
nemergut_operon <- read_delim("~/nemergut_operon.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#change time to factor
nemergut_operon$T<-as.factor(nemergut_operon$T)

#plot the data
library(ggplot2)
nemer_copy<-ggplot(nemergut_operon, aes(T, WeightedCopyNo))+
  geom_boxplot()+
   theme_bw()+
  xlab("Time Since Start of Exp")+
  ylab("Weighted Mean Ribosomal Copy Number")+
  coord_cartesian(ylim=c(0,10))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#read in Ferrenburg et al. (2013) data & plot it
ferrenburg_data <- read_delim("~/ferrenburg_data.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
ferrenburg_copy_no<-  ggplot(ferrenburg_data, aes(reorder(Type, Order), copynum, fill=Type2, colour=Type2))+
  geom_jitter()+
  geom_boxplot()+
  xlab("")+
  ylab("Weighted Mean Ribosomal Copy Number")+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+
  theme_bw()+
  coord_cartesian(ylim=c(0,10))+
  scale_fill_manual(values=c("grey", "white"))+scale_colour_manual(values=c("black", "black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# read in DeAngelis et al. (2015) data and plot it
library(readr)
deangelis_copy <- read_delim("~/deangelis_copy.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
deangelis_copy_no<-  ggplot(deangelis_copy, aes(reorder(warming_treatment_s, Order), Weighted_mean, colour=warming_treatment_s, fill=warming_treatment_s))+
  geom_boxplot(position=dodge2)+
  xlab("")+
  geom_jitter()+
  ylab("Weighted Mean Ribosomal Copy Number")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme_bw()+
  coord_cartesian(ylim=c(0,10))+
  scale_colour_manual(values=c(""))+scale_colour_manual(values=c("black", "black"))+
  scale_fill_manual(values=c("grey", "white"))


#make a 3 panel figure
library(gridExtra)
grid.arrange(nemer_copy, ferrenburg_copy_no, deangelis_copy_no, ncol=3)
#export as PDF with 6 by 24 resolution


#############maths

#variance for Nemergut data
bartlett.test(WeightedCopyNo ~ TimeFull2, data=nemergut_operon)

#Bartlett test of homogeneity of variances

#data:  WeightedCopyNo by TimeFull2
#Bartlett's K-squared = 21.699, df = 6, p-value = 0.001373

#por queeeee, why can't anything have equal variance
nemergut_operon$TimeFull2<-as.factor(nemergut_operon$TimeFull2)
kruskal.test(WeightedCopyNo ~ TimeFull2, data=nemergut_operon)

#Kruskal-Wallis rank sum test

#data:  WeightedCopyNo by TimeFull2
#Kruskal-Wallis chi-squared = 55.572, df = 6, p-value = 3.552e-10
#math works!

#post-hoc with a Dunn test, correct with a BH
library(dunn.test)
dunn.test(nemergut_operon$WeightedCopyNo, nemergut_operon$TimeFull2, method="bh")


###variance for Ferrenburg data
bartlett.test(copynum ~ Type, data=ferrenburg_data)

#Bartlett test of homogeneity of variances

#data:  copynum by Type
#Bartlett's K-squared = 439.33, df = 7, p-value < 2.2e-16

#who like's KW's? Pat does :)
ferrenburg_data$Type<-as.factor(ferrenburg_data$Type)
kruskal.test(copynum ~ Type, data=ferrenburg_data)

#Kruskal-Wallis rank sum test

#data:  copynum by Type
#Kruskal-Wallis chi-squared = 108.97, df = 7, p-value < 2.2e-16

#pariwise comparisons
library(dunn.test)
ferrenburg_data$copynum<-as.numeric(ferrenburg_data$copynum)
dunn.test(ferrenburg_data$copynum ~ ferrenburg_data$Type, method="bh")

###variance for DeAngelis data
bartlett.test(deangelis_copy$Weighted_mean, deangelis_copy$warming_treatment_s)

#surprise surpise, unequal variance
deangelis_copy$warming_treatment_s<-as.factor(deangelis_copy$warming_treatment_s)
kruskal.test(Weighted_mean ~ warming_treatment_s, data=deangelis_copy)

#Kruskal-Wallis rank sum test

#data:  Weighted_mean by warming_treatment_s
#Kruskal-Wallis chi-squared = 19.377, df = 1, p-value = 1.073e-05

#boom, significance

