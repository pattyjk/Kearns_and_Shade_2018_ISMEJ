####################################################################
##############################FIGURE 2##############################
####################################################################
setwd("/Users/patty/Dropbox/R/centralia/")
##Nemergut data

#read in data
library(readr)
nemergut_dorm <- read_delim("~/nemergut_dorm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#make it ggplto compatible 
library(reshape2)
nem_dorm_m<-melt(nemergut_dorm)
#plot it
library(ggplot2)

dodge2<-position_dodge(1)
nemer_dorm<-    ggplot(nem_dorm_m, aes(TimeFull, value, fill=TimeFull))+
  geom_boxplot(position=position_dodge(1))+
  xlab("")+
  ylab("Log Abundance")+
  scale_y_continuous(trans='log10')+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  theme_bw()+
  geom_jitter(position=position_dodge(1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("white", "white", "white", "white", "white", "white", "white", "white", "white"))+
  theme(legend.position="none")+
  facet_wrap(~variable)


##Ferrenburg data
# read in Ferrenburg dormancy data
library(readr)
ferrenburg_dorm <- read_delim("~/ferrenburg_dorm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#make it ggplot friendly
library(reshape2)
ferrenburg_dorm_m<-melt(ferrenburg_dorm)


#plot it
ferrenburg_dormy<-  ggplot(ferrenburg_dorm_m, aes(variable, value, fill=Time_Since_Burn))+
  geom_boxplot(position=dodge2)+
  xlab("")+
  ylab("Log Abundance")+
  scale_y_continuous(trans='log10')+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  theme_bw()+
  facet_wrap(~variable, scales="free")+
  geom_jitter(position=dodge2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white","grey", "white", "grey"))+
    theme(legend.position="none")


ferrenburg_dormy<-ggplot(ferrenburg_dorm_m, aes(variable, value, fill=Time_Since_Burn))+
  geom_boxplot(position=dodge2)+
  xlab("")+
  ylab("Log Abundance")+
    theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  theme_bw()+
  facet_wrap(~variable, scales="free")+
  geom_jitter(position=dodge2)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white","grey", "white", "grey"))+
  theme(legend.position="none")


##DeAngelis data
#read in data
library(readr)
deangelis_dorm <- read_delim("~/deangelis_dorm.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#make ggplot friendly
library(reshape2)
deangelis_dorm_m<-melt(deangelis_dorm)

#plot it
deangelis_dormancy<-  ggplot(deangelis_dorm_m, aes(variable, value, fill=warming_treatment_s))+
  geom_boxplot(position=dodge2)+
  geom_jitter(position=dodge2)+
  xlab("")+
  ylab("Log Abundance")+
  scale_y_continuous(trans='log10')+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("white", "grey"))+
  theme(legend.position="none")

layout_matrix = rbind(c(1,1,1),c(2,3,4))

#make figure
library(gridExtra)
grid.arrange(nemer_dorm, deangelis_dormancy, ferrenburg_dormy,  nrow=2, layout_matrix = rbind(c(1,2),c(3,3)))
#export as PDF with 8 by 14 resolution



####deangelis data
bartlett.test(value ~ variable, data=deangelis_dorm_m)
#Bartlett's K-squared = 1221.5, df = 2, p-value < 2.2e-16

kruskal.test(value ~ variable, data=deangelis_dorm_m)
#Kruskal-Wallis chi-squared = 198.02, df = 2, p-value < 2.2e-16

#####Ferrenburg data
bartlett.test(value ~ Burned, data=ferrenburg_dorm_m)
#Bartlett's K-squared = 22.561, df = 1, p-value = 2.036e-06

ferrenburg_dorm_m$Burned<-as.factor(ferrenburg_dorm_m$Burned)
kruskal.test(value ~ Burned, data=ferrenburg_dorm_m)
#Kruskal-Wallis chi-squared = 8.2258, df = 1, p-value = 0.00413
