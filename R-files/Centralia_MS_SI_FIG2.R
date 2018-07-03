####################################################################
##############################SI FIGURE 2###########################
####################################################################
#analyses performed in R version 3.4.1

#rrnDB and dormancy 
library(readr)

#read data
rrndb_copytxt <- read_delim("~/rrndb_copytxt.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#make a histogram of the distribution of the data
hist(rrndb_copytxt$Monkey, xlab="16S rRNA Operon No", main="")

#read in BLAST data
rrndb_dorm <- read_delim("~/rrndb_dorm.txt", "\t", escape_double = FALSE, col_types = cols(Resuscitation = col_integer()), na = "NA", trim_ws = TRUE)

#make the data ggplot ready, and plot it
library(reshape2)
library(ggplot2)
rrndb_m<-melt(rrndb_dorm)
str(rrndb_m)
ggplot(rrndb_m, aes(variable, value))+
  geom_boxplot()+
  ylab("Ribosomal  Operon Count")+
  xlab("")+
  theme_bw()+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#export as PDF with 5 by 7" resolution

##maths

bartlett.test(rrndb_dorm)

#Bartlett test of homogeneity of variances

#data:  rrndb_dorm
#Bartlett's K-squared = 257.29, df = 3, p-value < 2.2e-16

#variances aren't homogenous, need a non-parametric test

kruskal.test(value ~ variable,  data=rrndb_m)

#Kruskal-Wallis rank sum test

#data:  value by variable
#Kruskal-Wallis chi-squared = 1326.6, df = 3, p-value <
# 2.2e-16

#there's significance here, so we'll do a multiple comparisons with a Dunn test

library(dunn.test)
dunn.test(rrndb_m$value, rrndb_m$variable, method="bh")

#Kruskal-Wallis rank sum test

#data: x and group
#Kruskal-Wallis chi-squared = 1326.6029, df = 3, p-value = 0


#Comparison of x by group                            
#(Benjamini-Hochberg)                              
#Col Mean-|
# Row Mean |       None   Resuscit   Sporulat
---------+---------------------------------
  #Resuscit |  -33.63221
  #|    0.0000*
  #|
  #Sporulat |  -20.84430   8.874946
  #|    0.0000*    0.0000*
  #|
  #Toxin |  -7.992161   23.99443   12.83431
  #|    0.0000*    0.0000*    0.0000*
  