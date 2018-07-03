####################################################################
####################Global Correlation##############################
###################Dorm genes and operons###########################
####################################################################
#analyses done using R 3.4.1.


library(readr)
operon_dom_full <- read_delim("~/operon_dom_full.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

cor.test(operon_dom_full$Weighted_Copy, operon_dom_full$Sporulation, method='spearman')
#p<2e-16, rho=0.44
cor.test(operon_dom_full$Weighted_Copy, operon_dom_full$Resuscitation, method='spearman')
#p<2e-16, rho=0.73
cor.test(operon_dom_full$Weighted_Copy, operon_dom_full$Toxins, method='spearman')
#p<2e-16, rho=0.78

