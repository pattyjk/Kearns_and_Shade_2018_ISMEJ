#calculating weighted mean copy number

#read in data for the original OTU table. Make sure to remove the first line (#contructed from OTU table) and the # from the first column.
library(readr)
ferrenburg_otu_table <- read_delim("~/ferrenburg-otu_table.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)

#read in data for the copy number correctd OTU table. Make sure to remove the first line (#contructed from OTU table) and the # from the first column.
library(readr)
Ferrenburg_copy_table <- read_delim("~/Ferrenburg_copy_table.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#need to remove the OTU ID's from the first column
Ferrenburg_copy_table<-Ferrenburg_copy_table[,-1]
ferrenburg_otu_table<-ferrenburg_otu_table[,-1]

#Divide the original OTU table by the copy number OTU table and write a new table
copy_number_otus<-ferrenburg_otu_table/Ferrenburg_copy_table

#need to get rid of the NaN's so we can math things
copy_number_otus[is.na(copy_number_otus)] <-0

#calculate the relative abundance of each OTU in the original OTU table (can also get this from QIIME)
ferrenburg_normalized<-sweep(ferrenburg_normalized, 2, colSums(ferrenburg_normalized), '/')

#Multiple the relativized abundance OTU table by the copy number OTU table
weighted_mean_copy<-ferrenburg_normalized*copy_number_otus

#take the column sums for this new file and appened it to the column names for the table
weighted_mean_col_sum<-data.frame(colSums(weighted_mean_copy))

write.table(weighted_mean_col_sum, "copy_number.txt", sep=',')