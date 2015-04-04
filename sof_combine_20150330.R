##########################################################################
# Purpose: Combine SOF gwas and phenotype data
# Date: 12/06/2013
##########################################################################
setwd("C:/Users/JNettiksimmons/Documents/DATA")

sof_main = read.csv("sof_slopes.csv", header = TRUE)
names(sof_main) = tolower(names(sof_main))
load("sof_impC_20150313.rdata")
#---------------------------format data---------------------------------------
# get rid of SOF people without gwas data
imp_id = sof_imp_clean$id
#sof_main2 = subset(sof_main, id %in% imp_id)

#imputed
sof_ci = merge(sof_main, sof_imp_clean, by.all = 'id')
save(sof_ci, file = "sof_ci_20150330.rdata")

#baseline imputed
sof_cib = subset(sof_ci, visit == 1)
save(sof_cib, file = "sof_cib_20150330.rdata")
