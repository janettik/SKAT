##########################################################################
# Purpose: perform skat on sof and mros os (no birmingham)(imputed only)
# Date: 07/23/2014
##########################################################################
library('doHPC')
registerDoHPC('HOMER',Rbin="C:/R/R-2.15.2/bin/i386/Rterm.exe")

setwd("C:/Users/JNettiksimmons/Documents/DATA")
load("mros_cib_noBirm_20150330.rdata")
load("sof_cib_20150330.rdata")
load("snpedia_20150330.rdata")
allsnps = snpedia

source("C:/Users/JNettiksimmons/Documents/indiv_skatFunctions_20140722.r") 

#####################################################################
# gene list

geneList = c('abca7', 'picalm', 'clu', 'cr1' , 'ms4a6a', 'bin1', 'cd33', 
             'epha1', 'cd2ap', 'hla', 'ptk2b', 'sorl1', 'slc24a4',
			 'inpp5d', 'mef2c', 'nme8', 'zcwpw1', 'celf1', 'fermt2', 'cass4', 'ms4a6e')

###########################################################
# Imputed


sof_out_i = list()
for(j in 1:length(geneList)){
  
  g = geneList[j]
  sof_out_i[[j]] = 
        unlist( foreach(j=1:1,.packages="SKAT") %dopar% {  
      runAll(sof_cib, gene = g , outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
    }, recursive = FALSE)
  
  print(j)
}



sof_ires = formatResult(sof_out_i)
names(sof_ires) = geneList
save(sof_ires, file = "sof_impResult_20150330.rdata")

###########################################################
# Imputed


mros_out_i = list()
for(j in 1:length(geneList)){

g = geneList[j]
mros_out_i[[j]]= unlist( foreach(j=1:1,.packages="SKAT") %dopar% {  
  runAll(mros_cib_nb, gene = g , outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
}, recursive = FALSE)

print(j)

}

i.res = formatResult(mros_out_i)
names(i.res) = geneList
mros_ires = i.res  

save(mros_ires, file = "mros_imputedResult_noBirm_20150330.rdata")



 