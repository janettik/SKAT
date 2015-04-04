####################################################################
# Purpose: Format SKAT results
# Date: 11/7/2013
####################################################################
setwd("C:/Users/JNettiksimmons/Documents/DATA")
load("sof_impResult_20150330.rdata")


sof_res = list(sof_ires)

sof_resf1 = do.call('cbind', 
  lapply(sof_res, function(x){
    
    L = length(x)
    do.call('rbind', lapply(1:L, function(j){
      temp = as.data.frame(x[[j]])
      temp$gene = names(x)[j]
      temp
    }))

}))

write.csv(sof_resf1,file ="C:/Users/JNettiksimmons/Documents/output/sof_skatResults_20150330.csv",
            row.names = FALSE)

#-----------------------------------------------------------------------------------------------
# combine mros results (noBirm)


load("mros_imputedResult_noBirm_20150330.rdata")


mros_res = list( mros_ires)

mros_resf1 = do.call('cbind', 
  lapply(mros_res, function(x){
    
    L = length(x)
    do.call('rbind', lapply(1:L, function(j){
      temp = as.data.frame(x[[j]])
      temp$gene = names(x)[j]
      temp
    }))

}))

write.csv(mros_resf1,file ="C:/Users/JNettiksimmons/Documents/output/mros_skatResults_20150330.csv",
            row.names = FALSE)
