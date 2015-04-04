##########################################################################
# Purpose: SOF Subset
# Date: 12/06/2013
##########################################################################
library(nlme)
setwd("C:/Users/JNettiksimmons/Documents/DATA")

sof_main0 = read.csv("mergedSof.csv", header = TRUE)
names(sof_main0) = tolower(names(sof_main0))
#-----------------------------------------------------------------------

sof_main0[] = lapply(sof_main0, function(x){replace(x, x %in% c("A", "M", "", "."), NA)})
sof_main0$smmse = as.numeric(as.character(sof_main0$smmse))

#calculate age at baseline
u = unique(sof_main0$id)
sof_main0$ageb = as.integer(as.character(sof_main0$v1age))
sof_main0$agec70 = sof_main0$ageb-70

sof_main0$race = sof_main0$v1race

#make education categories - given in years
#making same cats as sof (less than hs, hs/some college, 4-year degree+)
sof_main0$v1educ = as.integer(as.character(sof_main0$v1educ))
ed = rep(NA, nrow(sof_main0))
ed[sof_main0$v1educ < 12] = 1
ed[sof_main0$v1educ %in% c(12,13,14,15)] = 2
ed[sof_main0$v1educ >= 16] = 3
sof_main0$edcat = factor(ed)
contrasts(sof_main0$edcat) = rbind( c(1,0), c(0,0), c(0,1))

#--------------------------------------------------------------------------
# subset down to smaller data


# NOTE: Missing depression, didn't get v2 quality of life merged before sas
#       license expiration. 

sof_select = c("id", "time", "visit", "smmse","ageb", "agec70", "race", "edcat", 
               "v1age", "v1race","v1educ",
               "v1diabcl","v1diabyr", "v1ediab", "v1skchol",
                "v1bmi" ,"v1htcm25", "v1kgs25",
               "v1hyten", "v1smkevr", "v1smklve","v1smknow","v1smksme",
               "v1smoke", "v1dr30", "v1drink", "v1drwk30", "v1pactwk",
               "v1ttkcal")
sof_main = subset(sof_main0, select = sof_select)

write.csv(sof_main,"sof_longSubset20131206.csv", na = "")
