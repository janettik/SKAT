##########################################################################
# Purpose: Calculate slopes, merge together MrOS gWAS and phenotype data
#  EXCLUDE BIRMINGHAM SUBJECTS
# Date: 11/08/2013
##########################################################################
library(nlme)
setwd("C:/Users/JNettiksimmons/Documents/DATA")

mros_main = read.csv("mergedMrOS_noB.csv", header = TRUE)
load("mros_impC_20150330.rdata")
#---------------------------format data---------------------------------------

names(mros_main) = tolower(names(mros_main))
mros_main[] = lapply(mros_main, function(x){replace(x, x %in% c("A", "M", "", "."), NA)})
mros_main$tmmscore = as.numeric(as.character(mros_main$tmmscore))

#calculate age at baseline
u = unique(mros_main$id)
ageb = as.integer(
tapply(mros_main$age, mros_main$id, function(x){
  
  a = min(x, na.rm=T)
}))

#filldown race
ethrace = unlist(lapply(u, function(x){
          temp = subset(mros_main, id == x, select = c('id', 'gierace'))
          n = nrow(temp)
          race = na.omit(temp$gierace)[1]
}))

mros_main$ethrace = ethrace[match(mros_main$id, u)]
mros_main$ageb = ageb[match(mros_main$id, u)]
mros_main$agec70 = mros_main$ageb-70 

#make education categories
ed = rep(NA, nrow(mros_main))
ed[mros_main$gieduc %in% c(1 ,2, 3)] = 1
ed[mros_main$gieduc %in% c(4, 5)] = 2
ed[mros_main$gieduc %in% c(6, 7, 8)] = 3
mros_main$edcat = factor(ed)
contrasts(mros_main$edcat) = rbind( c(1,0), c(0,0), c(0,1))

#identify people with fewer than two cognitive tests
num3ms = as.integer(
tapply(mros_main$tmmscore, mros_main$id, function(x){
  length(na.omit(x))
}))

mros_main$cognum = num3ms[match(mros_main$id, u)]
mros_cog = subset(mros_main, cognum > 1)


#-----------------------------------------------------------------
# calculate random slopes

lme_age = lme(tmmscore ~ time + agec70 + agec70*time, 
              random = ~ 1 + time| id, 
              data = mros_cog,
              na.action = na.omit)

rslopes = ranef(lme_age)

mros_slopes = as.data.frame(rslopes)
mros_slopes[,1] = NULL
names(mros_slopes) = 'slope_age'
mros_slopes$id = attr(rslopes,'row.names')
mros_slopes$zslope_age = scale(mros_slopes$slope_age)

#calculate education + age adjusted
lme_age_ed = lme(tmmscore ~ time + agec70 + agec70*time + edcat + edcat*time, 
              random = ~ 1 + time| id, 
              data = mros_cog,
              na.action = na.omit)

rslopes_ed = ranef(lme_age_ed)
mros_slopes$slope_ageEd = rslopes_ed$time 
mros_slopes$zslope_ageEd = scale(mros_slopes$slope_ageEd)

#---------------------------------------------------------------------------
# combine slopes back into mros_cog

mros_cog2 = merge(mros_cog, mros_slopes, by.all = "id")

#-------------------------------------------------------------------------
# make reduced cog data set containing only variables we need

mros_keep = c('id', 'time', 'visit_order', 'site', 'ageb', 'agec70', 'ethrace', 'edcat', 'gierace',
              'tmmscore',  'cognum',
              'slope_age', 'zslope_age', 'slope_ageEd', 'zslope_ageEd',
              'mhdiab', 'mhdiabt', 'd1fgluc', 'ast_chol', 
              'hwbmi', 'mhhgtcm', 'mhwgtkg', 'bpaaicat', 'mhbp', 'mhbpt',
              'mhdepr', 'mhdeprt', 'dpgds15', 'lfsm6qt1', 'tursmoke',
              'tudramt', 'tudravg', 'tudrinwk', 'tudrprwk')

mros_cog3 = subset(mros_cog2, select = mros_keep)

mros_cog3$diab_self = mros_cog3$mhdiab
mros_cog3$diabtx_self = ifelse(mros_cog3$diab_self == 1 & mros_cog3$mhdiabt == 1, 1, 0 )
mros_cog3$diabtx_self[is.na(mros_cog3$diab_self)] = NA

mros_cog3$depr_self = mros_cog3$mhdepr
mros_cog3$deprtx_self = ifelse(mros_cog3$depr_self == 1 & mros_cog3$mhdeprt == 1, 1, 0)
mros_cog3$deprtx_self[is.na(mros_cog3$depr_self)] = NA

mros_cog3$highbp_self = mros_cog3$mhbp
mros_cog3$highbptx_self = ifelse(mros_cog3$highbp_self ==1 & mros_cog3$mhbpt ==1, 1, 0)
mros_cog3$highbptx_self[is.na(mros_cog3$highbp_self)] = NA

mros_cog3$fastgluc = mros_cog3$d1fgluc
mros_cog3$v1cholest = mros_cog3$ast_chol
mros_cog3$gds15 = mros_cog3$dpgds15

mros_cog3 = mros_cog3[order(mros_cog3$id, mros_cog3$visit_order),]
#--------------------------------------------------------------------------
# combine cog data with gwas data
# (make imputed version and genotyped version)

#genotyped
# mros_cg_nb = merge(mros_cog3, mros_geno_clean, by.all = 'id')
# save(mros_cg_nb, file = "mros_cg_noBirm_20131122.rdata")
# 
# #baseline genotyped
# mros_cgb_nb = subset(mros_cg_nb, visit_order == 1)
# save(mros_cgb_nb, file = "mros_cgb_noBirm_20131122.rdata")

#imputed
mros_ci_nb = merge(mros_cog3, mros_imp_clean, by.all = 'id')
save(mros_ci_nb, file = "mros_ci_noBirm_20150330.rdata")

#baseline imputed
mros_cib_nb = subset(mros_ci_nb, visit_order == 1)
save(mros_cib_nb, file = "mros_cib_noBirm_20150330.rdata")

