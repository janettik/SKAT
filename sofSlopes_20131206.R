##########################################################################
# Purpose: Calculate slopes, reduce to needed variables (SOF)
# Date: 10/18/2013
##########################################################################
library(nlme)
setwd("C:/Users/JNettiksimmons/Documents/DATA")

sof_main0 = read.csv("T120613.csv", header = TRUE)
sof_main = subset(sof_main0, visit <=6)


#---------------------------format data---------------------------------------

names(sof_main) = tolower(names(sof_main))
sof_main[] = lapply(sof_main, function(x){replace(x, x %in% c("A", "M", "", "."), NA)})
sof_main$smmse = as.numeric(as.character(sof_main$smmse))

#identify people with fewer than two cognitive tests
num_mmse = as.integer(
tapply(sof_main$smmse, sof_main$id, function(x){
  length(na.omit(x))
}))

u = unique(sof_main$id)
sof_main$cognum = num_mmse[match(sof_main$id, u)]
sof_cog = subset(sof_main, cognum > 1 )


#-----------------------------------------------------------------
# calculate random slopes
ctrl = lmeControl(opt = 'optim')

lme_age = lme(smmse ~ time + agec70 + agec70*time, 
              random = ~ 1 + time| factor(id), 
              data = sof_cog,
              na.action = na.omit,
              control = ctrl)

rslopes = ranef(lme_age)

sof_slopes = as.data.frame(rslopes)
sof_slopes[,1] = NULL
names(sof_slopes) = 'slope_age'
sof_slopes$id = attr(rslopes,'row.names')
sof_slopes$zslope_age = scale(sof_slopes$slope_age)

#calculate education + age adjusted
lme_age_ed = lme(smmse ~ time + agec70 + agec70*time + edcat + edcat*time, 
              random = ~ 1 + time| factor(id), 
              data = sof_cog,
              na.action = na.omit,
              control = ctrl)

rslopes_ed = ranef(lme_age_ed)
edID = attr(rslopes_ed, 'row.names')
sof_slopes$slope_ageEd = rslopes_ed$time[match(sof_slopes$id, edID)] 
sof_slopes$zslope_ageEd = scale(sof_slopes$slope_ageEd)

#---------------------------------------------------------------------------
# combine slopes back into sof_cog

sof_cog2 = merge(sof_cog, sof_slopes, by.all = "id")
sof_cog3 = subset(sof_cog2, visit == 1)

#-------------------------------------------------------------------------

write.table(sof_cog3, file = "sof_slopes.csv", sep = ",", row.names = FALSE)

