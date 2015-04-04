##########################################################################
# Purpose: SOF/MROS ID'd snp table
# Date:06/19/2014
##########################################################################
library('doHPC')
registerDoHPC('HOMER',Rbin="C:/R/R-2.15.2/bin/i386/Rterm.exe")

setwd("C:/Users/JNettiksimmons/Documents/DATA")
load("mros_cgb_NoBirm_20131122.rdata")
load("mros_cib_NoBirm_20131122.rdata")
load("sof_cgb_20131206.rdata")
load("sof_cib_20131206.rdata")
load("nineGSnps_20131031.rdata")
allsnps = nineGSnps

source("C:/Users/JNettiksimmons/Documents/indiv_skatFunctions_20131206.r")
#-------------------------------------------------------------------------

#CR1 rs6701713

#BIN1 rs7561528

#CD2AP rs9349407

#EPHA1 rs11767557

#CLU rs1532278

#MS4A4A rs4938933 #not in our data

#PICALM rs561655

#ABCA7 rs3752246

#CD33 rs3865444


snps = c("rs3764650","rs744373","rs9349407","rs9296559","rs3865444",
		"rs11136000","rs3818361","rs11767557","rs610932","rs670139",
		"rs3851179", 
		"rs6701713", "rs7561528", "rs1532278",
		 "rs561655", "rs3752246" )

snps_m = na.omit(match(snps, allsnps$rsid))




#-------------------------------------------------------------------------
# Get frequency for Mros and sof

mros_snps = match(snps, names(mros_cib_nb))


mros_freq = lapply(snps, function(x){

		m = match(x, names(mros_cgb_nb))
		if(is.na(m)){
			m2 = match(x, names(mros_cib_nb))
			cutgeno = cut(mros_cib_nb[,m2], right = FALSE, breaks = c(0,0.5, 1.5, 2.5))
			tab = table(cutgeno)
			tabsum = sum(tab)
			tab/tabsum

		}else{
			
			tab = table(mros_cgb_nb[,m])
			tabsum = sum(tab)
			tab/tabsum

		}

		

})
names(mros_freq) = snps


sof_freq = lapply(snps, function(x){

		m = match(x, names(sof_cgb))
		if(is.na(m)){
			m2 = match(x, names(sof_cib))
			cutgeno = cut(sof_cib[,m2], right = FALSE, breaks = c(0,0.5, 1.5, 2.5))
			tab = table(cutgeno)
			tabsum = sum(tab)
			tab/tabsum

		}else{
			
			tab = table(sof_cgb[,m])
			tabsum = sum(tab)
			tab/tabsum

		}

		

})

names(sof_freq) = snps


#-------------------------------------------------------------------------
# check actual (not cut) frequency for Mros and sof

mros_snps = match(snps, names(mros_cib))


mros_freqactual = lapply(mros_snps, function(j){

		tab = table(mros_cib[,j])

})
names(mros_freq) = snps


sof_snps = match(snps, names(sof_cib))
sof_freqactual = lapply(sof_snps, function(j){

		tab = table(sof_cib[,j])

})
names(sof_freq) = snps


#--------------------------------------------------------------------------
# double check linear regression


littlelm= function(snp, data){
	subdata = subset(data, select = c('zslope_age', snp, 'EV1', 'EV2', 'EV3', 'EV4'))
	lm.out = lm(zslope_age ~ .,  data = subdata)
	summary(lm.out)
}


mros_reg = lapply(snps, function(x){

		littlelm(x, mros_cib_nb)
	})


sof_reg = lapply(snps, function(x){

		littlelm(x, sof_cib)
	})

