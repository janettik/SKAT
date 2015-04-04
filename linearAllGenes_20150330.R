##########################################################################
# Purpose: SOF and MrOS SNPs linear analysis all genes
# Date: 11/14/2013
##########################################################################
library('doHPC')
registerDoHPC('HOMER',Rbin="C:/R/R-2.15.2/bin/i386/Rterm.exe")

setwd("C:/Users/JNettiksimmons/Documents/DATA")
load("mros_cib_noBirm_20150330.rdata")
load("sof_cib_20150330.rdata")
load("snpedia_20150330.rdata")
allsnps = snpedia

source("C:/Users/JNettiksimmons/Documents/indiv_skatFunctions_20140722.r")
#########################################################################

geneList = c('abca7', 'picalm', 'clu', 'cr1' , 'ms4a6a', 'bin1', 'cd33', 
             'epha1', 'cd2ap', 'hla', 'ptk2b', 'sorl1', 'slc24a4',
			 'inpp5d', 'mef2c', 'nme8', 'zcwpw1', 'celf1', 'fermt2', 'cass4', 'ms4a6e')

getSigLin = function(data){
 do.call('rbind', lapply(geneList, function(g){
  
   lin = geneLin(gene = g, outcome = 'zslope_age', data = data, rarity = 0.02)
   temp_dat = getSnpDataMAF(gene = g, data = data, rarity = 0.02)
   rsCol = grep('rs', names(temp_dat))
   rs_dat = subset(temp_dat, select = rsCol)
   rs_dat2 = apply(rs_dat, 2, function(x){
          crs = cut(x, c(-.5, .5, 1.5, 2.5))
          levels(crs) = c(0,1,2)
          crs
   })
   rs_tab = t(apply(rs_dat2, 2, function(x){
            x2 = factor(x, levels = c(0, 1, 2))
            t = table(x2)
            t/sum(t)
   }))

   gtp_maf = (rs_tab[,2]+2*rs_tab[,3])/2
   orig_maf = apply(rs_dat, 2, function(x){
			sum(x)/nrow(rs_dat)/2
	})
   
   D = cbind(lin[[2]], round(rs_tab,3), gtp_maf, orig_maf)
   names(D) = c('est', 'p', 'snp', 'adjustP', '0', '1', '2', 'gtp_maf', 'imp_maf')
   D$gene = g
   D$impTot = ncol(rs_dat)
   D
  
}))
}

mros = getSigLin(data = mros_cib_nb)
sof = getSigLin(data = sof_cib)

m_sof = match(dimnames(mros)[[1]], dimnames(sof)[[1]])

comb = cbind(mros, sof[m_sof,])

names(comb) = c(paste('m', names(mros), sep = '_'), paste('s', names(sof), sep = '_'))

save(comb, file = "C:/Users/JNettiksimmons/Documents/output/completeLinear_20150330.rdata")


comb2 = do.call('rbind', lapply(1:nrow(comb), function(j){
  
     	temp = comb[j,]

	if(is.na(temp$m_p) == FALSE & is.na(temp$s_p) == FALSE){
    
        	if(temp$m_p <0.05 | temp$s_p <0.05){
          		temp
		}
        	
	}else if(is.na(temp$m_p)){
		if(temp$s_p<0.05){
			temp
		}
	}else if(is.na(temp$s_p)){
		if(temp$m_p<0.05){
			temp
		}
	}
  
}))



write.csv(comb2, file = "C:/Users/JNettiksimmons/Documents/output/sofMros_allSigLinear_20150330.csv")

#------------------------------------------------------------------------
# look specifically at GWAS snps
----------------------------------

hits = subset(allsnps, gwas.hit == 1)
hitnames = hits$rsid
#--------------------------------------------------------------
checkAssoc = function(cohort, hit){
  
  if(cohort == 'mros'){
      use_data = mros_cib_nb
    
  }else if(cohort == 'sof'){
      use_data = sof_cib
  }
  use_data$snp = subset(use_data, select = hit)[,1]
  lm.out = lm(zslope_age ~ snp + EV1 + EV2 + EV3 + EV4, data = use_data)
  summary(lm.out)
}

mros_out = lapply(hitnames, function(x){
            modcoef = checkAssoc(cohort = 'mros', hit = x)
            if(dim(modcoef$coefficients)[1] == 6){
              modcoef$coefficients[2,]
            }
})
            
names(mros_out) = hitnames

sof_out = lapply(hitnames, function(x){
            modcoef = checkAssoc(cohort = 'sof', hit = x)
            if(dim(modcoef$coefficients)[1] == 6){
              modcoef$coefficients[2,]
            }
})
names(sof_out) = hitnames

save(mros_out, file = "C:/Users/JNettiksimmons/Documents/output/mrosNB_hits_20150330.rdata")
save(sof_out,file = "C:/Users/JNettiksimmons/Documents/output/sof_hits_20150330.rdata") 


#-------------------------------------------------------------------
# get snp count with minor allele homozygous threshold of 2%


getCountMrOS =
lapply(geneList, function(g){
  
   temp_dat = getSnpData(gene = g, data = mros_cib_nb, rarity = .02)
	length(grep("^rs", names(temp_dat)))
   
  
})

getCountSOF =
lapply(geneList, function(g){
  
   temp_dat = getSnpData(gene = g, data = sof_cib, rarity = .02)
	length(grep("^rs", names(temp_dat)))
   
  
})

count = as.data.frame(cbind(getCountMrOS, getCountSOF))
count$names = geneList


getfreq = function(data){
lapply(geneList, function(g){
temp_dat = getSnpData(gene = g, data = mros_cib_nb, rarity = .02)
  rsCol = grep('rs', names(temp_dat))
   rs_dat = subset(temp_dat, select = rsCol)
   rs_dat2 = apply(rs_dat, 2, function(x){
          crs = cut(x, c(-.5, .5, 1.5, 2.5))
          levels(crs) = c(0,1,2)
          crs
})
}
 
#-------------------------------------------------------------------
# get snp count with MAF threshold of 2%


getCountMrOS =
lapply(geneList, function(g){
  
   temp_dat = getSnpDataMAF(gene = g, data = mros_cib_nb, rarity = .02)
	length(grep("^rs", names(temp_dat)))
   
  
})

getCountSOF =
lapply(geneList, function(g){
  
   temp_dat = getSnpDataMAF(gene = g, data = sof_cib, rarity = .02)
	length(grep("^rs", names(temp_dat)))
   
  
})

count = as.data.frame(cbind(getCountMrOS, getCountSOF))
count$names = geneList
count2 = as.matrix(count)
       
write.csv(count2, file = "C:/Users/JNettiksimmons/Documents/output/sofMrosSnpCount_20150330.csv")
