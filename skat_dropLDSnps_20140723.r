##########################################################################
# Purpose: do skat on skat hits dropping ld sections
# Date:06/19/2014
##########################################################################
library('doHPC')
registerDoHPC('HOMER',Rbin="C:/R/R-2.15.2/bin/i386/Rterm.exe")

setwd("C:/Users/JNettiksimmons/Documents/DATA")
load("mros_cib_noBirm_20131122.rdata")
load("sof_cib_20131206.rdata")
load("nineGSnps_20131031.rdata")
allsnps =  nineGSnps

source("C:/Users/JNettiksimmons/Documents/indiv_skatFunctions_20140722.r")
blocks = read.csv("C:/Users/JNettiksimmons/Documents/ldBlocksData.csv", header = F)
#-----------------------------------------------------------------------------

# MROS PICALM
# using R^2 = 0.8 as threshhold

picalmLD = as.character(blocks[blocks[,2] == 'picalm',1])

abca7LD = as.character(blocks[blocks[,2] == 'abca7',1])


m_mros = match(c(picalmLD, abca7LD), names(mros_cib_nb))
mros_newcib = mros_cib_nb[, -na.omit(m_mros)]

is_inm = is.na(m_mros) == FALSE
mros_gene = c(rep('picalm', length(picalmLD)), rep('abca7', length(abca7LD)))
table(is_inm, mros_gene)


#-------------------------------------------------------------------------
# SOF

# using R^2 = 0.8 as threshhold


cd33LD = as.character(blocks[blocks[,2] == 'cd33',1])
bin1LD = as.character(blocks[blocks[,2] == 'bin1',1])
cr1LD = as.character(blocks[blocks[,2] == 'cr1',1])

m_sof = match(c(cd33LD, bin1LD, cr1LD), names(sof_cib))
sof_newcib = sof_cib[, -na.omit(m_sof)]

is_in = is.na(m_sof) == FALSE
sof_gene = c(rep('cd33', length(cd33LD)), rep('bin1', length(bin1LD)), 
              rep('cr1', length(cr1LD)))
table(is_in, sof_gene)

#-----------run and save-------------------------


picalm_out=  unlist(foreach(j=1:1,.packages="SKAT") %dopar% {  
  runAll(mros_newcib, gene = 'picalm', outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
}, recursive = FALSE)  

abca7_out=  unlist(foreach(j=1:1,.packages="SKAT") %dopar% {  
  runAll(mros_newcib, gene = 'abca7', outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
}, recursive = FALSE) 
 

cd33_out =  unlist(foreach(j=1:1,.packages="SKAT") %dopar% {  
  runAll(sof_newcib, gene = 'cd33', outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
}, recursive = FALSE)

cr1_out =  unlist(foreach(j=1:1,.packages="SKAT") %dopar% {  
  runAll(sof_newcib, gene = 'cr1', outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
}, recursive = FALSE)

bin1_out =  unlist(foreach(j=1:1,.packages="SKAT") %dopar% {  
  runAll(sof_newcib, gene = 'bin1', outcome = 'zslope_age', covlist = c('EV1', 'EV2', 'EV3', 'EV4'))
}, recursive = FALSE)



dropLD = list(picalm_out, abca7_out, cd33_out, bin1_out, cr1_out)
names(dropLD) = c('picalm', 'abca7', 'cd33', 'bin1', 'cr1')

save(dropLD, file = "C:/Users/JNettiksimmons/Documents/output/sofMros_dropLD20140819.rdata")

