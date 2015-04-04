###################################################################
# Purpose: Functions for individual cohort skat analysis

#----functions/basics--------------------------------------------------

#function to select names of SNPs in specified gene
getsnps = function(gene, drophit, restrictB){
        temp = subset(allsnps, chunk == gene & (ingene == 1 | distbound < restrictB ))
        
        if(drophit == TRUE){
          hit = c(1:nrow(temp))[temp$gwas.hit == 1]
          temp = temp[, -hit]
        }
           
          
          
        temp$rsid
}


#function to select columns of genotyped data based on 
# names of snps. Specify a rarity level for minor allele to be excluded
# ie. exclude snps with < 2% in minor allele homozygous
getData = function(data, snps, rarity){
  
  #keep non-genetic info
  g = grep('rs', names(data))
  nongene = c(1:ncol(data))[-g]
  
  #figure out which snps in gene are in data
  m = match(snps, names(data))
  
  #keep nongene + gene specific genotype/imputed info
  dat = subset(data, select =  na.omit(m))
  
   #recode snps so that 2 is homozygous for the minor allele
  dat2 = apply(dat, 2, function(x){
    
    c = cut(x, c(-.5, .5, 1.5, 2.5))
    tab = table(c)
    if(tab[1] < tab[3]){
      2-x
    }else{
      x
    }
  })
 
  if(is.null(rarity)){
   dat3 =  dat2
  }else{
  #drop snps with less than x% in homozygous 2
    rarityCheck = apply(dat2, 2, function(x){
      c = cut(x, c(-.5, .5, 1.5, 2.5))
      tab = table(c)
      high = tab[3]/sum(tab)
      high > rarity        
    })
    dat3 = dat2[,rarityCheck]
  }
  
   #drop SNPs with more than 25% missing
  dat_check = apply(dat3, 2, function(x){ sum(is.na(x))})
  dat4 = dat3[,dat_check<.25*nrow(dat3)]
  dat5 = cbind(data[,nongene], dat4)
 
  dat5
}


#function to select columns of genotyped data based on 
# names of snps. Specify a rarity level for minor allele to be excluded
# ie. exclude snps with < 2% in minor allele frequency
getDataMAF = function(data, snps, rarity){
  
  #keep non-genetic info
  g = grep('rs', names(data))
  nongene = c(1:ncol(data))[-g]
  
  #figure out which snps in gene are in data
  m = match(snps, names(data))
  
  #keep nongene + gene specific genotype/imputed info
  dat = subset(data, select =  na.omit(m))
  
  #recode snps so that 2 is homozygous for the minor allele
  dat2 = apply(dat, 2, function(x){
    
    c = cut(x, c(-.5, .5, 1.5, 2.5))
    tab = table(c)
    if(tab[1] < tab[3]){
      2-x
    }else{
      x
    }
  })
 
  if(is.null(rarity)){
   dat3 =  dat2
  }else{
  #drop snps with less than x% MAF
    rarityCheck = apply(dat2, 2, function(x){
      sum(x)/(2*nrow(dat2)) >0.02    
    })
    dat3 = dat2[,rarityCheck]
  }
  
   #drop SNPs with more than 25% missing
  dat_check = apply(dat3, 2, function(x){ sum(is.na(x))})
  dat4 = dat3[,dat_check<.25*nrow(dat3)]
  dat5 = cbind(data[,nongene], dat4)
 
  dat5
}


getSnpData = function(gene, data, rarity, drophit = FALSE, restrictB= 100000){
  
   snps = getsnps(gene = gene, drophit = FALSE, restrictB = restrictB)
   use_data = getData(data,  snps = snps, rarity = rarity)
  
}

getSnpDataMAF = function(gene, data, rarity, drophit = FALSE, restrictB= 100000){
  
   snps = getsnps(gene = gene, drophit = FALSE, restrictB = restrictB)
   use_data = getDataMAF(data,  snps = snps, rarity = rarity)
  
}


#Function to actually run skat
# out = outcome
# cov = vector of covariates
runSKATGene = function(data, out, cov, kern="linear.weighted"){
  
  g = grep('rs', names(data))
  
  if(is.null(cov)){
    skatdat = list(y = subset(data, select = out)[,1],
                   Z = as.matrix(subset(data, select = g)))
    obj = with(skatdat,  SKAT_Null_Model(y ~ 1, out_type = "C"))
    
  }else{
    skatdat = list(y = subset(data, select = out)[,1],
                   X = as.matrix(subset(data, select = cov)),
                   Z = as.matrix(subset(data, select = g)))    
    obj = with(skatdat,  SKAT_Null_Model(y ~ X, out_type = "C"))
  }
  
  p = with(skatdat, SKAT(Z, obj, kernel = kern, is_dosage = TRUE))$p.value      
  list(p, kern, out, cov)
}
    
      
runAll = function(data, gene, outcome, drophit = FALSE, restrictBound = 100000, covlist = NULL, rarity = .02){
    gene_snps = getsnps(gene, drophit = drophit, restrictB = restrictBound)
    rel_data = getDataMAF(data, snps = gene_snps, rarity = rarity)
      
    lin = runSKATGene(data = rel_data, out = outcome, cov = covlist, kern = "linear") 
    linw = runSKATGene(data = rel_data, out = outcome, cov = covlist, kern = "linear.weighted")
    ibs = runSKATGene(data = rel_data, out = outcome, cov = covlist, kern = "IBS")
    ibsw = runSKATGene(data = rel_data, out = outcome, cov = covlist, kern = "IBS.weighted")
    x2 = runSKATGene(data = rel_data, out = outcome,  cov = covlist, kern = "2wayIX")
    
    list(lin, linw, ibs, ibsw, x2)
}
    
#############################################################################
# function to format results from runall
        
formatResult = function(reslist){
  
  out =  
  lapply(reslist, function(x){
    do.call('rbind', lapply(x, function(y){
      data.frame(kernel = y[[2]], p = y[[1]], outcome = y[[3]]) 
    }))
  })
      
  out
}    
#############################################################################
# Linear models for individual snps in a given gene
                            
geneLin = function(gene, outcome, data, rarity){
  #get data for specified gene
  snps = getsnps(gene = gene, drophit = FALSE, restrictB = 100000)
  use_data = getDataMAF(data,  snps = snps, rarity = rarity)
  #select outcome
  use_data$out = subset(use_data, select = outcome)[,1]
  #separate snp data
  rs = grep('rs', names(use_data))
  gendata = use_data[,rs]
  
  lm.out = apply(gendata, 2, function(x){
      summary(lm(out ~ x + EV1 + EV2 + EV3 + EV4, data = use_data))
  })
      
  lm.processed = lapply(lm.out, function(x){
    if(nrow(x$coefficient) == 6){
      x$coefficient[2,c(1,4)]
    }else{
      c(NA, NA)
    }
  })
  
    
  lm.processed2 = as.data.frame(do.call('rbind', lm.processed))
  lm.processed2$snp = names(lm.processed)
  adjustP = p.adjust(lm.processed2[,2], method = "fdr")
  lm.processed2$adjustP = adjustP
  list(lm.out, lm.processed2)
}
                            
#############################################################################
# Linear models for snpxsnp interaction individual snps in a given gene
                            
# geneInt = function(gene, outcome, data){
#   #get data for specified gene
#   snps = getsnps(gene = gene, drophit = FALSE, restrictB = 100000)
#   use_data = getData(data,  snps = snps)
#   #select outcome
#   use_data$out = subset(use_data, select = outcome)[,1]
#   #separate snp data
#   rs = grep('rs', names(use_data))
#   gendata = use_data[,rs]
#   
#   x2 = do.call('rbind', lapply(1:ncol(gendata), function(i){
#     
#         temp = lapply(1:ncol(gendata), function(j){
#           
#           s = summary(lm(use_data$out ~ gendata[,i]*gendata[,j]))
#           s$coefficient
#         })
#         
#         temp2 = unlist(lapply(temp, function(x){
#               nr = nrow(x)
#               if(nr >3){
#                 x[4,4]
#               }else{
#                 NA
#               }
#         }))
#   }))
#   
#  dimnames(x2)[[1]] = names(gendata)
#  dimnames(x2)[[2]] = names(gendata)
#   x2
# }

