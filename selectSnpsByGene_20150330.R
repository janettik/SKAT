#Purpose: Get subset of SNPs by gene
# Date: 10/25/2013


setwd("C:/Users/JNettiksimmons/Documents/code_formatData")
snpTab = read.csv('newSnps_20150330.csv', header = T)

setwd("C:/Users/JNettiksimmons/Documents/DATA")
load('gwasChromSnpPos.rdata')


snpPos$order = 1:nrow(snpPos)
###############################################################
# add chunk name to snp table

snpTab$chunk = tolower(snpTab$gene)
snpTab$chunk[grep('hla', snpTab$chunk)] = 'hla'

###############################################################
# Define Function
#############################################################

pullsnps = function( geneLow, geneHigh, upstream,
                     downstream = upstream,
                    chunkName, chr, hit = NULL){
  #get expanded range
  low = geneLow - downstream*1000
  high = geneHigh + upstream*1000
  
  #subset down to required range
  temp2 = subset(snpPos, chrom == chr & pos > low & pos < high)
  
  #see what is actually in gene boundaries
  ingene = c(1:nrow(temp2))[temp2$pos > geneLow & temp2$pos <geneHigh]
  temp2$ingene= 0
  temp2$ingene[ingene] = 1
  
  #mark gwas hit, if necessary
  temp2$gwas.hit = 0
  if(is.null(hit) == FALSE){
    hit2 = as.numeric(gsub('rs', '', hit))
    m = match(hit2, temp2$name )
    temp2$gwas.hit[m] = 1
  }
  
  #add additional info
  temp2$rsid = paste('rs', temp2$name, sep = '')
  temp2$chunk = chunkName
  temp2$geneLBound = geneLow
  temp2$geneUBound = geneHigh
  bounddistlow = with(temp2, geneLow - pos)
  bounddisthigh = with(temp2, geneHigh - pos)
  distbound = sapply(1:nrow(bounddistlow), function(j){
    
            if(temp2$ingene[j] == 1){
              NA
            }else if(bounddisthigh[j] < 0){
              abs(bounddisthigh[j])
            }else if(bounddistlow[j] >0){
              abs(bounddistlow[j])
            }
  })
  temp2$distbound = distbound
    
  temp2
}

#-------------------------------------------------------------
# run function

geneSnps = 
lapply(1:nrow(snpTab), function(j){
	
	pullsnps(geneLow = snpTab$start[j], geneHigh = snpTab$end[j],
		downstream = snpTab$downstream[j], upstream = snpTab$upstream[j],
		chr = snpTab$chrom[j], hit = snpTab$snp[j],
		chunkName = snpTab$chunk[j])

})


names(geneSnps) = snpTab$chunk


#combine
snpedia = do.call('rbind', geneSnps)

#-----------------------------------------------------------------------------
# add index information

#first, use 64bit, read huge SNP annot file, get index numbers for SNPs
snpAnnot<-read.table("Z:/MROS_SOF_201109/SNPs_related/SNP_ANNOT_09132011.txt",
                      header=T,sep="\t",stringsAsFactors=F,comment.char="")
#imputed SNP list file
snpImp<-read.csv(file="Z:/MROS_hapmap_imp/imputeSnpList.csv",
                      header=T,stringsAsFactors=F)

#find ncdf index associated with my SNPs
snpedia$index = snpAnnot$int.id[match(snpedia$rsid,snpAnnot$rs.id_b36)]
snpedia$index_imp = c(1:nrow(snpImp))[match(snpedia$rsid, snpImp[,1])]

#---------------------------------------------------------------------------
# output

save(snpedia, file = "snpedia_20150330.rdata")
