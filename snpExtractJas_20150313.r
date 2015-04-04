#SNP extract from netCDF files
#parallel SNP extraction using ncdf and foreach

library('doHPC')
library('ncdf')
registerDoHPC('HOMER',Rbin="C:/R/R-2.15.2/bin/i386/Rterm.exe")
setwd("C:/Users/JNettiksimmons/Documents/DATA")

#--------------------------------------------------------------------------------
# Jasmine's SNPs with index information
#--------------------------------------------------------------------------------
load("C:/Users/JNettiksimmons/Documents/DATA/snpedia_20150313.rdata")
snpSelect = snpedia

#--------------------------------------------------------------------------------
# function to extract multiple SNPs on one node


getGeno = function(snpIndices, impute = FALSE, mycohort = 'both'){
  result <- foreach(j=1:1,.packages="ncdf") %dopar% {
    
    #open up the right file depending on cohort and impute/direct
    #get the number of subjects/samples 
    if(impute == TRUE & mycohort == "MrOS"){
      nc <- open.ncdf('//HOMER/MROS_SOF/MROS_hapmap_imp/mrosAllChr.nc',readunlim=FALSE)
      numsub =  nc$dim$subject$len
    }else if (impute == TRUE & mycohort == "SOF"){
      nc <- open.ncdf('//HOMER/MROS_SOF/SOF_hapmap_imp/sofAllChr.nc',readunlim=FALSE)
      numsub =  nc$dim$subject$len
    }else if(impute == FALSE){
      nc <- open.ncdf('//HOMER/MROS_SOF/MROS_SOF_201109/NetCDF_files_Sept2011/MrOS_SOF_SW.OMNI.GENO.09062011.nc', 
                                                                    readunlim = FALSE)
      numsub =  nc$dim$sample$len
    }
    
    #make blank matrix to fill
    mydataset<-matrix(NaN,nrow=numsub,ncol=length(snpIndices))
    
    #pull out each snp in requested list
    if(impute == TRUE){
        for(i in 1:length(snpIndices)){
          mydat<-get.var.ncdf(nc, "genotypes", start=c(1,snpIndices[i]), count=c(-1,1))
          mydataset[,i]<-mydat
        }
        close.ncdf(nc)
      return(mydataset) 
    }else{
        for(i in 1:length(snpIndices)){
          mydat<-get.var.ncdf(nc, "genotype", start=c(snpIndices[i],1), count=c(1,-1))
          mydataset[,i]<-mydat
        }
        close.ncdf(nc)
      return(mydataset)
    }
  }
}


#-----------------------------------------------------------------------------
#           clean up function
#-----------------------------------------------------------------------------

cleanUp= function(mydat, mycohort, mypop, impute = FALSE){
  
  #read in annotation data
  samp0<-read.table("Z:/MROS_SOF_201109/samples_related/sampleAnnotDan2011-10-10.txt",
                                            header=T,stringsAsFactors=F,sep="\t")
  #process differently depending on whether imputed or not
  if(impute == FALSE){
    if(mycohort=="MrOS"){
      pc<-read.table("Z:/MROS_SOF_201109/PCA_EVs_files/MrOS_Whites_PCs.txt",
                                              header=T,stringsAsFactors=F)
    }else if(mycohort=="SOF" ) {
      pc<-read.table("Z:/MROS_SOF_201109/PCA_EVs_files/SOF_Whites_PCs.txt",
                                              header=T,stringsAsFactors=F)
    }
                     
    pc<-pc[,c("individual.id","EV1","EV2","EV3","EV4")]
    samp<-merge(samp0,pc,by="individual.id",all.x=T)
    samp<-samp[order(samp$sample.num),]
    mydat$id = samp$individual.id
    mydat2 = cbind(mydat, samp)
    rm(pc)
    gc()
    
  }else if(impute == TRUE & mycohort == "MrOS"){
    pc<-read.table("//HOMER/MROS_SOF/MROS_SOF_201109/PCA_EVs_files/MrOS_Whites_PCs.txt",
                                                                header=T,stringsAsFactors=F)
    gwasID<-read.csv("//HOMER/MROS_SOF/MROS_hapmap_imp/ids.csv",header=T,stringsAsFactors=F)
    gwasID = gwasID[order(gwasID$machIndex),]
    samp<-samp0
    samp<-samp[samp$worseDup==0,]
    samp$siteFact<-factor(samp$site)
    samp<-merge(gwasID,samp,by.x="ID",by.y="individual.id",all.x=T)
    mydat$id = gwasID$ID
    mydat$machIndex = gwasID$machIndex
 
    pc<-pc[,c("individual.id","EV1","EV2","EV3","EV4")]
    samp<-merge(samp,pc,by.x="ID", by.y="individual.id",all.x=T)
    samp<-samp[order(samp$machIndex),]
    mydat2 = cbind(mydat, samp)
    rm(pc)
    gc()
    
  }else if(impute == TRUE & mycohort == "SOF"){
    pc<-read.table("//HOMER/MROS_SOF/MROS_SOF_201109/PCA_EVs_files/SOF_Whites_PCs.txt",
                                                      header=T,stringsAsFactors=F)
    gwasID<-read.csv("//HOMER/MROS_SOF/SOF_hapmap_imp/ids.csv",header=T,stringsAsFactors=F)
    gwasID = gwasID[order(gwasID$machIndex),]
    mydat$id = gwasID$ID
    mydat$machIndex = gwasID$machIndex
    
    
    samp<-samp0
    samp<-samp[samp$worseDup==0,]
    samp$siteFact<-factor(samp$site)
    samp<-merge(gwasID,samp,by.x="ID",by.y="individual.id",all.x=T)
  
    pc<-pc[,c("individual.id","EV1","EV2","EV3","EV4")]
    samp<-merge(samp,pc,by.x="ID", by.y="individual.id",all.x=T)
    samp<-samp[order(samp$machIndex),]
    mydat2 = cbind(mydat, samp)
    
    rm(pc)
    gc()
  }
  
  #replace -1 with NA
  mydat2[] = lapply(mydat2, function(x){ replace(x, x== -1, NA)})
  #select correct group of individuals
  mydat3<-subset(mydat2, worseDup==0 & collection==mycohort & Ethnicity==mypop)

}

############################################################################
#                           run
############################################################################
#get rid of NA indexes in directly genotyped and imputed snp lists
snpSelectGeno = subset(snpSelect, is.na(index) == FALSE)
snpSelectImp= subset(snpSelect, is.na(index_imp) == FALSE)

#get directly genotyped snps from both studies combined
# both_geno0 = getGeno(snpSelectGeno$index)
# both_geno<-as.data.frame(both_geno0[[1]]) #reassign first element of list to new object names(both_geno) = snpSelectGeno$rsid
# save(both_geno, file = "both_geno_20150313.rdata")
# 
#get imputed snps from MrOS
mros_imp0 = getGeno(snpSelectImp$index_imp, impute = TRUE, mycohort = 'MrOS')
mros_imp<-as.data.frame(mros_imp0[[1]]) 
names(mros_imp) = snpSelectImp$rsid
save(mros_imp, file = "mros_imp_20150313.rdata")

#get imputed snps from SOF
sof_imp0 = getGeno(snpSelectImp$index_imp, impute = TRUE, mycohort = 'SOF')
sof_imp<-as.data.frame(sof_imp0[[1]]) #reassign first element of list to new object
names(sof_imp) = snpSelectImp$rsid
save(sof_imp, file = "sof_imp_20150313.rdata")

#clean and split up directly genotyped snps
# mros_geno_clean = cleanUp(both_geno, mycohort = "MrOS", mypop = "White")
# save(mros_geno_clean, file = "mros_genoC_20150313.rdata")
# 
# sof_geno_clean = cleanUp(both_geno, mycohort = "SOF", mypop = "White")
# save(sof_geno_clean, file = "sof_genoC_20150313.rdata")

#clean imputed snps
mros_imp_clean = cleanUp(mros_imp, mycohort = "MrOS", mypop = 'White', impute = TRUE)
save(mros_imp_clean, file = "mros_impC_20150313.rdata")

sof_imp_clean = cleanUp(sof_imp, mycohort = "SOF", mypop = 'White', impute = TRUE)
save(sof_imp_clean, file = "sof_impC_20150313.rdata")



