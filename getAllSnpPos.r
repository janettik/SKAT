library('doHPC')
library('ncdf')
registerDoHPC('HOMER',Rbin="C:/R/R-2.15.2/bin/i386/Rterm.exe")


result <- foreach(j=1:1,.packages="ncdf") %dopar% {
nc <- open.ncdf('//HOMER/MROS_SOF/MROS_hapmap_imp/mrosAllChr.nc',readunlim=FALSE)
pos<-get.var.ncdf(nc, "pos", start=1, count=3020515)
name <-get.var.ncdf(nc, "rsid", start=1, count=3020515)
chrom = get.var.ncdf(nc, "chr", start=1, count=3020515)
dat = data.frame(pos, name, chrom)
close.ncdf(nc)

dat
}

snpPos = result[[1]]
save(snpPos, file = "gwasChromSnpPos.rdata")


#check for chromosome specific
#tapply(snpPos$pos, snpPos$chrom, summary, na.rm=T)
#yep