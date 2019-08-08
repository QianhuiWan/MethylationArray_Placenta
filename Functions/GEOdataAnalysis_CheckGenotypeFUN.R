
# The function for using ewastools to check sex for placenta samples:
# let argument RGsetChara to be a charactor

CheckGenotype <- function(RGsetChara){
  registerDoParallel(cores=10)
  
  # get the object from character
  RGset <- get(RGsetChara)
  
  # get meta data
  phenoData <- pData(RGset)
  # phenoData to phenoData.dt
  phenoData.dt <- phenoData %>% as.data.frame() %>% rownames_to_column() %>% setDT()
  
  # read idat files and get beta
  sampleDirs <- phenoData$filenames
  meth = read_idats(idat_files = sampleDirs,quiet=FALSE)
  beta = dont_normalize(meth)
  
  # get detection P values
  detP <- minfi::detectionP(RGset)
  
  # get snps
  snps = meth$manifest[probe_type=='rs',index]
  snps = beta[snps,]
  genotypes = call_genotypes(snps,learn=TRUE)
  
  if(ncol(RGset)==1){
    if (!"outliers" %in% names(genotypes)) 
      stop("Invalid argument")
    log_odds = genotypes$outliers/(1 - genotypes$outliers)
    log_odds = mean(log2(log_odds), na.rm = TRUE)
    # log_odds
    phenoData.dt$log_odds = log_odds
  }else{
    
    phenoData.dt$log_odds = snp_outliers(genotypes)
  }
  
  table(phenoData.dt$log_odds > -4)
  
  phenoData.dt[phenoData.dt$log_odds > -4, ] %>% arrange(desc(log_odds))
  
  # assumed distribution of SNP-probe intensities. 
  # if probes are outoff distribution, there could be some technical variance (could to be used to indicate tissue contanmination)
  mxm_(genotypes)
  
  # plot genotype checking results
  par(mfrow=c(2,1))
  barplot(colMeans(detP), col=pal[6], las=2,
          cex.names=0.8, ylim=c(0,0.01),
          xlab = "",
          ylab = "Mean detection p-values",
          axisnames = FALSE)
  abline(h=0.01,col="red")
  
  barplot(phenoData.dt$log_odds,
          xlab = "placental tissue samples(n=53)",
          ylab = "Average log odds")
  abline(h=-4,col="red")
  
  rm(meth)
  
  # return phenoData_withPredict, the df with predicted fetal sex
  return(list(detP=detP, phenoData.dt=as_tibble(phenoData.dt)))
}

