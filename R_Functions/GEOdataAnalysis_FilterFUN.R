
# This function is used to filter failed probes from data sets (408 samples)
# TissueTypes: Amnoin, Chorion, Decidua, maternal whole blood, and umbilical cord blood.

# the parameter: ArrayType: can be 450K or EPIC


# get extended RG channel sets===========================================================


getRGset_extended <- function(ArrayType){
  # base directory
  GEO_baseDir <- file.path("/home/qianhui/DNAme/Process_decidual/GEO_dataSets/IDAT_all")
  
  if(ArrayType=="450K"){
    # get 450k meta data
    kept450k <- GEO_phenotypes$ArrayType=='450K'& 
      !(GEO_phenotypes$Sample_name%in%c("GSM1617002_7668610068_R01C01","GSM1616993_7668610068_R06C02"))
    GEO_phenotypes_sub <- GEO_phenotypes[kept450k,]
    
  }else if(ArrayType=="EPIC"){
    # get EPIC meta data
    GEO_phenotypes_sub <- GEO_phenotypes[GEO_phenotypes$ArrayType=='EPIC',]
  }
  
  # get the RGset for this tissue
  GEO_RGset_extended <- read.metharray.exp(base = GEO_baseDir, 
                                           targets = GEO_phenotypes_sub, 
                                           extended = TRUE, force = TRUE)
  
  return(GEO_RGset_extended)
  
}



# filter probes =========================================================================
# the parameters:
# RGsetChara: "RGset" that you are interested in;
# crossreactiveProbes: data frame that containing crossreactiveProbes from 450K or EPIC array;
# annotation: ann450k or annEPIC, a DataFrame
# DropXreactive: drop crossreactive probes, FALSE or TRUE;
# dropSNP: dropSNP=TRUE means drop all probes related with SNPs;
# dropXY: dropXY=TRUE means drop all probes on X or Y chromosomes.


FilterFun <- function(RGsetChara, crossreactiveProbes, annotation, 
                      DropXreactive, DropSNP, DropXY){
  
  registerDoParallel(cores = 10)
  # get the object from character
  RGset <- RGsetChara
  
  # get methylset
  methylset_raw <- preprocessRaw(RGset) # 866238
  
  # filtering of probes
  ## drop P
  detP <- minfi::detectionP(RGset)
  Failed <- detP > 0.01
  methylset_P <- methylset_raw[rowSums(Failed)==0,] #846991
  
  ## drop probes with <3 beads
  ## If a probe has less than 3 beads in 95% of the samples, we need to remove it because the intensity of this probe is not reliable.
  
  # if(identical(annotation, ann450k)){
  #   
  # EpicRGset_extended <- getRGset_extended(ArrayType = "450K")
  # 
  # }else if(identical(annotation, annEPIC)){
  #   
  # EpicRGset_extended <- getRGset_extended(ArrayType = "EPIC")
  # }
  
  EpicRGset_extended <- RGset
  BeadCount <- getNBeads(EpicRGset_extended)
  RemainProbe <- rowSums(BeadCount<3) < 0.05 * (ncol(BeadCount))
  EpicRGset_extended_Bead <- EpicRGset_extended[RemainProbe,]
  methylset_rawBead <- preprocessRaw(EpicRGset_extended_Bead) # bead number filter
  
  methylset_PBead <- methylset_P[rownames(methylset_P)%in%rownames(methylset_rawBead),] # 834351
  
  ## drop cross-reactive probes
  if(DropXreactive==TRUE){
    keep <- !(featureNames(methylset_PBead) %in% crossreactiveProbes$ProbeNames)
    methylset_PBeadX <- methylset_PBead[keep, ] # 792070
  }else{
    methylset_PBeadX <- methylset_PBead
  }
  
  ## drop SNPs
  ### if I use DropSNP argument, I drop all probes with SNPs
  if(DropSNP==TRUE){
    GRaSet_meth_PBeadX <- methylset_PBeadX %>% ratioConvert() %>% mapToGenome()
    GRaSet_meth_PBeadXAddSNPinfo <- addSnpInfo(GRaSet_meth_PBeadX)
    GRaSet_meth_PBeadXsnp <- dropLociWithSnps(GRaSet_meth_PBeadXAddSNPinfo,
                                              snps = c("CpG", "SBE", "Probe"), maf = 0)
    
    meth_PBeadXsnp <- methylset_PBeadX[rownames(methylset_PBeadX)%in% rownames(GRaSet_meth_PBeadXsnp), ]
    
    
  }else{
    ### else, I only drop unwanted SNP related probes, that is SNPs at CpG and SBE sites
    GRaSet_meth_PBeadX <- methylset_PBeadX %>% ratioConvert() %>% mapToGenome()
    GRaSet_meth_PBeadXAddSNPinfo <- addSnpInfo(GRaSet_meth_PBeadX)
    GRaSet_meth_PBeadXsnp <- dropLociWithSnps(GRaSet_meth_PBeadXAddSNPinfo, 
                                              snps = c("CpG", "SBE"), maf = 0)
    
    meth_PBeadXsnp <- methylset_PBeadX[rownames(methylset_PBeadX)%in% rownames(GRaSet_meth_PBeadXsnp), ]
    
  }
  
  ## drop probes on X, Y chr
  if(DropXY==TRUE){
    
    keep <- !(featureNames(meth_PBeadXsnp) %in% annotation$Name[annotation$chr %in% c("chrX", "chrY")])
    meth_PBeadXsnpXY <- meth_PBeadXsnp[keep, ] # 600k
    
  }else{
    
    meth_PBeadXsnpXY <- meth_PBeadXsnp
  }
  
  # the filtered CpGs
  # filteredCpGs <- setdiff(featureNames(methylset_raw), featureNames(meth_PBeadXsnpXY))
  
  # remove files not need
  rm(methylset_raw, EpicRGset_extended,     
     GRaSet_meth_PBeadX, GRaSet_meth_PBeadXAddSNPinfo, GRaSet_meth_PBeadXsnp)
  
  return(meth_PBeadXsnpXY)
}
