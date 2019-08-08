
# The function for using ewastools to check sex for placenta samples:
# let argument RGsetChara to be a charactor

CheckFetalSex <- function(RGsetChara){
  registerDoParallel(cores=10)
  # get the object from character
  RGset <- get(RGsetChara)
  # get meta data
  phenoData <- pData(RGset)
  
  # read idat files
  sampleDirs <- phenoData$filenames
  meth = read_idats(idat_files = sampleDirs,quiet=FALSE)
  
  # sex check
  ## dye bias corrected before sex check
  meth %<>% correct_dye_bias
  
  # phenoData to phenoData.dt
  phenoData.dt <- phenoData %>% as.data.frame() %>% rownames_to_column() %>% setDT()
  # add average total beta values for probes on X and Y chr
  phenoData.dt[ ,c("X","Y") := check_sex(meth)]
  
  # predict fetal sex
  phenoData.dt$predicted_sex = predict_sex(phenoData.dt$X,phenoData.dt$Y
                                           # which(phenoData.dt$Fetal_Sex%in%'Male'),
                                           # which(phenoData.dt$Fetal_Sex%in%'Female')
  )
  
  phenoData.dt_v1 <- phenoData.dt %>%
    mutate(predicted_sex = str_replace_all(predicted_sex, c("m"="Male","f" ="Female")))
  
  phenoData.dt_v2 <- setDT(phenoData.dt_v1)
  
  # plot the X/Y chr intensities
  # it is better to put pch into 3 groups (red color marked the misslabelled samples)
  plot(Y ~ X,data=phenoData.dt_v2,
       pch=ifelse(phenoData.dt_v2$Fetal_Sex=="Female",1,
                  ifelse(phenoData.dt_v2$Fetal_Sex=="Male", 4,13)),
       asp=1,xlab="Normalized X chromosome intensities",ylab="Normalized Y chromosome intensities")
  points(Y ~ X,data=phenoData.dt_v2[Fetal_Sex !=predicted_sex],
         pch=ifelse(phenoData.dt_v2[Fetal_Sex!=predicted_sex]$Fetal_Sex=="Female",1,
                    ifelse(phenoData.dt_v2[Fetal_Sex!=predicted_sex]$Fetal_Sex=="Male",4,13)),
         col=2)
  legend("topright",pch=c(1,4,13),legend=c("Female","Male", "NA"))
  title(main = RGsetChara)
  
  # samples with fetal sex in meta data not same as predicted fetal sex
  phenoData.dt_v2[Fetal_Sex !=predicted_sex] %>% as_tibble %>% 
    dplyr::select(rowname:Study, X:predicted_sex) %>% print(n=Inf)
  
  # remove intermediat data
  rm(meth)
  
  # return phenoData_withPredict, the df with predicted fetal sex
  return(as_tibble(phenoData.dt_v2))
}