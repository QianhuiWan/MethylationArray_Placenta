
# This function is used to manipulate meta data from each GEO dataset we used in technical report_placenta
# the parameter: GSEnum can be a GSE number (e.g. GSE66210) or "all", all means manipulating meta data from all datasets.

TrimGEOmetaDataFUN <- function(GSEnum){
  
  registerDoParallel(cores = 10)
  load(file = '/home/qianhui/DNAme/Process_decidual/RDSfiles/pdataOrMetadata_GEO_15datasets.Rdata')
  
  if(GSEnum%in%c("GSE66210", "all")){
    # establish a proper samplesheet for this dataset
    GSE66210.meta_v1 <- GSE66210.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample=str_replace_all(`tissue:ch1`, c('maternal whole blood'= 'Maternal whole blood', 'chorionic villus'='Placenta'))) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`pregnancy:ch1`, 'normal pregnancy', 'Uncomplicated')) %>% 
      inset('Gestation', value=rep('NA', nrow(GSE66210.meta))) %>% 
      inset('Trimester', value=rep('First', nrow(GSE66210.meta))) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE66210.meta))) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE66210.meta))) %>% 
      mutate(Fetal_Sex=str_replace_all(`individual:ch1`, c('female fetus'='Female', 'male fetus'='Male'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE66210.meta))) %>% 
      inset('Study', value=rep('GSE66210', nrow(GSE66210.meta))) %>% 
      mutate('Array'=str_replace_all(description.1, "Sa.+_", "")) %>% 
      mutate('Slide'=str_replace_all(description.1, c("Sa.+ "="", "_.+"=""))) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE66210.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE66210/GSE66210.meta_v1.csv")
    
  }
  
  ## `GEO_num[2]` from GEO database. 
  if(GSEnum%in%c("GSE74738", "all")){
    # establish a proper samplesheet for 
    GA <- GSE74738.meta$`fetal gestational age (weeks)/trimester:ch1`
    GA_trimester <- str_replace_all(GA, c('na'='NA', '1st_trimester'='First', '2nd_trimester'='Second', '3rd_trimester'='Third',
                                          '<10w'='10', '12w 35-38d'='12', '15w by dates 9w size 53d'='15',
                                          '8w 37d DA'='8', '8w 48 days DA'='8', '9w 26-30d DA'='9'))
    FSnames <- GSE74738.meta$`fetal sex/sex:ch1` %>% table() %>% names()
    GSE74738.meta$`fetal sex/sex:ch1`[GSE74738.meta$`fetal sex/sex:ch1`==FSnames[7]] <- "M"
    
    GSE74738.meta_v1 <- GSE74738.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample=str_replace_all(`sample tissue:ch1`, 
                                    c('amnioniotic membrane'='Amnion',
                                      'amniotic membrane'='Amnion',
                                      'chorionic membrane'='Chorion',
                                      'cord blood'='Umbilical cord blood', 'placental chorionic villi'='Placenta', 'placental decidua'='Decidua', 
                                      'whole maternal blood'='Maternal whole blood'))) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`status/group:ch1`, c('Ab. MSS '= 'Ab.MSS', 'Chromosomically Normal Miscarriage'='Miscarriage',
                                                             'control'='Uncomplicated'))) %>% 
      mutate(Gestation= str_replace_all(GA, c('na'='NA', '1st_trimester'='NA', '2nd_trimester'='NA', '3rd_trimester'='NA',
                                              '<10w'='10', '12w 35-38d'='12', '15w by dates 9w size 53d'='15',
                                              '8w 37d DA'='8', '8w 48 days DA'='8', '9w 26-30d DA'='9'))) %>% 
      inset('Trimester', value=ifelse(GA_trimester=='NA', 'NA', 
                                      ifelse(GA_trimester=='First', 'First', 
                                             ifelse(GA_trimester=='Second', 'Second', 
                                                    ifelse(GA_trimester=='Third', 'Third', 
                                                           ifelse(as.numeric(GA_trimester)>=37, 'Term', 'Third')))))) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE74738.meta))) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE74738.meta))) %>% 
      mutate(Fetal_Sex= str_replace_all(`fetal sex/sex:ch1`, c('F'='Female', 'M'='Male', ' '='', 'pool_na'='NA'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE74738.meta))) %>% 
      inset('Study', value=rep('GSE74738', nrow(GSE74738.meta))) %>% 
      mutate('Array'= `450k sentrix_position:ch1`) %>% 
      mutate('Slide'= `450k sentrix_id:ch1`) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    ControlSamples <- GSE74738.meta_v1 %>% filter(str_detect(Sample, "Placenta") & Disease=='Uncomplicated')
    ControlSamples$Gestation[ControlSamples$Gestation!='NA'] %>% range()
    
    write_csv(GSE74738.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE74738/GSE74738.meta_v1.csv")
  }
  
  
  
  ## `GEO_num[3]` from GEO database. 
  if(GSEnum%in%c("GSE69502", "all")){
    
    # establish a proper samplesheet for 
    GSE69502.meta_v1 <- GSE69502.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`sample tissue:ch1`, 'chorionic villi', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`ntd status:ch1`, 'control', 'Uncomplicated')) %>% 
      mutate(Gestation =`fetal gestational age (weeks):ch1`) %>% 
      inset('Trimester', value=rep('Second', nrow(GSE69502.meta))) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE69502.meta))) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE69502.meta))) %>% 
      mutate(Fetal_Sex= str_replace_all(`fetal sex:ch1`, c('F'='Female', 'M'='Male'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE69502.meta))) %>% 
      inset('Study', value=rep('GSE69502', nrow(GSE69502.meta))) %>% 
      mutate(Array= str_replace_all(description.3, ".+ ", "")) %>% 
      mutate(Slide= str_replace_all(description.2, ".+ ", "")) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE69502.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE69502/GSE69502.meta_v1.csv")
    ##getwd()
    # get RGset for GSE66210
    # GSE69502.RGset <- getRGset('GSE69502')
  }
  
  
  ## `GEO_num[4]` from GEO database. 
  if(GSEnum%in%c("GSE75196", "all")){
    # establish a proper samplesheet for 
    GSE75196.meta_v1 <- GSE75196.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample=str_replace_all(`tissue:ch1`, 'placenta', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`disease:ch1`, c('normal/healthy'='Uncomplicated', 'preeclampsia'='PE'))) %>% 
      mutate(Gestation =`gestation (wk):ch1`) %>% 
      inset('Trimester', value=ifelse(GSE75196.meta$`gestation (wk):ch1` >= 37, 'Term', 'Third')) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE75196.meta))) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE75196.meta))) %>% 
      mutate(Fetal_Sex= `Sex:ch1`) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE75196.meta))) %>% 
      inset('Study', value=rep('GSE75196', nrow(GSE75196.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE75196.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE75196/GSE75196.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE75196.RGset <- getRGset('GSE75196')
  }
  
  
  ## `GEO_num[5]` from GEO database. 
  if(GSEnum%in% c("GSE75248", "all")){
    # establish a proper samplesheet for 
    GSE75248.meta_v1 <- GSE75248.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(description, 'normal infant placenta tissue', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`birth weight group:ch1`, 'AGA', 'Uncomplicated')) %>% 
      inset('Gestation', value=rep('NA', nrow(GSE75248.meta))) %>% 
      inset('Trimester', value=rep('Term', nrow(GSE75248.meta))) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE75248.meta))) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE75248.meta))) %>% 
      mutate(Fetal_Sex= `gender:ch1`) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE75248.meta))) %>% 
      inset('Study', value=rep('GSE75248', nrow(GSE75248.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE75248.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE75248/GSE75248.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE75248.RGset <- getRGset('GSE75248')
  }
  
  ## `GEO_num[6]` from GEO database. GSE93429 has no IDAT files, so it was removed.
  
  ## `GEO_num[7]` GSE93208
  if(GSEnum%in% c("GSE93208", "all")){
    # establish a proper samplesheet for 
    GSE93208.meta_v1 <- GSE93208.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`tissue:ch1`, 'chorionic villus', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      inset('Disease', value = rep('Uncomplicated', nrow(GSE93208.meta))) %>% 
      # mutate(Gestation= str_replace_all(`gestational age (weeks):ch1`, c('8 to 10'='9', '12 to 14'='13'))) %>%
      mutate(Gestation= `gestational age (weeks):ch1`) %>%
      inset('Trimester', value=ifelse(GSE93208.meta$`gestational age (weeks):ch1`=='8 to 10', 'First', 'Second')) %>% 
      mutate(Additional_Info = title) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE93208.meta))) %>%  
      inset('Fetal_Sex', value=rep('NA', nrow(GSE93208.meta))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE93208.meta))) %>% 
      inset('Study', value=rep('GSE93208', nrow(GSE93208.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE93208.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE93208/GSE93208.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE93208.RGset <- getRGset('GSE93208')
  }
  
  
  ## `GEO_num[8]` GSE115508
  # additional info is based on pathology groups.
  if(GSEnum%in% c("GSE115508", "all")){
    
    # establish a proper samplesheet for 
    GSE115508.meta_v1 <- GSE115508.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`tissue:ch1`, 'Chorionic villi', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`pathology group:ch1`, 'non_chorioamnionitis', 'Uncomplicated')) %>% 
      mutate(Gestation= `gestational age:ch1`) %>% 
      inset('Trimester', value=ifelse(GSE115508.meta$`gestational age:ch1` >= 37, 'Term', 'Third')) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE115508.meta))) %>%  
      inset('Additional_Info_2', value=rep('NA', nrow(GSE115508.meta))) %>%  
      mutate(Fetal_Sex= str_replace_all(`fetal sex:ch1`, c('F'='Female', 'M'='Male'))) %>% 
      inset('ArrayType', value=rep('EPIC', nrow(GSE115508.meta))) %>% 
      inset('Study', value=rep('GSE115508', nrow(GSE115508.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE115508.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE115508/GSE115508.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE115508.RGset <- getRGset('GSE115508')
  }
  
  
  ## `GEO_num[9]` GSE120250
  #outlier information provided in Additional_Info column
  
  if(GSEnum%in% c("GSE120250", "all")){
    
    # establish a proper samplesheet for 
    GSE120250.meta_v1 <- GSE120250.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`title`, c('ART_..+'='Placenta', 'Cont..+'='Placenta'))) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`characteristics_ch1.2`, 'art treatment: NA', 'Uncomplicated')) %>% 
      inset('Gestation', value=rep('NA', nrow(GSE120250.meta))) %>% 
      inset('Trimester', value=rep('Term', nrow(GSE120250.meta))) %>% 
      mutate(Additional_Info= str_replace_all(`outlier:ch1`, c('N'='Standard', 'Y'='Outlier'))) %>%  #outliers
      inset('Additional_Info_2', value=rep('NA', nrow(GSE120250.meta))) %>%  
      mutate(Fetal_Sex= str_replace_all(`gender:ch1`, c('F'='Female', 'M'='Male'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE120250.meta))) %>% 
      inset('Study', value=rep('GSE120250', nrow(GSE120250.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE120250.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE120250/GSE120250.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE120250.RGset <- getRGset('GSE120250')
  }
  
  
  ## `GEO_num[11]` GSE100197 
  #Additional_Info: pathology group
  #additional_info_2: sample tissue
  if(GSEnum%in% c("GSE100197", "all")){
    
    # establish a proper samplesheet for 
    GSE100197.meta_v1 <- GSE100197.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`sample tissue:ch1`, 'Chorionic..+', 'Placenta')) %>%  
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`pathology group:ch1`, c('REPLICATE'='TechnicalReplicate', 'Term'='Uncomplicated'))) %>% 
      mutate(Gestation= `gestational age:ch1`) %>% 
      inset('Trimester', value=ifelse(GSE100197.meta$`gestational age:ch1` >= 37, 'Term', 'Third')) %>% 
      mutate(Additional_Info = `sample tissue:ch1`) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE100197.meta))) %>%  
      mutate(Fetal_Sex= str_replace_all(`fetal sex:ch1`, c('FEMALE'='Female', 'MALE'='Male'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE100197.meta))) %>% 
      inset('Study', value=rep('GSE100197', nrow(GSE100197.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "[^_]+_[^_]+")) %>% 
      # mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE100197.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE100197/GSE100197.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE100197.RGset <- getRGset('GSE100197')
  }
  
  
  
  ## `GEO_num[12]` GSE98224, 48 samples are methylation samples, other 48 are expression arrary data 
  #Additional_Info: previous hypertensive pregnancy
  #additional_info_2: previous miscarriage
  if(GSEnum%in% c("GSE98224", "all")){
    
    # establish a proper samplesheet for 
    GSE98224.meta_v1 <- GSE98224.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`tissue:ch1`, 'placenta', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>%
      mutate(Disease = str_replace_all(`diagnosis:ch1`, 'non-PE','Uncomplicated')) %>% 
      mutate(Gestation= `ga week:ch1`) %>% 
      inset('Trimester', value=ifelse(GSE98224.meta$`ga week:ch1` >= 37, 'Term', 'Third')) %>% 
      mutate(Additional_Info = `previous hypertensive pregnancy:ch1`) %>% 
      mutate(Additional_Info_2= `previous miscarriage:ch1`) %>%  
      mutate(Fetal_Sex= str_replace_all(`gender:ch1`, c('F'='Female', 'M'='Male'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE98224.meta))) %>% 
      inset('Study', value=rep('GSE98224', nrow(GSE98224.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE98224.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE98224/GSE98224.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE98224.RGset <- getRGset('GSE98224')
  }
  
  ## 13. GSE106089 
  #Additional_Info: antibiotic use
  #additional_info_2: NA
  
  if(GSEnum%in% c("GSE106089", "all")){
    
    # establish a proper samplesheet for 
    GSE106089.meta_v1 <- GSE106089.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(description, 'normal placental sample', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      inset('Disease', value = ifelse(GSE106089.meta$`bacteria present:ch1`=='None', 'Uncomplicated', 'BacteriaDetected')) %>% 
      inset('Gestation', value=rep('NA', nrow(GSE106089.meta))) %>% 
      inset('Trimester', value=rep('Third', nrow(GSE106089.meta))) %>% 
      mutate(Additional_Info = `characteristics_ch1.4`) %>% 
      inset('Additional_Info_2', value=rep('NA', nrow(GSE106089.meta))) %>%  
      mutate(Fetal_Sex= str_replace_all(`fetal sex:ch1`, c('^male'='Male','^female'='Female'))) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE106089.meta))) %>% 
      inset('Study', value=rep('GSE106089', nrow(GSE106089.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE106089.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE106089/GSE106089.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE106089.RGset <- getRGset('GSE106089')
  }
  
  ## 14. GSE71678, `GEO_num[13]` the meta data is wrong in `hyb_protocol` coplumn, it should be 450k not 27k.
  # Additional_Info: NA
  # additional_info_2: NA
  if(GSEnum%in% c("GSE71678", "all")){
    
    # establish a proper samplesheet for 
    GSE71678.meta_v1 <- GSE71678.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(description, 'normal infant placenta tissue', 'Placenta')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      inset('Disease', value = ifelse(GSE71678.meta$`placental as levels:ch1`=='NA', 'Uncomplicated', 'ArsenicDetected')) %>% 
      mutate(Gestation = `gestational age:ch1`) %>% 
      inset('Trimester', value=ifelse(GSE71678.meta$`gestational age:ch1` >= 37, 'Term', 'Third')) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE71678.meta))) %>%  
      inset('Additional_Info_2', value=rep('NA', nrow(GSE71678.meta))) %>%  
      mutate(Fetal_Sex= `infant gender:ch1`) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE71678.meta))) %>% 
      inset('Study', value=rep('GSE71678', nrow(GSE71678.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "R..C..")) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE71678.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE71678/GSE71678.meta_v1.csv")
    
    ##getwd()
    # get RGset for GSE66210
    # GSE71678.RGset <- getRGset('GSE71678')
  }
  
  
  ## 16. GSE113600 `GEO_num[14]`
  # Additional_Info: NA
  # additional_info_2: NA
  if(GSEnum%in% c("GSE113600", "all")){
    
    # establish a proper samplesheet for 
    GSE113600.meta_v1 <- GSE113600.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= str_replace_all(`tissue:ch1`, 'decidua', 'Decidua')) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      mutate(Disease = str_replace_all(`condition:ch1`, 'induced abortion', 'Uncomplicated')) %>% 
      inset('Gestation', value=rep('NA', nrow(GSE113600.meta))) %>% 
      inset('Trimester', value=rep('First', nrow(GSE113600.meta))) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE113600.meta))) %>%  
      inset('Additional_Info_2', value=rep('NA', nrow(GSE113600.meta))) %>%  
      inset('Fetal_Sex', value=rep('NA', nrow(GSE113600.meta))) %>% 
      inset('ArrayType', value=rep('EPIC', nrow(GSE113600.meta))) %>% 
      inset('Study', value=rep('GSE113600', nrow(GSE113600.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = "decidua.")) %>% 
      mutate(Slide= str_extract(Sample_name, pattern = "[^_]+")) %>%
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE113600.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE113600/GSE113600.meta_v1.csv")
    
    # GSE113600.RGset <- getRGset('GSE113600')
    # phenotypes <- read.metharray.sheet('/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE113600', pattern = "_v1.csv$")
    
  }
  
  
  ## 18. GSE66459, `GEO_num[16]`
  # Additional_Info: NA
  # additional_info_2: NA
  
  if(GSEnum%in% c("GSE66459", "all")){
    
    # establish a proper samplesheet for 
    GSE66459.meta_v1 <- GSE66459.meta %>% 
      mutate(GEO_accession=geo_accession) %>% 
      mutate(Sample= `tissue:ch1`) %>% 
      mutate(Sample_name = str_replace_all(supplementary_file, c("_(Red|Grn).idat.gz"="", ".+/"=""))) %>% 
      mutate(Sample_nameToRowname=Sample_name) %>% 
      inset('Disease', value = rep('Uncomplicated', nrow(GSE66459.meta))) %>% 
      mutate(Gestation = as.numeric(`gestational_age (days):ch1`)/7) %>% 
      inset('Trimester', value=ifelse(.$Gestation >= 37, 'Term', 'Third')) %>% 
      inset('Additional_Info', value=rep('NA', nrow(GSE66459.meta))) %>%  
      inset('Additional_Info_2', value=rep('NA', nrow(GSE66459.meta))) %>%  
      mutate(Fetal_Sex = `gender:ch1`) %>% 
      inset('ArrayType', value=rep('450K', nrow(GSE66459.meta))) %>% 
      inset('Study', value=rep('GSE66459', nrow(GSE66459.meta))) %>% 
      mutate(Array= str_extract(Sample_name, pattern = '(_.|_..)$')) %>% 
      mutate(Array= str_replace_all(Array, '_', '')) %>% 
      mutate(Slide=  str_extract(Sample_name, pattern= "_(.*?)_")) %>% 
      mutate(Slide=  str_replace_all(Slide, '_', '')) %>% 
      dplyr::select(GEO_accession:Slide) %>% 
      column_to_rownames(var='Sample_nameToRowname')
    
    write_csv(GSE66459.meta_v1, path = "/home/qianhui/DNAme/Process_decidual/GEO_dataSets/GSE66459/GSE66459.meta_v1.csv")
    
    # GSE66459.RGset <- getRGset('GSE66459')
  }
  
  save(GSE66210.meta_v1, GSE74738.meta_v1, GSE69502.meta_v1, GSE75196.meta_v1, 
       GSE75248.meta_v1, GSE93208.meta_v1, GSE115508.meta_v1, GSE120250.meta_v1, 
       GSE100197.meta_v1, GSE98224.meta_v1, GSE106089.meta_v1, GSE71678.meta_v1, 
       GSE113600.meta_v1, GSE66459.meta_v1,
       file = paste0("/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_allMetaData_v1_", GSEnum, ".RData"))
}
