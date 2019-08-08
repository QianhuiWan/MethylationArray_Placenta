# write the ggplot PCA function:================================
ggplotPCA <- function(PCA_Placenta_M, Matched_pd, value){
  
  library(ggpubr)
  
  new_plotDF <- tibble(x=PCA_Placenta_M$ind$coord[,1], 
                       y=PCA_Placenta_M$ind$coord[,2],
                       # GA = Matched_pd$Gestation,
                       FetalSex= Matched_pd$Fetal_Sex,
                       SampleNames= Matched_pd$Sample_name,
                       Tissue0= Matched_pd$Sample,
                       Tissue=str_replace_all(Matched_pd$Sample, c(
                         ".lacenta.+"="Placenta",
                         "umbilical"="Umbilical")),
                       Trimester=Matched_pd$Trimester,
                       Study=Matched_pd$Study,
                       Outlier=ifelse(Matched_pd$Additional_Info%in%"Outlier", "Outlier", "Pure"),
                       ArrayType=Matched_pd$ArrayType) %>% 
    mutate(Outlier2 = case_when(SampleNames %in% c("GSM1702177_7970368142_R05C02", "GSM1947213_6008581028_R01C01") ~ "Outlier",
                                !(SampleNames %in% c("GSM1702177_7970368142_R05C02", "GSM1947213_6008581028_R01C01")) ~ Outlier))
  
  colnames(new_plotDF) <- c("x", "y", #"Gestational age", 
                            "Sex", "SampleNames",
                            "Tissue0","Tissue", 
                            "Trimester", "Study", "Outlier", "ArrayType", "Outlier2")
  
  new_plotDF$Tissue[which(new_plotDF$Outlier2=="Outlier")] <- rep("Impure placenta", 11)
  
  # Plot the new dataframe
  ## color to use
  
  PCAplot_norm408 <- new_plotDF %>%
    ggplot(aes(x, y)) + 
    geom_point(aes(shape = Tissue, col=Trimester), alpha = 1, size = 5)+
    scale_shape_manual(values=1:11)+
    labs(x="PC1", y="PC2")+
    theme(text = element_text(size=20))

  # save the MDS plot as pdf:
  pdf(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", value, "_allGEOsamples408.pdf"),
      width = 13, height = 7)
  
  theme_set(theme_pubr())
  PCAplot_norm408
  dev.off()
  
  # show the plot
  theme_set(theme_pubr())
  PCAplot_norm408 %>% print()
  
  # return new_plotDF
  return(new_plotDF)
}
# end of ggplotPCA function=======================================

