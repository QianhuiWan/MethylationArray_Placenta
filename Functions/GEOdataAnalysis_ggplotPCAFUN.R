# write the ggplot PCA function:================================
ggplotPCA <- function(PCA_Placenta_M, Matched_pd, value){
  
  new_plotDF <- data.frame(x=PCA_Placenta_M$ind$coord[,1], 
                           y=PCA_Placenta_M$ind$coord[,2],
                           # GA = Matched_pd$Gestation,
                           Sex= Matched_pd$Fetal_Sex,
                           Tissue=str_replace_all(Matched_pd$Sample, c(
                             "placental"="Placental",
                             "umbilical"="Umbilical")),
                           Trimester=Matched_pd$Trimester,
                           Study=Matched_pd$Study,
                           Outlier=ifelse(Matched_pd$Additional_Info%in%"Outlier", "Outlier", "Standard"),
                           ArrayType=Matched_pd$ArrayType
  )
  
  colnames(new_plotDF) <- c("x", "y", #"Gestational age", 
                            "Sex", "Tissue", 
                            "Trimester", "Study", "Outlier", "ArrayType")
  
  # Plot the new dataframe
  PCAplot_norm551 <- new_plotDF %>%
    ggplot(aes(x, y)) + 
    geom_point(aes(shape = Tissue, col=Outlier), alpha = 1, size = 5)+ 
    # scale_color_manual(values = namedPal2)+
    scale_shape_manual(values=1:11)+
    labs(x="PC1", y="PC2")+
    theme(text = element_text(size=20))
  # axis.text.x = element_text(angle=90, hjust=1)) 
  
  # save the PCA plot as pdf:
  pdf(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/PCAplot-Normalised", value, "_GEOsamples.pdf"),
      width = 13, height = 7)
  
  theme_set(theme_bw())
  PCAplot_norm551
  dev.off()
  
  # show the plot
  theme_set(theme_bw())
  PCAplot_norm551 %>% print()
  
  # return new_plotDF
  return(new_plotDF)
}
# end of ggplotPCA function=======================================
