
# 1. Write the function to  divide traning and test data:
## the argument datMeth_placenta is NormBMlist$NormMAll

MclustTrainTestFUN <- function(datMeth_placenta){
  
  library(gridExtra)
  # Deviding data:
  ind <- sample.int(n=ncol(datMeth_placenta), size = floor(0.7*ncol(datMeth_placenta)), replace = FALSE) # floor() is used to return an integer instead of a decimal.
  
  datMethTraining <- datMeth_placenta[,ind]
  
  datMethTesting <- datMeth_placenta[,-ind]
  
  # # plot numbers of outliers in each data set:
  # TrainingData=table(colnames(datMethTraining)%in%OutlierDetail[OutlierDetail$Outlier%in%"Outlier",]$SampleNames) %>% 
  #     as.data.frame() %>% dplyr::select(Freq)
  # 
  # TestData=table(colnames(datMethTesting)%in%OutlierDetail[OutlierDetail$Outlier%in%"Outlier",]$SampleNames) %>% 
  #   as.data.frame() %>% dplyr::select(Freq) 
  # 
  # tibble(Samples=c("Pure placenta", "Outlier"), TrainingData,TestData) %>% 
  #   grid.table(rows=NULL)
  
  # Training data:
  PCA_traning <- PCA(t(datMethTraining))
  
  df_training <- PCA_traning$ind$coord %>% as.data.frame() %>% rownames_to_column(var="Sample_name") %>% 
    left_join(NormBMlist$phenoDataAll, by="Sample_name") %>% column_to_rownames(var="Sample_name")
  
  mclust_training <- mclust::Mclust(df_training[,c("Dim.1","Dim.2")],2) #for 2 groups
  
  table(mclust_training$classification,df_training$Sample)
  with(df_training,plot(Dim.2 ~ Dim.1,pch=as.numeric(as.factor(Sample)),col=mclust_training$classification))
  with(df_training,legend("topright",levels(as.factor(Sample)),pch=1:6))
  
  df_training$class <- mclust_training$classification
  
  # test data
  
  PCA_test <- predict.PCA(object = PCA_traning, newdata = t(datMethTesting))
  
  df_test <- PCA_test$coord %>% as.data.frame() %>% rownames_to_column(var="Sample_name") %>% 
    left_join(NormBMlist$phenoDataAll, by="Sample_name") %>% column_to_rownames(var="Sample_name")
  
  mclust_test <- predict.Mclust(object = mclust_training, newdata = df_test[,c("Dim.1","Dim.2")])
  
  df_test$class <- mclust_test$classification
  
  return(list(PCA_traning=PCA_traning, mclust_training=mclust_training, mclust_test=mclust_test,
              df_training=df_training, df_test=df_test))
  
}


# 2. Write the function to plot traning and test result of classification:
## df == the df_training or df_test
## PCA_training == the PCA result from taining data

PlotClassFUN <- function(df, PCA_traning){
  if(is.null(PCA_traning)==TRUE){
    PCAclassPlot <- df %>%
      ggplot(aes(Dim.1, Dim.2)) + 
      geom_point(aes(shape = Sample, col=as.factor(class)), alpha = 1, size = 2)+ 
      scale_shape_manual(values=c("Amnion"=0, "Chorion"=1, "Decidua"=2, "Maternal whole blood"=3, 
                                  "Placenta"=4, "Umbilical cord blood"=5))+
      scale_color_manual(name="Clusters", values = c("1"="#000000","2"="#D95F02"), 
                         labels=c("1"="Cluster 1", "2"="Cluster 2"))+
      labs(x="PC1", y="PC2")+
      scale_y_continuous(breaks=c(-600,-300,0, 300, 600), limits = c(-600, 800))+
      scale_x_continuous(breaks=c(-400,0, 400, 800), limits = c(-410, 1150))+
      # ylim(-600,600)+
      # xlim(-410, 1150)+
      guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
      theme(text = element_text(size=9),
            legend.title=element_text(size=9),
            legend.text=element_text(size=9),
            axis.text.x =element_text(size=9, colour = "black"),
            axis.text.y =element_text(size=9, colour = "black"),
            legend.position="top")
  }else{
    PCAclassPlot <- df %>%
      ggplot(aes(Dim.1, Dim.2)) + 
      geom_point(aes(shape = Sample, col=as.factor(class)), alpha = 1, size = 2)+ 
      scale_shape_manual(values=c("Amnion"=0, "Chorion"=1, "Decidua"=2, "Maternal whole blood"=3, 
                                   "Placenta"=4, "Umbilical cord blood"=5))+
      scale_color_manual(name="Clusters", values = c("1"="#000000","2"="#D95F02"), 
                         labels=c("1"="Cluster 1", "2"="Cluster 2"))+
      labs(x=paste0("PC1 (", PCA_traning$eig[, 2][1], ")"), 
           y=paste0("PC2 (", PCA_traning$eig[, 2][2], ")"))+
      scale_y_continuous(breaks=c(-600,-300,0, 300, 600), limits = c(-600, 800))+
      scale_x_continuous(breaks=c(-400,0, 400, 800), limits = c(-410, 1150))+
      # ylim(-600,600)+
      # xlim(-410, 1150)+
      guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
      theme(text = element_text(size=9),
            legend.title=element_text(size=9),
            legend.text=element_text(size=9),
            axis.text.x =element_text(size=9, colour = "black"),
            axis.text.y =element_text(size=9, colour = "black"),
            legend.position="top")
  }
  PCAclassPlot
}


# 3. Write the function to generate tables for results from traning and test data:

## df == the df_training or df_test

TableTrainTestFUN <- function(df){
  
  Placenta <- df$class %>% 
    table() %>% as.data.frame() %>% 
    dplyr::select(PlacentaSamples=Freq) 
  
  if(unique(df[!(df$Sample%in%"Placenta"),]$class)%in%"1"){
    NotPlacenta <- df[!(df$Sample%in%"Placenta"),]$class %>% 
      table() %>% as.data.frame() %>% 
      dplyr::select(NotPlacentaSamples=Freq) %>% rbind("2"=0)
  }
  
  cbind(Placenta,NotPlacenta) %>% 
    rownames_to_column(var="PredictedResults") %>% 
    mutate(PredictedResults=str_replace_all(PredictedResults, c("1"="Not placenta", "2"="Predicted placenta")))  %>% 
    grid.table(rows=NULL)
}

