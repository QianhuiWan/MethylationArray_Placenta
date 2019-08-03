

# 1. Write the function to plot traning and test result of classification:
## df == the df_training or df_test
## PCA_training == the PCA result from taining data

PlotClassFUN <- function(df, PCA_traning){
  if(is.null(PCA_traning)==TRUE){
    PCAclassPlot <- df %>%
      ggplot(aes(Dim.1, Dim.2)) + 
      geom_point(aes(shape = Sample, col=as.factor(class)), alpha = 1, size = 2)+ 
      scale_shape_manual(values=c(0:5))+
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
      scale_shape_manual(values=c(0:5))+
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

# 2. Write the function to generate tables for results from traning and test data:

## df == the df_training or df_test

TableTrainTestFUN <- function(df){
  
  Placenta <- df[df$Sample%in%"Placenta",]$class %>% 
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

