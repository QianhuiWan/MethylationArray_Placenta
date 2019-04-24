
library(tidyverse)
library(reshape2)
library(minfi)
library(magrittr)
library(doParallel)
# load data 
NormBMlist <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/NormBMlist.rds")

registerDoParallel(cores = 5)
table(colnames(NormBMlist$NormBAll) %in% NormBMlist$phenoDataAll$Sample_name)
# Matched_pd <- NormBMlist$phenoDataAll
# PhenodataPlacenta <- Matched_pd[Matched_pd$Trimester%in%"First",]
PhenodataPlacenta <- NormBMlist$phenoDataAll

PhenodataPlacenta$Outlier <- ifelse(PhenodataPlacenta$Additional_Info%in%"Outlier", "Outlier", "Standard")
PhenodataPlacenta$Sample2 <- str_replace_all(PhenodataPlacenta$Sample, c(".+lacent.+"="Pure placenta", "umbilical cord"="Umbilical cord"))
PhenodataPlacenta$Sample2[which(PhenodataPlacenta$Outlier=="Outlier")] <- rep("Impure placenta",9)

# mean beta values for each group
Placenta_betaNorm_SampleGroupMean <- NormBMlist$NormBAll %>% t() %>% dplyr::as_data_frame() %>% 
                                     inset('Sample2', value=PhenodataPlacenta$Sample2) %>% 
                                     reshape2::melt() %>% 
                                     group_by(Sample2, variable) %>% 
                                     summarise(SampleGroupMean=mean(value)) %>% 
                                     reshape2::dcast(variable~Sample2) %>% 
                                     column_to_rownames(var="variable") %>% 
                                     as.matrix()

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/Beta_density_plot_551.tiff", 
     width = 17, height = 12, units="cm",  res = 300)
par(cex=0.8)
densityPlot(Placenta_betaNorm_SampleGroupMean, 
            sampGroups = PhenodataPlacenta$Sample2 %>% unique() %>% sort(),
            legend=TRUE,
            xlab = "DNA methylation proportion")
dev.off()

pdf("/home/qianhui/DNAme/Process_decidual/figures/Beta_density_plot_551.pdf", 
    width = 13, height = 7);
par(cex=1)
densityPlot(Placenta_betaNorm_SampleGroupMean, 
            sampGroups = PhenodataPlacenta$Sample2 %>% unique() %>% sort(),
            legend=TRUE,
            xlab = "DNA methylation proportion")
dev.off()



# or use ggplot to plot density plot
Placenta_betaNorm_SampleGroupMeanDF <- NormBMlist$NormBAll %>% t() %>% as.data.frame(stringsAsFactors=FALSE) %>% 
                                       inset('Sample2', value=PhenodataPlacenta$Sample2) %>% 
                                       reshape2::melt() %>% 
                                       group_by(Sample2, variable) %>% 
                                       summarise(SampleGroupMean=mean(value))

density551 <- Placenta_betaNorm_SampleGroupMeanDF %>% 
  ggplot(aes(x=SampleGroupMean, color=Sample2))+geom_density(size=0.5)+ #aes(linetype=Sample2)
  scale_color_manual(name="Tissue", values = c("#000000", "#BD0026", "#984EA3", "#1B9E77", "#999999",  "#D95F02", "#F781BF"))+
  labs(x="Beta value", y="Density", linetype="Tissue")+
  # guides(colour=guide_legend(override.aes = list(size=0.5)))+
  theme(axis.text.y = element_text(size=9, colour ="black"),
        axis.text.x = element_text(size=9, colour ="black"),
        axis.title=element_text(size=9, colour ="black"),#,face="bold"
        legend.position = "right",
        legend.key.size = unit(0.9,"line"),
        text = element_text(size=9, colour ="black"))

# save density551
tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/Beta_density_ggplot_551_v1.tiff", 
     width = 17, height = 12, units="cm",  res = 300)
theme_set(theme_bw())
density551
dev.off()

# pdf("/home/qianhui/DNAme/Process_decidual/figures/Beta_density_ggplot_551.pdf", 
#     width = 13, height = 7);
# theme_set(theme_bw())
# density551
# dev.off()


saveRDS(density551, file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/density551Plot_ggplot.rds")


