library(tidyverse)
library(magrittr)
library(tibble)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(cluster)
library(RColorBrewer)
theme_set(theme_pubr())

registerDoParallel(cores=10)

# load data 
## GEOmeta, GEO_phenotypes
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_phenotypes_add10.RData")

NormBMlist <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/NormBMlist.rds")


# calculate impure outlier percentage for each trimester
table(colnames(NormBMlist$NormBAll) %in% NormBMlist$phenoDataAll$Sample_name)
Matched_pd <- NormBMlist$phenoDataAll
Matched_pd$Outlier <- ifelse(Matched_pd$Additional_Info%in%"Outlier",
                             # | Matched_pd$Sample_name%in%c("GSM1702177_7970368142_R05C02", "GSM1947213_6008581028_R01C01"), 
                             "Outlier", "Standard")


OutlierPerPlot <- Matched_pd %>% dplyr::select(Sample_name, Trimester, Outlier) %>% melt() %>% 
  ggplot(aes(Trimester, fill=Outlier))+
  geom_bar(position = "fill") +
  scale_fill_manual(name="Placenta", 
                    values = c("Standard"="#D95F02", "Outlier"="#1B9E77"),
                    # breaks = c(),
                    labels = c("Pure", "Impure"))+
  labs(y="Percentage")+
  theme(text = element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        axis.title.x = element_text(size=9),
        axis.text.x =element_text(size=9, colour = "black"),
        axis.text.y =element_text(size=9, colour = "black"),
        legend.position="right")

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/OutlierPerPlot.tiff",
     width = 17, height = 8, units = "cm", res = 300)
OutlierPerPlot 
dev.off()

saveRDS(OutlierPerPlot, file = "/home/qianhui/DNAme/Process_decidual/figures/OutlierPerPlot.rds")
