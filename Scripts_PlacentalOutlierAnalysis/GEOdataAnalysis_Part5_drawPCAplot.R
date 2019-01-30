
library(tidyverse)
library(magrittr)
library(tibble)
library(FactoMineR)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
registerDoParallel(cores=10)

# load data 
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/allSamplesResiduals_ControlCorrected.rds")
## GEOmeta, GEO_phenotypes
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_phenotypes_add10.RData")

# NormBMlist <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/NormBMlist.rds")

# BMall_combat, MDS_normCorrect551.Beta, MDS_normCorrect551.M
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_joinedSamples_NormBatchCorrect551_BetaM_MDS_OutlierConsidered.RData")

## annotation needed for using functions
ann450k <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/ann450k_GEO15sets.rds")
annEPIC <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/annEPIC_GEO15sets.rds")

# all maples PCA
table(colnames(Residuals) %in% GEO_phenotypes$Sample_name)

Matched_pd <- GEO_phenotypes[match(colnames(Residuals), GEO_phenotypes$Sample_name),]

# PCA_Placenta_M <- FactoMineR::PCA(base::t(Residuals))

# T1 PCA
PhenodataPlacentaT1 <- Matched_pd[Matched_pd$Sample%in%"Placenta" & Matched_pd$Trimester%in%"First",] 
Placenta_Mnorm_T1 <-  Residuals[,colnames(Residuals)%in%PhenodataPlacentaT1$Sample_name]
# Placenta_betaNorm <- 2^Placenta_Mnorm_v1 / (1+2^Placenta_Mnorm_v1)

table(PhenodataPlacentaT1$Sample_name==colnames(Placenta_Mnorm_T1))
# PCA_PlacentaT1_M <- FactoMineR::PCA(base::t(Placenta_Mnorm_T1))

# Term PCA
PhenodataPlacentaTerm <- Matched_pd[Matched_pd$Sample%in%"Placenta" &Matched_pd$Trimester%in%"Term",] 
Placenta_Mnorm_Term <-  Residuals[,colnames(Residuals)%in%PhenodataPlacentaTerm$Sample_name]

table(PhenodataPlacentaTerm$Sample_name==colnames(Placenta_Mnorm_Term))
# PCA_PlacentaTerm_M <- FactoMineR::PCA(base::t(Placenta_Mnorm_Term))

# save(PCA_Placenta_M, PCA_PlacentaT1_M, PCA_PlacentaTerm_M, file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term.RData")
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term.RData")


# write the ggplot PCA function:================================
ggplotPCA <- function(PCA_Placenta_M, Matched_pd, value){
  
new_plotDF <- tibble(x=PCA_Placenta_M$ind$coord[,1], 
                     y=PCA_Placenta_M$ind$coord[,2],
                     GA = Matched_pd$Gestation,
                     FetalSex= Matched_pd$Fetal_Sex,
                     Tissue0= Matched_pd$Sample,
                     Tissue=str_replace_all(Matched_pd$Sample, c(
                            ".lacenta.+"="Placenta",
                            "umbilical"="Umbilical")),
                     Trimester=Matched_pd$Trimester,
                     Study=Matched_pd$Study,
                     Outlier=ifelse(Matched_pd$Additional_Info%in%"Outlier", "Outlier", "Standard"),
                     ArrayType=Matched_pd$ArrayType
)

colnames(new_plotDF) <- c("x", "y", "Gestational age", "FetalSex", 
                          "Tissue0","Tissue", 
                          "Trimester", "Study", "Outlier", "ArrayType")

new_plotDF$Tissue[which(new_plotDF$Outlier=="Outlier")] <- rep("Placental outliers",9)

# Plot the new dataframe
## color to use

PCAplot_norm551 <- new_plotDF %>%
  ggplot(aes(x, -y)) + 
  geom_point(aes(shape = Tissue, col=Trimester), alpha = 1, size = 5)+
  # geom_point(aes(shape = Trimester, col=Tissue), alpha = 1, size = 5)+ 
  # scale_color_manual(values = namedPal2)+
  scale_shape_manual(values=1:11)+
  # scale_color_manual(name="Tissue", values = c("#000000", "#BD0026", "#984EA3", "#999999", "#D95F02", "#1B9E77", "#F781BF", "#A65628"))+
  labs(x="PC1", y="PC2")+
  theme(text = element_text(size=20))
# axis.text.x = element_text(angle=90, hjust=1)) 

# save the MDS plot as pdf:
pdf(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", value, "_joinedGEOsamples551.pdf"),
    width = 13, height = 7)

theme_set(theme_bw())
PCAplot_norm551
dev.off()

# save as tiff
tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", value, "_joinedGEOsamples551.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)

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



# use ggplotPCA function
## all sampoes
PCA551 <- ggplotPCA(PCA_Placenta_M = PCA_Placenta_M, Matched_pd = Matched_pd, value = "all551residuals")

PCA551plot <- PCA551 %>%
  ggplot(aes(x, -y)) + 
  # geom_point(aes(shape = Tissue, col=Trimester), alpha = 1, size = 5)+ 
  geom_point(aes(shape = Trimester, col=Tissue), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("#000000", "#BD0026", "#984EA3", "#999999", "#D95F02", "#1B9E77", "#F781BF", "#A65628"))+
  scale_shape_manual(values=c(0,1,5,2))+
  labs(x="PC1", y="PC2")+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        legend.position="top")

tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", "allPCA551", "_joinedGEOsamples551.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)
theme_set(theme_bw())
PCA551plot
dev.off()

# ## T1 sampels
# PCAT1 <- ggplotPCA(PCA_Placenta_M = PCA_PlacentaT1_M, Matched_pd = PhenodataPlacentaT1, value = "all_T1_residuals")
# ## Term sampels
# PCATerm <- ggplotPCA(PCA_Placenta_M = PCA_PlacentaTerm_M, Matched_pd = PhenodataPlacentaTerm, value = "all_Term_residuals")


# 1st trimester & term samples
PCA551_First <- PCA551[PCA551$Trimester%in%"First",]
PCA551_Term <- PCA551[PCA551$Trimester%in%"Term",]

PCAT1plot <- PCA551_First %>%
  ggplot(aes(x, -y)) + 
  geom_point(aes(shape = Trimester, col=Tissue), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("#984EA3", "#999999", "#D95F02", "#1B9E77"))+
  scale_shape_manual(values=0)+
  labs(x="PC1", y="PC2")+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        legend.position="top")

PCATermplot <- PCA551_Term %>%
  ggplot(aes(x, -y)) + 
  geom_point(aes(shape = Trimester, col=Tissue), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("#000000", "#BD0026", "#984EA3","#D95F02", "#1B9E77", "#A65628"))+
  scale_shape_manual(values=5)+
  labs(x="PC1", y="PC2")+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        legend.position="top")

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/T1Samples_PCA551_First.tiff",
     width = 6, height = 4, units = "cm", res = 300)
theme_set(theme_pubr())
PCAT1plot
dev.off()


tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/TermSamples_PCA551_First.tiff",
     width = 20, height = 12, units = "cm", res = 300)
theme_set(theme_bw())
PCATermplot
dev.off()

# use ggarange to put plots together
PCAs <- ggarrange(PCA551plot, PCAT1plot, PCATermplot, 
                  labels = c("A", "B", "C"),
                  ncol = 3, nrow = 1, 
                  font.label=list(size = 12, color = "black", face = "bold", family = "Arial"),
                  legend="top", common.legend = TRUE)


tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/T1Samples_PCA551_3plots.tiff",
     width = 17, height = 8, units = "cm", res = 300)
# theme_set(theme_pubr())
PCAs 
dev.off()


saveRDS(PCAs, file = "/home/qianhui/DNAme/Process_decidual/figures/T1TermSamples_PCA551_3plots.rds")
