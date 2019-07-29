
library(tidyverse)
library(magrittr)
library(tibble)
library(FactoMineR)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(mclust)
library(cluster)
library(RColorBrewer)
theme_set(theme_pubr())

registerDoParallel(cores=10)

# load data 
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/allSamplesResiduals_ControlCorrected.rds")
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/allSamplesResiduals_ControlCorrected_residuals.rds")
## GEOmeta, GEO_phenotypes
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_phenotypes_add10.RData")
# # PCA_Placenta_M, PCA_PlacentaT1_M, PCA_PlacentaTerm_M
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term.RData")
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term_BMval.RData")

NormBMlist <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/NormBMlist.rds")

# BMall_combat, MDS_normCorrect551.Beta, MDS_normCorrect551.M
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_joinedSamples_NormBatchCorrect551_BetaM_MDS_OutlierConsidered.RData")

## annotation needed for using functions
ann450k <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/ann450k_GEO15sets.rds")
annEPIC <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/annEPIC_GEO15sets.rds")

# all saples PCA
Residuals <- NormBMlist$NormMAll
table(colnames(Residuals) %in% GEO_phenotypes$Sample_name)

Matched_pd <- GEO_phenotypes[match(colnames(Residuals), GEO_phenotypes$Sample_name),]

df <- PCA_Placenta_M$ind$coord %>% as.data.frame() %>% rownames_to_column(var="Sample_name") %>% 
  left_join(Matched_pd, by="Sample_name")

a1 <- mclust::Mclust(df[,c("Dim.1","Dim.2")],2) #for 2 groups
table(a1$classification,df$Sample)
with(df,plot(Dim.2 ~ Dim.1,pch=as.numeric(as.factor(Sample)),col=a1$classification))
with(df,legend("topright",levels(as.factor(Sample)),pch=1:6))

# plot mclust result
df$class <- a1$classification

PCAclassPlot <- df %>%
  ggplot(aes(Dim.1, Dim.2)) + 
  geom_point(aes(shape = Sample, col=as.factor(class)), alpha = 1, size = 2)+ 
  scale_shape_manual(values=c(0:5))+
  scale_color_manual(name="Clusters", values = c("1"="#000000","2"="#D95F02"), 
                     labels=c("1"="Cluster 1", "2"="Cluster 2"))+
  labs(x="PC1 (25.35%)", y="PC2 (10.20%)")+
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

tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/PCAmclustClassPlot.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)
theme_set(theme_pubr())
PCAclassPlot
dev.off()

# PCA_Placenta_M <- FactoMineR::PCA(base::t(Residuals))
# saveRDS(PCA_Placenta_M, file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_all408_Mval.RData")
# PCA_Placenta_M <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta379_Mval.RData")

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

# T1,T2, Term PCA
# PhenodataNoT3 <- Matched_pd[!(Matched_pd$Trimester%in%"Third"),] 
# Mnorm_NoT3 <-  Residuals[,colnames(Residuals)%in%PhenodataNoT3$Sample_name]
# 
# table(PhenodataNoT3$Sample_name==colnames(Mnorm_NoT3))
# PCA_NoT3 <- FactoMineR::PCA(base::t(Mnorm_NoT3))

# save(PCA_Placenta_M, PCA_PlacentaT1_M, PCA_PlacentaTerm_M, file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term_BMval.RData")
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term.RData")

# write the ggplot PCA function:================================
ggplotPCA <- function(PCA_Placenta_M, Matched_pd, value){
  
new_plotDF <- tibble(x=PCA_Placenta_M$ind$coord[,1], 
                     y=PCA_Placenta_M$ind$coord[,2],
                     GA = Matched_pd$Gestation,
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
  
colnames(new_plotDF) <- c("x", "y", "Gestational age", "FetalSex", "SampleNames",
                          "Tissue0","Tissue", 
                          "Trimester", "Study", "Outlier", "ArrayType", "Outlier2")

# new_plotDF$Tissue[which(new_plotDF$Outlier2=="Outlier")] <- rep("Placental outliers", 11)
new_plotDF$Tissue[which(new_plotDF$Outlier2=="Outlier")] <- rep("Impure placenta", 11)

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

# theme_set(theme_bw())
PCAplot_norm551
dev.off()

# save as tiff
tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", value, "_joinedGEOsamples551.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)

# theme_set(theme_bw())
PCAplot_norm551
dev.off()

# show the plot
# theme_set(theme_bw())
PCAplot_norm551 %>% print()

# return new_plotDF
return(new_plotDF)
}
# end of ggplotPCA function=======================================



# use ggplotPCA function
## all sampoes
PCA551 <- ggplotPCA(PCA_Placenta_M = PCA_Placenta_M, Matched_pd = Matched_pd, value = "all551residuals")
# PCA435 <- ggplotPCA(PCA_Placenta_M = PCA_NoT3, Matched_pd = PhenodataNoT3, value = "NoT3residuals")

PCA551$TissueNoGreen <- str_replace_all(PCA551$Tissue, c("Impure placenta"="Placenta"))

PCA551plot <- PCA551 %>%
  ggplot(aes(x, y)) + 
  # geom_point(aes(shape = Tissue, col=Trimester), alpha = 1, size = 5)+ 
  geom_point(aes(shape = Trimester, col=TissueNoGreen), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("Amnion"="#000000", "Chorion"="#BD0026", 
                                               "Decidua"="#984EA3", "Maternal whole blood"="#999999", 
                                               "Placenta"="#D95F02", #"Impure placenta"="#1B9E77", 
                                               "Umbilical cord blood"="#F781BF"))+
  scale_shape_manual(values=c(0,1,5,2))+
  labs(x="PC1 (25.35%)", y="PC2 (10.20%)")+
  scale_y_continuous(breaks=c(-600,-300,0, 300, 600), limits = c(-600, 800))+
  scale_x_continuous(breaks=c(-400,0, 400, 800), limits = c(-410, 1150))+
  # ylim(-600,600)+
  # xlim(-410, 1150)+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  # stat_ellipse(aes(x, y, col=Tissue0), type = "t") +
  # stat_ellipse(type = "t")+ #type = "norm", linetype = 2
  geom_vline(aes(xintercept=quantile(PCA_Placenta_M$ind$coord[, 1], 0.75)+ 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])),
             linetype=4, colour="black")+
  # geom_vline(aes(xintercept=quantile(PCA_Placenta_M$ind$coord[, 1], 0.25)- 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])),
  #            linetype=4, colour="black")+
  geom_text(aes(x=155 , y=-600, label="Q3+1.5*IQR"), vjust=-0.4, hjust=0, size=3) + #angle=90,, "Q3+3*IQR")
  # geom_text(aes(x=-270 , y=-600, label="Q1-1.5*IQR"), vjust=-0.4, hjust=0, size=2) + #angle=90,, "Q3+3*IQR")
  theme(text = element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        axis.text.x =element_text(size=9, colour = "black"),
        axis.text.y =element_text(size=9, colour = "black"),
        legend.position="top")

# calculte Q3+1.5IQR
# quantile(PCA_Placenta_M$ind$coord[, 1], 0.75)+ 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])



tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", "allPCA6clusters", "_joinedGEOsamples408.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)
autoplot(pam(PCA551[,1:2],6), frame =TRUE, frame.type='norm') + labs(x = "PC1 (25.35%)", y="PC2 (10.20%)")
dev.off()

tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised", "allPCA551", "_joinedGEOsamples408.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)
# theme_set(theme_bw())
PCA551plot
dev.off()

# ## T1 sampels
# PCAT1 <- ggplotPCA(PCA_Placenta_M = PCA_PlacentaT1_M, Matched_pd = PhenodataPlacentaT1, value = "all_T1_residuals")
# ## Term sampels
# PCATerm <- ggplotPCA(PCA_Placenta_M = PCA_PlacentaTerm_M, Matched_pd = PhenodataPlacentaTerm, value = "all_Term_residuals")


# 1st 2nd trimester & term samples
PCA551_First <- PCA551[PCA551$Trimester%in%"First",]
PCA551_Second <- PCA551[PCA551$Trimester%in%"Second",]
PCA551_Term <- PCA551[PCA551$Trimester%in%"Term",]

## T1 

PCAT1plot <- PCA551_First %>%
  ggplot(aes(x, y)) + 
  geom_point(aes(shape = Trimester, col=TissueNoGreen), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("Decidua"="#984EA3", "Maternal whole blood"="#999999", 
                                               "Placenta"="#D95F02", "Impure placenta"="#1B9E77"))+ #, "Impure placenta"="#1B9E77"
  scale_shape_manual(values=0)+
  labs(x="PC1", y="PC2")+
  ylim(-600,600)+xlim(-410, 1150)+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  stat_ellipse(aes(x, y,col=Tissue), type = "norm", level = 0.9) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        legend.position="top")

PCAT2plot <- PCA551_Second %>%
  ggplot(aes(x, y)) + 
  geom_point(aes(shape = Trimester, col=TissueNoGreen), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("Placenta"="#D95F02", "Impure placenta"="#1B9E77"))+ #, "Impure placenta"="#1B9E77"
  scale_shape_manual(values=1)+
  labs(x="PC1", y="PC2")+
  ylim(-600,600)+xlim(-410, 1150)+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  stat_ellipse(aes(x, y,col=Tissue), type = "norm", level = 0.9) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        legend.position="top")

PCA551_Term$Tissue1 <- str_replace_all(PCA551_Term$Tissue, c("Impure placenta"="Placenta", 
                                                             "Amnion|Chorion|Decidua"="Impure placenta"))

# PCA551_Term$Tissue1[PCA551_Term$SampleNames%in%"GSM1947213_6008581028_R01C01"] <- "Impure placenta"

PCATermplot <- PCA551_Term %>%
  ggplot(aes(x, y)) + 
  geom_point(aes(shape = Trimester, col=TissueNoGreen), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("Amnion"="#000000", "Chorion"="#BD0026", 
                                               "Decidua"="#984EA3",  
                                               "Placenta"="#D95F02", "Impure placenta"="#1B9E77", 
                                               "Umbilical cord blood"="#F781BF"))+
  scale_shape_manual(values=5)+
  labs(x="PC1", y="PC2")+
  ylim(-600,600)+xlim(-410, 1150)+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  stat_ellipse(data =PCA551_Term, aes(x, y,col=Tissue1), 
               type = "norm", show.legend =FALSE, inherit.aes = FALSE, level = 0.9) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        legend.position="top")

autoplot(pam(PCA551_Term[,1:2],6), frame =TRUE, frame.type='norm')

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/T1Samples_PCA551_First.tiff",
     width = 6, height = 4, units = "cm", res = 300)
# theme_set(theme_pubr())
PCAT1plot
dev.off()

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/T1Samples_PCA551_First.tiff",
     width = 6, height = 4, units = "cm", res = 300)
# theme_set(theme_pubr())
PCAT2plot
dev.off()

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/TermSamples_PCA551_First.tiff",
     width = 20, height = 12, units = "cm", res = 300)
# theme_set(theme_bw())
PCATermplot
dev.off()

# use ggarange to put plots together
PCAs <- ggarrange(PCA551plot, PCAT1plot, PCAT2plot, PCATermplot, 
                  labels = c("A", "B", "C", "D"),
                  ncol = 4, nrow = 1, 
                  font.label=list(size = 12, color = "black", face = "bold", family = "Arial"),
                  legend="top", common.legend = TRUE)


tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/T1Samples_PCA551_4plots.tiff",
     width = 17, height = 8, units = "cm", res = 300)
# theme_set(theme_pubr())
PCAs 
dev.off()


# rbind df and plot them together by facet_wrap
PCA551_First$Tissue1 <- PCA551_First$Tissue
PCA551_Second$Tissue1 <- PCA551_Second$Tissue

All551rbind <- rbind(PCA551_First, PCA551_Second, PCA551_Term)

PCAallplot <- All551rbind %>%
  ggplot(aes(x, y)) + 
  geom_point(aes(shape = Trimester, col=TissueNoGreen), alpha = 1, size = 2)+ 
  scale_color_manual(name="Tissue", values = c("Amnion"="#000000", "Chorion"="#BD0026", 
                                               "Decidua"="#984EA3",  
                                               "Placenta"="#D95F02", "Impure placenta"="#1B9E77", 
                                               "Umbilical cord blood"="#F781BF","Maternal whole blood"="#999999"))+
  scale_shape_manual(values=c(0,1,5))+
  labs(x="PC1", y="PC2")+
  scale_y_continuous(breaks=c(-600,-300,0, 300, 600), limits = c(-600, 800))+
  scale_x_continuous(breaks=c(-400,0, 400, 800), limits = c(-410, 1150))+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  stat_ellipse(data =All551rbind, aes(x, y,col=Tissue1), 
               type = "norm", show.legend =FALSE, inherit.aes = FALSE, level = 0.9) +
  facet_wrap(~Trimester)+
  theme(text = element_text(size=9),
        axis.text.x =element_text(size=9, colour = "black"),
        legend.title=element_text(size=6),
        legend.text=element_text(size=6),
        legend.position="top")

saveRDS(PCA551, file = "/home/qianhui/DNAme/Process_decidual/figures/PCA408df.rds")

saveRDS(PCA551plot, file = "/home/qianhui/DNAme/Process_decidual/figures/T1TermSamples_PCA408samples.rds")
saveRDS(PCAT1plot, file = "/home/qianhui/DNAme/Process_decidual/figures/T1TermSamples_PCAT1plot.rds")
saveRDS(PCAT2plot, file = "/home/qianhui/DNAme/Process_decidual/figures/T1TermSamples_PCAT2plot.rds")
saveRDS(PCATermplot, file = "/home/qianhui/DNAme/Process_decidual/figures/T1TermSamples_PCATermPlot.rds")

saveRDS(PCAs, file = "/home/qianhui/DNAme/Process_decidual/figures/T1TermSamples_PCA408_3plots.rds")
saveRDS(PCAallplot, file = "/home/qianhui/DNAme/Process_decidual/figures/PCAallplot.rds")

