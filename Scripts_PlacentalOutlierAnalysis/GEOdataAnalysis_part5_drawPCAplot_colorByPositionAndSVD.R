
library(tidyverse)
library(magrittr)
library(tibble)
library(FactoMineR)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(cluster)
library(ChAMP)
library(RColorBrewer)
theme_set(theme_pubr())
registerDoParallel(cores=10)

# load data 
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/allSamplesResiduals_ControlCorrected.rds")
# load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/allSamplesResiduals_ControlCorrected_residuals.rds")
# PCA_Placenta_M <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta379_Mval.RData")

## GEOmeta, GEO_phenotypes
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_phenotypes_add10.RData")
# # PCA_Placenta_M, PCA_PlacentaT1_M, PCA_PlacentaTerm_M
load(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta_deciduals_allT1Term_BMval.RData")

NormBMlist <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/NormBMlist.rds")

## annotation needed for using functions
ann450k <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/ann450k_GEO15sets.rds")
annEPIC <- readRDS(file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/annEPIC_GEO15sets.rds")

# all maples PCA
table(colnames(NormBMlist$NormMAll)==NormBMlist$phenoDataAll$Sample_name)
M_all <- NormBMlist$NormMAll#[,NormBMlist$phenoDataAll$Sample%in%"Placenta"]

table(colnames(M_all) %in% GEO_phenotypes$Sample_name)

Matched_pd <- GEO_phenotypes[match(colnames(M_all), GEO_phenotypes$Sample_name),]

# PCA_Placenta_M <- FactoMineR::PCA(base::t(M_placenta))

# saveRDS(PCA_Placenta_M, file = "/home/qianhui/DNAme/Process_decidual/RDSfiles/PCA_placenta379_Mval.RData")

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

PCAplot_all <- new_plotDF %>%
  ggplot(aes(x, y)) + 
  geom_point(aes(shape = Tissue, col=Trimester), alpha = 1, size = 5)+
  # geom_point(aes(shape = Trimester, col=Tissue), alpha = 1, size = 5)+ 
  # scale_color_manual(values = namedPal2)+
  scale_shape_manual(values=1:11)+
  # scale_color_manual(name="Tissue", values = c("#000000", "#BD0026", "#984EA3", "#999999", "#D95F02", "#1B9E77", "#F781BF", "#A65628"))+
  labs(x="PC1", y="PC2")+
  theme(text = element_text(size=20))
# axis.text.x = element_text(angle=90, hjust=1)) 

# save as tiff
tiff(file = paste0("/home/qianhui/DNAme/Process_decidual/figures/MDSplot-Normalised_all", value, "_joinedGEOsamples.tiff"), 
     width = 20,height = 15, units="cm",  res = 300)

# theme_set(theme_bw())
PCAplot_all
dev.off()

# show the plot
# theme_set(theme_bw())
PCAplot_all %>% print()

# return new_plotDF
return(new_plotDF)
}
# end of ggplotPCA function=======================================



# use ggplotPCA function
## all sampoes
PCA_all <- ggplotPCA(PCA_Placenta_M = PCA_Placenta_M, Matched_pd = Matched_pd, value = "allSample")

PCA_placenta_addPostion <- PCA_all %>% 
  dplyr::mutate(Position=str_replace_all(Study, c(
    "GSE69502|GSE100197|GSE115508|GSE98224|ourStudy"="Villous",
    "GSE71678"="Adjacent to cord insertion", #adjacent to cord insertion, #full thickness placenta
    "GSE75248"="2cm from umbilical cord", #2cm below umbilical cord
    "GSE120250"="1cm below chorionic plate", #1cm below chorionic plate
    "GSE75196"="Fetal side of placenta",
    "GSE66459|GSE113600"="Not Placenta"))) %>% #fetal side of placenta 
  dplyr::mutate(Position2 = case_when(((Study=="GSE74738" | Study=="GSE66210") & Tissue=="Placenta") ~ "Villous",
                ((Study=="GSE74738" | Study=="GSE66210") & Tissue!="Placenta") ~ "Not Placenta",
                (Study!="GSE74738" | Study!="GSE66210") ~ Position))

PCA_placentaPlot <- PCA_placenta_addPostion %>%
  ggplot(aes(x, y)) + 
  # geom_point(aes(shape = Tissue, col=Trimester), alpha = 1, size = 5)+ 
  geom_point(aes(shape = Trimester, col=Position2), alpha = 1, size = 2)+ 
  # scale_color_manual(name="Tissue", values = c("#000000", "#BD0026", "#984EA3", "#999999", "#D95F02", "#1B9E77", "#F781BF"))+
  # scale_shape_manual(values=c(0,1,5,2))+
  labs(x="PC1 (25.35%)", y="PC2 (10.20%)")+
  guides(color=guide_legend(nrow = 3), shape=guide_legend(nrow = 2)) +
  # stat_ellipse(aes(x, y, col=Tissue0), type = "t") +
  # stat_ellipse(type = "t")+ #type = "norm", linetype = 2
  geom_vline(aes(xintercept=quantile(PCA_Placenta_M$ind$coord[, 1], 0.75)+ 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])),
             linetype=4, colour="black")+
  geom_vline(aes(xintercept=quantile(PCA_Placenta_M$ind$coord[, 1], 0.25)- 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])),
             linetype=4, colour="black")+
  geom_text(aes(x=155 , y=-600, label="Q3+1.5*IQR"), vjust=-0.4, hjust=0, size=3) + #angle=90,, "Q3+3*IQR")
  geom_text(aes(x=-270 , y=-600, label="Q1-1.5*IQR"), vjust=-0.4, hjust=0, size=3) + #angle=90,, "Q3+3*IQR")
  stat_ellipse(aes(x, y,col=Position2), type = "norm", level = 0.9) +
  theme(text = element_text(size=9),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        axis.text.x =element_text(size=9, colour = "black"),
        axis.text.y =element_text(size=9, colour = "black"),
        legend.position="top")

# calculte Q3+1.5IQR
quantile(PCA_Placenta_M$ind$coord[, 1], 0.75)+ 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])
quantile(PCA_Placenta_M$ind$coord[, 1], 0.25)- 1.5*IQR(PCA_Placenta_M$ind$coord[, 1])

tiff(file = "/home/qianhui/DNAme/Process_decidual/figures/PCA_placentaPlot_position.tiff", 
     width = 20, height = 15, units="cm",  res = 300)
theme_set(theme_pubr())
PCA_placentaPlot
dev.off()

saveRDS(PCA_placentaPlot, file = "/home/qianhui/DNAme/Process_decidual/figures/PCA_placentaPlot.rds")


# svd plot using M values

Mall <- NormBMlist$NormMAll
table(colnames(Mall) %in% NormBMlist$phenoDataAll$Sample_name)

NameWithoutGSE71678 <- NormBMlist$phenoDataAll[!(NormBMlist$phenoDataAll$Study%in%"GSE71678"),]$Sample_name
Mall <- Mall[,colnames(Mall)%in%NameWithoutGSE71678]
  
Matched_pd <- NormBMlist$phenoDataAll[match(colnames(Mall), NormBMlist$phenoDataAll$Sample_name),]

pd <- Matched_pd[, c("Sample_name", "Sample","Trimester", "Fetal_Sex", "ArrayType", "Study")]

# chage to right class for each column

cols2 = c(1:6)
pd[,cols2] <-  lapply(pd[,cols2], function(x){as.character(x)})


# SVD plots_all 408 samples
champ.SVD(beta = Mall, pd = pd,resultsDir = "/home/qianhui/DNAme/Process_decidual/figures/MPlot_allM408_SVD")


# Combat to correct position effects
# table(PCA_placenta_addPostion$SampleNames==NormBMlist$phenoDataAll$Sample_name)
# NormBMlist$phenoDataAll$Position2 <- PCA_placenta_addPostion$Position2

TrimesterAdjusted_Mall <- champ.runCombat(beta = Mall, 
                pd = pd,
                variablename = "Sample",
                batchname = "Trimester",
                logitTrans = FALSE)

champ.SVD(beta = TrimesterAdjusted_Mall, pd = pd, resultsDir = "/home/qianhui/DNAme/Process_decidual/figures/MPlot_allM408_Tadjusted_SVD")


TrimesterAdjusted_MallPCA <- FactoMineR::PCA(t(TrimesterAdjusted_Mall))

source(file = "~/DNAme/Process_decidual/ggplotPCAFUN.R")

TrimesterAdjusted_Mall_pcaRES <- ggplotPCA(PCA_Placenta_M = TrimesterAdjusted_MallPCA, 
                                           Matched_pd = Matched_pd, 
                                           value = "Mval_Tadjusted")

TrimesterAdjusted_Mall_pcaRES %>% 
  ggplot(aes(x, -y)) + 
  geom_point(aes(shape = Tissue, col=Study), alpha = 1, size = 5)+ 
  # scale_color_manual(values = namedPal2)+
  scale_shape_manual(values=1:11)+
  labs(x="PC1", y="PC2")+
  theme(text = element_text(size=20))


df <- TrimesterAdjusted_MallPCA$ind$coord %>% as.data.frame() %>% rownames_to_column(var="Sample_name") %>% 
  left_join(Matched_pd, by="Sample_name") %>% column_to_rownames(var="Sample_name")

a1 <- mclust::Mclust(df[,c("Dim.1","Dim.2")],2) #for 2 groups
table(a1$classification,df$Sample)
with(df,plot(Dim.2 ~ Dim.1,pch=as.numeric(as.factor(Sample)),col=a1$classification))
with(df,legend("topright",levels(as.factor(Sample)),pch=1:6))
