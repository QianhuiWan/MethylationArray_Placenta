
# MethylationArray_Placenta

The scripts in this repository were used to identify potential mixed placenta samples (placenta samples mixed with other tissue types) by processing placental DNA methylation data from Infinium Human Methylation 450K/EPIC BeadChip.

# For analysing placental outliers (potential mixed placenta samples), GEO data were used.

I use the publicly available data and 10 samples from our own study to find placental outliers (mixed placenta tissue samples) based on DNA methylation data.

# Steps of processing

## Part 1

Download the data sets from GEO and read raw IDAT files to RGchannel sets.

* Get raw IDAT files and meta data files for 13 GEO data sets (GSE66210, GSE74738, GSE69502, GSE75196, GSE75248, GSE120250, GSE100197, GSE98224, GSE71678, GSE66459, GSE115508, GSE113600 and GSE131945). There are 408 samples in total, including 379 placenta, 11 maternal whole blood, 11 umbilical cord blood, 3 decidua, 2 amnion and 2 chorion. Details of these samples was listed in meta data files (GEOmeta/GEO_phenotypes or NormBMlist$phenoDataAll).

* Modified all the meta data files, so they all have same column names (GEO_accession:Slide). I used `rbind` to join all the meta data files into one file.

* Subset data for uncomplicated placenta, amnion, chorion, decidua, maternal whole blood and umbilical cord blood and read these data separately using `minfi`. 

* Saved RGchannel sets

## Part 2

I did quality control for all 408 samples to look at whether there are any failed samples, so they can be removed.

* plot QC plots:

    + barplot for mean detection P-values across all samples to identify any failed samples
    + QC plot from `getQC` function from `minfi` package
    + Strip plots for control probes
  
* Then plot QC plots for all 408 samples, all samples of first and second trimester are from terminated pregnancies and all term samples are from uncomplicated pregnancies.

## Part 3

To have an overview of batch effects using those control probes on 450K and EPIC arrays, PCA were applied to control probes from RGchannel sets.

* Get control probe intensities from RGchannel sets, calculated beta and M values for control probes and did PCA (`FactoMineR`). 

* Saved PCA results for 408 samples.

## Part 4

Filtration of unwanted probes followed with background and dye bias correction and normalisation.

* Different types of unwanted probes (including detection P>0.01, BeadCount<3 in 95% samples, cross-reactive probes, probes related with SNPs, and probes on X, Y chromosomes) were filtered 

* Background noise and dye bias were corrected using `preprocessENmix` package.

* Normalisation and batch correction for data sets

    I tried 4 ways to normalise and correct data:
    
    * 1. First used `BMIQ_Horvath.R` which used BMIQ and then calibrated beta values to a gold standard to control batch effects. 
    * 2. QN1 (quantile normalisation, between array normalisation) before BMIQ normalisation.
    * 3. Normalising with BMIQ and then correct batch with Combat
    * 4. Normalising with BMIQ and then regress out batch (PC1 and PC2 from control probes), and use residuals for further analyses.
    
    Comparing these 4 ways of data correction:
    
    * I need to find a good gold standard for Horvath's calibration method. I used GSE74738, because it has 6 tissue types. However, the maternal blood samples in this data set (GSE74738) didn't have information on gestational age or trimester, so I have filtered samples out in previous step. GSE74738 doesn't contain samples from 6 tissue types after previous filtering, so we can't use Horvath's calibration method.
    
    * QN1 doesn't eliminate batch effects (if normalised all samples together, haven't tried to use QN1 separately according to tissue types because it's time consuming and we prefer to pre-process all samples together). 
    
    * Combat relies on a good experimental design, because if the unwanted variables and wanted variables were in the same group, this method will remove all the variables.
    
    * PC1 and PC2 from control probes (not pre-processed) could represent unwanted variables, but we can't adjust unwanted variables in processed data using unprocessed control probes. This method could also over process the data (see Vegard Nygaard's paper, Biostatistics. 2016).
  
    In summary, we used the normalised data (without adding extra step for batch correction) for the following analyses, since the statistical results (ANOVA-F test) showed that the different clusters in PCA were more likely to reflect tissue types rather than study batches.

## Part 5

  * Plotted PCA plot for normalised data
  
  We observed 5 outliers from first trimester, 4 from term and 2 potential outliers (one from term and one from second trimester), so we checked samples from first, second trimester and term respectively.
  
  * Plot supplementary figure (result of sample clustering)
  
  * All data (408 samples) were separated into training data (n=285) and test data (n=123) to test the clustering methods used for identifying outliers/potential mixed placenta samples.

## Part 6

* boxplot (mixed and pure placenta samples) for top 2% probes that significantly contributed to PC1 (delta beta > 0.2)

* generate Figure 1 (containing PCA plots and boxplots) 

* test correlation between PC1/PC2 and factors (including Tissue types, Trimester, sex, study) using ANOVA F test.

## Part 7

Show DNA methylation of mixed and pure placenta samples in PMDs (Partially methylated domain).

* Subset probes for PMD, calculate mean beta values in 10kb bins.

* Gviz plot to show DNA methylation of outliers (first, second trimester and term) in PMDs (Chr21)

## Part 8

Show DNA methylation of mixed and pure placenta samples in placenta specific ICRs (imprinting control regions)

* Gviz plot to show DNA methylation of outliers (first, second trimester and term) in ICRs

* Plot supplementary figure (difference of DNA methylation between mixed and pure placenta samples at all placenta-specific ICRs)

## Part 9

Additional quality control (including the check of fetal sex & genotypes) for 408 samples using `ewastools` R package.

* Identify miss labelled fetal sex for all 408 samples. Fetal sex predicated based on intensities from X/Y chromosomes.

* Check the homogeneity of genotypes (using 65 SNP probes to check genotype) for each sample. The suggested cutoff (log odds < -4) for samples contaminated with foreign DNA was used (Heiss JA & Just AC. 2018). 

* The result from `ewastools` were used to compare with the result from our method to see whether mixed placenta samples were identified.

* We probably need to run this part before part5, because the `GEO_phenotypes_ewastoolPurityCheck_output.csv` file generated from part9 were used in part5.

**If you require any further information, feel free to contact Qianhui (Email: qianhui.wan@adelaide.edu.au or wanqianhui@outlook.com).**
