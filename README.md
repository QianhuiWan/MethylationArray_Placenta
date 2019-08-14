# MethylationArray_Placenta

Identify potential mixed placenta samples (placenta samples mixed with other tissue types) by processing placental DNA methylation data from Illumina methylation array.


# For analysing placental outliers (potential mixed placenta samples), GEO data analysis were used.
I use the publica data and 10 samples from our lab to find placental outliers based on DNA methylation data.

# Steps of processing

## Part 1

Downloaded the data sets from GEO and read raw IDAT files to RGchannel sets.

* Get raw IDAT files and meta data files for 13 GEO data sets (DNA methylation arrays). All these data sets contain placenta samples, some have other tissue types as well. This information was listed in meta data files.

* Modified all the meta data files, so they have same column names (GEO_accession:Slide). I used `rbind` to join all the meta data files into one file.

* Subseted data for uncomplicated placenta, Amnion, Chorion, Decidua, maternal whole blood and umbilical cord blood and read these data separately using `minfi`. 

* Saved RGchannel sets

## Part 2

I did quality control for all 408 samples to look at whether there are any failed samples, so I can remove them.

* plot QC plots:

    + barplot for mean detection P-values across all samples to identify any failed samples
    + QC plot from `getQC` function from `minfi` package
    + Strip plots for control probes
  
* Then plot QC plots for all 408 samples, all samples are from uncomplicated pregnancies.

## Part 3

To have a overview of batch effects using those control probes on 450k an EPIC arrays, PCA were applied to control probes from RGchannel sets.

* Get control probe intensities from RGchannel sets, calculated beta and M values for control probes and did PCA (`FactoMineR`). 

* Saved PCA results for 408 samples.

## Part 4

Filtration of unwanted probes followed with background and dye bias correction and normalisation.

* Different types of unwanted probes (including detectionP>0.01, BeadCount<3 in 95% samples, cross-reactive probes, probes related with SNPs, and probes on X, Y chromosomes) were filtered 

* Background noise and dye bias were corrected using `preprocessENmix` package.

* Normalisation and batch correction for data sets

    I tried 4 ways to normalise and correct data:
    
    * 1. First used `BMIQ_Horvath.R` which used BMIQ and then calibrated beta values to a gold standard to control batch effects. 
    * 2. QN1 (quantile normalisation, between array normalisation) before BMIQ normalisation.
    * 3. Normalising with BMIQ and then correct batch with Combat
    * 4. Normalising with BMIQ and then regress out batch (PC1 and PC2 from control probes), and use residuals for further analyses.
    
    Comparing these 4 ways of data correction:
    
    * I need to find a good gold standard for Horvath's calibration method. I used GSE74738, because it has all these 6 tissue types. However, the maternal blood samples in this data set (GSE74738) didn't have information of gestational age or trimester, so I have filtered samples out in previous step. GSE74738 doesn't contain samples from 6 tissue types after previous filtering, so we can't use Horvath's calibration method.
    
    * QN1 doesn't eliminate batch effects (if normalised all samples together, haven't try to use QN1 separately according to tissue types because it's time consuming and we prefer to pre-process all samples together). 
    
    * Combat is more rely on a good experimental design, because if the unwanted variables and wanted variables were in the same group, this method will remove all the variables.
    
    * PC1 and PC2 from control probes (not pre-processed) could represent unwanted variables, but we can't adjust unwanted variables in processed data using unprocessed control probes. This method could also over processed the data (see Vegard Nygaard's paper, Biostatistics. 2016).
  
    In summary, we used the normalised data (without adding extra step for batch correction) for the following analyses.

## Part 5

  * Plotted PCA plot for normalised data
  
  We observed 5 outliers from 1st trimester, 4 from term and 2 potential outliers (one from term and one from second trimester), so we check samples from first, second trimester and term respectively.
  
  * Plot supplementary figure (result of sample clustering)

## Part 6

* boxplot (mixed and pure placenta samples) for top 2% probes that significantly contributing to PC1 (delta beta > 0.2)

* generate Figure 1 (containing PCA plots and boxplots)


## Part 7

Show DNA methylation of mixed and pure placenta samples in PMDs (Partially methylated domain).

* Subset probes for PMD, calculate mean beta values in 10kb bins.

* Gviz plot to show DNA methylation of outliers (first, second trimester and term) in PMDs (Chr21)

## Part 8

Show DNA methylation of mixed and pure placenta samples in placenta specific ICRs (Imprinting control regions)

* Gviz plot to show DNA methylation of outliers (first second trimester and term) in ICRs

* Plot supplementary figure (difference of DNA methylation between mixed and pure placenta samples at all placenta-specific ICRs)

If you require any further information, feel free to contact Qianhui (Email: qianhui.wan@adelaide.edu.au or wanqianhui@outlook.com).
