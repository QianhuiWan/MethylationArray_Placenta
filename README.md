# MethylationArray_Placenta
Processing placental DNA methylation data from Illumina methylation array.


# For analysing placental outliers, GEO data analysis were used.
I use the publica data and 10 samples from our lab to find placental outliers based on DNA methylation.

# Steps of processing

## Part 1

I mainly dowloaded the data sets from GEO and read raw IDAT files to RGchannel sets.

* I get raw IDAT files and meta data files for 15 GEO data sets (DNA methylation arrays). All these data sets contain placenta samples, some have other tissue types as well. These information were listed in meta data files.

* I modified all the meta data files, so they have same colum names (GEO_accession:Slide). I used `rbind` to join all the meta data files into one file.

* I subseted data for uncomplicated placenta, Amnoin, Chorion, Decidua, maternal whole blood, placental mesenchyme, placental trophoblast, umbilical cord and umbilical cord blood and read these data seperatly using `minfi`. 

* Saved RGchannel sets

## Part 2

I did quality contorl to look at whether there are any failed smaples, so I can remove them.

* I wrote the function to plot QC plots:
    + barplot for mean detection P-values across all samples to identify any failed samples
    + QC plot from `getQC` function from `minfi` package
    + Strip plots for contorl probes
  
* Then plot QC plots for all RGchannel sets using for loop (RGceannel sets by tissue type). 551 samples in total, all sampels are uncomplicated samples.

## Part 3

To look at batch effects using those contorl probes on 450k an EPIC arrays, I get the contorl probes from RGchannel sets and did PCA.

* I wrote a functon for getting control probe intensities from RGchannel sets, calculated beta and M values for control probes and did PCA (`FactoMineR`). Then I used this function for doing PCA for each tissue type.

* Row bind all control probes (for all 551 samples) to plot one PCA plot.

* Saved PCA results for 408 saples.

## Part 4

Filter unwanted probes.

* Wrote a function to filter diiffernt types of unwanted probes (including detectionP>0.01, BeadCount<3 in 95% samples, cross-reactive probes, porbes related with SNPs, and probes on X, Y chr) and correct bacground use `preprocessENmix`.

* Use function for samples of 6 tissue types.

* Plot 1 PCA plot for 408 smaples all together.

* Saved all filtered data sets: MethlSets, and all beta values for 408 samples in one matrix.

## Part 5

I did normalisation and batch correction for data sets

* I tried 4 ways to normalie and correct data:
    * 1. First used `BMIQ_Horvath.R` which used BMIQ and then calibrated beta values to a gold standard to control batch effects. 
    * 2. Also tried to use QN1 (quantile normalisation, between array normalisaiton) before BMIQ normalisation.
    * 3. Then tried normalising with BMIQ and then correct batch with Combat
    * 4. Also tried normalising with BMIQ and then regress out batch (PC1 and PC2 from control probes) to get residuls.
    
    1 ,2 and 3 are not correcting batch effect very well: 
    
    1--there is not one study that contain all tissue types
    2--QN1 doesn't eliminate batch effects
    3--Combat can remove "Trimester" variable, but can't remove "Study" batch ("Study" is contained/confounded with "tissue type")
    4--this method use the linear model to remove batch, and use residuals generate PCA plots.

* Ploted PCA plot for normalised and corrected data

* Compraring these 3 ways of data correction, 
    * I need to find a good gold standard for Horvath's calibration method. I used GSE74738, because it has all these 9 tissue types, but the maternal blood samples in this data set didn't have information of gestational age or trimester. I have filtered samples without trimester infomations in previous step, so I'm not sure I should add these maternal blood samples back or not.
    * QN1 doesn't eliminate batch effects (if normalised all samples together, haven't try to use QN1 seperately according to tissue types). 
    * Combat is more rely on a good expeimental design, because if the unwanted variables and wanted variables were in the same group, this method will remove them all.
    * Since PC1 and PC2 from control probes are unwanted variables, the result form this method for fine for PCA. However, this method could also over processed the data (see Vegard Nygaard's paper, Biostatistics. 2016).
    
## Part 6

We observed 5 outliers from 1st trimester, 4 from term and 2 potential outliers (one from term and one from second trimester), so we check samples from frist, second trimester and term respectively.

Since we saw outliers were seperated with other samples on PC1 of PCA plots, we want to look at the probes contributing to PC1. DNA methylaiton of these probes should be most different between outliers and other pure placenta smaples.

* PCA for samples from frist, second trimester and term

* boxplot (outliers and other samples) for top 2% probes that significantly contributing to PC1


## Part 7

Show DNA methylation of outliers in PMDs (Parcially methylated domain).

* Subset probes for PMD, calculate mean beta values in 10kb bins.

* Gviz plot to show DNA methylation of outliers (frist, second trimester and term) in PMDs (Chr21)

## Part 8

Show DNA methylation of outliers in placenta specific ICRs (Imprinting contorl regions)

* Gviz plot to show DNA methylation of outliers (first second trimester and term) in ICRs

