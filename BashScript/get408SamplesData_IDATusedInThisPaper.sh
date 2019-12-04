#!/bin/bash/

# read each line in ourStudy10samples.txt
for i in $(cat /Users/qianhuiwan/R_projects/MethylationArray_Placenta/BashScript/all408sampleNames.txt);
do

# change to each our array data directory
cd /Volumes/QianhuiDrive/Freya/DNAme/ProcessWF/Process_decidual/GEO_dataSets/IDAT_all

# copy all the 10 IDAT files in the directory to IDAT_all folder
cp ${i}* /Users/qianhuiwan/R_projects/MethylationArray_Placenta/RawDataFromGEO_408/IDAT_408

done

# Then zip this folder
cd /Users/qianhuiwan/R_projects/MethylationArray_Placenta/RawDataFromGEO_408
zip -r -X IDAT_408.zip IDAT_408

# Copy the csv file (meta data) to RawDataFromGEO_408 direcroty as well 
cd /Volumes/QianhuiDrive/Freya/DNAme/ProcessWF/Process_decidual/GEO_dataSets/IDAT_all
cp *.csv /Users/qianhuiwan/R_projects/MethylationArray_Placenta/RawDataFromGEO_408
