#!/bin/bash/
cd /home/qianhui/DNAme/Process_decidual/GEO_dataSets

# read each line in GEO_GSEids.txt
for i in $(cat /home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_GSEids.txt);
do

# change to each GSExx directory
cd /home/qianhui/DNAme/Process_decidual/GEO_dataSets/${i}

# copy all the .tar files in the directory to IDAT_all folder
cp *.idat /home/qianhui/DNAme/Process_decidual/GEO_dataSets/IDAT_all

done

