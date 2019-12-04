#!/bin/bash/

# read each line in ourStudy10samples.txt
for i in $(cat /home/qianhui/DNAme/Process_decidual/RDSfiles/ourStudy10samples.txt);
do

# change to each our EPIC array directory
cd /home/qianhui/DNAme/EPICarrayData/idat_files

# copy all the 10 IDAT files in the directory to IDAT_all folder
cp ${i}* /home/qianhui/DNAme/Process_decidual/GEO_dataSets/IDAT_all

done

