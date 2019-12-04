#!/bin/bash/
cd /home/qianhui/DNAme/Process_decidual/GEO_dataSets

# read each line in GEO_GSEids.txt
for i in $(cat /home/qianhui/DNAme/Process_decidual/RDSfiles/GEO_GSEids.txt);do

 # change to each GSExx directory
 cd /home/qianhui/DNAme/Process_decidual/GEO_dataSets/${i}
 
 # count the unmber of tar files in the folder
 TarFileNum=$(ls *.tar -l | wc -l)

 # untar files if there are >=1 tar file in the folder
 if [ $TarFileNum -ge 1 ]; then 
 # untar all the .tar files in the directory
    ls *.tar |xargs -n1 tar -xvf
 else
    echo 'no tar files in this folder'
 fi

done

