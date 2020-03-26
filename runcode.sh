#!/bin/bash

module load gcc/7.3.0 r/3.6.0
module load udunits
module load gdal

for (( i=0; i<100; ++i)); do
    Rscript /home/jpyanez/mader/mb_lightweight/code_noinstall.r
    sleep 1
done
