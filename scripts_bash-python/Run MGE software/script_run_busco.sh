#!/bin/bash

lineage=enterococcus_odb12

cd data/genomasP/

for genoma in *fna
do
busco -m genome -i "$genoma" -o data/genomas_busco/"${genoma}_busco"  -l "$lineage"
done


