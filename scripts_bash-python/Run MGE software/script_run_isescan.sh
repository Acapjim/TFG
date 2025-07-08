#!/bin/bash

cd /home/adria_tfg/tfg/data/copia_genomas/isescan_batch3

for genoma in *fna
do
isescan.py --seqfile "$genoma" --output /home/adria_tfg/tfg/data/isescan_resultsP/"${genoma}_isescan" --nthread 6
done

