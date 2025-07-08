#!/bin/bash

directorio=/home/adria_tfg/tfg/data/Mobile_elements/IS/isescan_tsv_kept/

for file in "$directorio"/*.tsv
do
python3 /home/adria_tfg/tfg/scripts/script_isescan_extract.py "$file"
done
