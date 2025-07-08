#!/bin/bash

directorio=/home/adria_tfg/tfg/data/Mobile_elements/icescreen_results/icescreen_adapted_tsv/

for file in "$directorio"/*.tsv
do
python3 /home/adria_tfg/tfg/scripts/script_icescreen_extract.py "$file"
done
