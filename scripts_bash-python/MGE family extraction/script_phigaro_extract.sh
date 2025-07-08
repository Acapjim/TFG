!/bin/bash

directorio=/home/adria_tfg/tfg/data/Mobile_elements/phigaro_results/phigaro_tsv_kept

for file in "$directorio"/*.tsv
do
python3 /home/adria_tfg/tfg/scripts/script_phigaro_extract.py "$file"
done
