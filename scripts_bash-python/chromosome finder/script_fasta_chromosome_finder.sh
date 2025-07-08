#!/bin/bash


directorio=/home/adria_tfg/tfg/data/genomas_sin_duplicaciones

for file in "$directorio"/*.fna
do
python3 /home/adria_tfg/tfg/scripts/script_fasta_chromosome_finder.py "$file"
done
