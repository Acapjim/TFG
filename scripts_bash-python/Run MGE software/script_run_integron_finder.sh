#!/bin/bash
cd /home/adria_tfg/tfg/data/genomas_taken
contador=0

for genoma in *fna
do
integron_finder "$genoma" --local-max --func-annot --outdir /home/adria_tfg/tfg/data/Mobile_elements/integron_finder_results_trial
((contador++))

echo "Archivos procesados: $contador"

done

echo "Fin"
