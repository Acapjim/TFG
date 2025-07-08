#!/bin/bash
cd /home/adria_tfg/tfg/data/genomasP

for genoma in *fna
do 
phigaro -f "$genoma" -o /home/adria_tfg/tfg/data/phigaro_result/"$genoma" -t 6
done

echo "Fin"


