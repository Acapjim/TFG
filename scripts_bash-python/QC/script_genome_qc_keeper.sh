# este script toma un txt venido de R con los nombres de los genomas y una carpeta con genomas y copia a otra carpeta aquellos que esten en el tsv
#!/bin/bash

# Definir rutas de los archivos y carpetas
txt_file="/home/adria_tfg/tfg/data/QC/QC_Tables/kept_genomes.txt"  # Ruta al archivo .txt
genomas_dir="/home/adria_tfg/tfg/data/genomas_sin_duplicaciones"  # Carpeta con archivos .fna
dest_dir="/home/adria_tfg/tfg/data/genomas_taken"  # Carpeta destino para los archivos copiados

# Crear la carpeta de destino si no existe
mkdir -p "$dest_dir"

# Leer el archivo línea por línea
while IFS= read -r line; do
    # Limpiar la línea para eliminar posibles espacios o saltos de línea extra
    line=$(echo "$line" | tr -d '"')  # Elimina comillas
    line_part=$(echo "$line" | cut -d'_' -f1,2,3)  # Tomar solo la parte clave (antes del primer '_')
    
    # Mostrar la línea clave que se está procesando
    echo "Procesando línea clave: $line_part"
    
    # Buscar archivos .fna que contengan esta parte del nombre en su nombre
    found=false
    for fna_file in "$genomas_dir"/*.fna; do
        # Mostrar el archivo .fna actual
        echo "Verificando archivo: $fna_file"
        
        # Si el nombre del archivo contiene la línea clave, copiarlo a la carpeta destino
        if [[ "$fna_file" == *"$line_part"* ]]; then
            cp "$fna_file" "$dest_dir"
            echo "Copiado: $fna_file"
            found=true
        fi
    done
    
    # Si no se encuentra el archivo, mostrar un mensaje
    if [ "$found" = false ]; then
        echo "No se encontró ningún archivo que coincida con: $line_part"
    fi
done < "$txt_file"

echo "Proceso completado."
