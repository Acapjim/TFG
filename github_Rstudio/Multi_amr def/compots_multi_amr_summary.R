########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")


library(stringr)
library(dplyr)
library(gtsummary)
library(tidyverse)
library(data.table)
library(R.utils)
library(writexl)
library(ape)
library(rjson)
library(readxl)
library(ggplot2)

#loading data 
dir = "W://AMR results/"

input_amr_dataset = paste(dir, "amr_composis_encontrados.tsv", sep = "")

amr_composis_found_dataset = read.delim(input_amr_dataset, sep = "\t")

# Paso 1: agrupar por seqID, Start.x y Stop.x
# Paso 2: quedarnos solo con los grupos que tienen más de una fila
# Paso 3: concatenar los valores de item separados por comas

# Crear la tabla Stop.x
resultado <- amr_composis_found_dataset %>%
  group_by(seqID, Start.x, Stop.x) %>%
  filter(n() > 1) %>%  # Solo grupos con más de una fila
  ungroup() %>%
  distinct(seqID, Start.x, Stop.x, nombre)  # Eliminar duplicados

# Asegurarse que los seqIDs de columnas estén como texto
resultado$seqID <- as.character(resultado$seqID)
resultado$Start.x <- as.character(resultado$Start.x)
resultado$Stop.x <- as.character(resultado$Stop.x)

amr_composis_found_dataset$seqID <- as.character(amr_composis_found_dataset$seqID)
amr_composis_found_dataset$Start.x <- as.character(amr_composis_found_dataset$Start.x)
amr_composis_found_dataset$Stop.x <- as.character(amr_composis_found_dataset$Stop.x)
amr_composis_found_dataset$Gene.symbol <- as.character(amr_composis_found_dataset$Gene.symbol)

# Crear la columna "elementos" en resultado
resultado$elementos <- NA

# Recorrer cada fila de resultado
for (i in 1:nrow(resultado)) {
  # Filtrar amr_composis_found_dataset con coincidencias exactas en seqID, Start.x y Stop.x
  coincidencias <- amr_composis_found_dataset %>%
    filter(seqID == resultado$seqID[i],
           Start.x == resultado$Start.x[i],
           Stop.x == resultado$Stop.x[i])
  
  # Unir los símbolos separados por coma
  simbolos_unidos <- paste(coincidencias$Gene.symbol, collapse = ",")
  
  # Asignar el resultado a la columna "elementos"
  resultado$elementos[i] <- simbolos_unidos
}


# resumimos
resultado_ordenado <- resultado %>%
  mutate(
    elementos_ordenados = sapply(strsplit(elementos, ","), function(x) {
      x <- trimws(x)          # Quitar espacios
      x <- sort(x)            # Ordenar
      paste(x, collapse = ",") # Volver a unir en string
    })
  )

resultado_final <- resultado_ordenado %>%
  group_by(nombre, elementos_ordenados) %>%
  summarise(n = n()) %>%  
  ungroup() %>%
  rename(elementos = elementos_ordenados)

write.table(resultado_final, file = "mdr_compots_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
