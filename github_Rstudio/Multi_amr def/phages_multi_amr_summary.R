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

input_amr_dataset = paste(dir, "amr_phage_encontrados.tsv", sep = "")

amr_phage_found_dataset = read.delim(input_amr_dataset, sep = "\t")



# Paso 1: agrupar por scaffold, begin y end
# Paso 2: quedarnos solo con los grupos que tienen más de una fila
# Paso 3: concatenar los valores de item separados por comas

# Crear la tabla end
resultado <- amr_phage_found_dataset %>%
  group_by(scaffold, begin, end) %>%
  filter(n() > 1) %>%  # Solo grupos con más de una fila
  ungroup() %>%
  distinct(scaffold, begin, end, taxonomy)  # Eliminar duplicados


# Asegurarse que los scaffolds de columnas estén como texto
resultado$scaffold <- as.character(resultado$scaffold)
resultado$begin <- as.character(resultado$begin)
resultado$end <- as.character(resultado$end)

amr_phage_found_dataset$scaffold <- as.character(amr_phage_found_dataset$scaffold)
amr_phage_found_dataset$begin <- as.character(amr_phage_found_dataset$begin)
amr_phage_found_dataset$end <- as.character(amr_phage_found_dataset$end)
amr_phage_found_dataset$Gene.symbol <- as.character(amr_phage_found_dataset$Gene.symbol)

# Crear la columna "elementos" en resultado
resultado$elementos <- NA

# Recorrer cada fila de resultado
for (i in 1:nrow(resultado)) {
  # Filtrar amr_phage_found_dataset con coincidencias exactas en scaffold, begin y end
  coincidencias <- amr_phage_found_dataset %>%
    filter(scaffold == resultado$scaffold[i],
           begin == resultado$begin[i],
           end == resultado$end[i])
  
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
  group_by(taxonomy, elementos_ordenados) %>%
  summarise(n = n()) %>%  
  ungroup() %>%
  rename(elementos = elementos_ordenados)

write.table(resultado_final, file = "mdr_phages_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

