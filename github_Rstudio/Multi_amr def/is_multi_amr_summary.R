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

input_amr_dataset = paste(dir, "Amr_is_validados.tsv", sep = "")

amr_is_found_dataset = read.delim(input_amr_dataset, sep = "\t")



# Paso 1: agrupar por seqID, isBegin y isEnd
# Paso 2: quedarnos solo con los grupos que tienen más de una fila
# Paso 3: concatenar los valores de item separados por comas

# Crear la tabla isEnd
resultado <- amr_is_found_dataset %>%
  group_by(seqID, isBegin, isEnd) %>%
  filter(n() > 1) %>%  # Solo grupos con más de una fila
  ungroup() %>%
  distinct(seqID, isBegin, isEnd, family)  # Eliminar duplicados

# Asegurarse que los seqIDs de columnas estén como texto
resultado$seqID <- as.character(resultado$seqID)
resultado$isBegin <- as.character(resultado$isBegin)
resultado$isEnd <- as.character(resultado$isEnd)

amr_is_found_dataset$seqID <- as.character(amr_is_found_dataset$seqID)
amr_is_found_dataset$isBegin <- as.character(amr_is_found_dataset$isBegin)
amr_is_found_dataset$isEnd <- as.character(amr_is_found_dataset$isEnd)
amr_is_found_dataset$Gene.symbol <- as.character(amr_is_found_dataset$Gene.symbol)

# Crear la columna "elementos" en resultado
resultado$elementos <- NA

# Recorrer cada fila de resultado
for (i in 1:nrow(resultado)) {
  # Filtrar amr_is_found_dataset con coincidencias exactas en seqID, isBegin y isEnd
  coincidencias <- amr_is_found_dataset %>%
    filter(seqID == resultado$seqID[i],
           isBegin == resultado$isBegin[i],
           isEnd == resultado$isEnd[i])
  
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
  group_by(family, elementos_ordenados) %>%
  summarise(n = n()) %>%  
  ungroup() %>%
  rename(elementos = elementos_ordenados)

write.table(resultado_final, file = "mdr_is_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
