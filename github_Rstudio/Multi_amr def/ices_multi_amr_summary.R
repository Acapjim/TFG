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

input_amr_dataset = paste(dir, "amr_ices_encontrados.tsv", sep = "")

amr_ices_found_dataset = read.delim(input_amr_dataset, sep = "\t")



# Paso 1: agrupar por contig, iceBegin y Stop_of_most_downstream_SP
# Paso 2: quedarnos solo con los grupos que tienen más de una fila
# Paso 3: concatenar los valores de item separados por comas

# Crear la tabla Stop_of_most_downstream_SP
resultado <- amr_ices_found_dataset %>%
  group_by(contig, Start_of_most_upstream_SP, Stop_of_most_downstream_SP) %>%
  filter(n() > 1) %>%  # Solo grupos con más de una fila
  ungroup() %>%
  distinct(contig, Start_of_most_upstream_SP, Stop_of_most_downstream_SP, ICE_consensus_superfamily_SP_conj_module)  # Eliminar duplicados

# Asegurarse que los contigs de columnas estén como texto
resultado$contig <- as.character(resultado$contig)
resultado$Start_of_most_upstream_SP <- as.character(resultado$Start_of_most_upstream_SP)
resultado$Stop_of_most_downstream_SP <- as.character(resultado$Stop_of_most_downstream_SP)

amr_ices_found_dataset$contig <- as.character(amr_ices_found_dataset$contig)
amr_ices_found_dataset$Start_of_most_upstream_SP <- as.character(amr_ices_found_dataset$Start_of_most_upstream_SP)
amr_ices_found_dataset$Stop_of_most_downstream_SP <- as.character(amr_ices_found_dataset$Stop_of_most_downstream_SP)
amr_ices_found_dataset$Gene.symbol <- as.character(amr_ices_found_dataset$Gene.symbol)

# Crear la columna "elementos" en resultado
resultado$elementos <- NA

# Recorrer cada fila de resultado
for (i in 1:nrow(resultado)) {
  # Filtrar amr_ices_found_dataset con coincidencias exactas en contig, Start_of_most_upstream_SP y Stop_of_most_downstream_SP
  coincidencias <- amr_ices_found_dataset %>%
    filter(contig == resultado$contig[i],
           Start_of_most_upstream_SP == resultado$Start_of_most_upstream_SP[i],
           Stop_of_most_downstream_SP == resultado$Stop_of_most_downstream_SP[i])
  
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
  group_by(ICE_consensus_superfamily_SP_conj_module, elementos_ordenados) %>%
  summarise(n = n()) %>%  
  ungroup() %>%
  rename(elementos = elementos_ordenados)
write.table(resultado_final, file = "mdr_ices_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

