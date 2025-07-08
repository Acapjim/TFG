########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")


library(dplyr)
library(gtsummary)
library(tidyverse)
library(plyr)
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


# Seleccionar columnas por nombre directamente con select()
IS_amr_gene_dataset <- amr_is_found_dataset %>% select(Gene.symbol, Subclass,family,seqID)

# Agregar los nombres en listas
family_agrupados <- aggregate(family ~ Gene.symbol + Subclass, data = IS_amr_gene_dataset, FUN = function(x) paste(x, collapse = ", "))

# Contar los contigs para cada grupo
seqID_agrupados <- aggregate(seqID ~ Gene.symbol + Subclass, data = IS_amr_gene_dataset, FUN = function(x) length(unique(x)))

# Unir los dos dataframes por "nombre" y "subclase"
gene_amr_dataset <- merge(family_agrupados, seqID_agrupados, by = c("Gene.symbol", "Subclass"))

# Crear la columna final con el formato deseado
gene_amr_dataset$Resultado <- paste(gene_amr_dataset$Gene.symbol, gene_amr_dataset$Subclass, gene_amr_dataset$family, gene_amr_dataset$seqID, sep = " ")


# Eliminar la columna 'Resultado' utilizando subset()
gene_amr_dataset <- subset(gene_amr_dataset, select = -Resultado)


# Eliminar duplicados 
gene_amr_dataset$family <- sapply(gene_amr_dataset$family, function(x) {
  # Dividir la cadena en una lista
  family_list <- strsplit(x, ", ")[[1]]
  
  # Eliminar los duplicados
  unique_family <- unique(family_list)
  
  # Volver a unir la lista sin duplicados en una cadena separada por comas
  paste(unique_family, collapse = ", ")
})

# Mostrar el dataframe actualizado

write.table(gene_amr_dataset, file = "IS_amr_gene_summary.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


