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

input_amr_dataset = paste(dir, "amr_phage_encontrados.tsv", sep = "")

amr_phage_found_dataset = read.delim(input_amr_dataset, sep = "\t")


# Seleccionar columnas por nombre directamente con select()
phage_amr_gene_dataset <- amr_phage_found_dataset %>% select(Gene.symbol, Subclass,taxonomy,scaffold)


# Agregar los nombres en listas
taxonomy_agrupados <- aggregate(taxonomy ~ Gene.symbol + Subclass, data = phage_amr_gene_dataset, FUN = function(x) paste(x, collapse = ", "))

# Contar los contigs para cada grupo
scaffold_agrupados <- aggregate(scaffold ~ Gene.symbol + Subclass, data = phage_amr_gene_dataset, FUN = function(x) length(unique(x)))

# Unir los dos dataframes por "nombre" y "subclase"
gene_amr_dataset <- merge(taxonomy_agrupados, scaffold_agrupados, by = c("Gene.symbol", "Subclass"))

# Crear la columna final con el formato deseado
gene_amr_dataset$Resultado <- paste(gene_amr_dataset$Gene.symbol, gene_amr_dataset$Subclass, gene_amr_dataset$taxonomy, gene_amr_dataset$scaffold, sep = " ")


# Eliminar la columna 'Resultado' utilizando subset()
gene_amr_dataset <- subset(gene_amr_dataset, select = -Resultado)

# Eliminar duplicados
gene_amr_dataset$taxonomy <- sapply(gene_amr_dataset$taxonomy, function(x) {
  # Dividir la cadena en una lista
  taxonomy_list <- strsplit(x, ", ")[[1]]
  
  # Eliminar los duplicados
  unique_taxonomy <- unique(taxonomy_list)
  
  # Volver a unir la lista sin duplicados en una cadena separada por comas
  paste(unique_taxonomy, collapse = ", ")
})

# Mostrar el dataframe actualizado
write.table(gene_amr_dataset, file = "phage_amr_gene_summary.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


