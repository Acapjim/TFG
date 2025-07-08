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

input_amr_dataset = paste(dir, "amr_ices_encontrados.tsv", sep = "")

amr_ices_found_dataset = read.delim(input_amr_dataset, sep = "\t")


# Seleccionar columnas por nombre directamente con select()
ices_amr_gene_dataset <- amr_ices_found_dataset %>% select(Gene.symbol, Subclass,ICE_consensus_superfamily_SP_conj_module,contig)

# Agregar los nombres en listas
family_agrupados <- aggregate(ICE_consensus_superfamily_SP_conj_module ~ Gene.symbol + Subclass, data = ices_amr_gene_dataset, FUN = function(x) paste(x, collapse = ", "))

# Contar los contigs para cada grupo
seqID_agrupados <- aggregate(contig ~ Gene.symbol + Subclass, data = ices_amr_gene_dataset, FUN = function(x) length(unique(x)))

# Unir los dos dataframes por "nombre" y "subclase"
gene_amr_dataset <- merge(family_agrupados, seqID_agrupados, by = c("Gene.symbol", "Subclass"))

# Crear la columna final con el formato deseado
gene_amr_dataset$Resultado <- paste(gene_amr_dataset$Gene.symbol, gene_amr_dataset$Subclass, gene_amr_dataset$family, gene_amr_dataset$seqID, sep = " ")


# Eliminar la columna 'Resultado' utilizando subset()
gene_amr_dataset <- subset(gene_amr_dataset, select = -Resultado)

# Eliminar duplicados 
gene_amr_dataset$family <- sapply(gene_amr_dataset$ICE_consensus_superfamily_SP_conj_module, function(x) {
  # Dividir la cadena en una lista
  family_list <- strsplit(x, ", ")[[1]]
  
  # Eliminar los duplicados
  unique_family <- unique(family_list)
  
  # Volver a unir la lista sin duplicados en una cadena separada por comas
  paste(unique_family, collapse = ", ")
})

# Mostrar el dataframe actualizado
gene_amr_dataset <- subset(gene_amr_dataset, select = -ICE_consensus_superfamily_SP_conj_module)


write.table(gene_amr_dataset, file = "ICES_amr_gene_summary.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


