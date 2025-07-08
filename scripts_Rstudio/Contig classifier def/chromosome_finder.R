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


# Cargar los archivos
# Leer archivo de texto
file1_lines <- readLines("accesion_contig_chromosome_list.txt")

# Filtrar solo las líneas que terminan con "_chromosome"
complete_lines <- file1_lines[grepl("_chromosome$", file1_lines)]

# Quitar el sufijo "_chromosome"
complete_lines_clean <- sub("_chromosome$", "", complete_lines)

# Extraer la primera palabra (el código) de cada línea
codes_to_complete <- sub("^(\\S+).*", "\\1", complete_lines_clean)

# Leer el archivo TSV
file2_df <- read.delim("compots_family.tsv", stringsAsFactors = FALSE)

# Verificar que la columna "contig" existe
if (!"seqID" %in% colnames(file2_df)) {
  stop("La columna 'contig' no se encuentra en el archivo TSV.")
}

# Añadir "_complete" a los nombres que están en la lista
file2_df$seqID <- ifelse(file2_df$seqID %in% codes_to_complete,
                           paste0(file2_df$seqID, "_chromosome"),
                           file2_df$seqID)

# Guardar el resultado en un nuevo archivo
write.table(file2_df, "compots_family_chromosomed.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

