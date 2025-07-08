########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")

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

# Leer el archivo .txt
lines <- readLines("accesion_contig_chromosome_list.txt")
# Leer el archivo TSV con las etiquetas
compois_dataset <- read.delim("pre_composis_family.tsv", stringsAsFactors = FALSE)

# Obtener valores únicos de la columna "nombre"
etiquetas_unicas <- unique(compois_dataset$nombre)

# Crear nuevas columnas vacías (NA) con los nombres de cada etiqueta única
for (etiqueta in etiquetas_unicas) {
  compois_dataset[[etiqueta]] <- NA  
}
compois_dataset$Acc <- NULL

# Este script sumará las ocurrencias de cada grupo de familia por ID

# Primero, agrupamos los datos por 'ID' y contamos la cantidad de veces que aparece cada familia
resultado <- compois_dataset %>%
  group_by(seqID, nombre) %>%
  summarise(cantidad = n(), .groups = "drop")

# Ahora, pivotamos los datos para tener una columna por cada valor de 'familia'
resultado_pivotado <- resultado %>%
  pivot_wider(names_from = nombre, values_from = cantidad, values_fill = list(cantidad = 0))

# El dataframe 'resultado_pivotado' debería tener una columna por cada grupo de familia
# y la cantidad de veces que cada familia aparece por ID
compois_dataset <- resultado_pivotado

etiqueta_a_acc <- list()
acc_actual <- NA

for (line in lines) {
  if (grepl("^GCF|^GCA", line)) {
    acc_actual <- strsplit(line, "\\s+")[[1]][1]
  } else if (nzchar(line)) {
    etiqueta <- strsplit(line, "\\s+")[[1]][1]
    etiqueta_a_acc[[etiqueta]] <- acc_actual
  }
}

# Convertir la lista a data.frame
df_map <- data.frame(
  Etiqueta = names(etiqueta_a_acc),
  Acc = unlist(etiqueta_a_acc),
  stringsAsFactors = FALSE
)

compois_dataset <- merge(compois_dataset, df_map, by.x = "seqID", by.y = "Etiqueta", all.x = TRUE)


# Ahora queremos agregar una fila por cada valor distinto de 'Acc'

# Filtramos las columnas numéricas (exceptuando 'Acc' y 'ID')
numerical_columns <- compois_dataset %>% select(-c(Acc, seqID)) %>% select(where(is.numeric))

# Agrupamos por 'Acc' y sumamos las columnas numéricas
resultado <- compois_dataset %>%
  group_by(Acc) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE))
# Renombramos la columna 'Acc
resultado <- resultado %>%
  rename(seqID = Acc)
compois_dataset$Acc <- NULL 

# Aseguramos que las columnas de ambos dataframes sean exactamente las mismas
compo_columns <- colnames(compois_dataset)
resultado <- resultado[, compo_columns]

# Ahora añadimos las filas sin duplicar columnas
compois_dataset <- bind_rows(compois_dataset, resultado)

# Cargar librerías necesarias
library(dplyr)

# Leer el archivo .txt
txt_lines <- readLines("accesion_contig_chromosome_list.txt")

ids_txt <- sapply(txt_lines, function(x) strsplit(x, " ")[[1]][1])


# Reorganizar el dataset según el orden de los IDs en el archivo .txt
compois_dataset <- compois_dataset %>%
  filter(seqID %in% ids_txt) %>%   
  arrange(match(seqID, ids_txt)) 

write.table(compois_dataset, "compots_family.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

