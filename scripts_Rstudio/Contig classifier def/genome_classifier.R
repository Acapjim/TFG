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




# Leer la tabla desde el archivo Tareas.tsv
tareas <- read.delim("compots_family_chromosomed.tsv", stringsAsFactors = FALSE)

tareas <- tareas %>%
  rename(contig = seqID)

# Modificar la columna 'contig' según tus condiciones:
tareas$contig <- ifelse(
  grepl("^GCF|^GCA", tareas$contig) | grepl("_chromosome$", tareas$contig),
  tareas$contig,
  paste0(tareas$contig, "_plasmid")
)

# Identificar las filas que son GCF o GCA
es_gcf_gca <- grepl("^GCF|^GCA", tareas$contig)

# Crear una columna de bloque: cada vez que aparece GCF/GCA, se incrementa un contador
tareas$bloque_id <- cumsum(es_gcf_gca)

# Crear una lista para almacenar resultados
resultado <- list()

# Procesar por bloque
for (bloque in unique(tareas$bloque_id)) {
  datos_bloque <- tareas[tareas$bloque_id == bloque, ]
  
  # GCF/GCA al inicio del bloque
  referencia <- datos_bloque$contig[1]
  
  # 1. Añadir línea GCF/GCA tal cual
  resultado[[length(resultado) + 1]] <- datos_bloque[1, ]
  
  # 2. Línea con "_chromosome"
  fila_complete <- datos_bloque[grepl("_chromosome$", datos_bloque$contig), ]
  if (nrow(fila_complete) == 1) {
    fila_complete$contig <- paste0(referencia, "_chromosome")
    resultado[[length(resultado) + 1]] <- fila_complete
  }
  
  # 3. Filas "_incomplete" sumar
  filas_incompletas <- datos_bloque[grepl("_plasmid$", datos_bloque$contig), ]
  if (nrow(filas_incompletas) > 0) {
    fila_suma <- filas_incompletas[1, ]
    fila_suma$contig <- paste0(referencia, "_plasmid")
    
    # Sumar columnas numéricas (tratando NA como 0)
    for (col in names(fila_suma)) {
      if (col != "contig" && col != "bloque_id") {
        valores <- suppressWarnings(as.numeric(filas_incompletas[[col]]))
        valores[is.na(valores)] <- 0
        fila_suma[[col]] <- sum(valores)
      }
    }
    
    resultado[[length(resultado) + 1]] <- fila_suma
  }
}

# Combinar todo
resultado_final <- do.call(rbind, resultado)

# Quitar la columna de bloque_id antes de guardar
resultado_final$bloque_id <- NULL

# Leer el archivo final
final <- resultado_final

# Reemplazar los 0 por NA en todas las columnas excepto 'contig'
for (col in names(final)) {
  if (col != "contig") {
    valores <- suppressWarnings(as.numeric(final[[col]]))
    final[[col]] <- ifelse(!is.na(valores) & valores == 0, NA, final[[col]])
  }
}

# Guardar el resultado con los ceros convertidos a NA
write.table(final, "compots_family_chromobile.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

