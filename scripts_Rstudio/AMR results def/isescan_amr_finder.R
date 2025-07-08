########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")

library(purrr)
library(stringr)
library(tidyr)
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

input_amr_dataset = paste(dir, "amr_results_combinado.tsv", sep = "")
input_isescan_dataset = paste(dir, "documento_combinado_isescan.tsv", sep = "")

amr_dataset = read.delim(input_amr_dataset, sep = "\t")
isescan_dataset = read.delim(input_isescan_dataset, sep = "\t")


amr_dataset_elements <- amr_dataset[amr_dataset$Element.subtype == "AMR", ]


# Primero: renombramos para tener una clave comun
amr_dataset_elements$seqID <- amr_dataset_elements$Contig.id

# Segundo: hacemos merge interno por nombre (esto genera todas las combinaciones posibles)
combinado <- merge(isescan_dataset, amr_dataset_elements, by = "seqID", allow.cartesian = TRUE)

#lógica para la detección 
combinado$Amr_IS <- ifelse(combinado$isBegin <= combinado$Start & 
                             combinado$Stop <= combinado$isEnd, 
                           "Yes", "No")

combinado$Amr_IS_proof <- ifelse(combinado$end1 <= combinado$Start & 
                             combinado$Stop <= combinado$start2, 
                           "Yes", "No")

Amr_is_encontrados <- subset(combinado, Amr_IS == "Yes")
Amr_is_validados <- subset(Amr_is_encontrados, Amr_IS_proof == "Yes")

#evaluamos si la Tpasa está completa
Amr_is_validados$Tpase_solapada <- ifelse(
  (Amr_is_validados$orfBegin <= Amr_is_validados$Start & Amr_is_validados$Start <= Amr_is_validados$orfEnd) |
    (Amr_is_validados$orfBegin <= Amr_is_validados$Stop & Amr_is_validados$Stop <= Amr_is_validados$orfEnd),
  "Yes", "No"
)
Amr_is_validados$Tpase_complete <- "No"

Amr_is_validados$Tpase_complete <- ifelse(
  Amr_is_validados$Tpase_solapada == "Yes" & Amr_is_validados$Strand == Amr_is_validados$strand, 
  "No", 
  "Yes"
)


write.table(Amr_is_encontrados, file = "Amr_is_encontrados.tsv", 
             sep = "\t", row.names = FALSE, quote = FALSE)


write.table(Amr_is_validados, file = "Amr_is_validados.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

#COMPOSITE TRANSPOSON FINDER
# hacemos el dataset
compo_trans_dataset <- isescan_dataset %>% select(seqID,isBegin, isEnd)

# metemos la ubicacion de cada IS en un bloque
compo_trans_dataset <- compo_trans_dataset %>% 
  mutate(IS = paste(isBegin, isEnd, sep = "-")) %>%
  select(seqID, IS)

# para cada seq las coordenadas de los IS que tiene
# Paso 1: Anyadir un numero de fila dentro de cada grupo para numerar los bloques
compo_trans_dataset <- compo_trans_dataset %>%
  group_by(seqID) %>%
  mutate(IS_id = paste0("IS", row_number())) %>%
  ungroup()

# Paso 2: Reordenar para que cada bloque este en una columna
compo_trans_dataset <- compo_trans_dataset %>%
  select(seqID, IS_id, IS) %>%
  pivot_wider(names_from = IS_id, values_from = IS)

VALOR_UMBRAL_sup <- 10001  
VALOR_UMBRAL_inf <- 2499  


compo_trans_dataset$donde <- NA  # Inicializamos la nueva columna

for (fila in 1:nrow(compo_trans_dataset)) {
  fila_actual <- compo_trans_dataset[fila, -1]  # Excluimos la columna "nombre"
  nombres_columnas <- names(fila_actual)
  
  valores_validos <- which(!is.na(fila_actual))
  
  # Solo seguir si hay al menos dos columnas validas
  if (length(valores_validos) >= 2) {
    coincidencias <- c()  # Guardar las parejas que cumplen
    
    for (i in 1:(length(valores_validos) - 1)) {
      idx1 <- valores_validos[i]
      idx2 <- valores_validos[i + 1]
      
      if (!is.na(fila_actual[[idx1]]) && !is.na(fila_actual[[idx2]])) {
        val1 <- strsplit(fila_actual[[idx1]], "-")[[1]]
        val2 <- strsplit(fila_actual[[idx2]], "-")[[1]]
        
        if (length(val1) == 2 && length(val2) == 2) {
          fin_1 <- as.numeric(val1[2])
          ini_2 <- as.numeric(val2[1])
          
          if ((ini_2 - fin_1) < VALOR_UMBRAL_sup & (ini_2 - fin_1) > VALOR_UMBRAL_inf){
            coincidencias <- c(coincidencias, paste0(nombres_columnas[idx1], "-", nombres_columnas[idx2]))
          }
        }
      }
    }
    
    # Si se encontraron coincidencias, unirlas con comas
    if (length(coincidencias) > 0) {
      compo_trans_dataset$donde[fila] <- paste(coincidencias, collapse = ", ")
    }
  }
}


procesar_compo <- function(df, col_nombre = "seqID", col_donde = "donde") {
  # Asegurarse de nombres limpios
  colnames(df) <- str_trim(colnames(df))
  df[[col_donde]] <- str_trim(df[[col_donde]])
  
  result <- map_dfr(seq_len(nrow(df)), function(i) {
    donde_str <- df[[col_donde]][i]
    nombre <- df[[col_nombre]][i]
    
    if (is.na(donde_str) || donde_str == "") return(NULL)
    
    donde_list <- str_split(donde_str, ",")[[1]]
    
    map_dfr(donde_list, function(rango) {
      partes <- str_split(rango, "-")[[1]]
      col_start <- str_trim(partes[1])
      col_stop <- str_trim(partes[2])
      
      val_start <- if (col_start %in% names(df)) df[i, col_start, drop = TRUE] else {
        warning(paste("Columna no encontrada:", col_start, "en fila", i))
        NA
      }
      
      val_stop <- if (col_stop %in% names(df)) df[i, col_stop, drop = TRUE] else {
        warning(paste("Columna no encontrada:", col_stop, "en fila", i))
        NA
      }
      
      tibble(
        !!col_nombre := nombre,
        Start = val_start,
        Stop = val_stop,
        !!col_donde := rango
      )
    })
  })
  
  return(result)
}


compo_trans_dataset_1_per_row <- procesar_compo(compo_trans_dataset)

# Funcion para extraer el numero antes y despues del guion en la columna 'Start'
extraer_numeros <- function(start_string) {
  partes <- strsplit(start_string, "-")[[1]] # Divide la cadena por el guion
  if (length(partes) == 2) {
    return(c(partes[1], partes[2])) # Devuelve ambos numeros
  } else {
    return(c(NA, NA)) # Devuelve NA si no hay guion o hay mas de uno
  }
}

# Aplica la funcion a cada valor de la columna 'Start' y crea nuevas columnas 'Begin' y 'End'
resultados <- lapply(compo_trans_dataset_1_per_row$Start, extraer_numeros)
compo_trans_dataset_1_per_row$isBegin <- as.numeric(sapply(resultados, function(x) x[1])) 
compo_trans_dataset_1_per_row$isEnd <- as.numeric(sapply(resultados, function(x) x[2])) 



compo_trans_dataset_1_per_row$family1 <- ""

# Hacemos el join para traer la columna 'family'
compo_trans_dataset_1_per_row <- compo_trans_dataset_1_per_row %>%
  left_join(isescan_dataset %>% select(seqID, isBegin, isEnd, family), by = c("seqID", "isBegin", "isEnd"))

# Copiamos el contenido de 'family' a 'family1'
compo_trans_dataset_1_per_row$family1 <- compo_trans_dataset_1_per_row$family

# Eliminamos la columna original 'family'
compo_trans_dataset_1_per_row$family <- NULL
compo_trans_dataset_1_per_row$isBegin <- NULL
compo_trans_dataset_1_per_row$isEnd <- NULL
compo_trans_dataset_1_per_row$Begin1 <- NULL
compo_trans_dataset_1_per_row$End1 <- NULL


# encontramos family2
# Funcion para extraer el numero antes y despues del guion en la columna 'Start'
extraer_numeros <- function(start_string) {
  partes <- strsplit(start_string, "-")[[1]] # Divide la cadena por el guion
  if (length(partes) == 2) {
    return(c(partes[1], partes[2])) # Devuelve ambos numeros
  } else {
    return(c(NA, NA)) # Devuelve NA si no hay guion o hay mas de uno
  }
}

# Aplica la funcion a cada valor de la columna 'Start' y crea nuevas columnas 'Begin' y 'End'
resultados <- lapply(compo_trans_dataset_1_per_row$Stop, extraer_numeros)
compo_trans_dataset_1_per_row$isBegin <- as.numeric(sapply(resultados, function(x) x[1])) 
compo_trans_dataset_1_per_row$isEnd <- as.numeric(sapply(resultados, function(x) x[2])) 


compo_trans_dataset_1_per_row$family2 <- ""

# Hacemos el join para traer la columna 'family'
compo_trans_dataset_1_per_row <- compo_trans_dataset_1_per_row %>%
  left_join(isescan_dataset %>% select(seqID, isBegin, isEnd, family), by = c("seqID", "isBegin", "isEnd"))

# Copiamos el contenido de 'family' a 'family2'
compo_trans_dataset_1_per_row$family2 <- compo_trans_dataset_1_per_row$family

# Eliminamos la columna original 'family' si no la necesitas
compo_trans_dataset_1_per_row$family <- NULL
compo_trans_dataset_1_per_row$isBegin <- NULL
compo_trans_dataset_1_per_row$isEnd <- NULL
compo_trans_dataset_1_per_row$Begin1 <- NULL
compo_trans_dataset_1_per_row$End1 <- NULL

#unificamos los valores

compo_trans_dataset_1_per_row <- compo_trans_dataset_1_per_row %>%
  mutate(
    nombre = paste(family1, family2, sep = "-"),
    Coincidencia = if_else(family1 == family2, "Yes", "No")
  )


# guardamos el resultado
compo_trans_dataset_1_per_row$donde <- NULL
compois_dataset <- compo_trans_dataset_1_per_row %>% filter(Coincidencia == "Yes")
compois_dataset$Coincidencia <- NULL
compois_dataset$family1 <- NULL
compois_dataset$family2 <- NULL


write.table(compois_dataset, file = "composis_family.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

# buscamos AMR en los posibles composite transposons encontrados.

# Primero: renombramos para tener una clave comun

# Segundo: hacemos merge interno por nombre (esto genera todas las combinaciones posibles)
combinado_compo_amr <- merge(compo_trans_dataset_1_per_row, amr_dataset_elements, by = "seqID", allow.cartesian = TRUE, suffixes = c(".x", ".y"))

# sacamos el final del primer IS y el principio del ultimo 

# Extraer los valores numericos desde los strings Start y Stop
combinado_compo_amr$Start_val <- as.numeric(sub(".*-", "", combinado_compo_amr$Start.x))  
combinado_compo_amr$Stop_val <- as.numeric(sub("-.*", "", combinado_compo_amr$Stop.x))    

# Aplicar la logica condicional
  combinado_compo_amr$AMR_found <- ifelse(combinado_compo_amr$Start_val <= combinado_compo_amr$Start.y & combinado_compo_amr$Stop.y <= combinado_compo_amr$Stop_val, "Yes", "No")
  combinado_compo_amr <- combinado_compo_amr %>% filter(Coincidencia == "Yes")
  combinado_compo_amr <- combinado_compo_amr %>% filter(AMR_found == "Yes")
  
  
write.table(combinado_compo_amr, file = "amr_composis_encontrados.tsv", 
              sep = "\t", row.names = FALSE, quote = FALSE)
