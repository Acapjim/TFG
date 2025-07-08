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
dir = "W://QC_metrics_R/"

input_checkm_dataset = paste(dir, "bin_stats_ext_efaecium_tsv_results.tsv", sep = "")

input_busco_dataset = paste(dir, "busco.tsv", sep = "")
input_sylph_dataset = paste(dir, "sylph_results.tsv", sep = "")

checkm_dataset = read.delim(input_checkm_dataset, sep = "\t")
colnames(checkm_dataset)[-1] <- paste("checkm_", colnames(checkm_dataset)[-1], sep = "")

busco_dataset = read.delim(input_busco_dataset, sep = "\t")
colnames(busco_dataset)[-1] <- paste("busco_", colnames(busco_dataset)[-1], sep = "")


sylph_dataset = read.delim(input_sylph_dataset, sep = "\t")
colnames(sylph_dataset)[-1] <- paste("sylph_", colnames(sylph_dataset)[-1], sep = "")


# 2. Modificar la columna de Busco_df (eliminar el sufijo ".txt")
busco_dataset$accesion <- sub("\\.fna$", "", busco_dataset$accesion)

# 3. Extraer el último término de la columna Simple de Sylph_df (eliminar el sufijo ".exe")
sylph_dataset$Sample_file <- sub("\\.fna_1.fastq$", "", basename(sylph_dataset$Sample_file))

# 4. Unir los tres dataframes en un solo dataframe master_df
master_df <- checkm_dataset %>%
  inner_join(busco_dataset, by = c("acc" = "accesion")) %>%
  inner_join(sylph_dataset, by = c("acc" = "Sample_file"))


# Seleccionar solo los genomas sin duplicaciones del listado .txt

txt_file <- "lista_genomas_sin_duplicaciones.txt"
df_txt <- data.frame(
  acc = gsub("\\.fna$", "", trimws(readLines(txt_file))),
  stringsAsFactors = FALSE
)

# Limpiar TSV también
master_df_sindu <- master_df

master_df_sindu$acc = trimws(master_df_sindu$acc)

# Filtrar usando semi_join para conservar coincidencias
master_df_sindu <- semi_join(master_df_sindu, df_txt, by = "acc")

write.table(master_df_sindu,file = "Master_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

