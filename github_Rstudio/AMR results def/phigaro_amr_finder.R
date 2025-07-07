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

input_amr_dataset = paste(dir, "amr_results_combinado.tsv", sep = "")
input_phigaro_dataset = paste(dir, "documento_combinado_phigaro.tsv", sep = "")

amr_dataset = read.delim(input_amr_dataset, sep = "\t")
phigaro_dataset = read.delim(input_phigaro_dataset, sep = "\t")


amr_dataset_elements <- amr_dataset[amr_dataset$Element.subtype == "AMR", ]


# Primero: renombramos para tener una clave comun
amr_dataset_elements$scaffold <- amr_dataset_elements$Contig.id

# Segundo: hacemos merge interno por nombre (esto genera todas las combinaciones posibles)
amr_phage_dataset <- merge(phigaro_dataset, amr_dataset_elements, by = "scaffold", allow.cartesian = TRUE)

#linea logica par encontrar AMR
amr_phage_dataset$Amr_phage <- ifelse(amr_phage_dataset$Start >= amr_phage_dataset$begin & 
                             amr_phage_dataset$end >= amr_phage_dataset$Stop, 
                           "Yes", "No")

amr_phage_dataset_found <- subset(amr_phage_dataset, Amr_phage == "Yes")


write.table(amr_phage_dataset_found, file = "amr_phage_encontrados.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)



