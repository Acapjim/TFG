########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")


library(dplyr)
library(plyr)
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
input_icescreen_dataset = paste(dir, "documento_combinado_icescreen.tsv", sep = "")

amr_dataset = read.delim(input_amr_dataset, sep = "\t")
icescreen_dataset = read.delim(input_icescreen_dataset, sep = "\t")


amr_dataset_elements <- amr_dataset[amr_dataset$Element.subtype == "AMR", ]


# Primero: renombramos para tener una clave comun
icescreen_dataset <- icescreen_dataset %>% dplyr::rename(contig = Genome_accession)
amr_dataset_elements$contig <- amr_dataset_elements$Contig.id


# Segundo: hacemos merge interno por nombre (esto genera todas las combinaciones posibles)
amr_ices_dataset <- merge(icescreen_dataset, amr_dataset_elements, by = "contig", allow.cartesian = TRUE)

#combinado$Amr_phage <- ifelse((combinado$Stop - combinado$Start) + combinado$begin >= combinado$end, "Yes", "No")
amr_ices_dataset$Amr_ice <- ifelse(amr_ices_dataset$Start_of_most_upstream_SP <= amr_ices_dataset$Start & 
                             amr_ices_dataset$Stop <= amr_ices_dataset$Stop_of_most_downstream_SP, 
                           "Yes", "No")

amr_ices_dataset_found <- subset(amr_ices_dataset, Amr_ice == "Yes")


write.table(amr_ices_dataset_found, file = "amr_ices_encontrados.tsv", 
           sep = "\t", row.names = FALSE, quote = FALSE)



