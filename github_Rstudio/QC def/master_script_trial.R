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

#input_extracted_dataset = paste(dir, "extracted_dataset_updated.csv", sep = "")
input_master_dataset = paste(dir, "Master_table_QC2.tsv", sep = "")
#extracted_dataset = read.delim(input_extracted_dataset, sep = ",")
master_dataset = read.delim(input_master_dataset, sep = "\t")
###################################################################################################################
#multi-panel for QC metrics 
###################################################################################################################


## ordenar per acc
master_dataset = master_dataset %>% arrange(acc)


### fer el dataset de amb les variables de seleccio

select_dataset = bind_cols(master_dataset %>% select(acc, checkm_Genome.size, checkm_N50..contigs., checkm_X..contigs,checkm_Completeness,checkm_Contamination, sylph_Contig_name, busco_X.C, busco_X.S, busco_X.D, busco_X.F, busco_X.M,busco_X.gap, Assembly.Sequencing.Tech))

### triar els cutoff per cada variable
#### checkm_Genome.size
select_dataset$good_lenght = "no"
# Calcular la media
media_lenght <- mean(select_dataset$checkm_Genome.size, na.rm = TRUE) 

# Calcular el intervalor permitido
intervalo = 10
lmax <- media_lenght * (1+intervalo/100)
lmin <- media_lenght * (1-intervalo/100)

lenght_cutoff = 3000000
select_dataset$good_lenght[which(as.numeric(select_dataset$checkm_Genome.size) > lmin & as.numeric(select_dataset$checkm_Genome.size) < lmax)] <- "yes"
select_dataset$good_lenght= factor(select_dataset$good_lenght)

# Paso 2: Mover 'good_lenght' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "checkm_Genome.size")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_lenght",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]
# borramos la duplicada
select_dataset$good_lenght.1 <- NULL

#### checkm_N50..contig

select_dataset$good_n50co = "no"
n50co_cutoff = 2416265
select_dataset$good_n50co[which(as.numeric(select_dataset$checkm_N50..contigs.) > n50co_cutoff)] = "yes"
select_dataset$good_n50co= factor(select_dataset$good_n50co)

# Paso 2: Mover 'good_n50co' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "checkm_N50..contigs.")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_n50co",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]


# borramos la duplicada
select_dataset$good_n50co.1 <- NULL


#### checkm_contig

select_dataset$good_ncon = "no"
ncon_cutoff = 18
select_dataset$good_ncon[which(as.numeric(select_dataset$checkm_X..contigs) < ncon_cutoff)] = "yes"
select_dataset$good_ncon= factor(select_dataset$good_ncon)

# Paso 2: Mover 'good_n50co' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "checkm_X..contigs")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_ncon",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]


# borramos la duplicada
select_dataset$good_ncon.1 <- NULL
#### checkm_Completeness

select_dataset$good_comp = "no"
comp_cutoff = 98
select_dataset$good_comp[which(as.numeric(select_dataset$checkm_Completeness) > comp_cutoff)] = "yes"
select_dataset$good_comp= factor(select_dataset$good_comp)

# Paso 2: Mover 'low_lenght' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "checkm_Completeness")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_comp",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]


# borramos la duplicada
select_dataset$good_comp.1 <- NULL


#### checkm_Contamination

select_dataset$good_cont = "no"
cont_cutoff = 1
select_dataset$good_cont[which(as.numeric(select_dataset$checkm_Contamination) < cont_cutoff)] = "yes"
select_dataset$good_cont= factor(select_dataset$good_cont)

# Paso 2: Mover 'good_cont' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "checkm_Contamination")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_cont",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]

# borramos la duplicada
select_dataset$good_cont.1 <- NULL


### checkm_Strain.heterogeneity

#select_dataset$good_hete = "no"
#hete_cutoff = 50
#select_dataset$good_hete[which(as.numeric(select_dataset$checkm_Strain.heterogeneity) > hete_cutoff)] = "yes"
#select_dataset$good_hete= factor(select_dataset$good_hete)

# Paso 2: Mover 'good_hete' justo después de 'lenght'
# Obtener nombres de columnas
#colnames_all <- colnames(select_dataset)

# Encontrar posición de la columna 'lenght'
#pos_lenght <- which(colnames_all == "checkm_Strain.heterogeneity")

# Reorganizar columnas
#select_dataset <- select_dataset[, c(
#  colnames_all[1:pos_lenght],
#  "good_hete",
#  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
#)]

# borramos la duplicada
#select_dataset$good_hete.1 <- NULL
# esto li usaramos iría dentro de nombres personalizados. good_hete="Heterogeneidad"

### BUsco_C

select_dataset$good_C = "no"
C_cutoff = 91
select_dataset$good_C[which(as.numeric(select_dataset$busco_X.C) > C_cutoff)] = "yes"
select_dataset$good_C= factor(select_dataset$good_C)

# Paso 2: Mover 'good_hete' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "busco_X.C")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_C",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]

# borramos la duplicada
select_dataset$good_C.1 <- NULL

### BUsco_S

select_dataset$good_S = "no"
S_cutoff = 91
select_dataset$good_S[which(as.numeric(select_dataset$busco_X.S) > S_cutoff)] = "yes"
select_dataset$good_S = factor(select_dataset$good_S)

# Paso 2: Mover 'good_hete' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "busco_X.S")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_S",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]

# borramos la duplicada
select_dataset$good_S.1 <- NULL

### BUsco_D

select_dataset$good_D = "no"
D_cutoff = 1
select_dataset$good_D[which(as.numeric(select_dataset$busco_X.D) < D_cutoff)] = "yes"
select_dataset$good_D = factor(select_dataset$good_D)

# Paso 2: Mover 'good_hete' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "busco_X.D")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_D",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]

# borramos la duplicada
select_dataset$good_D.1 <- NULL


### BUsco_F

select_dataset$good_F = "no"
F_cutoff = 1
select_dataset$good_F[which(as.numeric(select_dataset$busco_X.F) < F_cutoff)] = "yes"
select_dataset$good_F = factor(select_dataset$good_F)

# Paso 2: Mover 'good_hete' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "busco_X.F")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_F",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]

# borramos la duplicada
select_dataset$good_F.1 <- NULL

select_dataset$good_M = "no"
M_cutoff = 4.1
select_dataset$good_M[which(as.numeric(select_dataset$busco_X.M) < M_cutoff)] = "yes"
select_dataset$good_M = factor(select_dataset$good_M)

# Paso 2: Mover 'good_hete' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "busco_X.M")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_M",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]

# borramos la duplicada
select_dataset$good_M.1 <- NULL


#gap%
select_dataset$busco_X.gap <- as.numeric(gsub("%", "", select_dataset$busco_X.gap))

select_dataset$good_gap = "no"
gap_cutoff = 0.2
select_dataset$good_gap[which(as.numeric(select_dataset$busco_X.gap) < gap_cutoff)] = "yes"
select_dataset$good_gap = factor(select_dataset$good_gap)

# Paso 2: Mover 'good_hete' justo despues de 'lenght'
# Obtener nombres de columnas
colnames_all <- colnames(select_dataset)

# Encontrar posicion de la columna 'lenght'
pos_lenght <- which(colnames_all == "busco_X.gap")

# Reorganizar columnas
select_dataset <- select_dataset[, c(
  colnames_all[1:pos_lenght],
  "good_gap",
  colnames_all[(pos_lenght + 1):length(colnames_all)][colnames_all[(pos_lenght + 1):length(colnames_all)] != "low_lenght"]
)]


# borramos la duplicada
select_dataset$good_gap.1 <- NULL

# listamos los criterios de seleccion no superados
nombres_personalizados <- c( good_n50co = "n50 contig",good_ncon = "Ncont", good_comp = "Completeness", good_cont = "Contamination", good_C = "Busco C", good_S = "Busco S", good_D = "Busco D", good_F = "Busco F", good_M = "Busco M", good_gap = "gaps")

select_dataset$Exclusion <- apply(select_dataset[, names(nombres_personalizados)], 1, function(fila) {
  ## Selecciona los nombres personalizados donde el valor es "NO"
seleccion <- nombres_personalizados[fila == "no"]
  paste(seleccion, collapse = ", ") })
good_lenght = "Longitud"
# Indicar seq_tech no fiables
select_dataset$Exclusion <- apply(select_dataset, 1, function(fila) {
  # Identifica exclusiones personalizadas (columnas con valor "no")
  exclusiones <- nombres_personalizados[fila[names(nombres_personalizados)] == "no"]
  
  # Anyadir exclusiones segun la tecnologia
 # tech_val <- fila["Assembly.Sequencing.Tech"]
 # if (is.na(tech_val) || tech_val == "") {exclusiones <- c(exclusiones, "No-tech")}
 # else if (tech_val == "Oxford Nanopore MinION") {exclusiones <- c(exclusiones, "Oxford Nanopore MinION")}
#  else if (tech_val == "Oxford Nanopore MiniION") {exclusiones <- c(exclusiones, "Oxford Nanopore MiniION")}
  
  # Combina todo con coma y espacio
  paste(exclusiones, collapse = ", ")})

# contamos variables defectuosas
select_dataset$N_Exclusiones <- sapply(strsplit(select_dataset$Exclusion, ",\\s*"), length)
# mes de una variable defectuosa la eliminem
select_dataset$kept <- ifelse(select_dataset$N_Exclusiones < 2, "Yes", "No")

# Filtrar las filas donde 'kept' es igual a "Yes" y guarar acc
kept_dataset <- select_dataset[select_dataset$kept == "Yes", "acc", drop = FALSE]

write.table(select_dataset,file = "Taula_valors_decisi�.tsv", sep = "\t", row.names = FALSE)
write.table(kept_dataset,file = "kept_genomes.tsv", sep = "\t", row.names = FALSE)

