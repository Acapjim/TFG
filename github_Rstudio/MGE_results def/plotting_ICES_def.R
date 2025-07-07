########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")
install.packages("svglite")
library(svglite)
library(gridExtra)

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

# configuració del gràfic
# Definir el tema base con los ajustes globales de texto
theme_base <- theme(
  plot.title = element_text(size = 30, face = "bold"),   
  axis.title.x = element_text(size = 20),                 
  axis.title.y = element_text(size = 20),                 
  axis.text.x = element_text(angle = 90, hjust = 1, size = 18), 
  axis.text.y = element_text(size = 18),                  
  legend.title = element_text(size = 20),                 
  legend.text = element_text(size = 18)                   
)


#loading data 
dir = "W://MGE_results/"

input_ICE_dataset = paste(dir, "icescreen_family_chromobile.tsv", sep = "")

ICE_dataset = read.delim(input_ICE_dataset, sep = "\t")


ICE_dataset_genomas <- ICE_dataset %>% filter(grepl("gen$", contig))

# Asignar la columna "contig" como los nombres de fila
rownames(ICE_dataset_genomas) <- ICE_dataset_genomas$contig


# Eliminar la columna "contig" ya que ahora es parte de los rownames
ICE_dataset_genomas$contig <- NULL


# Calcular el porcentaje de valores distintos de "NaN" por columna
porcentaje_validos <- sapply(ICE_dataset_genomas, function(col) {
  sum(!is.na(col)) / length(col) * 100
})


# Redondeamos a 2 decimales
porcentaje_validos <- round(porcentaje_validos, 2)

# Anyadirlo como una nueva fila al final del data.frame
ICE_dataset_genomas <- rbind(ICE_dataset_genomas, as.list(porcentaje_validos))

# Opcional: ponerle un nombre a la fila anyadida
rownames(ICE_dataset_genomas)[nrow(ICE_dataset_genomas)] <- "Porcentaje_abundancia"

#graficar
# Extraer la fila de porcentaje_abundancia como un vector numerico
valores <- as.numeric(ICE_dataset_genomas["Porcentaje_abundancia", ])

# Obtener los nombres de las columnas como etiquetas
etiquetas <- colnames(ICE_dataset_genomas)

# Crear un data frame para ggplot2
df_plot <- data.frame(
  categoria = etiquetas,
  porcentaje = valores
)

# Usar ggplot2 para graficar
library(ggplot2)

genpcent_p <- ggplot(df_plot, aes(x = categoria, y = porcentaje)) +
  geom_bar(stat = "identity", fill = "#00BA38") +
  labs(title = "A",
       x = "Família",
       y = "ICE trobats (%)") +
   theme_minimal() +  
  theme_base + theme(axis.text.x = element_text(angle = 45, hjust = 1))

genpcent_p

# Guardar el grafico como un archivo SVG con un tamanyo adecuado
ggsave("ICE_genomic_pc.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300)   


# Crear la columna "total" sumando solo las columnas numericas
ICE_dataset_genomas <- ICE_dataset %>% filter(grepl("gen$", contig))
ICE_dataset_genomas$contig <- NULL

# Seleccionamos las columnas cuyo % son mayores a un umbral
selected_cols <- colnames(ICE_dataset_genomas)[valores > 0]

# Seleccionamos las columnas de ICE_dataset segun los nombres de las columnas seleccionadas
# Y agregamos la columna "contig" al dataframe final
general_grafic_df <- ICE_dataset[, c("contig", selected_cols), drop = FALSE]

general_grafic_df$contig <- ICE_dataset$contig

# Crear un nuevo dataframe con solo esas columnas
df_selected <- ICE_dataset[, selected_cols, drop = FALSE]

# Convertimos columnas a numericas si es necesario
# y calculamos media y sd (ignorando "NaN" si es string)

summary_stats <- data.frame(
  columna = character(),
  media = numeric(),
  sd = numeric(),
  stringsAsFactors = FALSE
)

for (col_name in names(df_selected)) {
  col_data <- df_selected[[col_name]]
  
  # Filtrar valores validos
  numeric_values <- as.numeric(col_data[!is.na(col_data)])
  
  # Calcular estadisticas
  media <- mean(numeric_values, na.rm = TRUE)
  sd <- sd(numeric_values, na.rm = TRUE)
  
  summary_stats <- rbind(summary_stats, data.frame(
    columna = col_name,
    media = media,
    sd = sd
  ))
}


gencount_p <- ggplot(summary_stats, aes(x = columna, y = media)) +
  geom_bar(stat = "identity", fill = "#00BA38") +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd), width = 0.2) +
  labs(title = "A",
       x = "Família", y = "ICE per mostra") + 
   theme_minimal() + theme_base + theme(axis.text.x = element_text(angle = 45, hjust = 1))

gencount_p

ggsave("ICE_genomic_average.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300) 

# PARA CROMOSOMA

ICE_dataset_cromosoma <- ICE_dataset %>% filter(grepl("chromosome$", contig))

# Asignar la columna "contig" como los nombres de fila
rownames(ICE_dataset_cromosoma) <- ICE_dataset_cromosoma$contig


# Eliminar la columna "contig" ya que ahora es parte de los rownames
ICE_dataset_cromosoma$contig <- NULL


# Calcular el porcentaje de valores distintos de "NaN" por columna
porcentaje_validos <- sapply(ICE_dataset_cromosoma, function(col) {
  sum(!is.na(col)) / length(col) * 100
})


# Redondeamos a 2 decimales
porcentaje_validos <- round(porcentaje_validos, 2)

# AÃ±adirlo como una nueva fila al final del data.frame
ICE_dataset_cromosoma <- rbind(ICE_dataset_cromosoma, as.list(porcentaje_validos))

# Ponerle un nombre a la fila anyadida
rownames(ICE_dataset_cromosoma)[nrow(ICE_dataset_cromosoma)] <- "Porcentaje_abundancia"

#graficar
# Extraer la fila de porcentaje_abundancia como un vector numerico
valores <- as.numeric(ICE_dataset_cromosoma["Porcentaje_abundancia", ])

# Obtener los nombres de las columnas como etiquetas
etiquetas <- colnames(ICE_dataset_cromosoma)

# Crear un data frame para ggplot2
df_plot <- data.frame(
  categoria = etiquetas,
  porcentaje = valores
)

# Usar ggplot2 para graficar
library(ggplot2)

crompcent_p <- ggplot(df_plot, aes(x = categoria, y = porcentaje)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  labs(title = "B",
       x = "Famíla",
       y = "ICE trobats (%)") +
   theme_minimal() + theme_base + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  theme_base
crompcent_p

ggsave("ICE_cromosomic_pc.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300) 

# Crear la columna "total" sumando solo las columnas numericas
ICE_dataset_cromosoma <- ICE_dataset %>% filter(grepl("chromosome$", contig))
ICE_dataset_cromosoma$contig <- NULL

# Filtrar columnas con mas de X valores validos
selected_cols <- colnames(ICE_dataset_cromosoma)[valores > 0]

# Crear un nuevo dataframe con solo esas columnas
df_selected <- ICE_dataset_cromosoma[, selected_cols, drop = FALSE]

# Convertimos columnas a numÃ©ricas si es necesario
# y calculamos media y sd (ignorando "NaN" si es string)

summary_stats <- data.frame(
  columna = character(),
  media = numeric(),
  sd = numeric(),
  stringsAsFactors = FALSE
)

for (col_name in names(df_selected)) {
  col_data <- df_selected[[col_name]]
  
  # Filtrar valores validos
  numeric_values <- as.numeric(col_data[!is.na(col_data)])
  
  # Calcular estadisticas
  media <- mean(numeric_values, na.rm = TRUE)
  sd <- sd(numeric_values, na.rm = TRUE)
  
  summary_stats <- rbind(summary_stats, data.frame(
    columna = col_name,
    media = media,
    sd = sd
  ))
}


cromcount_p <- ggplot(summary_stats, aes(x = columna, y = media)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd), width = 0.2) +
  labs(title = "B",
       x = "Família", y = "ICE per mostra") +
    theme_minimal() + theme_base + theme(axis.text.x = element_text(angle = 45, hjust = 1))
cromcount_p
ggsave("ICE_cromosomic_average.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300) 


# PARA PLASMIDO
ICE_dataset_plasmid <- ICE_dataset %>% filter(grepl("plasmid$", contig))
# Asignar la columna "contig" como los nombres de fila
rownames(ICE_dataset_plasmid) <- ICE_dataset_plasmid$contig


# Eliminar la columna "contig" ya que ahora es parte de los rownames
ICE_dataset_plasmid$contig <- NULL


# Calcular el porcentaje de valores distintos de "NaN" por columna
porcentaje_validos <- sapply(ICE_dataset_plasmid, function(col) {
  sum(!is.na(col)) / length(col) * 100
})


# Redondeamos a 2 decimales
porcentaje_validos <- round(porcentaje_validos, 2)

# Anyaadirlo como una nueva fila al final del data.frame
ICE_dataset_plasmid <- rbind(ICE_dataset_plasmid, as.list(porcentaje_validos))

# Ponerle un nombre a la fila anyadida
rownames(ICE_dataset_plasmid)[nrow(ICE_dataset_plasmid)] <- "Porcentaje_abundancia"

#graficar
# Extraer la fila de porcentaje_abundancia como un vector numerico
valores <- as.numeric(ICE_dataset_plasmid["Porcentaje_abundancia", ])

# Obtener los nombres de las columnas como etiquetas
etiquetas <- colnames(ICE_dataset_plasmid)

# Crear un data frame para ggplot2
df_plot <- data.frame(
  categoria = etiquetas,
  porcentaje = valores
)

# Usar ggplot2 para graficar
library(ggplot2)

plaspcent_p <- ggplot(df_plot, aes(x = categoria, y = porcentaje)) +
  geom_bar(stat = "identity", fill = "#619CFF") +
  labs(title = "C",
       x = "Família",
       y = "ICE trobats (%)") +
    theme_minimal() + theme_base +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

plaspcent_p
ggsave("ICE_plasmid_pc.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300) 

# Crear la columna "total" sumando solo las columnas numericas
ICE_dataset_plasmid <- ICE_dataset %>% filter(grepl("plasmid$", contig))
ICE_dataset_plasmid$contig <- NULL

# Filtrar columnas con mas de X valores validos
selected_cols <- colnames(ICE_dataset_plasmid)[valores > 0]

# Crear un nuevo dataframe con solo esas columnas
df_selected <- ICE_dataset_plasmid[, selected_cols, drop = FALSE]

# Convertimos columnas a numericas si es necesario
# y calculamos media y sd (ignorando "NaN" si es string)

summary_stats <- data.frame(
  columna = character(),
  media = numeric(),
  sd = numeric(),
  stringsAsFactors = FALSE
)

for (col_name in names(df_selected)) {
  col_data <- df_selected[[col_name]]
  
  # Filtrar valores validos
  numeric_values <- as.numeric(col_data[!is.na(col_data)])
  
  # Calcular estadisticas
  media <- mean(numeric_values, na.rm = TRUE)
  sd <- sd(numeric_values, na.rm = TRUE)
  
  summary_stats <- rbind(summary_stats, data.frame(
    columna = col_name,
    media = media,
    sd = sd
  ))
}


plascount_p <- ggplot(summary_stats, aes(x = columna, y = media)) +
  geom_bar(stat = "identity", fill = "#619CFF") +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd), width = 0.2) +
  labs(title = "C",
       x = "Família", y = "ICE per mostra") +
   theme_minimal() + theme_base + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plascount_p
ggsave("ICE_plasmid_average.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300) 

# calcular el total per a cada dataset
ICE_dataset_genomas$total <- rowSums(ICE_dataset_genomas[, sapply(ICE_dataset_genomas, is.numeric)], na.rm = TRUE)
ICE_dataset_cromosoma$total <- rowSums(ICE_dataset_cromosoma[, sapply(ICE_dataset_cromosoma, is.numeric)], na.rm = TRUE)
ICE_dataset_plasmid$total <- rowSums(ICE_dataset_plasmid[, sapply(ICE_dataset_plasmid, is.numeric)], na.rm = TRUE)



# Calcular la media y desviacion estandar para la columna "total" de cada dataset
media_pri <- mean(ICE_dataset_genomas$total, na.rm = TRUE)
desviacion_pri <- sd(ICE_dataset_genomas$total, na.rm = TRUE)

media_seg <- mean(ICE_dataset_cromosoma$total, na.rm = TRUE)
desviacion_seg <- sd(ICE_dataset_cromosoma$total, na.rm = TRUE)

media_ter <- mean(ICE_dataset_plasmid$total, na.rm = TRUE)
desviacion_ter <- sd(ICE_dataset_plasmid$total, na.rm = TRUE)


# Crear un dataframe con los valores que necesitas para graficar
data_grafico <- data.frame(
  Dataset = c("Genoma", "Cromosoma", "Plasmidis"),
  Media = c(media_pri, media_seg, media_ter),
  Desviacion = c(desviacion_pri, desviacion_seg, desviacion_ter)
)

# Graficar los resultados con ggplot

totalcount_p <- ggplot(data_grafico, aes(x = Dataset, y = Media, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  geom_errorbar(aes(ymin = Media - Desviacion, ymax = Media + Desviacion), width = 0.2) +
  labs(title = "D",
       y = "ICE trobats totals", 
       x = "") +
  scale_fill_manual(values = c("Cromosoma" = "#F8766D", 
                               "Plasmidis" = "#619CFF", 
                               "Genoma" = "#00BA38")) +
   theme_minimal() + theme_base + theme(axis.text.x = element_text(angle = 45, hjust = 1))
totalcount_p

ggsave("ICE_total_average.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300) 

# Extraer el sufijo del nombre
datos <- general_grafic_df %>%
  mutate(grupo = case_when(
    grepl("_gen$", contig) ~ "Genoma",
    grepl("_chromosome$", contig) ~ "Cromosoma",
    grepl("_plasmid$", contig) ~ "Plasmidis"
  ))

# Eliminar la columna 'contig' si ya no la necesitas
# Convertir a formato largo (long format)
datos_largos <- datos %>%
  pivot_longer(cols = -c(contig, grupo), names_to = "variable", values_to = "valor")


resumen <- datos_largos %>%
  group_by(grupo, variable) %>%
  summarise(
    media = mean(valor, na.rm = TRUE),
    sd = sd(valor, na.rm = TRUE),
    .groups = "drop"
  )

# Calcular el maximo de forma segura, ignorando valores NA
max_y <- max(resumen$media + resumen$sd, na.rm = TRUE)

# Crear el grafico con la cuadricula y el eje Y mejorados
grafico_general <- ggplot(resumen, aes(x = variable, y = media, fill = grupo)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd),
                width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Família", y = "ICE per mostra", fill = "Grupo") +
   theme_minimal() + theme_base +
  scale_fill_manual(values = c("Cromosoma" = "#F8766D", 
                               "Plasmidis" = "#619CFF", 
                               "Genoma" = "#00BA38")) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),  
    axis.title.y = element_text(size = 14), 
    panel.grid.major.y = element_line(linewidth = 0.5, color = "gray"), 
    panel.grid.minor.y = element_line(linewidth = 0.2, color = "gray", linetype = "dashed"),  
    panel.grid.major.x = element_line(linewidth = 0.2, color = "lightgray") 
  ) +
  scale_y_continuous(breaks = seq(0, max_y, by = 5))  

grafico_general

# Guardar el grafico como un archivo SVG con un tamaño adecuado
ggsave("ICE_count_complete.svg", 
       plot = last_plot(),       
       device = "svg",           
       width = 12,               
       height = 8,               
       dpi = 300)   

svg("ICE_pct_multipanel.svg", width = 24, height = 16)  
grid.arrange(genpcent_p, crompcent_p, plaspcent_p, ncol = 2, nrow = 2)
dev.off()

svg("ICE_count_multipanel.svg", width = 24, height = 16)  
grid.arrange(gencount_p, cromcount_p, plascount_p, totalcount_p, ncol = 2, nrow = 2)
dev.off()

