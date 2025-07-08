########################################################################################################
#SETUP
########################################################################################################
#packages and libraries 
install.packages("gtsummary")

library(gridExtra)
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


# Definir el tema base con los ajustes globales de texto
theme_base <- theme(
  plot.title = element_text(size = 28, face = "bold"),        
  plot.subtitle = element_text(size = 24),                    
  axis.title.x = element_text(size = 26),                     
  axis.title.y = element_text(size = 26),                     
  axis.text.x = element_text(size = 22),                      
  axis.text.y = element_text(size = 22),                      
  axis.line = element_line(colour = "black", linewidth = 0.8),
  axis.ticks = element_line(linewidth = 0.8),
  
  legend.title = element_text(size = 24),                     
  legend.text = element_text(size = 22),                      
  
  strip.text = element_text(size = 24),                       
  strip.background = element_rect(fill = "grey90", colour = NA)
)



#loading data 
dir = "W://QC_metrics_R/"

input_QC_dataset = paste(dir, "busco.tsv", sep = "")

busco_dataset = read.delim(input_QC_dataset, sep = "\t")
genomas_sindu <- readLines( "lista_genomas_sin_duplicaciones.txt")

#filtramos sin duplicaciones
busco_dataset <- subset(busco_dataset, accesion %in% genomas_sindu)

###################################################################################################################
#multi-panel plot for QC metrics 
###################################################################################################################
x_axis_size <- 12
size_axis_lines <-12
##busco

### C-%
#C pruebA

C_pcent = busco_dataset %>%
  select(accesion, X.C) %>%
  mutate(X.C = as.numeric(X.C)) %>%
  arrange(X.C)

C_cutoff = 91
C_pcent$good_C = ifelse(C_pcent$X.C > C_cutoff, "yes", "no")
C_pcent$good_C = factor(C_pcent$good_C)
C_pcent$index = seq_along(C_pcent$X.C)

Cpcent_p = ggplot(C_pcent, aes(x = index, y = X.C)) +
  geom_point(aes(color = good_C)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = C_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Complete Buscos %",
    subtitle = "A"
  )

Cpcent_p

svg("busco_Cpecent.svg", width = 8, height = 6)  
print(Cpcent_p)  
dev.off()  


### S-%

S_pcent = busco_dataset %>% select(accesion, X.S)

# S 
S_pcent = busco_dataset %>%
  select(accesion, X.S) %>%
  mutate(X.S = as.numeric(X.S)) %>%
  arrange(X.S)

S_cutoff = 91  # Ajusta el umbral segun tu criterio
S_pcent$good_S = ifelse(S_pcent$X.S > S_cutoff, "yes", "no")
S_pcent$good_S = factor(S_pcent$good_S)
S_pcent$index = seq_along(S_pcent$X.S)

Spcent_p = ggplot(S_pcent, aes(x = index, y = X.S)) +
  geom_point(aes(color = good_S)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = S_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Single-copy Buscos %",
    subtitle = "B"
  )

Spcent_p

svg("busco_Spcent_p.svg", width = 8, height = 6)  
print(Spcent_p)  
dev.off()  


# D-%
D_pcent = busco_dataset %>%
  select(accesion, X.D) %>%
  mutate(X.D = as.numeric(X.D)) %>%
  arrange(X.D)

D_cutoff = 1  # Bajo es mejor
D_pcent$good_D = ifelse(D_pcent$X.D < D_cutoff, "yes", "no")
D_pcent$good_D = factor(D_pcent$good_D)
D_pcent$index = seq_along(D_pcent$X.D)

Dpcent_p = ggplot(D_pcent, aes(x = index, y = X.D)) +
  geom_point(aes(color = good_D)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = D_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Duplicated Buscos %",
    subtitle = "C"
  )

Dpcent_p


svg("busco_Dpcent_p.svg", width = 8, height = 6)  
print(Dpcent_p)  
dev.off()  


# F%

F_pcent = busco_dataset %>%
  select(accesion, X.F) %>%
  mutate(X.F = as.numeric(X.F)) %>%
  arrange(X.F)

F_cutoff = 1
F_pcent$good_F = ifelse(F_pcent$X.F > F_cutoff, "yes", "no")
F_pcent$good_F = factor(F_pcent$good_F)
F_pcent$index = seq_along(F_pcent$X.F)

Fpcent_p = ggplot(F_pcent, aes(x = index, y = X.F)) +
  geom_point(aes(color = good_F)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = F_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Fragmented Buscos %",
    subtitle = "D"
  )

Fpcent_p

svg("busco_Fpcent_p.svg", width = 8, height = 6)  
print(Fpcent_p)  
dev.off()  


#M%
M_pcent = busco_dataset %>%
  select(accesion, X.M) %>%
  mutate(X.M = as.numeric(X.M)) %>%
  arrange(X.M)

M_cutoff = 4.1 # Por ejemplo, bajo es bueno
M_pcent$good_M = ifelse(M_pcent$X.M < M_cutoff, "yes", "no")
M_pcent$good_M = factor(M_pcent$good_M)
M_pcent$index = seq_along(M_pcent$X.M)

Mpcent_p = ggplot(M_pcent, aes(x = index, y = X.M)) +
  geom_point(aes(color = good_M)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = M_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Missing Buscos %",
    subtitle = "E"
  )

Mpcent_p


svg("busco_Mpcent_p.svg", width = 8, height = 6)  
print(Mpcent_p)  
dev.off()  

# %gap
gap_pcent = busco_dataset %>% select(accesion, X.gap)

gap_pcent$X.gap <- as.numeric(gsub("%", "", gap_pcent$X.gap))  
gap_pcent = gap_pcent %>% arrange(X.gap)

gap_cutoff = 0.2  # Umbral
gap_pcent$good_gap = ifelse(gap_pcent$X.gap > gap_cutoff, "yes", "no")
gap_pcent$good_gap = factor(gap_pcent$good_gap)

gap_pcent = gap_pcent %>%
  mutate(index = seq_along(X.gap))

gappcent_p = ggplot(gap_pcent, aes(x = index, y = X.gap)) +
  geom_point(aes(color = good_gap)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = gap_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Gap %",
    subtitle = "F"
  )

gappcent_p

svg("busco_gappcent_p.svg", width = 8, height = 6)  
print(gappcent_p)  
dev.off()  


#guardem multipanel

svg("busco_multipanel.svg", width = 20, height = 20)  


grid.arrange(Cpcent_p, Spcent_p, Dpcent_p, Fpcent_p, Mpcent_p, gappcent_p, ncol = 2, nrow = 3)
dev.off()

