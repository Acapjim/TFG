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


input_QC_dataset1 = paste(dir, "bin_stats_ext_efaecium_tsv_results.tsv", sep = "")
checkm_dataset = read.delim(input_QC_dataset1, sep = "\t")

genomas_sindu <- readLines("lista_genomas_sin_duplicaciones.txt")
genomas_sindu <- sub("\\.fna$", "", genomas_sindu)

#filtramos sin duplicaciones
checkm_dataset <- subset(checkm_dataset, acc %in% genomas_sindu)

###################################################################################################################
#multi-panel plot for QC metrics 
###################################################################################################################
x_axis_size <- 1
size_axis_lines <-1 
##checkm

### Genomic complete

Geno_comp = checkm_dataset %>% select(acc,Completeness)

#### Graficar genomic completeness
Geno_comp = Geno_comp %>% arrange(Completeness)
Geno_comp$good_comp = "no"
comp_cutoff = 98
Geno_comp$good_comp[which(as.numeric(Geno_comp$Completeness) > comp_cutoff)] = "yes"
Geno_comp$good_comp= factor(Geno_comp$good_comp)

Geno_comp <- Geno_comp %>%
  mutate(index = seq_along(Geno_comp[[1]]))

Geno_comp_p = ggplot(Geno_comp, aes(x = index, y = Completeness)) +
  geom_point(aes(color = good_comp)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = comp_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Completeness %",
    subtitle = "A"
  )
Geno_comp_p
svg("checkm_Completeness.svg", width = 8, height = 6)  
print(Geno_comp_p)  
dev.off()  


### Genomic contamination

Geno_cont = checkm_dataset %>% select(acc,Contamination)

#### Graficar genomic contamination
Geno_cont = Geno_cont %>% arrange(Contamination)
Geno_cont$good_cont = "no"
cont_cutoff = 1
Geno_cont$good_cont[which(as.numeric(Geno_cont$Contamination) < cont_cutoff)] = "yes"
Geno_cont$good_cont= factor(Geno_cont$good_cont)


Geno_cont <- Geno_cont %>%
  mutate(index = seq_along(Geno_cont[[1]]))

Geno_cont_p = ggplot(Geno_cont, aes(x = index, y = Contamination)) +
  geom_point(aes(color = good_cont)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = cont_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Contamination %",
    subtitle = "B"
  )
Geno_cont_p

svg("checkm_Contamination.svg", width = 8, height = 6)  
print(Geno_cont_p)  
dev.off()  

### N50 contig
N50_cont = checkm_dataset %>% select(acc,N50..contigs.)

#### Graficar genomic lenght

N50_cont = N50_cont %>% arrange(N50..contigs.)
n50_cutoff = 2416265
N50_cont$good_n50 = "no"
N50_cont$good_n50[which(as.numeric(N50_cont$N50..contigs.) > n50_cutoff)] = "yes"
N50_cont$good_n50 = factor(N50_cont$good_n50)


N50_cont <- N50_cont %>%
  mutate(index = seq_along(N50_cont[[1]]))

N50_cont_p = ggplot(N50_cont, aes(x = index, y =N50..contigs. )) +
  geom_point(aes(color = good_n50)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = n50_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
  theme_base +
  labs(
    x = "Isolate Index",
    y = "n50 contig %",
    subtitle = "C"
  )
N50_cont_p

svg("checkm_n50.svg", width = 8, height = 6)  
print(N50_cont_p)  
dev.off()  


### Genome size
Geno_len = checkm_dataset %>% select(acc,Genome.size)

#### Graficar genomic lenght
#### checkm_Genome.size
Geno_len = Geno_len %>% arrange(Genome.size)

Geno_len$good_lenght = "no"
# Calcular la media
media_lenght <- mean(checkm_dataset$Genome.size, na.rm = TRUE)  

# Calcular el intervalor permitido
intervalo = 10
lmax <- media_lenght * (1+intervalo/100)
lmin <- media_lenght * (1-intervalo/100)

Geno_len$good_lenght[which(as.numeric(Geno_len$Genome.size) > lmin & as.numeric(Geno_len$Genome.size) < lmax)] <- "yes"
Geno_len$good_lenght= factor(Geno_len$good_lenght)


Geno_len <- Geno_len %>%
  mutate(index = seq_along(Geno_len[[1]]))

Geno_len_p = ggplot(Geno_len, aes(x = index, y =Genome.size )) +
  geom_point(aes(color = good_lenght)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = lmax, linetype = "dashed", color = "red", linewidth = 1) + geom_hline(yintercept = lmin, linetype = "dashed", color = "red", linewidth = 1)
  theme_base +
  labs(
    x = "Isolate Index",
    y = "Genome size",
    subtitle = "Checkm"
  )
Geno_len_p


svg("checkm_Genomic_lenght.svg", width = 8, height = 6)  
print(Geno_len_p)  
dev.off()



### Nombre contig
contig_n = checkm_dataset %>% select(X..contigs,acc)

#### Graficar 
contig_n = contig_n %>% arrange(X..contigs)
ncon_cutoff = 18
contig_n$good_ncon = "no"
contig_n$good_ncon[which(as.numeric(contig_n$X..contigs) < ncon_cutoff)] = "yes"
contig_n$good_ncon = factor(contig_n$good_ncon)


contig_n <- contig_n %>%
  mutate(index = seq_along(contig_n[[1]]))


contig_n_p = ggplot(contig_n, aes(x = index, y =X..contigs)) +
  geom_point(aes(color = good_ncon)) +
  scale_color_manual(breaks = c("no", "yes"),
                     values = c("black", "red")) +
  geom_hline(yintercept = ncon_cutoff, linetype = "dashed", color = "red", linewidth = 1) +
theme_base +
  labs(
    x = "Isolate Index",
    y = "Nº of contigs",
    subtitle = "D"
  )
contig_n_p

svg("checkm_n_contigs.svg", width = 8, height = 6)  
print(contig_n_p)  
dev.off()  



svg("checkm_multipanel.svg", width = 20, height = 16)  
grid.arrange(Geno_comp_p, Geno_cont_p, N50_cont_p, contig_n_p, ncol = 2, nrow = 2)
dev.off()

