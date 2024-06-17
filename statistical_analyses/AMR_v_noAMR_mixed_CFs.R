#AMR vs no AMR mixed CFs comparisons

# Loading libraries that will be used
library(tidyverse)

#Loading plasmid data
#The "Heatmap_annotations.tsv" file contains data used to produce Fig. 4 from the paper (See PAE_plasmids-Figures.R file in GitHub)
Data_heatmap<-read.csv("Heatmap_annotations.tsv", sep = "\t")

#Size
Size_AMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_AMR") %>% 
  select(size_mean) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(size_mean)

Size_noAMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_NoAMR") %>% 
  select(size_mean) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(size_mean)

test_Size.AMR_noAMR<-ks.test(Size_AMR, Size_noAMR)
test_Size.AMR_noAMR$p.value


#Virulence genes
Vir_AMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_AMR") %>% 
  select(Vir_genes_mean) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(Vir_genes_mean)

Vir_noAMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_NoAMR") %>% 
  select(Vir_genes_mean) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(Vir_genes_mean)

test_Vir.AMR_noAMR<-ks.test(Vir_AMR, Vir_noAMR)
test_Vir.AMR_noAMR$p.value


#Transposases
Transposase_AMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_AMR") %>% 
  select(Transposases_mean) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(Transposases_mean)

Transposase_NoAMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_NoAMR") %>% 
  select(Transposases_mean) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(Transposases_mean)

test_Tran.AMR_noAMR<-ks.test(Transposase_AMR, Transposase_NoAMR)
test_Tran.AMR_noAMR$p.value


#Core genome
Core_genome_AMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_AMR") %>% 
  select(Core_genome_frc) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(Core_genome_frc)

Core_genome_NoAMR<-Data_heatmap %>% 
  filter(CF_type == "PAEModern_NoAMR") %>% 
  select(Core_genome_frc) %>% 
  slice_sample(n=100, replace = T) %>% 
  pull(Core_genome_frc)

test_Core.AMR_noAMR<-ks.test(Core_genome_AMR, Core_genome_NoAMR)
test_Core.AMR_noAMR$p.value
