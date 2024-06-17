# Loading libraries that will be used
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggthemes)

#Loading plamids data
#Note that the file "mt.tsv" is equivalent to Table S1 from the paper
mt<-read.csv("mt.tsv", sep="\t")
mt_comparison<-mt[which(mt$Lineage_type=="PAEModern_NoAMR" | mt$Lineage_type=="PAEModern_AMR"),] #Choosing the lineages that have both PAE and modern plasmids

Categories.shared<-unique(as.factor(mt_comparison$Lineagef)) #Vector of lineages that contain both PAE and modern plasmids

#Create an empty vectors/dataframes
Vir_dif.PAE_i_t<-NULL
Vir_dif.mod_i_t<-NULL
Vir_dif.PAE.mod_i_t<-NULL
Vir_DF_i<-data.frame() #Create an empty dataframe to store results of PAE vs modern plasmid
Vir_DF_all<-data.frame()
STAT.vir<-data.frame()
STAT.vir_all<-data.frame()

#Comparison per lineage

for (i in 1:88) { # Number of lineages shared between PAE and modern plasmids
  Vir_PAE_i.t1 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "PAE") %>% #Extract only PAE plasmids from lineage i
    select(Vir_genes) %>% #Select the column containing plasmids virulence gene number
    slice_sample(n=100, replace = T) %>% #randomly sample with replacement 100 elements
    pull(Vir_genes)
  
  Vir_PAE_i.t2 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "PAE") %>% 
    select(Vir_genes) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(Vir_genes)
  
  Vir_dif.PAE_i_t<-abs(Vir_PAE_i.t1-Vir_PAE_i.t2) #vector of pairwise virulence genes difference between PAE plasmids in lineage i
  
  
  #Sampling of modern plasmids
  Vir_mod_i.t1 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "Modern") %>%  #Extract only modern plasmids from lineage i
    select(Vir_genes) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(Vir_genes)
  
  Vir_mod_i.t2 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "Modern") %>% 
    select(Vir_genes) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(Vir_genes)
  
  Vir_dif.mod_i_t<-abs(Vir_mod_i.t1-Vir_mod_i.t2) #vector of pairwise virulence gene difference between modern plasmids in lineage i
  
  
  Vir_dif.PAE.mod_i_t<-abs(Vir_PAE_i.t1-Vir_mod_i.t2) #vector of pairwise virulence gene difference between PAE and modern plasmids in lineage i
  
  
  Vir_DF_i<-data.frame(Lineage=Categories.shared[i],
                        Vir_dif.PAE=Vir_dif.PAE_i_t,
                        Vir_dif.Mod=Vir_dif.mod_i_t,
                        Vir_dif.PAEvsmod=Vir_dif.PAE.mod_i_t)  
  
  Vir_DF_all<-bind_rows(Vir_DF_all, Vir_DF_i)
  
  Vir_testPAE_vs_mod<-ks.test(Vir_dif.PAE_i_t, Vir_dif.mod_i_t)
  Vir_testPAE_vs_PAE.mod<-ks.test(Vir_dif.PAE_i_t, Vir_dif.PAE.mod_i_t)
  Vir_testmod_vs_PAE.mod<-ks.test(Vir_dif.mod_i_t, Vir_dif.PAE.mod_i_t)
  
  STAT.vir<-data.frame(Lineage= Categories.shared[i],
                        Vir_P.value_PAEvsmod=Vir_testPAE_vs_mod$p.value,
                       Vir_P.value_PAEvsPAEvmod=Vir_testPAE_vs_PAE.mod$p.value,
                       Vir_P.value_modvsPAEvmod=Vir_testmod_vs_PAE.mod$p.value)
  STAT.vir_all<-bind_rows(STAT.vir_all, STAT.vir)
  
}

STAT.vir_all #Dataframe with p-values for comparisons of each lineage

write.table(STAT.vir_all,file="Virulence_difference_STATS", sep="\t", row.names = FALSE)
