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
Transp_dif.PAE_i_t<-NULL
Transp_dif.mod_i_t<-NULL
Transp_dif.PAE.mod_i_t<-NULL
Transp_DF_i<-data.frame() #Create an empty dataframe to store results of PAE vs modern plasmid
Transp_DF_all<-data.frame()
STAT.Transp<-data.frame()
STAT.Transp_all<-data.frame()

#Comparison per lineage

for (i in 1:88) { # Number of lineages shared between PAE and modern plasmids
  Transp_PAE_i.t1 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "PAE") %>% #Extract only PAE plasmids from lineage i
    select(Transposases) %>% #Select the column containing plasmids transposase gene number
    slice_sample(n=100, replace = T) %>% #randomly sample with replacement 100 elements
    pull(Transposases)
  
  Transp_PAE_i.t2 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "PAE") %>% 
    select(Transposases) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(Transposases)
  
  Transp_dif.PAE_i_t<-abs(Transp_PAE_i.t1-Transp_PAE_i.t2) #vector of pairwise transposase genes difference between PAE plasmids in lineage i
  
  
  #Sampling of modern plasmids
  Transp_mod_i.t1 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "Modern") %>%  #Extract only modern plasmids from lineage i
    select(Transposases) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(Transposases)
  
  Transp_mod_i.t2 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "Modern") %>% 
    select(Transposases) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(Transposases)
  
  Transp_dif.mod_i_t<-abs(Transp_mod_i.t1-Transp_mod_i.t2) #vector of pairwise transposase gene difference between modern plasmids in lineage i
  
  
  Transp_dif.PAE.mod_i_t<-abs(Transp_PAE_i.t1-Transp_mod_i.t2) #vector of pairwise transposase gene difference between PAE and modern plasmids in lineage i
  
  
  Transp_DF_i<-data.frame(Lineage=Categories.shared[i],
                       Transp_dif.PAE=Transp_dif.PAE_i_t,
                       Transp_dif.Mod=Transp_dif.mod_i_t,
                       Transp_dif.PAEvsmod=Transp_dif.PAE.mod_i_t)  
  
  Transp_DF_all<-bind_rows(Transp_DF_all, Transp_DF_i)
  
  Transp_testPAE_vs_mod<-ks.test(Transp_dif.PAE_i_t, Transp_dif.mod_i_t)
  Transp_testPAE_vs_PAE.mod<-ks.test(Transp_dif.PAE_i_t, Transp_dif.PAE.mod_i_t)
  Transp_testmod_vs_PAE.mod<-ks.test(Transp_dif.mod_i_t, Transp_dif.PAE.mod_i_t)
  
  STAT.Transp<-data.frame(Lineage= Categories.shared[i],
                       Transp_P.value_PAEvsmod=Transp_testPAE_vs_mod$p.value,
                       Transp_P.value_PAEvsPAEvmod=Transp_testPAE_vs_PAE.mod$p.value,
                       Transp_P.value_modvsPAEvmod=Transp_testmod_vs_PAE.mod$p.value)
  STAT.Transp_all<-bind_rows(STAT.Transp_all, STAT.Transp)
  
}

STAT.Transp_all #Dataframe with p-values for comparisons of each lineage

write.table(STAT.Transp_all,file="transposase_difference_STATS", sep="\t", row.names = FALSE)
