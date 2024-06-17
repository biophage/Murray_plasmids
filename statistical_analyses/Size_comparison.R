# Loading libraries that will be used
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggthemes)

#Loading plasmid data
#Note that the file "mt.tsv" is equivalent to Table S1 from the paper
mt<-read.csv("mt.tsv", sep="\t")
mt_comparison<-mt[which(mt$Lineage_type=="PAEModern_NoAMR" | mt$Lineage_type=="PAEModern_AMR"),] #Choosing the lineages that have both PAE and modern plasmids

Categories.shared<-unique(as.factor(mt_comparison$Lineagef)) #Vector of lineages that contain both PAE and modern plasmids

#Create empty vectors/dataframes
dif.PAE_i_t<-NULL
dif.mod_i_t<-NULL
dif.PAE.mod_i_t<-NULL
Size_DF_i<-data.frame() #Create an empty dataframe to store results of PAE vs modern plasmid
Size_DF_all<-data.frame()
STAT.size<-data.frame()
STAT.size_all<-data.frame()

#Comparison per lineage

for (i in 1:88) { # Number of lineages shared between PAE and modern plasmids
  PAE_i.t1 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "PAE") %>% #Extract only PAE plasmids from lineage i
    select(size) %>% #Select the column containing plasmids size
    slice_sample(n=100, replace = T) %>% #randomly sample with replacement 100 elements
    pull(size)
  
  PAE_i.t2 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "PAE") %>% 
    select(size) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(size)
  
  dif.PAE_i_t<-abs(PAE_i.t1-PAE_i.t2) #vector of pairwise size difference between PAE plasmids in lineage i
  
  
  #Sampling of modern plasmids
  mod_i.t1 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "Modern") %>%  #Extract only modern plasmids from lineage i
    select(size) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(size)
  
  mod_i.t2 <- mt_comparison %>% 
    filter(Lineagef == Categories.shared[i] & Category == "Modern") %>% 
    select(size) %>% 
    slice_sample(n=100, replace = T) %>% 
    pull(size)
  
  dif.mod_i_t<-abs(mod_i.t1-mod_i.t2) #vector of pairwise size difference between modern plasmids in lineage i
  
  
  dif.PAE.mod_i_t<-abs(PAE_i.t1-mod_i.t2) #vector of pairwise size difference between PAE and modern plasmids in lineage i
  
#Create a dataframe with the size differences
Size_DF_i<-data.frame(Lineage=Categories.shared[i],
                        dif.PAE=dif.PAE_i_t,
                        dif.Mod=dif.mod_i_t,
                        dif.PAEvsmod=dif.PAE.mod_i_t)  

#Accumulate the results for all the lineages in one dataframe
Size_DF_all<-bind_rows(Size_DF_all, Size_DF_i)

#Statistical tests
testPAE_vs_mod<-ks.test(dif.PAE_i_t, dif.mod_i_t)
testPAE_vs_PAE.mod<-ks.test(dif.PAE_i_t, dif.PAE.mod_i_t)
testmod_vs_PAE.mod<-ks.test(dif.mod_i_t, dif.PAE.mod_i_t)

#Stats results per lineage
STAT.size<-data.frame(Lineage= Categories.shared[i],
                      P.value_PAEvsmod=testPAE_vs_mod$p.value,
                      P.value_PAEvsPAEvmod=testPAE_vs_PAE.mod$p.value,
                      P.value_modvsPAEvmod=testmod_vs_PAE.mod$p.value)

#Stats results for all the lineages
STAT.size_all<-bind_rows(STAT.size_all, STAT.size)

}
write.table(STAT.size_all,file="Size_difference_STATS", sep="\t", row.names = FALSE)


#Switiching format of size differences dataframe to long version for plots
Size_DF_all_long <- pivot_longer(Size_DF_all, cols = "dif.PAE":"dif.PAEvsmod", names_to = "Comparison",
                             values_to = "Size.diference", values_drop_na = TRUE)



#All plots
#The plots can be found in the "Size_comparison_Allplots" directory in GitHub
Lineages<-unique(as.factor(Size_DF_all_long$Lineage))

for (i in 1:88) { #Number of lineages
  Size_for_plot <- Size_DF_all_long %>% 
    filter(Lineage == Lineages[i])
  
  name<-as.character(paste(Lineages[i],".svg"))
  
  
  p<-ggplot(Size_for_plot, aes(x=Comparison, y=(Size.diference+1), coulour=Comparison))+
    geom_violin(trim=FALSE,fill="grey")+
    geom_jitter(shape=16, position = position_jitter(0.2), alpha=0.3)+
    theme_classic()+
    labs(x="Plasmid origin", y="Size difference")+
    theme(axis.text=element_text(size=10, vjust=0.5))+
    theme(axis.title.x = element_text(vjust = 0.5, size = 11))+
    theme(axis.title.y = element_text(vjust = 0.5, size = 11))+
    theme(aspect.ratio = 0.5)+
    theme(legend.text = element_text(size=10))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(1, 10^7))+
    ggtitle(paste("CF", Lineages[i]))
  
  svg(filename = paste0("Plot_CF", Lineages[i], ".svg"), width = 7, height = 5)
  print(p)
  dev.off()
  
}


#Examples for Fig. S10: lineages 17 and 35, 107

ggplot(Size_DF_all_long[which(Size_DF_all_long$Lineage==17 | Size_DF_all_long$Lineage==35 |Size_DF_all_long$Lineage==107),], aes(x=Comparison, y=(Size.diference+1), coulour=Comparison))+
  geom_violin(trim=FALSE,fill="grey")+
  geom_jitter(shape=16, position = position_jitter(0.2), alpha=0.3)+
  theme_few()+
  labs(x="Comparison", y="Size difference")+
  theme(axis.text=element_text(size=8, vjust=0.5))+
  theme(axis.title.x = element_text(vjust = 0.5, size = 11))+
  theme(axis.title.y = element_text(vjust = 0.5, size = 11))+
  theme(aspect.ratio = 1.5)+
  theme(legend.text = element_text(size=10))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits=c(1, 10^7))+
  facet_wrap(.~Lineage)
  