### This file describes the R code used to produce main and supplementary figures for the paper "Pre and Post antibiotic epoch: insights into the historical spread of antimicrobial resistance"

# Figures list

# Main
# Fig. 3. Mobility, AMR gene carriage and size distribution across plasmid CFs.
# Fig. 4. Evolution of PAE plasmids that have maintained a clonal frame.

# Supplementary
# Fig. S1. Plasmids identified in genomes of the Murray collection.
# Fig. S2. Characteristics of PAE plasmids and their Modern relatives.
# Fig. S3. CFs' pangenomes.
# Fig. S5. Murray plasmids containment within Modern-CF plasmids.
# Fig. S6. Replicons identified in PAE plasmids and their modern close relatives.
# Fig. S8. Contribution of PAE-related plasmids to AMR.
# Fig. S9.  Distribution of plasmid mobility types across CFs.
# Fig. S10. Differences between PAE and Modern plasmids' average size within mixed CFs.
# Fig. S11. Relationship of Murray plasmids to all sequences from the integrated database.
# Fig. S12. Comparison of characteristics between plasmids related and unrelated to Murray plasmids.
# Fig. S14. Examples of putative Murray plasmid sequences identified as contained in bacterial chromosomes. <- Data from BLASTn online.

# The remaining figures from the paper, listed below, were not produced using R
# Fig. 1. Framework for the analysis of PAE plasmids and their modern relatives. <- Made with BioRender
# Fig. 2. Diversity of PAE and related modern plasmids and their contribution to AMR. <- Networks made with Cytoscape, integrated in Inkscape
# Fig. 5. Landscape of known plasmids and their links to Murray plasmids and AMR. <- Networks made with Cytoscape, integrated in Inkscape
# Fig. S4. Distribution of plasmids in close families (CFs). <- Made with BioRender
# Fig. S7. Pairwise comparisons between PAE plasmids and their modern relatives. <- Maps made with Easyfig
# Fig. S13. Characterisation and classification of plasmid sequences from Murray genomes. <- Made in BioRender



### Make output directories
dir.create("Figures/Main", recursive = T)
dir.create("Figures/Supplementary")



### Libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(scales)
library(ggpattern)
library(treemapify)
library(ggExtra)



### Figure 3
dir.create("Figures/Main/Figure3", recursive = T)

# Data used
# mt: Plasmids master table <- Table S1
# CFsExt_mt: CFs master table <- Table S2
# CFs_df: Table of CFs split by PAE/Modern plasmids and mobility. Obtained from parsing "mt"

# Mobility distribution
# Mobility frequency as proportion - bar plot
CFsMobdist <- ggplot(CFsExt_mt, aes(x=CF_type, fill=top_mobility)) +
  labs(x="CF type", y="Relative proportion of CFs", fill="Mobility") +
  theme_linedraw() + geom_bar(position = "fill", width = 0.55) +
  scale_fill_manual(values = c("#008B45","#00E5BF","#002147")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "top",
        text = element_text(family = "serif")) +
  scale_x_discrete(labels = c("Mixed AMR+", "Mixed AMR-",
                              "Modern AMR+", "Modern AMR-",
                              "PAE AMR-"))

# AMR genes distribution 
# Average number of AMR genes per lineage - scatter plot
fnames <- as_labeller(c(
  `PAEModern_AMR` = "Mixed-CFs AMR+ (n=39)",
  `PAEModern_NoAMR` = "Mixed-CFs AMR- (n=49)",
  `Divergent_modern_AMR` = "Modern-CFs AMR+ (n=203)",
  `Divergent_modern_NoAMR` = "Modern-CFs AMR- (n=546)",
  `Divergent_PAE_NoAMR` = "PAE-CFs AMR- (n=103)"
))
CFsAMRdist <- CFs_df %>% 
  filter(CF_type == "PAEModern_AMR" | CF_type == "Divergent_modern_AMR") %>% 
  ggplot(aes(x=predicted_mobility, y=mean_AMR_genes, colour = Category)) +
  labs(x="Mobility", y="Average number of AMR genes per CF", colour="Plasmid category") +
  theme_linedraw() + geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(aes(colour = Category), method = "pseudorandom", alpha = 0.5, dodge.width = 0.75,
                   groupOnX=T) +
  scale_colour_manual(values = c("steelblue", "indianred")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "top",
        text = element_text(family = "serif")) +
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames)

# Size distribution 
# Average plasmid size per lineage - scatter plot
CFsSizedist <- CFs_df %>% 
  filter(CF_type == "PAEModern_AMR" | CF_type == "Divergent_modern_AMR") %>% 
  ggplot(aes(x=predicted_mobility, y=mean_size, colour = Category)) +
  labs(x="Mobility", y="Average plasmid size per CF", colour="Plasmid category") +
  theme_linedraw() + geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(method = "quasirandom", alpha = 0.5, dodge.width = 0.75, groupOnX=T) +
  scale_colour_manual(values = c("steelblue", "indianred")) + theme_linedraw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "top",
        text = element_text(family = "serif")) +
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames)

# Combine and plot
plot_grid(CFsMobdist, CFsAMRdist, CFsSizedist, labels = "AUTO", 
          nrow = 1, align = "hv", axis = "b", rel_widths = c(1, 1.2, 1.2),
          hjust = -1)
ggsave2('Figures/Main/Figure3/Figure3_raw.svg', scale = 0.6, 
        width = 22, height = 7.5, device = "svg")
rm(CFsMobdist, CFsAMRdist, CFsSizedist)



### Figure 4
dir.create("Figures/Main/Figure4", recursive = T)

# Data used
# CFsExt_mt: CFs master table <- Table S2
# PAEModernCFs_sr: Table of the size ratio metrics for PAE/Modern plasmids clustered in the same CF
# PAEModernCFs_fd: Table of the cargo gene metrics for PAE/Modern plasmids clustered in the same CF

# Heatmap

# Annotation data
col_anno <- CFsExt_mt %>% 
  filter(CF_type == "PAEModern_AMR" | CF_type == "PAEModern_NoAMR" ) %>%
  select(CF,CF_type,top_mobility,top_host_range,AMR_genes_mean,Vir_genes_mean,
         Transposases_mean,Integrons_mean,size_mean,R_plasmidfreq,
         PAE_freq,Plas_total)
pangenome_anno <- CFs_pangenome %>% 
  filter(CF_type == "PAEModern_AMR" | CF_type == "PAEModern_NoAMR" ) %>%
  select(CF, Core_genome_frc)
col_anno <- left_join(col_anno, pangenome_anno, by='CF')
col_anno <- col_anno[order(-col_anno$AMR_genes_mean),]
Names <- CFsExt_mt %>% 
  filter(CF_type == "PAEModern_AMR" | CF_type == "PAEModern_NoAMR" ) %>%
  select(CF,AMR_genes_mean) %>% 
  arrange(desc(AMR_genes_mean)) %>% 
  select(CF)

# 1: Modern vs PAE size ratio
SizeR <- PAEModernCFs_sr %>% 
  select(CF,size_ratio_mean) %>% 
  remove_rownames %>% column_to_rownames(var="CF") %>% as.matrix %>% t
col_ha_sr <- HeatmapAnnotation(CF_type = col_anno$CF_type, mob = col_anno$top_mobility, 
                               col = list(
                                 CF_type = c("PAEModern_AMR" = "#FF6619","PAEModern_NoAMR" = "#762A83"),
                                 mob = c("conjugative" = "#002147","non-mobilizable" = "#00E5BF","mobilizable" = "#008B45")),
                               mean_size = anno_points(col_anno$size_mean,
                                                       height = unit(1.75, "cm"),
                                                       ylim = c(0, 3e+05),
                                                       axis_param = list(at = c(0,1e+05,2e+05,2.8e+05),
                                                                         labels = c("0","1e+05","2e+05","2.8e+05"))),
                               border = T,
                               gap = unit(c(1,1), "mm"),
                               simple_anno_size = unit(0.65, "cm"))
hm_sr <- Heatmap(SizeR,
                 col = colorRamp2(c(0,0.8,1,1.2,3), c("skyblue3","white","white","white","hotpink3")),
                 show_heatmap_legend = TRUE,
                 name = "Size ratio",
                 heatmap_legend_param = list(at = c(0.4,0.7,1,1.3,2,3), direction = "horizontal",
                                             labels_rot = 90),
                 top_annotation = col_ha_sr,
                 row_names_side = "left",
                 row_labels = c("Size ratio"),
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 column_labels = as.matrix(Names),
                 show_column_names = TRUE,
                 column_names_side="bottom",
                 column_names_gp=gpar(cex=0.95, fontface="plain"),
                 height = unit(0.65,"cm"))
rm(SizeR, col_ha_sr)

# 2: Modern vs PAE AMR genes mean difference
AMRdiff <- PAEModernCFs_fd %>% 
  select(CF,AMR_genes_mean_diff) %>% 
  remove_rownames %>% column_to_rownames(var="CF") %>% as.matrix %>% t
col_ha_amrdiff <- HeatmapAnnotation(Core_gen = anno_barplot(col_anno$Core_genome_frc, bar_width = 1, gp = gpar(fill = "#ADB6B6", lwd = 0.5)),
                                    R_freq = anno_barplot(col_anno$R_plasmidfreq, bar_width = 1, gp = gpar(fill = "#ADB6B6", lwd = 0.5)),
                                    AMR_mean = anno_barplot(col_anno$AMR_genes_mean, bar_width = 1, gp = gpar(fill = "#ADB6B6", lwd = 0.5)),
                                    Trans_mean = anno_barplot(col_anno$Transposases_mean, bar_width = 1, gp = gpar(fill = "#ADB6B6", lwd = 0.5)),
                                    Inte_mean = anno_barplot(col_anno$Integrons_mean, bar_width = 1, gp = gpar(fill = "#ADB6B6", lwd = 0.5)),
                                    Vir_mean = anno_barplot(col_anno$Vir_genes_mean, bar_width = 1, gp = gpar(fill = "#ADB6B6", lwd = 0.5)),
                                    border = F,
                                    gap = unit(c(1.5,1.5,1.5,1.5,1.5,1.5), "mm"))
hm_amrdiff <- Heatmap(AMRdiff,
                      col = colorRamp2(c(-0.6,0,1,7,14), c("#0074A6","white","gold2","indianred","darkred")),
                      show_heatmap_legend = TRUE,
                      name = "AMR mean diff",
                      heatmap_legend_param = list(at = c(-0.6,0,3.5,7,10.5,14), direction = "horizontal",
                                                  labels_rot = 90),
                      top_annotation = col_ha_amrdiff,
                      row_names_side = "left",
                      row_labels = c("AMR genes"),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      column_labels = as.matrix(Names),
                      show_column_names = TRUE,
                      column_names_side="bottom",
                      column_names_gp=gpar(cex=0.95, fontface="plain"),
                      column_split = data.frame(col_anno$CF_type),
                      cluster_column_slices = FALSE,
                      column_title = NULL,
                      height = unit(0.65,"cm"))
rm(AMRdiff,col_ha_amrdiff)

# 3: Modern vs PAE Vir genes mean difference
Virdiff <- PAEModernCFs_fd %>% 
  select(CF,Vir_genes_mean_diff) %>% 
  remove_rownames %>% column_to_rownames(var="CF") %>% as.matrix %>% t
hm_virdiff <- Heatmap(Virdiff,
                      col = colorRamp2(c(-2,0,1,12,24), c("#0074A6","white","gold2","indianred","darkred")),
                      show_heatmap_legend = TRUE,
                      name = "Vir mean diff",
                      heatmap_legend_param = list(at = c(-2,0,6,12,18,24), direction = "horizontal",
                                                  labels_rot = 90),
                      row_names_side = "left",
                      row_labels = c("Vir genes"),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      column_labels = as.matrix(Names),
                      show_column_names = TRUE,
                      column_names_side="bottom",
                      column_names_gp=gpar(cex=0.95, fontface="plain"),
                      height = unit(0.65,"cm"))
rm(Virdiff)

# 4: Modern vs PAE transposases mean difference
Trandiff <- PAEModernCFs_fd %>% 
  select(CF,Transposases_mean_diff) %>% 
  remove_rownames %>% column_to_rownames(var="CF") %>% as.matrix %>% t
hm_trandiff <- Heatmap(Trandiff,
                       col = colorRamp2(c(-13,0,1,28,57), c("#0074A6","white","gold2","indianred","darkred")),
                       show_heatmap_legend = TRUE,
                       name = "Trans mean diff",
                       heatmap_legend_param = list(at = c(-13,0,14,28,42,58), direction = "horizontal",
                                                   labels_rot = 90),
                       row_names_side = "left",
                       row_labels = c("Transposases"),
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       column_labels = as.matrix(Names),
                       show_column_names = TRUE,
                       column_names_side="bottom",
                       column_names_gp=gpar(cex=0.95, fontface="plain"),
                       height = unit(0.65,"cm"))
rm(Trandiff)

# 5: Modern vs PAE integrons mean difference
Intdiff <- PAEModernCFs_fd %>% 
  select(CF,Integrons_mean_diff) %>% 
  remove_rownames %>% column_to_rownames(var="CF") %>% as.matrix %>% t
col_ha_intdiff <- HeatmapAnnotation(PAE_frq = col_anno$PAE_freq,
                                    Total_plasmids = col_anno$Plas_total,
                                    host_range = col_anno$top_host_range,
                                    col = list(PAE_frq = circlize::colorRamp2(c(0,1), hcl_palette = "Mako", reverse = T),
                                               Total_plasmids = circlize::colorRamp2(c(2,732), hcl_palette = "Mako", reverse = T),
                                               host_range = c("multi-phylla"="#725663",
                                                              "class"="#8A8B79",
                                                              "order"="#79AF97",
                                                              "family"="#ADB17D",
                                                              "genus"="#D6D6CE",
                                                              "ND"="#D49464")),
                                    annotation_legend_param = list(
                                      PAE_frq = list(at = c(0.004,0.3,0.6,0.93), direction = "horizontal",
                                                     labels_rot = 90),
                                      Total_plasmids = list(at = c(2,150,300,450,600,732), direction = "horizontal",
                                                            labels_rot = 90),
                                      host_range = list(ncol = 3)),
                                    border = c(F,F,T),
                                    gap = unit(c(1,1,1), "mm"),
                                    simple_anno_size = unit(0.65, "cm"))
hm_intdiff <- Heatmap(Intdiff,
                      col = colorRamp2(c(0,1,1.3), c("white","gold2","indianred")),
                      show_heatmap_legend = TRUE,
                      name = "Int mean diff",
                      heatmap_legend_param = list(at = c(0,1,0.5,1.3), direction = "horizontal",
                                                  labels_rot = 90),
                      bottom_annotation = col_ha_intdiff,
                      row_names_side = "left",
                      row_labels = c("Integrons"),
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      column_labels = as.matrix(Names),
                      show_column_names = TRUE,
                      column_names_side="bottom",
                      column_names_gp=gpar(cex=0.95, fontface="plain"),
                      height = unit(0.65,"cm"))
rm(Intdiff, col_ha_intdiff)

# Combine heatmaps
hm_list <- hm_sr %v% hm_amrdiff %v% hm_virdiff %v% hm_trandiff %v% hm_intdiff

# Draw composed heatmap
svg(filename = 'Figures/Main/Figure4/Heatmap-raw.svg', width = 16, height = 9)
draw(hm_list, ht_gap = unit(1, "mm"), main_heatmap = 2)
dev.off()

# Write table for statistical analysis of heatmap comparisons
write.table(col_anno, "Supplementary_Tables/Heatmap_annotations.tsv", 
            sep = "\t", row.names = F, quote = F)

# Remove objects
rm(Names, col_anno, hm_sr, hm_amrdiff, hm_virdiff, hm_trandiff, hm_intdiff, hm_list)

# AMR vs NoAMR Mixed-CFs mean size - plot to be used as inset in the heatmap
CFsExt_mt %>% 
  filter(CF_type == "PAEModern_AMR" | CF_type == "PAEModern_NoAMR") %>% 
  ggplot(aes(x=CF_type, y=size_mean, fill=CF_type)) +
  geom_violin(scale = "width", adjust= 0.5, alpha = 0.9) +
  scale_fill_manual(values = c("#FF6619", "#762A83")) +
  theme_classic() + geom_boxplot(width=0.2, color="black", alpha=0.25, 
                                 position = position_dodge(width=0.9), lwd=0.75) +
  scale_y_continuous(breaks = c(0,100000,200000,280000), labels = scientific) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(),
        legend.position = "none")
ggsave('Figures/Main/Figure4/Avr_size_distribution.svg',  scale = 0.75, device = "svg")

# The svg figure files were integrated and edited in Inkscape to make fhe final version of the figure



### Figure S1
dir.create("Figures/Supplementary/FigureS1/", recursive = T)

# Data used
# mt: Plasmids master table <- Table S1

# Murray plasmids dataframe
murray <- mt %>% 
  filter(Source == "Murray") %>% 
  separate(ID, c("prefix","suffix"), "-", remove = F) %>% 
  separate(suffix, c("Strain","remove"), "_") %>% 
  select(-prefix, -remove, -Location, -Country, -Region)

# Cumulative plasmid size per strain
cumulative_size <- murray %>% 
  select(Strain,genus,size) %>% group_by(Strain,genus) %>% summarise(count=n(),sum(size)) %>% 
  mutate(genus=as.character(genus))

# Plasmid content per strain (Number of plasmids VS cumulative plasmid size)
CumSize <- ggplot(cumulative_size, aes(x= count, y=`sum(size)`, fill= genus)) +
  geom_point(aes(x= count, color= genus), alpha=1.5, size=0.5, stroke=1) +
  geom_violin(aes(x= count), width=1, alpha = 0.15, linewidth=0.25) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(breaks = c(1e3, 1e4, 2.5e4, 5e4, 1e5, 1.5e6, 2e5, 3.5e5), 
                     label = c("1", "10", "25", "50", "100", "150", "200", "350")) +
  scale_fill_manual(name = "Genus", values = c("Escherichia" = "#053061",
                                               "Shigella" = "#2166ac",
                                               "Salmonella" = "#8e0152",
                                               "Klebsiella" = "#ff0000",
                                               "Proteus" = "#cdcd00",
                                               "Raoultella" = "#A65628")) +
  scale_color_manual(values = c("Escherichia" = "#053061",
                                "Shigella" = "#2166ac",
                                "Salmonella" = "#8e0152",
                                "Klebsiella" = "#ff0000",
                                "Proteus" = "#cdcd00",
                                "Raoultella" = "#A65628")) +
  labs(x="Number of plasmids per strain", y="Cumulative plasmid size (kb)") +
  guides(colour = "none", fill =  guide_legend(nrow = 1)) +
  theme_linedraw() + theme(legend.position = "top", text = element_text(family = "serif")) +
  facet_wrap(factor(cumulative_size$genus, c("Escherichia", "Shigella", "Salmonella", "Klebsiella", "Proteus", "Raoultella")), nrow = 1)

# Murray plasmids size VS GC
SizeGC <- murray %>% select(ID, genus, predicted_mobility, size, gc) %>% 
  mutate(predicted_mobility=factor(predicted_mobility, levels = c("conjugative", "mobilizable","non-mobilizable"))) %>%
  ggplot(aes(x= genus, y= size, fill=genus)) +
  geom_jitter(aes(color= gc), alpha=1.5, size= 0.5, stroke=1) +
  geom_violin(aes(x=genus), width=1, alpha =0.15, size=0.25) +
  scale_y_continuous(trans='log10', breaks = c(1000, 2500, 5000, 10000, 25000, 50000, 100000, 150000, 200000), 
                     label = c("1", "2.5", "5", "10", "25", "50", "100", "150", "200")) +
  scale_x_discrete(breaks = NULL, labels = NULL) +
  scale_fill_manual(name = "Genus", values = c("Escherichia" = "#053061",
                                               "Shigella" = "#2166ac",
                                               "Salmonella" = "#8e0152",
                                               "Klebsiella" = "#ff0000",
                                               "Proteus" = "#cdcd00",
                                               "Raoultella" = "#A65628")) +
  scale_color_binned(name = "GC %", type = "viridis") +
  guides(fill = "none") +
  labs(x="", y="Plasmid size (kb)") +
  theme_linedraw() + theme(legend.position = c(0.89, 0.25), text = element_text(family = "serif")) +
  facet_wrap(predicted_mobility ~ factor(murray$genus, c("Escherichia", "Shigella", "Salmonella", "Klebsiella", "Proteus", "Raoultella")), scale = "free_x", ncol = 8)

# Combine and plot
plot_grid(CumSize, SizeGC, labels = "AUTO", 
          ncol = 1, align = "hv", axis = "b", hjust = -1, label_size = 20,
          label_fontfamily = "serif")
ggsave('Figures/Supplementary/FigureS1/FigureS1.png', scale = 0.5, 
       width = 22, height = 25)
rm(CumSize, SizeGC, murray)



### Figure S2
dir.create("Figures/Supplementary/FigureS2/", recursive = T)

# Data used
# mt: Plasmids master table <- Table S1

# Panel A: Year of isolation
p <- 
  ggplot(mt, aes(x=Isolation_year, y=size, colour=predicted_mobility)) + 
  geom_point(alpha=0.9) + labs(x="Year of isolation", y="Plasmid size (bp)", colour="Predicted mobility") +
  geom_smooth(method= "gam", linewidth=0.7, fill = "lightgrey") + scale_x_continuous(breaks = seq(1910,2022,10)) + 
  scale_colour_manual(values = c("#008B45","#00E5BF","#002147")) +
  scale_y_continuous(breaks = seq(0,600000,100000)) + theme_linedraw() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "bottom", text = element_text(family = "serif"))
A_Year <- ggMarginal(p, type="histogram", margins = "both", groupColour = T, groupFill = T, alpha = 0.9,
                   xparams = list(bins=70, size=0.01), yparams = list(bins=70, size=0.01))

# Panel B: Top genera (Set to NULL; added separately)
genera <- mt %>% 
  group_by(top_genera) %>% 
  summarise(count=n()) %>% 
  mutate(percent = count / length(mt$ID),
         Sample = "plasmid",
         top_genera = factor(top_genera, levels = c("Enterobacter","Shigella","Other","Salmonella","Klebsiella","Escherichia")))
ggplot(genera, aes(x=Sample, y=percent)) +
  geom_col_pattern(aes(pattern=top_genera, pattern_type=top_genera, pattern_shape=top_genera),
                   fill='white',
                   colour='black',
                   pattern_fill='darkgrey',
                   pattern_density=0.4,
                   pattern_angle = 45) + 
  scale_pattern_manual(values=c('pch', 'weave', 'pch','wave', 'pch', 'stripe')) +
  scale_pattern_type_manual(values=c('hex', NA, 'hex', NA, 'hex', NA)) +
  scale_pattern_shape_manual(values=c(1, NA, 15, NA, 3, NA)) +
  scale_pattern_spacing_discrete(range = c(0.01, 0.05)) + coord_flip() + 
  scale_fill_discrete(breaks = rev(levels(genera$top_genera))) + 
  theme_minimal() +
  theme(legend.direction = "horizontal", legend.position = "top",
        axis.title = element_blank(), axis.text.y = element_blank(), 
        text = element_text(family = "serif"),
        panel.grid = element_blank(), axis.ticks.x.bottom = element_line(linewidth = 0.75),
        aspect.ratio = 0.3) 
ggsave('Figures/Supplementary/FigureS2/Topgenera_sinbar.svg', scale = 0.75)
rm(genera)

# Panel C: Network coloured by PAE and Modern (Set to NULL; added separately)

# Panel D: Network coloured by mobility (Set to NULL; added separately)

# Panel E: Geographical distribution of plasmids from databases
country <- mt %>% 
  group_by(Country) %>% 
  summarise(count=n()) %>% 
  filter(Country != "NP" & Country != "Unknown")
# Treemap
E_Country <- ggplot(country, aes(area = count, fill = count, label = Country)) +
  geom_treemap() + theme(text = element_text(family = "serif")) +
  scale_fill_gradient(low = "darkcyan", high = "navyblue") +
  geom_treemap_text(family = "serif", fontface = "italic", colour = "white", 
                    place = "centre", grow = TRUE)
rm(country)

# Panel F: Virulence genes
library(readr)
vir_genes <- amr_genes <- read_tsv("Supplementary_Tables/TableS6-Vir_genes_detected.tsv")
vir_genes_count <- vir_genes %>% 
  group_by(GENE) %>% 
  summarise(count=n())
# Treemap
F_Virulence <- ggplot(vir_genes_count, aes(area = count, fill = count, label = GENE)) +
  geom_treemap() + theme(text = element_text(family = "serif")) +
  scale_fill_gradient(low = "darkcyan", high = "navyblue") +
  geom_treemap_text(family = "serif", fontface = "italic", colour = "white", 
                    place = "centre", grow = TRUE)
rm(vir_genes,vir_genes_count)

# Panel G: Size vs mobility per lineage type - facets
fnames <- as_labeller(c(
  `PAEModern_AMR` = "Mixed-CFs AMR+ (n=39)",
  `PAEModern_NoAMR` = "Mixed-CFs AMR- (n=49)",
  `Divergent_modern_AMR` = "Modern-CFs AMR+ (n=203)",
  `Divergent_modern_NoAMR` = "Modern-CFs AMR- (n=546)",
  `Divergent_PAE_NoAMR` = "PAE-CFs AMR- (n=103)"
))
G_SizeMob <- ggplot(CFs_df, aes(x=predicted_mobility, y=mean_size, fill = Category, colour = Category)) +
  labs(x="Mobility", y="Average plasmid size per CF (bp)") +
  theme_linedraw() + geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(method = "quasirandom", alpha = 0.5, dodge.width = 0.75) +
  scale_fill_manual(values = c("white", "white")) +
  scale_colour_manual(values = c("steelblue", "indianred")) + theme_linedraw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = c(0.89, 0.25), text = element_text(family = "serif")) + 
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames)

# Panel H: Size vs mobility per lineage type - facets
H_TransMob <- ggplot(CFs_df, aes(x=predicted_mobility, y=mean_Transposases, fill = Category, colour = Category)) +
  labs(x="Mobility", y="Average number of transposases per CF") +
  theme_linedraw() + geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(method = "quasirandom", alpha = 0.5, dodge.width = 0.75) +
  scale_fill_manual(values = c("white", "white")) +
  scale_colour_manual(values = c("steelblue", "indianred")) + theme_linedraw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = c(0.89, 0.25), text = element_text(family = "serif")) + 
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames)

# Merge plots and export
plot_grid(A_Year, NULL, NULL, NULL, E_Country, F_Virulence, G_SizeMob, H_TransMob,
          labels = c("A","B","C","D","E","F","G","H"), ncol = 2, 
          label_fontfamily = "serif", hjust = -1, label_size = 20)
ggsave('Figures/Supplementary/FigureS2/FigureS2.svg', scale = 0.5, 
       width = 25, height = 37)
rm(A_Year, E_Country, F_Virulence, G_SizeMob, H_TransMob)



### Figure S3
dir.create("Figures/Supplementary/FigureS3/", recursive = T)

# Data used
# CFs_pangenome: NN <- NN

# Panel A: CFs core genome
# Average proportion of core genes per plasmid per CF
left <- ggplot(CFs_pangenome, aes(x=CF_type, y=Core_genome_frc_per_plasmid, colour = CF_type)) +
  labs(x="CF type", y="Average proportion of core genes in plasmid per CF", colour="CF type") +
  theme_linedraw() + geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(aes(colour = CF_type), method = "pseudorandom", dodge.width = 0.75,
                   groupOnX=T) +
  scale_x_discrete(labels=c(
    PAEModern_AMR = "Mixed-CFs AMR+",
    PAEModern_NoAMR = "Mixed-CFs AMR-",
    Divergent_modern_AMR = "Modern-CFs AMR+",
    Divergent_modern_NoAMR = "Modern-CFs AMR-",
    Divergent_PAE_NoAMR = "PAE-CFs AMR-")) + 
  scale_colour_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none",
        text = element_text(family = "serif"))
# Proportion of core genes in pangenome per CF
right <- ggplot(CFs_pangenome, aes(x=CF_type, y=Core_genome_frc, colour = CF_type)) +
  labs(x="CF type", y="Proportion of core genes in pangenome per CF", colour="CF type") +
  theme_linedraw() + geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(aes(colour = CF_type), method = "pseudorandom", dodge.width = 0.75,
                   groupOnX=T) +
  scale_x_discrete(labels=c(
    PAEModern_AMR = "Mixed-CFs AMR+",
    PAEModern_NoAMR = "Mixed-CFs AMR-",
    Divergent_modern_AMR = "Modern-CFs AMR+",
    Divergent_modern_NoAMR = "Modern-CFs AMR-",
    Divergent_PAE_NoAMR = "PAE-CFs AMR-")) + 
  scale_colour_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4")) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none",
        text = element_text(family = "serif"))
# Merge plots
A_core <- plot_grid(left, right, labels = c("A",""), nrow = 1,
               label_fontfamily = "serif", hjust = -1, label_size = 20)
# Merge plots and export
plot_grid(A_core, NULL, labels = c("","B"), ncol = 1, 
          label_fontfamily = "serif", hjust = -1, label_size = 20)
ggsave('Figures/Supplementary/FigureS3/FigureS3.svg', scale = 0.5, 
       width = 25, height = 25)
rm(left, right, A_core)

# The pairwise comparisons shown in panel B figures were generated with clinker



### Figure S5
dir.create("Figures/Supplementary/FigureS5/", recursive = T)

# Data used
# DivModern_vs_pMUR_cat: NN <- NN

# Histogram of Modern-CF plasmids size
# Mean size label
meansize <- as.character(round(mean(DivModern_vs_pMUR_cat$Modern_size)))
# Plot
modernsizehist <- ggplot(DivModern_vs_pMUR_cat, aes(x=Modern_size, fill=Modern_CF_type)) + geom_histogram(bins = 100) +
  scale_x_continuous(limits = c(0,610000)) +
  labs(x="Modern-CF plasmids size (bp)", y="Count") + 
  scale_fill_manual(values = c("#FDAE61","#C2A5CF"),
                    name = "CF type",
                    labels = c("Modern-CF AMR+","Modern-CF AMR-")) +
  geom_vline(xintercept = mean(DivModern_vs_pMUR_cat$Modern_size), colour = "black") + 
  geom_text(x=mean((DivModern_vs_pMUR_cat$Modern_size)+50000), y=550, label=meansize, 
            stat = "unique", family = "serif", size=4, aes(fontface = 1)) +
  theme_linedraw() + theme(legend.position = "top", text = element_text(family = "serif"))

# Histogram of contained PAE plasmids size
# Mean size label
PAE_contained <- DivModern_vs_pMUR %>% 
  select(Murray, Murray_size, Murray_CF_type) %>% 
  distinct(.keep_all = T)
meansize <- as.character(round(mean(PAE_contained$Murray_size)))
# Plot
containedsizehist <- ggplot(PAE_contained, aes(x=Murray_size, fill=Murray_CF_type)) + geom_histogram(bins = 100) +
  scale_x_continuous(limits = c(0,610000)) +
  labs(x="Contained PAE plasmids size (bp)", y="Count") + 
  scale_fill_manual(values = c("#FF6619","#762A83","#41B6C4"),
                    name = "CF type",
                    labels = c("Mixed-CF AMR+","Mixed-CF AMR-","PAE-CF")) +
  geom_vline(xintercept = mean(PAE_contained$Murray_size), colour = "black") + 
  geom_text(x=mean((PAE_contained$Murray_size)+50000), y=110, label=meansize, 
            stat = "unique", family = "serif", size=4, aes(fontface = 1)) +
  theme_linedraw() + theme(legend.position = "top", text = element_text(family = "serif"))

# Panel A: Modern-CF and contained PAE plasmids size distribution
leftupper <- plot_grid(modernsizehist, containedsizehist, ncol = 1)
rm(PAE_contained, meansize, modernsizehist, containedsizehist)

# Panel B: Modern vs contained PAE plasmid size
rightupper <- ggplot(DivModern_vs_pMUR, aes(x=Murray_size, y=Modern_size)) + 
  geom_point(alpha=0.1, colour="darkcyan", size=1) + 
  labs(x="Contained PAE plasmid size (bp)", y="Modern-CF plasmid size (bp)") +
  scale_x_continuous(limits = c(0,610000)) +
  scale_y_continuous(limits = c(0,610000)) +
  geom_abline(slope = c(0.7,1,1.3), colour = c("grey70","black","grey70"),
              size = c(0.4,0.4,0.4), linetype = c(2,2,2)) + 
  theme_linedraw() + theme(legend.position = "top", text = element_text(family = "serif"))

# Panel C: Histogram of size ratio for Modern-CF vs contained PAE plasmids
# No outliers
# Mean size label
meansize <- as.character(round(mean(DivModern_vs_pMUR_cat$Size_ratio_mean)))
# Plot
leftbottom <- ggplot(DivModern_vs_pMUR_cat, aes(x=Size_ratio_mean, fill=Modern_CF_type)) + 
  geom_histogram(bins = 100) + theme_linedraw() +
  labs(x="Size ratio (Modern-CF vs contained PAE plasmids)", y="Count") +
  scale_fill_manual(values = c("#FDAE61","#C2A5CF"),
                    name = "CF type",
                    labels = c("Modern-CF AMR+","Modern-CF AMR-")) +
  scale_x_continuous(limits = c(0,66)) +
  geom_vline(xintercept = mean(DivModern_vs_pMUR_cat$Size_ratio_mean), colour = "black") +
  geom_text(x=mean((DivModern_vs_pMUR_cat$Size_ratio_mean)+5), y=470, label=meansize, 
            stat = "unique", family = "serif", size=4, aes(fontface = 1)) +
  theme_linedraw() + theme(legend.position = "none", text = element_text(family = "serif"))
# Including outliers
rightbottom <- ggplot(DivModern_vs_pMUR_cat, aes(x=Size_ratio_mean, fill=Modern_CF_type)) + 
  geom_histogram(bins = 100) + theme_linedraw() +
  labs(x="Size ratio (Modern-CF vs contained PAE plasmids)", y="Count") +
  scale_fill_manual(values = c("#FDAE61","#C2A5CF")) +
  geom_vline(xintercept = mean(DivModern_vs_pMUR_cat$Size_ratio_mean), colour = "black") +
  geom_text(x=mean((DivModern_vs_pMUR_cat$Size_ratio_mean)+25), y=1100, label=meansize, 
            stat = "unique", family = "serif", size=4, aes(fontface = 1)) +
  theme_linedraw() + theme(legend.position = "none", text = element_text(family = "serif"))

# Merge plots and export
upper <- plot_grid(leftupper, rightupper, ncol = 2, labels = "AUTO", label_fontfamily = "serif",
                   hjust = -1, label_size = 20, rel_widths = c(1,1.5))
bottom <- plot_grid(leftbottom, rightbottom, ncol = 2, label_fontfamily = "serif")
plot_grid(upper, bottom, labels = c("","C"), ncol = 1, label_fontfamily = "serif",
          hjust = -1, vjust = 0.5, label_size = 20, rel_heights = c(2,1))
ggsave('Figures/Supplementary/FigureS5/FigureS5.png', scale = 0.5, 
       width = 22, height = 22)
rm(leftupper, rightupper, upper, leftbottom, rightbottom, bottom)



### Figure S6
dir.create("Figures/Supplementary/FigureS6/", recursive = T)

# Data used
# DivModern_vs_pMUR_cat: NN <- NN

# Panel A: Average number of replicons per CF
A_RepAve <- ggplot(CFs_df, aes(x=predicted_mobility, y=mean_RepMT, fill = Category, colour = Category)) +
  labs(x="Mobility", y="Average number of replicons per CF") +
  geom_boxplot(outlier.shape = NA) + 
  geom_quasirandom(method = "quasirandom", alpha = 0.5, dodge.width = 0.75) +
  scale_fill_manual(values = c("white", "white")) +
  scale_colour_manual(values = c("steelblue", "indianred")) + theme_linedraw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = c(0.85, 0.25), text = element_text(family = "serif")) +
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames)

# Panel B: Number of replicons per CF
B_RepPerCF <- ggplot(mt, aes(x=reorder(CF, rep(-1,length(CF)), FUN=sum), y=Replicons_MobTyper, colour=CF_type)) +
  labs(x="CF", y="Number of replicons") +
  geom_hline(yintercept = 1, linetype =2) + theme_linedraw() + geom_boxplot() +
  scale_colour_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1),
        legend.position = "none", text = element_text(family = "serif")) +
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames, ncol = 1)

# Panel C: Proportion of single- vs multi-replicon plasmids per CF
# Extract data for the analysis
# Select required columns from the master data frame and add Replicon (no, single, multi) categories
replicons_df <- mt %>% 
  select(ID, Replicons_MobTyper, predicted_mobility, Category, CF, CF_type) %>% 
  mutate(Rep_category = 
           case_when(Replicons_MobTyper < 1 ~ "no-rep",
                     Replicons_MobTyper == 1 ~ "single-rep",
                     TRUE ~ "multi-rep"),
         Rep_category=factor(Rep_category, levels = c("no-rep", "single-rep","multi-rep")))
# Remove plasmids for which no replicons were detected - we don't need these for these plots
# Add a column to join outlier/singleton categories with the plasmid IDs - this creates unique IDs for plasmids in these lineage types
replicons_df <- replicons_df %>% 
  filter(Rep_category != "no-rep") %>% 
  mutate(CF_Ext = 
           case_when(CF == "Outlier" | CF == "Singleton" ~ paste(paste(CF, ID, sep = ".")),
                     TRUE ~ CF))
# Proportion of single- vs multi-replicon plasmids per mobility type
C_RepCat <- ggplot(replicons_df, aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(width = 0.55) + labs(x="Mobility", y="Count", fill="Rep category") + 
  theme_linedraw() + scale_fill_manual(values = c("darkgoldenrod","maroon")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1),
        legend.position = "right", text = element_text(family = "serif"))
  
# Panel D
# Proportion of single- vs multi-replicon plasmids per lineage - adjusted
D_RepCatCF <- ggplot(replicons_df, aes(x=CF, fill=Rep_category)) +
  geom_bar(position = "fill") + labs(y="Proportion of plasmids in CF") +
  theme_linedraw() + scale_fill_manual(values = c("darkgoldenrod","maroon")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1),
        legend.position = "none", text = element_text(family = "serif")) +
  facet_wrap( ~ CF_type, scales = "free", ncol = 1)

# Panel E: Number of distinct replicons identified per CF
# Make new data frame merging ourliers/singletons with ID
subclus <- mt %>% 
  filter(CF != "Outlier" & CF != "Singleton" ) %>% 
  select(CF, CF_type,`rep_type(s)`)
outsin <- mt %>%
  filter(CF == "Outlier" | CF == "Singleton" ) %>% 
  mutate(CF = as.character(CF),
         CF = paste(CF, ID, sep = "."),
         CF = as.factor(CF),) %>% 
  select(CF, CF_type,`rep_type(s)`)
linrep <- rbind(subclus,outsin)
rm(subclus,outsin)
# Extract total plasmids per replicon/CF and calculate frequency
CFsrepcount <- linrep %>% 
  separate_rows(`rep_type(s)`, sep = ",") %>% 
  group_by(CF, CF_type,`rep_type(s)`) %>% 
  summarise(rep_count = n(), .groups = "drop")
CFsrepcount <- CFsrepcount %>% 
  mutate(Plasmids = CFsExt_mt[match(CFsrepcount$CF, CFsExt_mt$CF),"Plas_total"],
         rep_freq = rep_count / Plasmids)
nam <- names(CFsrepcount)
CFsrepcount <- unnest(CFsrepcount)
colnames(CFsrepcount) <- nam
rm(nam, linrep)
# Number of different replicons (No ND & Outliers/Singletons) per CF
left <- CFsrepcount %>% 
  filter(!grepl("Outlier|Singleton", CF) & `rep_type(s)` != "ND") %>% 
  ggplot(aes(x=reorder(CF, rep(-1,length(CF)), FUN=sum), fill=CF_type)) +
  geom_bar() + labs(x="CF (No outliers/singletons, n=335)", y="Number of distinct replicons identified in CF") +
  scale_fill_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4"),
                    name = "CF type",
                    labels = c("Mixed-CF AMR+", "Mixed-CF AMR-", "Modern-CF AMR+",
                               "Modern-CF AMR-", "PAE-CF AMR-")) +
  scale_y_continuous(breaks = seq(0,18,4)) + theme_linedraw() +
  theme(axis.text.x=element_blank(), panel.grid.major.x = element_blank(), 
        text = element_text(family = "serif"), legend.position = c(0.75, 0.8),
        legend.background = element_blank())
# Number of different replicons (No ND) per CF
right <- CFsrepcount %>% 
  filter(`rep_type(s)` != "ND") %>% 
  ggplot(aes(x=reorder(CF, rep(-1,length(CF)), FUN=sum), fill=CF_type)) +
  geom_bar() + labs(x="CF (n=940)", y="Number of distinct replicons identified in CF") +
  scale_fill_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4")) +
  scale_y_continuous(breaks = seq(0,18,4)) + theme_linedraw() +
  theme(axis.text.x=element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "none", text = element_text(family = "serif"))
# Merge plots
E_DiffRepinCFs <-  plot_grid(left, right, ncol = 2, label_fontfamily = "serif",
                             hjust = -1, vjust = 0.5, label_size = 20)

# Panel F: Replicons occurrence per CF (No Outliers/Singletons)
left <- CFsrepcount %>% 
  filter(!grepl("Outlier|Singleton", CF)) %>% 
  ggplot(aes(x=reorder(`rep_type(s)`, rep(-1,length(`rep_type(s)`)), FUN=sum), fill=CF_type)) +
  geom_bar() + labs(x="Replicon (n=84)", y="Number of CFs with replicon (No outliers/singletons)") +
  scale_fill_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4"),
                    name = "CF type",
                    labels = c("Mixed-CF AMR+", "Mixed-CF AMR-", "Modern-CF AMR+",
                               "Modern-CF AMR-", "PAE-CF AMR-")) +
  scale_y_continuous(breaks = seq(0,150,20)) + theme_linedraw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1), panel.grid.major.x = element_blank(), 
        text = element_text(family = "serif"), legend.position = c(0.75, 0.8),
        legend.background = element_blank())
# Replicons occurrence per CF
right <- ggplot(CFsrepcount, aes(x=reorder(`rep_type(s)`, rep(-1,length(`rep_type(s)`)), FUN=sum), fill=CF_type)) +
  geom_bar() + labs(x="Replicon (n=84)", y="Number of CFs with replicon") +
  scale_fill_manual(values = c("#FF6619","#762A83","#FDAE61","#C2A5CF","#41B6C4")) +
  scale_y_continuous(breaks = seq(0,380,40)) + theme_linedraw() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1), panel.grid.major.x = element_blank(), 
        text = element_text(family = "serif"), legend.position = "none",
        legend.background = element_blank())
# Merge plots
F_RepOccinCFs <- plot_grid(left, right, ncol = 2, label_fontfamily = "serif",
                           hjust = -1, vjust = 0.5, label_size = 20)

# Merge all plots and export
top <- plot_grid(A_RepAve, B_RepPerCF, C_RepCat, D_RepCatCF, labels = "AUTO", 
                 ncol = 2, label_fontfamily = "serif", scale = c(.95,1,.75,1),
                 hjust = -1, vjust = 1, label_size = 20, rel_widths = c(1,1.5))
plot_grid(top, E_DiffRepinCFs, F_RepOccinCFs, labels = c("","E","F"), ncol = 1, 
          label_fontfamily = "serif", hjust = -1, vjust = 0.5, label_size = 20, 
          rel_heights = c(2.5,.8,1))
ggsave('Figures/Supplementary/FigureS6/FigureS6.png', scale = 0.5, 
       width = 30, height = 45)
rm(A_RepAve, B_RepPerCF, C_RepCat, D_RepCatCF, top, E_DiffRepinCFs, F_RepOccinCFs)



### Figure S8
dir.create("Figures/Supplementary/FigureS8/", recursive = T)

# Data used
# amr_genes_count: NN <- NN

# Panel A: Prevalence of AMR genes across the plasmids network (Set as NULL)

# Panel B: Diversity of AMR genes (AMR genes counts)
# Treemap
B_AMRgenesTreemap <- ggplot(amr_genes_count, aes(area = count, fill = count, label = GENE)) +
  geom_treemap() + theme(text = element_text(family = "serif")) +
  scale_fill_viridis_c(option = "rocket", direction = -1) +
  geom_treemap_text(family = "serif", fontface = "italic", colour = "black", 
                    place = "centre", grow = TRUE)

# Panel C: Average number of AMR genes VS R plasmids frequency per CF
p <- CFsExt_mt %>% 
  filter(AMR_genes_sum != 0) %>% 
  ggplot(aes(y=R_plasmidfreq, x=AMR_genes_mean, colour=CF_type)) +
  geom_point(alpha=0.85, aes(size=Plas_total)) + theme_linedraw() +
  labs(y="Frequency of R plasmids in CF", x="Average number of AMR genes in CF",
       size="Plasmids in CF") +
  scale_size_area(breaks = c(1,100,300,500,700), limits = c(1,732)) +
  scale_colour_manual(values = c("#FF6619","#FDAE61"),
                    name = "CF type",
                    labels = c("Mixed-CF AMR+", "Modern-CF AMR+")) +
  theme(text = element_text(family = "serif"), legend.position = "bottom",
        legend.box = "vertical")
left <- ggMarginal(p, type="boxplot", groupFill = T, alpha=0.85)
# No singletons/outliers
p <- CFsExt_mt %>% 
  filter(AMR_genes_sum != 0 & Plas_total > 1) %>% 
  ggplot(aes(y=R_plasmidfreq, x=AMR_genes_mean, colour=CF_type)) +
  geom_point(alpha=0.85, aes(size=Plas_total)) + theme_linedraw() +
  labs(y="Frequency of R plasmids in CF", x="Average number of AMR genes in CF",
       size="Plasmids in CF") +
  scale_size_area(breaks = c(2,100,300,500,700), limits = c(2,732)) +
  scale_colour_manual(values = c("#FF6619","#FDAE61"),
                      name = "CF type",
                      labels = c("Mixed-CF AMR+", "Modern-CF AMR+")) +
  theme(text = element_text(family = "serif"), legend.position = "bottom",
        legend.box = "vertical")
right <- ggMarginal(p, type="boxplot", groupFill = T, alpha=0.85, size=7)

# Merge all plots and export
top <- plot_grid(NULL, B_AMRgenesTreemap, labels = "AUTO", 
                 ncol = 2, label_fontfamily = "serif", scale = c(1,.85),
                 hjust = -1, vjust = 1, label_size = 20, rel_widths = c(1.25,1))
bottom <- plot_grid(left, right, ncol = 2, label_fontfamily = "serif",
                    align = "h")
plot_grid(top, bottom, labels = c("","C"), ncol = 1, 
          label_fontfamily = "serif", hjust = -1, vjust = 0.5, label_size = 20)
ggsave('Figures/Supplementary/FigureS8/FigureS8.svg', scale = 0.5, 
       width = 25, height = 25)
rm(B_AMRgenesTreemap, top, left, right, bottom)



### Figure S9
dir.create("Figures/Supplementary/FigureS9/", recursive = T)

# Data used
# mt: Plasmids master table <- Table S1

# Panel A: Frequency of mobility types per CF
A_MobFreq <- ggplot(mt, aes(x=reorder(CF, rep(-1,length(CF)), FUN=sum), fill=predicted_mobility)) +
  labs(x="CF", y="Plasmid count") +
  theme_linedraw() + geom_bar() +
  scale_fill_manual(values = c("#002147","#008B45","#00E5BF"),
                    breaks = c("conjugative","mobilizable","non-mobilizable"),
                    name = "Mobility") +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 90, hjust = 1),
        legend.position = "top") +
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames, ncol = 1)

# Panel B: Frequency of mobility types per CF - Adjusted to 1
B_MobFreqAdj <- ggplot(mt, aes(x=reorder(CF, rep(-1,length(CF)), FUN=sum), fill=predicted_mobility)) +
  labs(x="CF", y="Proportion") +
  theme_linedraw() + geom_bar(position = "fill") +
  scale_fill_manual(values = c("#002147","#008B45","#00E5BF"),
                    breaks = c("conjugative","mobilizable","non-mobilizable"),
                    name = "Mobility") +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 90, hjust = 1),
        legend.position = "top") +
  facet_wrap( ~ CF_type, scales = "free_x", labeller = fnames, ncol = 1)

# Merge plots and export
plot_grid(A_MobFreq, B_MobFreqAdj, labels = "AUTO", ncol = 1, 
          label_fontfamily = "serif", hjust = -1, vjust = 1, label_size = 20)
ggsave('Figures/Supplementary/FigureS9/FigureS9.png', scale = 0.5, 
       width = 25, height = 35)
rm(A_MobFreq, B_MobFreqAdj)



### Figure S10

# See code in "Size_comparison.R" file from the GitHub repository



### Figure S11
dir.create("Figures/Supplementary/FigureS11/", recursive = T)

# Data used
# iPlDB: NNN <- Table SN

# Panel A: Matches between plasmids from databases and Murray plasmids: pMUR coverage VS pMUR length
p <- iPlDB %>% 
  filter(pMUR_bestmatch_cov > 0) %>% 
  ggplot(aes(y=pMUR_bestmatch_len, x=pMUR_bestmatch_cov)) +
  geom_point(alpha=0.25, colour="steelblue") + labs(x="Murray plasmid coverage (percent)", y="Murray plasmid length (log10)") +
  scale_x_continuous(breaks = seq(0,100,20)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept = 70, linetype = "dashed", colour = "magenta") + theme_linedraw() +
  theme(text = element_text(family = "serif"))
A_pMURCovLen <- ggMarginal(p, type="histogram", margins = "both", alpha = 0.9,
                           xparams = list(bins=70, size=0.01), yparams = list(bins=70, size=0.01))

# Panel B: Matches between plasmids from databases and Murray plasmids: Size ratio VS pMUR length
p <- iPlDB %>% 
  filter(pMUR_bestmatch_cov > 0) %>% 
  ggplot(aes(y=Size_ratio, x=pMUR_bestmatch_cov)) +
  geom_point(alpha=0.25, colour="steelblue") + labs(x="pMUR coverage percent", y="Size ratio (log10)") +
  scale_x_continuous(breaks = seq(0,100,20)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  geom_vline(xintercept = 70, linetype = "dashed", colour = "magenta") + theme_linedraw() +
  theme(text = element_text(family = "serif"))
B_pMURCovSR <- ggMarginal(p, type="histogram", margins = "both", alpha = 0.9,
                          xparams = list(bins=70, size=0.01), yparams = list(bins=70, size=0.01))

# Panel C: Taxonomic composition (Class) of pMUR close relatives
taxclass <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  group_by(taxon_class_name) %>% 
  summarise(count=n())
C_TaxCloseRel <- ggplot(taxclass, aes(area = count, fill = count, label = taxon_class_name)) +
  geom_treemap() +
  scale_fill_distiller(palette = "Set3", direction = -1) +
  geom_treemap_text(family = "serif", fontface = "italic", colour = "grey30",
                    place = "centre", grow = TRUE)

# Panel D: Taxonomic composition (Class) of pMUR unrelated plasmids
taxclass <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  group_by(taxon_class_name) %>% 
  summarise(count=n())
D_TaxUnrel <- ggplot(taxclass, aes(area = count, fill = count, label = taxon_class_name)) +
  geom_treemap() +
  scale_fill_distiller(palette = "Set3", direction = -1) +
  geom_treemap_text(family = "serif", fontface = "italic", colour = "grey30", 
                    place = "centre", grow = TRUE)
rm(taxclass)

# Merge plots and export
plot_grid(A_pMURCovLen, B_pMURCovSR, C_TaxCloseRel, D_TaxUnrel, labels = "AUTO", 
          ncol = 2, label_fontfamily = "serif", hjust = -1, vjust = 1, 
          label_size = 20, scale = c(1,1,.9,.9), rel_heights = c(1,.8))
ggsave('Figures/Supplementary/FigureS11/FigureS11.png', scale = 0.5, 
       width = 25, height = 20)
rm(A_pMURCovLen, B_pMURCovSR, C_TaxCloseRel, D_TaxUnrel)



### Figure S12
dir.create("Figures/Supplementary/FigureS12/", recursive = T)

# Data used
# iPlDB: NNN <- Table SN

# Row 1: Murray-related vs unrelated distribution of R plasmids across mobility types (count)
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Mobility", y="Murray-related plasmids (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids: Gammaproteobacteria (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 2 - Enterococcus (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 3 - Staphylococcus (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# Merge plots
Row1 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Row 2: Murray-related vs unrelated distribution of R plasmids across mobility types (proportion)
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Mobility", y="Murray-related plasmids (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids: Gammaproteobacteria (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 2 - Enterococcus (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3) %>% 
  ggplot(aes(x=predicted_mobility, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 3 - Staphylococcus (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# Merge plots
Row2 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Row 3: Murray-related vs unrelated distribution of R plasmids across rep categories (count)
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Rep category", y="Murray-related plasmids (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids: Gammaproteobacteria (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids: Cluster 2 - Enterococcus (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids: Cluster 3 - Staphylococcus (count)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# Merge plots
Row3 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Row 4: Murray-related vs unrelated distribution of R plasmids across rep categories (proportion)
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Rep category", y="Murray-related plasmids (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids: Gammaproteobacteria (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids: Cluster 2 - Enterococcus (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3) %>% 
  ggplot(aes(x=Rep_category, fill=R_category)) +
  geom_bar(position = "fill", width = 0.55) + 
  labs(x="Rep category", y="Murray-unrelated plasmids: Cluster 3 - Staphylococcus (proportion)", fill="Resistance category") +
  scale_fill_manual(values = c("navyblue","darkred")) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# Merge plots
Row4 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Row 5: Murray-related vs unrelated distribution of AMR genes across rep categories
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70 & AMR_genes > 0) %>% 
  ggplot(aes(x=Rep_category, y=AMR_genes)) +
  labs(x="Rep category", y="AMR genes") +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(method = "pseudorandom", alpha = 0.25, dodge.width = 0.75,
                   colour="hotpink4", size = 0.75, groupOnX=T) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & AMR_genes > 0) %>% 
  ggplot(aes(x=Rep_category, y=AMR_genes)) +
  labs(x="Rep category", y="AMR genes") +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(method = "pseudorandom", alpha = 0.25, dodge.width = 0.75,
                   colour="hotpink4", size = 0.75, groupOnX=T) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & AMR_genes > 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=Rep_category, y=AMR_genes)) +
  labs(x="Rep category", y="AMR genes") +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(method = "pseudorandom", alpha = 0.25, dodge.width = 0.75,
                   colour="hotpink4", size = 0.75, groupOnX=T) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2 & AMR_genes > 0) %>% 
  ggplot(aes(x=Rep_category, y=AMR_genes)) +
  labs(x="Rep category", y="AMR genes") +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(method = "pseudorandom", alpha = 0.25, dodge.width = 0.75,
                   colour="hotpink4", size = 0.75, groupOnX=T) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3 & AMR_genes > 0) %>% 
  ggplot(aes(x=Rep_category, y=AMR_genes)) +
  labs(x="Rep category", y="AMR genes") +
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(method = "pseudorandom", alpha = 0.25, dodge.width = 0.75,
                   colour="hotpink4", size = 0.75, groupOnX=T) + theme_linedraw() +
  theme(text = element_text(family = "serif"),
        axis.text.x=element_text(angle = 45, hjust = 1),
        legend.position = "none")
# Merge plots
Row5 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Row 6: Murray-related vs unrelated plasmids distribution across rep categories & mobility types (count)
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-related plasmids (count)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids (count)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids: Gammaproteobacteria (count)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 2 - Enterococcus (count)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 3 - Staphylococcus (count)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# Merge plots
Row6 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Row 7: Murray-related vs unrelated plasmids distribution across rep categories & mobility types (proportion)
# A: Murray close relatives (>=70% pMUR cov)
A <- iPlDB %>% 
  filter(pMUR_bestmatch_cov >= 70) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(position = "fill", width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-related plasmids (proportion)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# B: Murray unrelated plasmids
B <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(position = "fill", width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids (proportion)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# C: Murray-unrelated (Gammaproteobacteria) plasmids
C <- iPlDB %>% 
  filter(pMUR_bestmatch_cov == 0 & taxon_class_name == "Gammaproteobacteria") %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(position = "fill", width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids: Gammaproteobacteria (proportion)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# D: Murray-unrelated (Cluster 2 - Enterococcus) plasmids
D <- iPlDB %>% 
  filter(`__ccCluster` == 2) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(position = "fill", width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 2 - Enterococcus (proportion)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# E: Murray-unrelated (Cluster 3 - Staphylococcus) plasmids
E <- iPlDB %>% 
  filter(`__ccCluster` == 3) %>% 
  ggplot(aes(x=predicted_mobility, fill=Rep_category)) +
  geom_bar(position = "fill", width = 0.55) + scale_fill_manual(values = c("grey","darkgoldenrod","maroon")) +
  labs(x="Mobility", y="Murray-unrelated plasmids: Cluster 3 - Staphylococcus (proportion)", fill="Rep category") +
  theme_linedraw() + theme(text = element_text(family = "serif"),
                           axis.text.x=element_text(angle = 45, hjust = 1),
                           legend.position = "none")
# Merge plots
Row7 <- plot_grid(A, B, C, D, E, nrow = 1, label_fontfamily = "serif")

# Merge all rows and export
plot_grid(Row1, Row2, Row3, Row4, Row5, Row6, Row7, 
          ncol = 1, label_fontfamily = "serif")
ggsave('Figures/Supplementary/FigureS12/FigureS12.svg', scale = 0.75, 
       width = 25, height = 32.5)
rm(A, B, C, D, E, Row1, Row2, Row3, Row4, Row5, Row6, Row7)