# 4 TB datasets - reanalysis

library(tidyverse)
library(phyloseq)
library(qiime2R)
library(ape)
library(gridExtra)
library(microbiome)
library(ggpubr)
library(pairwiseAdonis)
library(metagenomeSeq)
library(vegan)
library(ANCOMBC)
library(microeco)
library(file2meco)
suppressPackageStartupMessages(library(microViz))
library(fs)
library(magrittr)
library(scales)
library(Biostrings)


# Set seed
set.seed(380)

# Setting theme
theme_set(theme_bw())
pal = "Dark2"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

five_shades <- c("#519872", "#DAB707", "#533A71", "#6184D8", "#EF3054")
# five_shades_1 <- c("#06d6a0", "#ffc24b", "#118ab2", "#d08c60", "#f47068")
# five_shades_2 <- c("#519872", "#DAB707", "#533A71", "#EF3054", "#6184D8")
colors_30 <- c("#f47068", "#ffc24b", "#06d6a0", "#ffff00", "#8a2be2", "#00ff7f", "#ff8c00", "#00ff00", "#00bfff", "#0000ff",
               "#ff00ff", "#1e90ff", "#db7093", "#f0e68c", "#ff1493", "#ffa07a", "#ee82ee", "#696969",
               "#556b2f", "#483d8b", "#b22222", "#3cb371", "#000080", "#9acd32", "#8b008b", "#48d1cc",
               "#519872", "#DAB707", "#533A71", "#6184D8")
colors_30_mod <- c("#ffc24b", "#adb5bd", "#06d6a0", "#ffff00", "#8a2be2", "#00ff7f", "#ff8c00", "#00ff00", "#00bfff", "#0000ff",
                   "#ff00ff", "#1e90ff", "#db7093", "#f0e68c", "#ff1493", "#ffa07a", "#ee82ee", "#696969",
                   "#556b2f", "#483d8b", "#b22222", "#3cb371", "#000080", "#9acd32", "#8b008b", "#48d1cc",
                   "#519872", "#DAB707", "#533A71", "#6184D8")

custom_text <- theme(text=element_text(family = "Roboto", size = 15))

`%notin%` <- Negate(`%in%`)

# =================
# Add more columns to the metadata
metadata <- read.csv("metadata.tsv", header = T, sep = "\t")
metadata %<>% 
  group_by(Type) %>% 
  mutate(SampleName = paste(substr(Location, 1, 1), Type, row_number(), sep = "_")) %>% 
  ungroup() %>% 
  mutate(Ty_Loc = paste(substr(Location, 1, 1), Type, sep = "_"))

write_csv(metadata, "metadata_ext.tsv", quote = "none")

# Import qiime2 artifacts - using qiime2R
tb <- qza_to_phyloseq(features = "tb_table.qza", tree = "tb_rooted_tree.qza",
                      taxonomy = "tb_taxonomy.qza", metadata = "metadata_ext.tsv")

tb
# 7360 taxa
sample_variables(tb)
taxa_names(tb)[1:5]
range(taxa_sums(tb))
# 1 2116226

# Add representative sequences
rep_seqs <- readDNAStringSet("tb_rep_seqs/data/dna-sequences.fasta")
tb <- merge_phyloseq(tb, rep_seqs)
tb
# 7360 taxa

head(refseq(tb))
range(width(refseq(tb)))
# plot(width(refseq(tb)))
hist(width(refseq(tb)), labels = T, breaks = c(20, 130, 200, 300, 450), freq = T)
# There are 53 sequences which are <130 bases
short_seqs <- names(refseq(tb)[width(refseq(tb)) < 130])
# Number of occurrences of these short taxa
taxa_sums(tb)[short_seqs]
range(taxa_sums(tb)[short_seqs])
hist(taxa_sums(tb)[short_seqs], labels = T, breaks = c(2, 50, 1000, 119701), freq = T)

# Making a df
short_seqs_df <- data.frame(short_seqs, width(refseq(tb)[short_seqs]), taxa_sums(tb)[short_seqs])
# Test whether the df has taken the correct values
width(refseq(tb)["633c222e60993fe9d6c526590c9384d7"])
taxa_sums(tb)["633c222e60993fe9d6c526590c9384d7"]

# Shortest sequence has highest taxa_sums
# Check what it is actually classified into?
# ac35b2a1e68f71b43269b818d9f7f4b5 - 20 bp - 11051 freq - classified up to c__Bacteroidia
# a33157b8d8960c079576f0e4282b297b - 21 bp - 119707 freq - classified up to d_Bacteria
# After that the sequences up to length 60 were either d__Bacteria or Unassigned

# Filtering the object to exclude sequences with less than 60 bases length
tb <- prune_taxa(width(refseq(tb)) >= 60, tb)
tb
# 7353 taxa
range(taxa_sums(tb))
# 1 2116226
range(sample_sums(tb))
# 12603 760187
# The min sample_sum did not change even after filtering the sequences based on seq length

View(tax_table(tb))
all_seqs <- names(refseq(tb))
all_seqs_df <- data.frame(all_seqs, width(refseq(tb)[all_seqs]), taxa_sums(tb)[all_seqs])
# tax_table(tb)[, "Kingdom"] == "Unassigned"
# Filter unassigned, mitochondria, chloroplast
tb_noUn <- subset_taxa(tb, !is.na(Phylum) & 
                         Phylum %notin% c("", "uncharacterized") & 
                         Genus %notin% c("Chloroplast", "Mitochondria"))
tb_noUn
# 6311 taxa
View(tax_table(tb_noUn))
range(taxa_sums(tb_noUn))
# 1 2116226
range(sample_sums(tb_noUn))
# 12430 759613
# Rename uncultured classifications & clean tax_table
tb_clTaxa <- tb_noUn %>% 
  tax_fix(
    min_length = 4,
    unknowns = c("", "uncultured"),
    sep = " ", anon_unique = T,
    suffix_rank = "current"
  )
View(tax_table(tb_clTaxa))

# Rarefy - check the feature count, this might have changed after performing the above filtering steps
tb_clTaxa
range(taxa_sums(tb_clTaxa))
# 1 2116226
range(sample_sums(tb_clTaxa))
# 12430 759613

# Going with tax_glom followed by rarefaction - Check the RAREFACTION TEST section at the end of this document
tb_gen <- tax_glom(tb_clTaxa, "Genus")
tb_gen
# 402 taxa
range(taxa_sums(tb_gen))
# 3 9920692
range(sample_sums(tb_gen))
# 12430 759613

tb_final <- rarefy_even_depth(tb_gen, sample.size = min(sample_sums(tb_gen)), rngseed = T)
# 20 OTUs removed
tb_final
# 382 taxa
range(taxa_sums(tb_final))
# 1 608828
range(sample_sums(tb_final))
# 12430 12430

# -----------------
# Compute prevalence of each feature, store as data.frame
# Adapted from DADA2 R workflow paper - https://f1000research.com/articles/5-1492
prev_df <- apply(X = otu_table(tb_final),
                 MARGIN = ifelse(taxa_are_rows(tb_final), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prev_df <- data.frame(Prevalence = prev_df,
                      TotalAbundance = taxa_sums(tb_final),
                      tax_table(tb_final))

# Prevalence sum and Mean
prev_mean <- plyr::ddply(prev_df, "Phylum", function(df1){cbind(mean(df1$Prevalence), sum(df1$Prevalence))})
colnames(prev_mean) <- c("Phylum", "Mean", "Sum")
View(prev_mean)
# 16 phyla

# Prevalence vs Count plots (yintercept represents 10% of the samples)
p_v_c <- ggplot(prev_df, aes(TotalAbundance, Prevalence / nsamples(tb_final), color=Phylum)) +
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +
  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position = "none") + scale_color_discrete()
p_v_c

# -----------------
# Alpha diversity
alpha_2 <- c("Observed", "Shannon")
sample_variables(tb_final)
# Amplicon (3 observations)
View(unique(sample_data(tb_final)[, "Amplicon"]))
# V3, V3-V5, V4
amp_shades <- c("#533A71", "#6184D8", "#EF3054")
amp_comp <- list(c("V3", "V3-V5"), c("V3", "V4"), c("V3-V5", "V4"))
alpha_amp <- plot_richness(tb_final, x = "Amplicon", measures = alpha_2) +
  geom_jitter(aes(color = Amplicon), show.legend = F, width = 0.2, alpha = 0.8)
alpha_amp$layers <- alpha_amp$layers[-1]
alpha_amp + geom_boxplot(aes(fill = Amplicon, color = Amplicon), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = amp_shades) +
  scale_color_manual(values = amp_shades, guide = "none") + custom_text +
  stat_compare_means(comparisons = amp_comp, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1) +
  labs(x = "")


# Location (4 observations)
View(unique(sample_data(tb_final)[, "Location"]))
# India, S_Africa, Taiwan, W_Africa
loc_shades <- c("#519872", "#533A71", "#DAB707", "#6184D8")
loc_comp <- list(c("India", "S_Africa"), c("India", "Taiwan"), c("India", "W_Africa"),
                 c("S_Africa", "Taiwan"), c("S_Africa", "W_Africa"), c("Taiwan", "W_Africa"))
alpha_loc <- plot_richness(tb_final, x = "Location", measures = alpha_2) +
  geom_jitter(aes(color = Location), show.legend = F, width = 0.2, alpha = 0.8)
alpha_loc$layers <- alpha_loc$layers[-1]
alpha_loc + geom_boxplot(aes(fill = Location, color = Location), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = loc_shades) +
  scale_color_manual(values = loc_shades, guide = "none") + custom_text +
  stat_compare_means(comparisons = loc_comp, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1) +
  labs(x = "")


# Type (7 observations)
View(unique(sample_data(tb_final)[, "Type"]))
# HHC: healthy households contacts; HC: Healthy control; CC: Close contacts of PTB patients; TB: PTB patients
# before initiation of anti tubercular therapy; TBW: PTB patients at one week of anti tubercular therapy; TBM:
# PTB patients at one month of anti tubercular therapy; TB2M: PTB patients at two months of anti tubercular therapy
# CC, HC, HHC, TB, TB2M, TBM, TBW
t_shades <- c("#533A71", "#519872", "#6184D8", "#b22222", "#DAB707", "#ff8c00", "#EF3054")
t_comp <- list(c("CC", "HC"), c("CC", "HHC"), c("CC", "TB"), c("CC", "TB2M"), c("CC", "TBM"), c("CC", "TBW"),
               c("HC", "HHC"), c("HC", "TB"), c("HC", "TB2M"), c("HC", "TBM"), c("HC", "TBW"),
               c("HHC", "TB"), c("HHC", "TB2M"), c("HHC", "TBM"), c("HHC", "TBW"),
               c("TB", "TB2M"), c("TB", "TBM"), c("TB", "TBW"),
               c("TB2M", "TBM"), c("TB2M", "TBW"), 
               c("TBM", "TBW"))
alpha_typ <- plot_richness(tb_final, x = "Type", measures = alpha_2) +
  geom_jitter(aes(color = Type), show.legend = F, width = 0.2, alpha = 0.8)
alpha_typ$layers <- alpha_typ$layers[-1]
alpha_typ + geom_boxplot(aes(fill = Type, color = Type), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = t_shades) +
  scale_color_manual(values = t_shades, guide = "none") + custom_text +
  stat_compare_means(comparisons = t_comp, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1) +
  labs(x = "")


# Ty_Loc (11 observations)
View(unique(sample_data(tb_final)[, "Ty_Loc"]))
# I_HHC, I_TB, I_TBM, I_TBW, S_CC, S_TB, T_HC, T_TB, W_HC, W_TB, W_TB2M
tl_shades <- c("#6184D8", "#b22222", "#ff8c00", "#EF3054", 
               "#533A71", "#b22222", 
               "#519872", "#b22222", 
               "#519872", "#b22222", "#DAB707")
tl_comp <- list(c("I_HHC", "I_TB"), c("I_HHC", "I_TBW"), c("I_HHC", "I_TBM"), c("I_HHC", "T_HC"), c("I_HHC", "T_TB"), 
                c("I_HHC", "W_TB"), c("I_HHC", "W_TB2M"), c("I_HHC", "W_HC"), c("I_HHC", "S_CC"), c("I_HHC", "S_TB"),
                c("I_TB", "I_TBW"), c("I_TB", "I_TBM"), c("I_TB", "T_HC"), c("I_TB", "T_TB"), c("I_TB", "W_TB"), 
                c("I_TB", "W_TB2M"), c("I_TB", "W_HC"), c("I_TB", "S_CC"), c("I_TB", "S_TB"),
                c("I_TBW", "I_TBM"), c("I_TBW", "T_HC"), c("I_TBW", "T_TB"), c("I_TBW", "W_TB"), c("I_TBW", "W_TB2M"), 
                c("I_TBW", "W_HC"), c("I_TBW", "S_CC"), c("I_TBW", "S_TB"),
                c("I_TBM", "T_HC"), c("I_TBM", "T_TB"), c("I_TBM", "W_TB"), c("I_TBM", "W_TB2M"), c("I_TBM", "W_HC"), 
                c("I_TBM", "S_CC"), c("I_TBM", "S_TB"),
                c("T_HC", "T_TB"), c("T_HC", "W_TB"), c("T_HC", "W_TB2M"), c("T_HC", "W_HC"), c("T_HC", "S_CC"), 
                c("T_HC", "S_TB"),
                c("T_TB", "W_TB"), c("T_TB", "W_TB2M"), c("T_TB", "W_HC"), c("T_TB", "S_CC"), c("T_TB", "S_TB"),
                c("W_TB", "W_TB2M"), c("W_TB", "W_HC"), c("W_TB", "S_CC"), c("W_TB", "S_TB"),
                c("W_TB2M", "W_HC"), c("W_TB2M", "S_CC"), c("W_TB2M", "S_TB"),
                c("W_HC", "S_CC"), c("W_HC", "S_TB"),
                c("S_CC", "S_TB"))
alpha_tl <- plot_richness(tb_final, x = "Ty_Loc", measures = alpha_2) +
  geom_jitter(aes(color = Ty_Loc), show.legend = F, width = 0.2, alpha = 0.8)
alpha_tl$layers <- alpha_tl$layers[-1]
alpha_tl + geom_boxplot(aes(fill = Ty_Loc, color = Ty_Loc), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = tl_shades) +
  scale_color_manual(values = tl_shades, guide = "none") + custom_text +
  # stat_compare_means(comparisons = tl_comp, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1) +
  labs(x = "")

# Custom comparisons - Within Locations
tl_comp_mod <- list(c("I_HHC", "I_TB"), c("I_HHC", "I_TBW"), c("I_HHC", "I_TBM"),
                    c("I_TB", "I_TBW"), c("I_TB", "I_TBM"), 
                    c("I_TBW", "I_TBM"),
                    c("S_CC", "S_TB"),
                    c("T_HC", "T_TB"), 
                    c("W_TB", "W_TB2M"), c("W_TB", "W_HC"), 
                    c("W_TB2M", "W_HC"))

tiff("plots_new2/alpha_type.tiff", units="in", width=9.5, height=7.5, res=300)
alpha_tl + geom_boxplot(aes(fill = Ty_Loc, color = Ty_Loc), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = tl_shades, guide = "none") +
  scale_color_manual(values = tl_shades, guide = "none") + custom_text +
  stat_compare_means(comparisons = tl_comp_mod, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1) +
  labs(x = "")
dev.off()

# Custom comparisons - Actual variations - May not be required

# -----------------
# Beta diversity
# Weighted UniFrac
ord_uf <- ordinate(tb_final, "PCoA", "unifrac", weighted = T)
# Amplicon
pcoa_uf_amp <- plot_ordination(tb_final, ord_uf, color = "Amplicon")
pcoa_uf_amp$layers[[1]]$aes_params$size = 0.3
pcoa_uf_amp$layers[[1]]$aes_params$alpha = 0.5
pcoa_uf_amp + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = amp_shades) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Amplicon"))

# Location
pcoa_uf_loc <- plot_ordination(tb_final, ord_uf, color = "Location")
pcoa_uf_loc$layers[[1]]$aes_params$size = 0.3
pcoa_uf_loc$layers[[1]]$aes_params$alpha = 0.5
pcoa_uf_loc + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = loc_shades) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Location"))

# Type
shades_11 <- c("#533A71", "#519872", "#6184D8", "#b22222", "#DAB707", "#ff8c00", "#EF3054",
               "#519872", "#533A71", "#DAB707", "#6184D8")

shades_4 <- c("#519872", "#533A71", "#DAB707", "#6184D8")

tax_vars <- sort(unique(sample_data(tb_final)[,"Type"]$Type))

pcoa_uf_typ <- plot_ordination(tb_final, ord_uf, color = "Type", shape = "Location")
pcoa_uf_typ$layers[[1]]$aes_params$size = 0.1
pcoa_uf_typ$layers[[1]]$aes_params$alpha = 0.5
tiff("plots_new2/beta_wUF.tiff", units="in", width=8.5, height=7.5, res=300)
pcoa_uf_typ + geom_point(size = 3, alpha = 0.4) +
  scale_color_manual(values = shades_11, breaks = tax_vars) +
  scale_shape_manual(values = c(15, 19, 17, 1)) +
  stat_ellipse(aes(colour = Location, fill = Location), 
               geom = "polygon", type = "t", linewidth = 0, alpha = 0.10) +
  scale_fill_manual(values = shades_4) +
  custom_text
dev.off()


# Ty_Loc - Not required
pcoa_uf_tl <- plot_ordination(tb_final, ord_uf, color = "Ty_Loc")
pcoa_uf_tl$layers[[1]]$aes_params$size = 0.3
pcoa_uf_tl$layers[[1]]$aes_params$alpha = 0.5
pcoa_uf_tl + geom_point(size = 3, alpha = 0.4) +
  scale_color_manual(values = tl_shades) +
  # stat_ellipse(type = "t", linetype = 2, alpha = 0.5) +
  custom_text +
  guides(color=guide_legend(title="Type & Location"))

# PERMANOVA
dist_wUF <- UniFrac(tb_final, weighted = T, normalized = T, parallel = T)
perm_oa <- adonis2(dist_wUF ~ sample_data(tb_final)$Ty_Loc, permutations = 10000)
perm_pw <- pairwise.adonis(dist_wUF, sample_data(tb_final)$Ty_Loc, perm = 10000)

write.csv(perm_oa, "plots_new/beta_overall_tl.csv", quote = F)
write.csv(perm_pw, "plots_new/beta_pairwise_tl.csv", quote = F)


# -----------------
taxa_names(tb_final)[1:5]
# Change the taxa_names
table(duplicated(tax_table(tb_final)[, "Genus"]))
# 1 duplicate found
View(tax_table(tb_final))
tax_table(tb_final)[duplicated(tax_table(tb_final)[, "Genus"])]
# Rename the second instance
tax_table(tb_final)[duplicated(tax_table(tb_final)[, "Genus"])][,"Genus"] <- "Peptostreptococcales-Tissierellales Genus2"
table(duplicated(tax_table(tb_final)[, "Genus"]))
# No duplicates
taxa_names(tb_final) <- tax_table(tb_final)[, "Genus"]
taxa_names(tb_final)[1:5]
phy_tree(tb_final)
refseq(tb_final)

# Phylum level bar plot
tb_phyla <- tax_glom(tb_final, taxrank = "Phylum")
tb_phyla
# 16 phyla
taxa_names(tb_phyla) <- tax_table(tb_phyla)[, "Phylum"]
taxa_names(tb_phyla)
tb_phyla_rel <- transform(tb_phyla, "compositional")
# Per sample bar plot - not useful
bar_phy <- plot_composition(tb_phyla_rel, sample.sort = "Location", otu.sort = "abundance", x.label = "SampleName") +
  guides(fill = guide_legend(title = "Phylum")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative Abundance (%)")
bar_phy + scale_fill_manual(values = colors_30)
# Type+Location merged plot
bar_phy_tl <- plot_composition(tb_phyla_rel, average_by = "Ty_Loc", 
                               transform = "compositional", otu.sort = "abundance") +
  guides(fill = guide_legend(title = "Phylum")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative Abundance (%)")
tiff("plots_new2/bar_phyla.tiff", units="in", width=7, height=5.5, res=300)
bar_phy_tl + scale_fill_manual(values = colors_30) + custom_text
dev.off()

# Genus level bar plot
tb_rel <- transform(tb_final, "compositional")
# Aggregating at 50% prevalence and 0.01% detection threshold
tb_rel_agg <- aggregate_rare(tb_rel, level = "Genus", detection = 1/1000, prevalence = 0.5)
tb_rel_agg
# 25 taxa
bar_gen_tl <- plot_composition(tb_rel_agg, average_by = "Ty_Loc",
                               transform = "compositional", otu.sort = "abundance") +
  guides(fill = guide_legend(title = "Genus")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Relative Abundance (%)")
tiff("plots_new2/bar_genera.tiff", units="in", width=10.5, height=5.5, res=300)
bar_gen_tl + scale_fill_manual(values = colors_30_mod) + custom_text
dev.off()

# -----------------
# Subset the data
tb_ind <- subset_samples(tb_final, Location == "India")
tb_ind <- prune_taxa(taxa_sums(tb_ind) > 0, tb_ind)
tb_ind
# 156 taxa
tb_ind_agg <- aggregate_rare(tb_ind, level = "Genus", detection = 1/1000, prevalence = 0.75)
tb_ind_agg
# 22 taxa
# Remove "Other" category - Keep only dominant genera
tb_ind_dom <- subset_taxa(tb_ind_agg, Genus != "Other")
tb_ind_dom
# 21 taxa

tb_saf <- subset_samples(tb_final, Location == "S_Africa")
tb_saf <- prune_taxa(taxa_sums(tb_saf) > 0, tb_saf)
tb_saf
# 294 taxa
tb_saf_agg <- aggregate_rare(tb_saf, level = "Genus", detection = 1/1000, prevalence = 0.75)
tb_saf_agg
# 27 taxa
# Remove "Other" category - Keep only dominant genera
tb_saf_dom <- subset_taxa(tb_saf_agg, Genus != "Other")
tb_saf_dom
# 26 taxa

tb_tai <- subset_samples(tb_final, Location == "Taiwan")
tb_tai <- prune_taxa(taxa_sums(tb_tai) > 0, tb_tai)
tb_tai
# 253 taxa
tb_tai_agg <- aggregate_rare(tb_tai, level = "Genus", detection = 1/1000, prevalence = 0.75)
tb_tai_agg
# 13 taxa
# Remove "Other" category - Keep only dominant genera
tb_tai_dom <- subset_taxa(tb_tai_agg, Genus != "Other")
tb_tai_dom
# 12 taxa

tb_waf <- subset_samples(tb_final, Location == "W_Africa")
tb_waf <- prune_taxa(taxa_sums(tb_waf) > 0, tb_waf)
tb_waf
# 227 taxa
tb_waf_agg <- aggregate_rare(tb_waf, level = "Genus", detection = 1/1000, prevalence = 0.75)
tb_waf_agg
# 16 taxa
# Remove "Other" category - Keep only dominant genera
tb_waf_dom <- subset_taxa(tb_waf_agg, Genus != "Other")
tb_waf_dom
# 15 taxa

# -----------------
# heatmaps for dominant genera
hm_ind <- plot_heatmap(tb_ind_dom, "NMDS", "bray", taxa.label = "Genus", sample.order = "Type", 
                       sample.label = "SampleName",
                       low = "#fee08b", high = "#d53e4f", na.value = "grey95")
tiff("plots_new/tiff_plots/hm_ind.tiff", units="in", width=13, height=7, res=300)
hm_ind + geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(.~Type, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) +
  labs(x = "")
dev.off()

hm_saf <- plot_heatmap(tb_saf_dom, "NMDS", "bray", taxa.label = "Genus", sample.order = "Type", 
                       sample.label = "SampleName",
                       low = "#fee08b", high = "#d53e4f", na.value = "grey95")
tiff("plots_new/tiff_plots/hm_saf.tiff", units="in", width=23, height=8, res=300)
hm_saf + geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(.~Type, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) +
  labs(x = "")
dev.off()

hm_tai <- plot_heatmap(tb_tai_dom, "NMDS", "bray", taxa.label = "Genus", sample.order = "Type", 
                       sample.label = "SampleName",
                       low = "#fee08b", high = "#d53e4f", na.value = "grey95")
tiff("plots_new/tiff_plots/hm_tai.tiff", units="in", width=18, height=5, res=300)
hm_tai + geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(.~Type, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) +
  labs(x = "")
dev.off()

hm_waf <- plot_heatmap(tb_waf_dom, "NMDS", "bray", taxa.label = "Genus", sample.order = "Type", 
                       sample.label = "SampleName",
                       low = "#fee08b", high = "#d53e4f", na.value = "grey95")
tiff("plots_new/tiff_plots/hm_waf.tiff", units="in", width=18, height=6, res=300)
hm_waf + geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(.~Type, scales = "free_x", space = "free_x") +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) +
  labs(x = "")
dev.off()

# -----------------
# Convert to microeco tables
tb_ind_me <- phyloseq2meco(tb_ind)
tb_saf_me <- phyloseq2meco(tb_saf)
tb_tai_me <- phyloseq2meco(tb_tai)
tb_waf_me <- phyloseq2meco(tb_waf)

# LEfSe analysis
me_ind_type <- trans_diff$new(dataset = tb_ind_me, method = 'lefse', group = "Type", taxa_level = "Genus")
# No discriminant taxa in the India samples. Setting the p_adjust_method to none (Not recommended)
me_ind_type <- trans_diff$new(dataset = tb_ind_me, method = 'lefse', group = "Type", taxa_level = "Genus", 
                              p_adjust_method = "none")
View(me_ind_type$res_diff)
ind_lefse <- me_ind_type$plot_diff_bar()
ind_bar <- me_ind_type$plot_diff_abund(select_taxa = me_ind_type$plot_diff_bar_taxa,
                                       plot_type = "barerrorbar", errorbar_addpoint = FALSE, errorbar_color_black = TRUE)
ind_lefse <- ind_lefse + theme(legend.position = "none")
ind_bar <- ind_bar + theme(axis.ticks.y = element_blank(), panel.border = element_blank())
tiff("plots_new2/lefse_ind.tiff", units="in", width=13, height=3, res=300)
ind_lefse %>% aplot::insert_right(ind_bar)
dev.off()
write.csv(me_ind_type$res_diff, "lefse_ind.csv")

me_saf_type <- trans_diff$new(dataset = tb_saf_me, method = 'lefse', group = "Type", taxa_level = "Genus")
# No discriminant taxa in the S_Africa samples. Setting the p_adjust_method to none (Not recommended)
me_saf_type <- trans_diff$new(dataset = tb_saf_me, method = 'lefse', group = "Type", taxa_level = "Genus",
                              p_adjust_method = "none")
View(me_saf_type$res_diff)
saf_lefse <- me_saf_type$plot_diff_bar(threshold = 2)
saf_bar <- me_saf_type$plot_diff_abund(select_taxa = me_saf_type$plot_diff_bar_taxa,
                                       plot_type = "barerrorbar", errorbar_addpoint = FALSE, errorbar_color_black = TRUE)
saf_lefse <- saf_lefse + theme(legend.position = "none")
saf_bar <- saf_bar + theme(axis.ticks.y = element_blank(), panel.border = element_blank())
tiff("plots_new2/lefse_saf.tiff", units="in", width=13, height=7, res=300)
saf_lefse %>% aplot::insert_right(saf_bar)
dev.off()
write.csv(me_saf_type$res_diff, "lefse_saf.csv")

me_tai_type <- trans_diff$new(dataset = tb_tai_me, method = 'lefse', group = "Type", taxa_level = "Genus")
# No discriminant taxa in the Taiwan samples. Setting the p_adjust_method to none (Not recommended)
me_tai_type <- trans_diff$new(dataset = tb_tai_me, method = 'lefse', group = "Type", taxa_level = "Genus",
                              p_adjust_method = "none")
View(me_tai_type$res_diff)
tai_lefse <- me_tai_type$plot_diff_bar(threshold = 2)
tai_bar <- me_tai_type$plot_diff_abund(select_taxa = me_tai_type$plot_diff_bar_taxa,
                                       plot_type = "barerrorbar", errorbar_addpoint = FALSE, errorbar_color_black = TRUE)
tai_lefse <- tai_lefse + theme(legend.position = "none")
tai_bar <- tai_bar + theme(axis.ticks.y = element_blank(), panel.border = element_blank())
tiff("plots_new2/lefse_tai.tiff", units="in", width=13, height=7, res=300)
tai_lefse %>% aplot::insert_right(tai_bar)
dev.off()
write.csv(me_tai_type$res_diff, "lefse_tai.csv")

me_waf_type <- trans_diff$new(dataset = tb_waf_me, method = 'lefse', group = "Type", taxa_level = "Genus")
# 3 discriminant taxa found in W_Africa samples. Setting the p_adjust_method to none (Not recommended)
me_waf_type <- trans_diff$new(dataset = tb_waf_me, method = 'lefse', group = "Type", taxa_level = "Genus",
                              p_adjust_method = "none")
View(me_waf_type$res_diff)
waf_lefse <- me_waf_type$plot_diff_bar(threshold = 2.6)
waf_bar <- me_waf_type$plot_diff_abund(select_taxa = me_waf_type$plot_diff_bar_taxa,
                                       plot_type = "barerrorbar", errorbar_addpoint = FALSE, errorbar_color_black = TRUE)
waf_lefse <- waf_lefse + theme(legend.position = "none")
waf_bar <- waf_bar + theme(axis.ticks.y = element_blank(), panel.border = element_blank())
tiff("plots_new2/lefse_waf.tiff", units="in", width=13, height=12, res=300)
waf_lefse %>% aplot::insert_right(waf_bar)
dev.off()
write.csv(me_waf_type$res_diff, "lefse_waf_unadj.csv")

# -----------------
# ANCOM-BC
# It requires raw counts instead of the rarefied table
# tb_gen is the agglomerated raw count table
# Make subsets
tb_ind_raw <- subset_samples(tb_gen, Location == "India")
tb_ind_raw <- prune_taxa(taxa_sums(tb_ind_raw) > 0, tb_ind_raw)
tb_ind_raw
# 188 taxa, 24 samples
tb_saf_raw <- subset_samples(tb_gen, Location == "S_Africa")
tb_saf_raw <- prune_taxa(taxa_sums(tb_saf_raw) > 0, tb_saf_raw)
tb_saf_raw
# 302 taxa, 88 samples
tb_tai_raw <- subset_samples(tb_gen, Location == "Taiwan")
tb_tai_raw <- prune_taxa(taxa_sums(tb_tai_raw) > 0, tb_tai_raw)
tb_tai_raw
# 263 taxa, 51 samples
tb_waf_raw <- subset_samples(tb_gen, Location == "W_Africa")
tb_waf_raw <- prune_taxa(taxa_sums(tb_waf_raw) > 0, tb_waf_raw)
tb_waf_raw
# 229 taxa, 50 samples

# Convert to microeco tables
me_ind_raw <- phyloseq2meco(tb_ind_raw)
me_saf_raw <- phyloseq2meco(tb_saf_raw)
me_tai_raw <- phyloseq2meco(tb_tai_raw)
me_waf_raw <- phyloseq2meco(tb_waf_raw)

# ANCOM objects
me_ind_ancom <- trans_diff$new(dataset = me_ind_raw, method = 'ancombc2', group = "Type", taxa_level = "Genus")
View(me_ind_ancom$res_diff)
me_ind_ancom$plot_diff_bar(keep_full_name = T, heatmap_cell = "P.adj", heatmap_sig = "Significance", heatmap_x = "Factors", heatmap_y = "Taxa")

me_saf_ancom <- trans_diff$new(dataset = me_saf_raw, method = 'ancombc2', group = "Type", taxa_level = "Genus")
View(me_saf_ancom$res_diff)

me_tai_ancom <- trans_diff$new(dataset = me_tai_raw, method = 'ancombc2', group = "Type", taxa_level = "Genus")
View(me_tai_ancom$res_diff)

me_waf_ancom <- trans_diff$new(dataset = me_waf_raw, method = 'ancombc2', group = "Type", taxa_level = "Genus")
View(me_waf_ancom$res_diff)


# =========== RAREFACTION TEST START ============
# Whether rarefying at this step is ok or can be done even after agglomeration at Genus level
# The sample sums should not change even after agglomeration
# The random sampling of reads may be done at Genus level rather than ASV level?!!!

# Testing by running tax_glom first, followed by rarefaction
test_ps <- tax_glom(tb_clTaxa, "Genus")
test_ps
# 402 taxa
range(taxa_sums(test_ps))
# 3 9920692
range(sample_sums(test_ps))
# 12430 759613
test_ps_rarefied <- rarefy_even_depth(test_ps, sample.size = min(sample_sums(test_ps)), rngseed = T)
test_ps_rarefied
# 382 taxa
range(taxa_sums(test_ps_rarefied))
# 1 608828
range(sample_sums(test_ps_rarefied))
# 12430 12430

# Testing by running rarefaction first, followed by tax_glom
test_ps2_rarefied <- rarefy_even_depth(tb_clTaxa, sample.size = min(sample_sums(tb_clTaxa)), rngseed = T) 
test_ps2_rarefied
# 5921 taxa
range(taxa_sums(test_ps2_rarefied))
# 1 158069
range(sample_sums(test_ps2_rarefied))
# 12430 12430
test_ps2 <- tax_glom(test_ps2_rarefied, "Genus")
test_ps2
# 378 taxa
range(taxa_sums(test_ps2))
# 1 609049
range(sample_sums(test_ps2))
# 12430 12430

# There are minor changes between the numbers
# Check this thread - https://github.com/joey711/phyloseq/issues/1077

# Test both the objects by plotting diversity indices
View(tax_table(test_ps_rarefied))
View(tax_table(test_ps2))

sample_variables(test_ps2)
unique(sample_data(test_ps2)[,"Ty_Loc"])
# 11 observations for Ty_Loc column
# 7 observations for Type column

# Alpha diversity
alpha_2 <- c("Observed", "Shannon")
tl_comp <- list(c("I_HHC", "I_TB"), c("I_HHC", "I_TBW"), c("I_HHC", "I_TBM"), c("I_HHC", "T_HC"), c("I_HHC", "T_TB"), 
                c("I_HHC", "W_TB"), c("I_HHC", "W_TB2M"), c("I_HHC", "W_HC"), c("I_HHC", "S_CC"), c("I_HHC", "S_TB"),
                c("I_TB", "I_TBW"), c("I_TB", "I_TBM"), c("I_TB", "T_HC"), c("I_TB", "T_TB"), c("I_TB", "W_TB"), 
                c("I_TB", "W_TB2M"), c("I_TB", "W_HC"), c("I_TB", "S_CC"), c("I_TB", "S_TB"),
                c("I_TBW", "I_TBM"), c("I_TBW", "T_HC"), c("I_TBW", "T_TB"), c("I_TBW", "W_TB"), c("I_TBW", "W_TB2M"), 
                c("I_TBW", "W_HC"), c("I_TBW", "S_CC"), c("I_TBW", "S_TB"),
                c("I_TBM", "T_HC"), c("I_TBM", "T_TB"), c("I_TBM", "W_TB"), c("I_TBM", "W_TB2M"), c("I_TBM", "W_HC"), 
                c("I_TBM", "S_CC"), c("I_TBM", "S_TB"),
                c("T_HC", "T_TB"), c("T_HC", "W_TB"), c("T_HC", "W_TB2M"), c("T_HC", "W_HC"), c("T_HC", "S_CC"), 
                c("T_HC", "S_TB"),
                c("T_TB", "W_TB"), c("T_TB", "W_TB2M"), c("T_TB", "W_HC"), c("T_TB", "S_CC"), c("T_TB", "S_TB"),
                c("W_TB", "W_TB2M"), c("W_TB", "W_HC"), c("W_TB", "S_CC"), c("W_TB", "S_TB"),
                c("W_TB2M", "W_HC"), c("W_TB2M", "S_CC"), c("W_TB2M", "S_TB"),
                c("W_HC", "S_CC"), c("W_HC", "S_TB"),
                c("S_CC", "S_TB"))

# test_ps_rarefied object
alpha_Grp <- plot_richness(test_ps_rarefied, x = "Ty_Loc", measures = alpha_2) +
  geom_jitter(aes(color = Location), show.legend = F, width = 0.2, alpha = 0.8)
alpha_Grp$layers <- alpha_Grp$layers[-1]

alpha_Grp + geom_boxplot(aes(fill = Location, color = Location), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = five_shades) + 
  scale_color_manual(values = five_shades, guide = "none") + custom_text +
  stat_compare_means(comparisons = tl_comp, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1)
  # labs(x = "", fill = "Ty_Loc")

# test_ps2 object
alpha_Grp2 <- plot_richness(test_ps2, x = "Ty_Loc", measures = alpha_2) +
  geom_jitter(aes(color = Location), show.legend = F, width = 0.2, alpha = 0.8)
alpha_Grp2$layers <- alpha_Grp2$layers[-1]

alpha_Grp2 + geom_boxplot(aes(fill = Location, color = Location), alpha = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = five_shades) + 
  scale_color_manual(values = five_shades, guide = "none") + custom_text +
  stat_compare_means(comparisons = tl_comp, label = "p.signif") +
  stat_compare_means(label.y.npc = 0, label.x.npc = 0.1)
  # labs(x = "", fill = "Ty_Loc")

# Beta diversity
# Bray-curtis
# test_ps_rarefied object
ord_bray <- ordinate(test_ps_rarefied, "PCoA", "bray")
pcoa_bray <- plot_ordination(test_ps_rarefied, ord_bray, color = "Location")
pcoa_bray + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = five_shades) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Location"))

pcoa_bray <- plot_ordination(test_ps_rarefied, ord_bray, color = "Type")
pcoa_bray + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = colors_30) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Type"))

# test_ps2 object
ord_bray2 <- ordinate(test_ps2, "PCoA", "bray")
pcoa_bray2 <- plot_ordination(test_ps2, ord_bray2, color = "Location")
pcoa_bray2 + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = five_shades) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Location"))

pcoa_bray2 <- plot_ordination(test_ps2, ord_bray2, color = "Type")
pcoa_bray2 + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = colors_30) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Type"))

# Weighted UniFrac
# test_ps_rarefied object
ord_uf <- ordinate(test_ps_rarefied, "PCoA", "unifrac", weighted = TRUE)
pcoa_uf <- plot_ordination(test_ps_rarefied, ord_uf, color = "Location")
pcoa_uf + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = five_shades) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Location"))

pcoa_uf <- plot_ordination(test_ps_rarefied, ord_uf, color = "Type")
pcoa_uf + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = colors_30) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Type"))

# test_ps2 object
ord_uf2 <- ordinate(test_ps2, "PCoA", "unifrac", weighted = TRUE)
pcoa_uf2 <- plot_ordination(test_ps2, ord_uf2, color = "Location")
pcoa_uf2 + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = five_shades) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Location"))

pcoa_uf2 <- plot_ordination(test_ps2, ord_uf2, color = "Type")
pcoa_uf2 + geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = colors_30) +
  stat_ellipse(type = "t", linetype = 2, alpha = 0.5) + custom_text +
  guides(color=guide_legend(title="Type"))

# taxa_plots
# Genus level
# test_ps_rarefied object
test_ps_rarefied_agg <- aggregate_rare(test_ps_rarefied, level = "Genus", detection = 1/1000, prevalence = 0.63)
test_ps_rarefied_agg
bar_genera <- plot_composition(test_ps_rarefied_agg, average_by = "Ty_Loc", 
                               transform = "compositional", otu.sort = "abundance") +
  guides(fill = guide_legend(title = "Genus")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Abundance")
bar_genera + scale_fill_manual(values = colors_30_mod)

# test_ps2 object
test_ps2_agg <- aggregate_rare(test_ps2, level = "Genus", detection = 1/1000, prevalence = 0.63)
test_ps2_agg
bar_genera2 <- plot_composition(test_ps2_agg, average_by = "Ty_Loc", 
                               transform = "compositional", otu.sort = "abundance") +
  guides(fill = guide_legend(title = "Genus")) +
  # scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Abundance")
bar_genera2 + scale_fill_manual(values = colors_30_mod)

# Both the objects are almost similar
# In alpha diversity, test_ps2 showed an extra significant variation in both Observed as well as Shannon
# In beta diversity, test_ps2 showed more variation across x-axis when using weighted UniFrac
# Genus level bar plots are also mostly similar
# Proceeding with agglomeration followed by rarefaction (test_ps_rarefied object)
# =========== RAREFACTION TEST END ==============