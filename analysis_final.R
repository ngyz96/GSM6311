##### loading libraries for general plotting and data wrangling #####
library(mia)
library(miaViz)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(vegan)
library(scater)
library(knitr)
library(kableExtra)
set.seed(100) # required for permutation based analysis

##### Helper functions to read in the files #####
# function to add taxid column to each NCBI ref number, collapse NCBI entries to the same Taxid
convert_data <- function(df, taxa){
    df$Taxid <- taxa$Taxid[match(df$Ref,taxa$Ref)]
    df <- aggregate(Count ~ Taxid, df, sum)
    df
}

# Function to read in alignment files with specified prefix in directory, default to current directory
# uses convert_data
read_data <- function(prefix, taxa, dir) {
    if(missing(dir))
        dir <- getwd()
    file_list <- list.files(dir)
    file_list <- file_list[startsWith(basename(file_list), prefix)]
    file_names <- sapply(strsplit(basename(file_list), "_"), `[`, 1)
    columns <- c("Taxid")
    main_df <- data.frame(matrix(nrow=0, ncol = length(columns))) 
    colnames(main_df) <- columns
    for (i in file_names){
        temp_df <- read.delim(file_list[which(file_names==i)])[,c(1,2)]
        temp_df <- convert_data(temp_df, taxa)
        main_df <- merge(main_df, temp_df, by="Taxid", all=T)
    }
    colnames(main_df) <- c("Taxid", file_names)
    rownames(main_df) <- main_df$Taxid
    main_df <- main_df[,!(colnames(main_df) %in% columns)]
    main_df <- as.data.frame(t(main_df))
    main_df[is.na(main_df)] <- 0
    main_df  
}


##### Loading in data #####
# Read in taxonomy required files
taxa <- read.delim("ref2taxid.targloci.tsv", col.names = c("Ref", "Taxid"))
# Read in all primary aligned files, starts with barcode by default
df <- read_data("barcode", taxa)
# remove control sample
control <- df[(row.names(df) %in% c("barcode03")),] 
df <- df[!(row.names(df) %in% c("barcode03")),] 
df <- df[,!colSums(df) == 0]

# DO NOT RUN unless needed
# Get taxonomy mapping (taxid -> taxonomy classes), only run when needed
# library(taxonomizr)
# taxids <- as.numeric(colnames(df))
# tmp <- sapply(taxids, getTaxonomy,'accessionTaxa.sql')
# tmp <- t(tmp)
# rownames(tmp) <- taxids
# colnames(tmp) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# tmp <- gsub(" ", "_", tmp)
# saveRDS(tmp, "collapse_taxid_map.rds")

# RUN FROM HERE
# load in taxonomy mapping files with only IDs present in sample
collapse_taxid_map <- readRDS("collapse_taxid_map.rds") 

# sanity check to ensure taxonomic information are all present
all(rownames(collapse_taxid_map) == colnames(df))

# load in clinical data 
metadata <- read.delim("metadata.tsv")
mapping <- read.delim("mapping.tsv")
metadata <- dplyr::filter(metadata,SampleID %in% mapping$Sample_ID)
rownames(metadata) <- mapping$barcode[match(metadata$SampleID,mapping$Sample_ID)]
metadata$sample <- rownames(metadata)

#create TSE object for mia
tax <- DataFrame(collapse_taxid_map)
counts <- t(as.matrix(df))
counts <- counts[rownames(tax), rownames(metadata)]
tse <- TreeSummarizedExperiment(assays =  SimpleList(counts = counts),
                                     colData = metadata,
                                     rowData = tax)
tse_nonfilt <- tse

# QC plotting
tse_nonfilt <- addPerCellQC(tse_nonfilt)
df <- as.data.frame(colData(tse_nonfilt))
raw_reads <- read.delim("track_reads.tsv")
df <- merge(df,raw_reads)
df <- dplyr::arrange(df,sum)
df <- df[,c("sample", "sum", "reads")]
df <- tidyr::pivot_longer(df, cols = c("sum", "reads"),
                          names_to = "type", values_to = "count")
read_histo <- ggplot(df, aes(fill = type, y = count, x = sample)) +
    geom_bar(position="stack", stat="identity", alpha=0.65) +  
    labs(x = "Samples index", y = "Read counts", fill = "Read mapping") + 
    theme_minimal() + scale_fill_manual(values=c('#FF9E4A','#729ECE', alpha), labels = c("Unmapped", "Mapped")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "#A2A2A2")) +
    scale_y_continuous(labels = ~ format(.x, scientific = FALSE))
rare_data <- rarecurve(t(assay(tse_nonfilt)), step=25, tidy = T)
rare_label <- dplyr::group_by(rare_data, Site)
rare_label <-  dplyr::summarize(rare_label, pos = which.max(Sample),
              x = Sample[pos],
              y = Species[pos])
rare_curve <- ggplot(rare_data, aes(x = Sample, y = Species, group=Site)) + 
    geom_line(linewidth=0.2) + 
    xlab("Number of Reads") +
    ylab("Number of Species Observed") +
    theme_minimal() +     
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "#A2A2A2")) 
wrap_plots(read_histo, rare_curve) + plot_annotation(tag_levels = list(c('B','C'))) & theme(plot.tag = element_text(size = 14))

# filter for at least taxa with at least 10 reads across all samples
tse <- tse[rowSums(assay(tse, "counts")) >= 10, ]


##### Code to get patient data for table 1 ######
# Age
median(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "Age"]) # 9.030055
median(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "Age"]) # 10.97777
iqr(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "Age"]) # 10.58082
iqr(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "Age"]) # 10.18659
wilcox.test(Age ~ RefractoryStatus, data = colData(tse)) #0.5063

# AgeGroup
table(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "AgeGroup"]) / 15
table(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "AgeGroup"]) / 8
fisher.test(table(colData(tse)[,c("RefractoryStatus", "AgeGroup")])) #0.369

# Gender
table(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "Gender"]) / 15
table(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "Gender"]) / 8
fisher.test(table(colData(tse)[,c("RefractoryStatus", "Gender")])) #0.685

# Race
table(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "Race"]) / 15
table(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "Race"]) / 8
fisher.test(table(colData(tse)[,c("RefractoryStatus", "Race")])) #0.1216

# SeizureType
table(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "SeizureType"]) / 15 
table(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "SeizureType"]) / 8 
fisher.test(table(colData(tse)[,c("RefractoryStatus", "SeizureType")])) #0.6815

# SeizureFrequency
table(colData(tse)[colData(tse)$RefractoryStatus == "Refractory", "SeizureFrequency"]) / 15 
table(colData(tse)[colData(tse)$RefractoryStatus == "Susceptible", "SeizureFrequency"]) / 8 
fisher.test(table(colData(tse)[,c("RefractoryStatus", "SeizureFrequency")])) #0.6815

#ActiveAEDs
drug_re <- as.data.frame(table(unlist(lapply(strsplit(colData(tse)[colData(tse)$RefractoryStatus=="Refractory", "ActiveAEDs"], split = ","), trimws))))
drug_sus <- as.data.frame(table(unlist(lapply(strsplit(colData(tse)[colData(tse)$RefractoryStatus=="Susceptible", "ActiveAEDs"], split = ","), trimws))))


##### Alpha diversity #####
tse <- estimateRichness(tse, assay.type = "counts", 
                             index = "chao1", name="Chao1")
tse <- estimateEvenness(tse, assay.type = "counts", 
                        index="simpson", name = "Simpson")
tse <- estimateDiversity(tse, assay.type = "counts",
                              index = "shannon", name = "Shannon")

plots_alpha <- lapply(c("Chao1", "Simpson", "Shannon"),
                plotColData, 
                object = tse,
                x = "RefractoryStatus",
                colour_by = "RefractoryStatus",
                shape="AgeGroup", point_size=2)

plots_alpha <- lapply(plots_alpha, "+", 
                theme(axis.text.x = element_blank(),
                      axis.title.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.y = element_blank()))

plots_alpha <- lapply(plots_alpha, "+", geom_signif(comparisons = list(c("Refractory", "Susceptible")),
                                                    tip_length=0,test="wilcox.test", 
                                                    map_signif_level = function(p) sprintf("p = %.3g", p)))

plots_alpha[[1]]<-plots_alpha[[1]]+ggtitle("Chao1") + theme(legend.position="none")
plots_alpha[[2]]<-plots_alpha[[2]]+ggtitle("Simpson") + theme(legend.position="none")
plots_alpha[[3]]<-plots_alpha[[3]]+ggtitle("Shannon") + theme(legend.position="none")
tse <- addPerCellQC(tse)
cor.test(colData(tse)$sum, colData(tse)$Chao1, method = "pearson")
cor.test(colData(tse)$sum, colData(tse)$Simpson, method = "pearson")
cor.test(colData(tse)$sum, colData(tse)$Shannon, method = "pearson")
wilcox.test(Chao1 ~ AgeGroup, data=colData(tse))



##### Beta Diversity #####
tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
tse <- transformAssay(tse, assay.type = "relabundance",
                      method = "clr", pseudocount = TRUE)

# Run NMDS on clr assay with Euclidean distances
tse <- runNMDS(tse,
               FUN = vegan::vegdist,
               method = "euclidean",
               assay.type = "clr",
               name = "NMDS_aitchison",ncomponents=3)

# stress testing NMDS fit, < 0.2 generally acceptable, 0.1-0.2 is fair
nmds_stress <- metaMDS(t(assay(tse, "clr")), distance = "euclidean") #0.14943
nmds_stress_3 <- metaMDS(assay(tse, "clr"), distance = "euclidean", k=3) #0.1016405 
nmds_stress_4 <- metaMDS(assay(tse, "clr"), distance = "euclidean", k=4) #0.08362942
nmds_stress_5 <- metaMDS(assay(tse, "clr"), distance = "euclidean", k=5) #0.07126484

# Create beta diversity plot based on aitchison distance with NMDS
plot_beta <- plotReducedDim(tse, "NMDS_aitchison", colour_by = "RefractoryStatus",
                            shape_by = "AgeGroup", point_size=3, text_by = "sample", text_size = 3) +
    stat_ellipse(aes(group = colour_by, color=colour_by))

# plot alpha + beta diversity
design <- "
123
444
"
# Figure 3 
png(filename = "diversity_plot.png", width = 3000, height = 2880, res=300)
plots_alpha[[1]] + plots_alpha[[2]] + plots_alpha[[3]] + 
    plot_beta + plot_layout(design = design, guides = "collect") +
plot_annotation(tag_levels = list(c('A','','','B'))) & theme(plot.tag = element_text(size = 14))
dev.off()

# Sanity check for PERMANOVA usage
anova(betadisper(vegdist(t(assay(tse, "counts")), method = "aitchison", pseudocount=4.96073577633035e-06), colData(tse)$Gender)) # p=0.7657
anova(betadisper(vegdist(t(assay(tse, "counts")), method = "aitchison", pseudocount=4.96073577633035e-06), colData(tse)$Race)) # p=2.823e-05
anova(betadisper(vegdist(t(assay(tse, "counts")), method = "aitchison", pseudocount=4.96073577633035e-06), colData(tse)$AgeGroup)) # p=0.003094 **
anova(betadisper(vegdist(t(assay(tse, "counts")), method = "aitchison", pseudocount=4.96073577633035e-06), colData(tse)$RefractoryStatus)) # p=0.08381
# anova(betadisper(vegdist(t(assay(tse, "counts")), method = "aitchison", pseudocount=4.96073577633035e-06), colData(tse)$FeedingMethod)) # p=0.000291
anova(betadisper(vegdist(t(assay(tse, "counts")), method = "aitchison", pseudocount=4.96073577633035e-06), colData(tse)$SeizureType)) # p=0.1184

# # Species
# tse_species <- mergeFeaturesByRank(tse, rank = "Species")
# tse_species <- transformAssay(tse_species, assay.type = "counts", method = "relabundance")
# tse_species <- transformAssay(tse_species, assay.type = "relabundance",method = "clr", pseudocount = TRUE)
# set.seed(100)
# permanova_species <- adonis2(t(assay(tse_species, "clr"))~Gender+Race+SeizureType+AgeGroup+RefractoryStatus,
#         by = "margin",
#         data = colData(tse_species),
#         method = "euclidean",
#         permutations = 100000)

# Genus
tse_genus <- mergeFeaturesByRank(tse, rank = "Genus")
tse_genus <- transformAssay(tse_genus, assay.type = "counts", method = "relabundance")
tse_genus <- transformAssay(tse_genus, assay.type = "relabundance",method = "clr", pseudocount = TRUE)
set.seed(100)
permanova_genus <- adonis2(t(assay(tse_genus, "clr"))~Gender+Race+SeizureType+AgeGroup+RefractoryStatus,
                     by = "margin",
                     data = colData(tse_genus),
                     method = "euclidean",
                     permutations = 100000)
kable(permanova_species, format="html",booktabs = TRUE,caption = "PERMANOVA at Genus level ~ Gender+Race+SeizureType+AgeGroup+RefractoryStatus") %>% 
    kable_styling("striped")# %>% save_kable("permanova.png") Table 3 1025x380

# Phylum
tse_phylum <- mergeFeaturesByRank(tse, rank = "Phylum")
tse_phylum <- transformAssay(tse_phylum, assay.type = "counts", method = "relabundance")
tse_phylum <- transformAssay(tse_phylum, assay.type = "relabundance",method = "clr", pseudocount = TRUE)
set.seed(100)
permanova_phylum <- adonis2(t(assay(tse_phylum, "clr"))~Gender+Race+SeizureType+AgeGroup+RefractoryStatus,
        by = "margin",
        data = colData(tse_phylum),
        method = "euclidean",
        permutations = 100000)
kable(permanova_phylum, format="html",booktabs = TRUE,caption = "PERMANOVA at Phylum level ~ Gender+Race+SeizureType+AgeGroup+RefractoryStatus") %>% 
    kable_styling("striped") # %>% save_kable("permanova.png") Table 4 1025x400


##### Differential Abundance testing #####
library(ALDEx2)
set.seed(100)
tse_aldex <- subsetByPrevalentFeatures(tse,
                                         rank = "Genus",
                                         prevalence = 10/ 100)

x <- aldex.clr(assay(tse_aldex, "counts"), tse_aldex$RefractoryStatus)     
x_tt <- aldex.ttest(x, paired.test = FALSE)
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
aldex_out <- data.frame(x_tt, x_effect)
par(mfrow = c(1, 3))

aldex.plot(aldex_out,
           type = "MA",
           test = "effect",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "MW",
           test = "effect",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "volcano",
           test = "welch",
           xlab = "Difference",
           cutoff = 0.05)

library(ANCOMBC)
tse_ancombc <- subsetByPrevalentFeatures(tse,
                                       rank = "Genus",
                                       prevalence = 10/ 100)
ancombc2_out <- ancombc2(data = tse,
                         assay.type = "counts",
                         fix_formula = "RefractoryStatus+sum",
                         p_adj_method = "fdr",
                         prv_cut = 0,
                         group = "RefractoryStatus",
                         struc_zero = TRUE,
                         neg_lb = TRUE,
                         # multi group comparison is deactivated automatically
                         global = TRUE)
ancombc2_out$res %>%
    dplyr::select(taxon, lfc_RefractoryStatusSusceptible, q_RefractoryStatusSusceptible) %>%
    dplyr::filter(q_RefractoryStatusSusceptible < 0.05)



##### Taxonomic composition plotting #####
tse_phylum <- mergeFeaturesByRank(tse, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopFeatures(tse_phylum,top = 6, assay.type = "relabundance")
phylum_renamed <- lapply(rowData(tse)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
tse_phylum <- tse
rowData(tse_phylum)$Phylum <- as.character(phylum_renamed)

plot_tax <- plotAbundance(tse_phylum, assay.type="relabundance", rank = "Phylum",
                          order_rank_by="abund",  order_sample_by = "Bacteroidota", 
                          features="AgeGroup", decreasing=F)
plot_tax[[1]] <- plot_tax[[1]] + theme(axis.title.x = element_blank())
plot_tax[[2]] <- plot_tax[[2]] + theme(axis.title=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.ticks.x=element_blank(),
                                       axis.line=element_blank())
plot_tax_age <- plot_tax
plot_tax <- plotAbundance(tse_phylum, assay.type="relabundance", rank = "Phylum",
                          order_rank_by="abund",  order_sample_by = "Bacteroidota",
                          features="RefractoryStatus",decreasing=F)
plot_tax[[1]] <- plot_tax[[1]] + theme(axis.title.x = element_blank())
plot_tax[[2]] <- plot_tax[[2]] + theme(axis.title=element_blank(),
                                       axis.text.x=element_blank(),
                                       axis.ticks.x=element_blank(),
                                       axis.line=element_blank())
plot_tax_ref <- plot_tax
plot_tax_all <- wrap_plots(list(plot_tax_ref[[1]],plot_tax_ref[[2]],plot_tax_age[[2]]),
                           ncol = 1, heights=c(0.8, 0.1, 0.1)) +
    plot_layout(guides = "collect") 
#Figure 4
png(filename = "taxa_plot.png", width = 3000, height = 2880, res=300)
plot_tax_all
dev.off()


##### Old group only analysis #####
tse_old <- tse_nonfilt[ ,tse_nonfilt$AgeGroup == "Old"] #17 samples
tse_old <- tse_old[rowSums(assay(tse_old, "counts")) >= 10, ]
tse_old <- tse_old[ ,colSums(assay(tse_old, "counts")) != 0]
tse_old <- transformAssay(tse_old, assay.type = "counts", method = "relabundance")
tse_old <- transformAssay(tse_old, assay.type = "relabundance",
                      method = "clr", pseudocount = TRUE)
# alpha diversity
tse_old <- estimateRichness(tse_old, assay.type = "counts", 
                        index = "chao1", name="Chao1")
tse_old <- estimateEvenness(tse_old, assay.type = "counts", 
                        index="simpson", name = "Simpson")
tse_old <- estimateDiversity(tse_old, assay.type = "counts",
                         index = "shannon", name = "Shannon")

plots_alpha_old <- lapply(c("Chao1", "Simpson", "Shannon"),
                      plotColData, 
                      object = tse_old,
                      x = "RefractoryStatus",
                      colour_by = "RefractoryStatus",
                      point_size=2)

plots_alpha_old <- lapply(plots_alpha_old, "+", 
                      theme(axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.title.y = element_blank()))

plots_alpha_old <- lapply(plots_alpha_old, "+", geom_signif(comparisons = list(c("Refractory", "Susceptible")),
                                                    tip_length=0,test="wilcox.test", 
                                                    map_signif_level = function(p) sprintf("p = %.3g", p)))

plots_alpha_old[[1]]<-plots_alpha_old[[1]]+ggtitle("Chao1") + theme(legend.position="none")
plots_alpha_old[[2]]<-plots_alpha_old[[2]]+ggtitle("Simpson") + theme(legend.position="none")
plots_alpha_old[[3]]<-plots_alpha_old[[3]]+ggtitle("Shannon") + theme(legend.position="none")

# Run NMDS on clr assay with Euclidean distances
tse_old <- runNMDS(tse_old,
               FUN = vegan::vegdist,
               method = "euclidean",
               assay.type = "clr",
               name = "NMDS_aitchison",ncomponents=3)

# stress testing NMDS fit, < 0.2 generally acceptable, 0.1-0.2 is fair
# nmds_stress_test <- metaMDS(t(assay(tse_old, "clr")), distance = "euclidean") #0.1546646 
# nmds_stress_t2 <- metaMDS(assay(tse_old, "clr"), distance = "euclidean", k=3) #0.1046511  


# Create beta diversity plot based on aitchison distance with NMDS
plot_beta_old <- plotReducedDim(tse_old, "NMDS_aitchison", colour_by = "RefractoryStatus",
                            point_size=3, text_by = "sample", text_size = 3) +
    stat_ellipse(aes(group = colour_by, color=colour_by))

# plot alpha + beta diversity
design <- "
123
444
"
# Figure 5
png(filename = "diversity_plot_old.png", width = 3000, height = 2880, res=300)
plots_alpha_old[[1]] + plots_alpha_old[[2]] + plots_alpha_old[[3]] + 
    plot_beta_old + plot_layout(design = design, guides = "collect") +
    plot_annotation(tag_levels = list(c('A','','','B'))) & theme(plot.tag = element_text(size = 14))
dev.off()

# Sanity check for PERMANOVA usage
anova(betadisper(vegdist(t(assay(tse_old, "counts")), method = "aitchison", pseudocount=4.96105571265565e-06), colData(tse_old)$Gender)) # p=0.03962
anova(betadisper(vegdist(t(assay(tse_old, "counts")), method = "aitchison", pseudocount=4.96105571265565e-06), colData(tse_old)$Race)) # p=9.299e-07 
anova(betadisper(vegdist(t(assay(tse_old, "counts")), method = "aitchison", pseudocount=4.96105571265565e-06), colData(tse_old)$RefractoryStatus)) # p=0.1478
# anova(betadisper(vegdist(t(assay(tse_old, "counts")), method = "aitchison", pseudocount=4.96105571265565e-06), colData(tse_old)$FeedingMethod)) # p=0.1362
anova(betadisper(vegdist(t(assay(tse_old, "counts")), method = "aitchison", pseudocount=4.96105571265565e-06), colData(tse_old)$SeizureType)) #5.835e-06 

# Genus
tse_old_genus <- mergeFeaturesByRank(tse_old, rank = "Genus")
tse_old_genus <- transformAssay(tse_old_genus, assay.type = "counts", method = "relabundance")
tse_old_genus <- transformAssay(tse_old_genus, assay.type = "relabundance",method = "clr", pseudocount = TRUE)
set.seed(100)
permanova_old_genus <- adonis2(t(assay(tse_old_genus, "clr"))~Gender+Race+SeizureType+RefractoryStatus,
                           by = "margin",
                           data = colData(tse_old_genus),
                           method = "euclidean",
                           permutations = 100000)
kable(permanova_old_genus, format="html",booktabs = TRUE,caption = "PERMANOVA for Old age group at Genus level ~ Gender+Race+SeizureType+RefractoryStatus") %>% 
    kable_styling("striped") # %>% save_kable("permanova.png") #Table 3 1025x340

# differential abundance
# set.seed(100)
# tse_aldex <- subsetByPrevalentFeatures(tse_old_genus,
#                                        rank = "Genus",
#                                        prevalence = 10/ 100)
# 
# x <- aldex.clr(assay(tse_aldex, "counts"), tse_aldex$RefractoryStatus)     
# x_tt <- aldex.ttest(x, paired.test = FALSE)
# x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
# aldex_out <- data.frame(x_tt, x_effect)
# par(mfrow = c(1, 2))
# 
# aldex.plot(aldex_out,
#            type = "MA",
#            test = "welch",
#            xlab = "Log-ratio abundance",
#            ylab = "Difference",
#            cutoff = 0.05)
# 
# aldex.plot(aldex_out,
#            type = "MW",
#            test = "welch",
#            xlab = "Dispersion",
#            ylab = "Difference",
#            cutoff = 0.05)
# 
# aldex.plot(aldex_out,
#            type = "volcano",
#            test = "welch",
#            xlab = "Difference",
#            cutoff = 0.05)
# 
# aldex_out %>%
#     tibble::rownames_to_column(var = "Genus") %>%
#     dplyr::filter(we.eBH <= 0.05)  %>% #authors recommend welch as it is not as sensitive to sample size
#     dplyr::select(Genus, we.eBH, wi.eBH, effect, effect.high, effect.low, overlap) %>%
#     knitr::kable(format="html",booktabs = TRUE) %>%
#     kable_styling("striped")


##### LEFSE #####
library(lefser)
lefse_out <- lefser(tse, groupCol = "RefractoryStatus", assay = "relabundance")
library(phyloseq)
phy <- makePhyloseqFromTreeSE(tse, assay.type="relabundance")
counts <- unclass(otu_table(phy))
colData <- as(sample_data(phy), "data.frame")
## create a SummarizedExperiment object
sme <- SummarizedExperiment(
    assays = list(counts = counts), colData = colData
)
lefse_out <- lefser(sme, groupCol = "RefractoryStatus")
lefserPlot(lefse_out)
# Lefse_input<-print(assay(tse, "counts"))
# Lefse_input<-cbind(rownames(Lefse_input), Lefse_input)
# Lefse_input<-rbind(c("Group", tse$RefractoryStatus), colnames(Lefse_input), Lefse_input)
# write.table(Lefse_input, "Lefse_input.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)

##### control stuff #####
control_tse <- TreeSummarizedExperiment(assays =  SimpleList(counts = t(as.matrix(control))),
colData=control_data, rowData = DataFrame(control_taxid))
control_tse <- mergeFeaturesByRank(control_tse, rank = "Genus")
control_tse <- transformAssay(control_tse, assay.type = "counts", method = "relabundance")
top_taxa <- getTopFeatures(control_tse,top = 10, assay.type = "relabundance")
top_taxa <- sub('......', '', top_taxa)
phylum_renamed <- lapply(rowData(control_tse)$Genus,
function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(control_tse)$Genus <- as.character(phylum_renamed)
plotAbundance(control_tse, assay.type="relabundance", rank = "Genus",
order_rank_by="abund")
control_df <- assay(control_tse, "relabundance")
control_df
                     
test_kingdom <- mergeFeaturesByRank(tse, rank = "Kingdom")
test_phylum <- mergeFeaturesByRank(tse, rank = "Phylum")
test_class <- mergeFeaturesByRank(tse, rank = "Class")
test_order <- mergeFeaturesByRank(tse, rank = "Order")
test_family <- mergeFeaturesByRank(tse, rank = "Family")
test_genus <- mergeFeaturesByRank(tse, rank = "Genus")
test_species <- mergeFeaturesByRank(tse, rank = "Species")


write.table(rbind(colData(tse)$RefractoryStatus, colData(tse)$AgeGroup, colData(tse)$sample, 
                  assay(test_kingdom, "counts"), assay(test_phylum, "counts"), assay(test_class, "counts"), 
                  assay(test_order, "counts"), assay(test_family, "counts"), assay(test_genus, "counts"),
                  assay(test_species, "counts")), file="lefse.tsv", quote=FALSE, sep='\t', col.names = NA)