###First, create a working directory to save scripts and data.
getwd()
setwd("/Users/raymondo.garcia")
#in Roar 
#setwd("/gpfs/scratch/[your PSU ID]/DAWG2021F-2nd/")

####Install and Load Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq","Biostrings","RDPutils","ggplot2","dplyr", "ggpubr", "vegan", "emmeans", "plyr"))

#In case you recieved the following message "is not available for R version 4.0.2", use the code below
install.packages("remotes")
remotes::install_github("jfq3/RDPutils")
devtools::install_github("tidyverse/ggplot2")

#Load Packages
library(phyloseq)
library(Biostrings)
library(RDPutils)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(emmeans)
library(plyr)

#Alternatively you can use p_load to load the packages installed, it can also try to install the packages if they are not install
#For the latter, include the: ..., install = TRUE
install.packages("pacman")
pacman::p_load(phyloseq, Biostrings, RDPutils, ggplot2, dplyr, ggpubr, vegan, emmeans)

####Import Metadata or Mapping File

#Read Metadata or Mapping file using read.table (). You could also use "Import Dataset"
Metadata <- read.table("pathname to your file", header = TRUE, row.names = 1, stringsAsFactors = F, na.strings = c(".", "-"))

#Verify Metadata or Mapping file structure. We want to make sure that R assigned the variables correctly
str(Metadata)

Meta$Sample_ID <- as.factor(Meta$Sample_ID) #changing to factor
Meta$Stream <- as.factor(Metadata$Stream)
Metadata$SampleType <- as.factor(Metadata$SampleType)

####Create a Phyloseq object. Here we merged the otu table, taxonomy table, and Metadata/Mapping File.
ps.bac <- merge_phyloseq(otu_table(yourOTUtable, taxa_are_rows = FALSE), tax_table(yourTaxtable), sample_data(yourMetadata))

#####DADA2 outputs a OTU/ASV table in the sample (rows) by OTU/ASV (columns) format. Many packages might required a OTU/ASV (rows) by sample (columns) by sample (columns) format, so to transpose use t(otu_table(nameoffile))
#Also, each column is a sequence, but we want to exchange the sequence for an 'ASV' label, to make things easier downstream
ps.bac@otu_table #to view your OTU table

#Duplicate the phyloseq object to preserve original
ps.bac1 <- ps.bac

#Change row names so that they are ASVs rather than individual nucleotide sequences
DNA <- Biostrings::DNAStringSet(taxa_names(ps.bac1))
names(DNA) <- taxa_names(ps.bac1)
ps.bac2 <- merge_phyloseq(ps.bac1, DNA)
taxa_names(ps.bac2) <- paste0("ASV", seq(ntaxa(ps.bac2)))

ps.bac2@otu_table #check to make sure it worked
ps.bac2@tax_table

####Order samples by the number of reads or counts or library size or depth of coverage
ordered(sample_sums(ps.bac2))

####Save RDS and Read RDS
saveRDS(ps.bac2, "PS_DAWG_11_5_2021_.rds")

#Read .RDS By using <-, we are restoring the RDS under the same name we used for our Phyloseq object.
#ps.bac <- readRDS("/Users/raymondo.garcia/Desktop/P. syringae/DAWG/PS_DAWG_11_5_2021_.rds")
ps.bac <- readRDS("/gpfs/scratch/[Your PSU ID]/DAWG2021F_2nd/ps_DAWG_2021_F.rds
######Now, let's start Filtering, Subsetting, and Performing Richnness, Ordination, and Multivariate analysis!

#First, show available taxa ranks in the dataset and create a table with number of features for each taxa.
rank_names(ps.bac)
table(tax_table(ps.bac)[, "Phylum"], exclude = NULL) #For Phylum

#Here we see that 2,636 features were annotated with a Phylum "NA", these features might be artifacts and should be filtered. 
#Other Phyla, including WS4, 10bav-F6, CK-2C2-2, and Schekmanbacteria, have only one features and should be remove too.

table(tax_table(ps.bac)[, "Kingdom"], exclude = NULL) #For Domain/Kingdom

#Again, here we see that "NA" features were annotated with a Kingdom "NA" and these features should be filtered.
#For the purpose of this workshop I will focus just on the bacterial community, thus I will filter the Archaea and Eukaryota (i.e. - Fungi)

####Remove Phyla with low number of features and "NA" and the Kingdom Archaea and Eukaryota. 
filterPhyla <- c("WS4", "10bav-F6", "CK-2C2-2", "Schekmanbacteria", NA, "Cyanobacteria") 
#Notice that we removed Cyanobacteria, as it might be host contaminants (i.e. - chloroplasts) released during rhizodeposition.

filterKingdom <- c("Archaea", "Eukaryota", NA)

#Filtering
ps.bac.filt <- subset_taxa(ps.bac, !Phylum %in% filterPhyla)
table(tax_table(ps.bac.filt)[, "Phylum"], exclude = NULL) #to verify if it worked

ps.bac.filt <- subset_taxa(ps.bac.filt, !Kingdom %in% filterKingdom)
table(tax_table(ps.bac.filt)[, "Kingdom"], exclude = NULL)

####Let's remove features with ambiguous Phylum annotation, such as "uncharacterized" or "NA".
ps.bac.filt@tax_table #to view our tax_table

#These were detected at the Order, Family and Genus levels.

ps.bac.filt <- subset_taxa(ps.bac.filt, !is.na(Order) & !Order %in% c("NA"))
ps.bac.filt <- subset_taxa(ps.bac.filt, !is.na(Family) & !Family %in% c("NA"))
ps.bac.filt <- subset_taxa(ps.bac.filt, !is.na(Genus) & !Genus %in% c("NA"))

ps.bac.filt@tax_table #to view our tax_table

####Let's verify the prevalence of each taxon, in other words, the number of samples in which a taxon appears at least once.

#Prevalence of each feature and store as data.frame
prevdf <- apply(X = otu_table(ps.bac.filt),
               MARGIN = ifelse(taxa_are_rows(ps.bac.filt), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

#Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps.bac.filt),
                    tax_table(ps.bac.filt))

#Average prevalence and total read counts of the features in Phyla
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Plot the Prevalence vs Total Abundance or Counts
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps.bac.filt),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + #Include a guess for parameter for prevalence treshold, here is 5%
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#Define a prevalence treshold based on the plot, here we used 5% of total samples
prevalenceThreshold <- 0.05*nsamples(ps.bac.filt)

#Finally, let's filter based on the prevalence threshold established
keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps.bac.filt.2 <- prune_taxa(keepTaxa, ps.bac.filt)

#Let's see how many taxa after filtering
length(get_taxa_unique(ps.bac.filt.2, taxonomic.rank = "Genus")) #382 Genera after filtering

###Plotting - Stacked Bar Plot

#In case that you would like to keep only the most abundant taxa, here the top 10 most abundant genera.
ps.bac.filt.top.genera <- tapply(taxa_sums(ps.bac.filt.2), tax_table(ps.bac.filt.2)[, "Genus"], sum, na.rm=TRUE)
ps.bac.filt.top.10.genera <- names(sort(ps.bac.filt.top.genera, TRUE)) [1:10]
ps.bac.filt.3 <- prune_taxa((tax_table(ps.bac.filt) [, "Genus"] %in% ps.bac.filt.top.10.genera), ps.bac.filt)

#Convert to simple proportions or relative abundance. 
ps.bac.filt.2.normalized <- transform_sample_counts(ps.bac.filt.2, function (x) {x/ sum(x)}) ###counts are divided by the total library size. 
ps.bac.filt.2.normalized@otu_table #to make sure it worked or
View(ps.bac.filt.2.normalized)

####Agglomerate taxa at the Genus level and convert to long format for plotting.
ps.bac.filt.agglo <- tax_glom(ps.bac.filt.2.normalized, taxrank = "Phylum")

ps.bac.filt.agglo.long <- ps.bac.filt.agglo %>% #this is a function of the dplyr package not R base
  tax_glom(taxrank = "Family") %>%                     
  psmelt() %>%                                         
  arrange(Family)

#Condense low abundance taxa into an "Other" category for the barplot
data <- psmelt(ps.bac.filt.agglo) #create a dataframe from the phyloseq object
data$Phylum <- as.character(data$Phylum) #convert to character
data$Phylum[data$Abundance < 0.01] <- "< 1% abund." #simple way to rename phyla with < 1% abundance
medians <- ddply(data, ~Phylum, function(x) c(median=median(x$Abundance))) #group dataframe by Genus, calculate median rel. abundance
remainder <- medians[medians$median <= 0.01,]$Phylum
data[data$Phylum %in% remainder,]$Phylum <- "Genera < 1% abundance"

#rename phyla with < 1% relative abundance
data$Phylum[data$Abundance < 0.01] <- "Phyla < 1% abundance"

#####Create a stacked bar plot to observe community composition at the Phylum Level.
colors <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
            "#AD6F3B", "#673770","#D14285", "#652926", "#C84248","#8569D5", 
            "#5E738F","#D1A33D", "#8A7C64", "#599861","#999999", "#E69F00",
            "#56B4E9", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", 
            "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", 
            "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
            "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", 
            "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", 
            "darkorange1", "cyan1", "darkgrey")

Phylum_Stacked_Bar_Plot <- ggplot(data = data, aes(x = Sample_ID, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position="stack") +
  theme_grey() + 
  theme () + 
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance") +
  xlab("Sample") + theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.5)) +
  ggtitle("Phylum Composition") + facet_grid (SampleType ~ Stream, scales = 'free')
print(Phylum_Stacked_Bar_Plot)

####Now, let's explore the alpha Diversity and have a better understanding of species richnness and evenness

#Rarefy to even depth. Other normalization approach besides utilizing simple proportions. 
ps.bac.filt.rare <- rarefy_even_depth(ps.bac.filt.2, sample.size = min(sample_sums(ps.bac.filt)),
                                      rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

#One option is to use the function plot_richness
plot_richness(ps.bac.filt.rare, x = "SampleType") #you can specify a sample variable, in this case "SampleType

#You can also choose some of the alpha diveristy indexes. I like ACE because it has a standard error bar; Chao1 estimate diversity from abundance data, highlights the importance of rare ASVs
plot_richness(ps.bac.filt.rare, measures = c("Shannon", "Chao1", "Observed", "ACE"), x = "SampleType", color = "SampleType", title = "") #you could add shape = "" for other variables

#Alternatively we can output estimates of alpha diveristy/richness and creating a data frame.
Bacteria_Richness_DF <- estimate_richness(ps.bac.filt.rare) 
Bacteria_Richness_DF$SampleType <- ps.bac.filt.rare@sam_data$SampleType

##Create a box plot for each of the richness estimates.
A <- ggplot(Bacteria_Richness_DF, aes(x=SampleType, y=Shannon, fill=SampleType)) + geom_boxplot() + theme(axis.title.x=element_blank()) +
  stat_summary(fun=mean, geom = "point", shape = 23, size = 4) +
  scale_fill_brewer(palette = "Dark2")
print(A)
B <- ggplot(Bacteria_Richness_DF, aes(x=SampleType, y=Chao1, fill=SampleType)) + geom_boxplot() + theme(axis.title.x=element_blank()) +
  stat_summary(fun=mean, geom = "point", shape = 23, size = 4) +
  scale_fill_brewer(palette = "Dark2")
print(B)
C <- ggplot(Bacteria_Richness_DF, aes(x=SampleType, y=Observed, fill=SampleType)) + geom_boxplot() + theme(axis.title.x=element_blank()) +
  stat_summary(fun=mean, geom = "point", shape = 23, size = 4) +
  scale_fill_brewer(palette = "Dark2")
print(C)
D <- ggplot(Bacteria_Richness_DF, aes(x=SampleType, y=Chao1, fill=SampleType)) + geom_boxplot() + theme(axis.title.x=element_blank()) +
  stat_summary(fun=mean, geom = "point", shape = 23, size = 4) +
  scale_fill_brewer(palette = "Dark2")
print(D)

#Install package ggpubr
install.packages("ggpubr") #arrange multiple ggplots over multiple pages
library(ggpubr)

Richness_Boxplot_Figure <- ggarrange(A, B, C, D, nrow = 2, ncol = 2, widths = 2, heights = 2, label.x = 0,
                                   labels = c("A","B", "C", "D"), common.legend = TRUE, legend = FALSE) + theme(axis.title.x=element_blank())
print(Richness_Boxplot_Figure)

#####Multivariate Community Analysis (PCoA, PERMANOVA, etc.)

####Now, let's normalize data. Abundance value transformation to relative abundance or simple proportions. 
####Count proportions outperformed rarefied counts in most simulations due to better sensitivity, 
####but also suffered from a higher rate of false positives at larger values of effect size.
ps.bac.filt.2.normalized <- transform_sample_counts(ps.bac.filt.2, function (x) {x/ sum(x)}) ###counts are divided by the total library size. 
ps.bac.filt.2.normalized@otu_table #to make sure it worked or
View(ps.bac.filt.2.normalized)

###First, we fix sample_ID in sample data and then we ordinate at the samples level using PCoA and Bray-Curtis distances.  
sample_data(ps.bac.filt.2.normalized)$SampleType <- factor (sample_data(ps.bac.filt.2.normalized)$SampleType, levels = c("Sediment","Water"))

###Unconstrained Ordination and plot ordination.
PCoA <- ordinate (ps.bac.filt.2.normalized, method =  "PCoA", distance = "bray") #method could be NMDS or other ordinations methods

PCoA_Plot <- plot_ordination(physeq = ps.bac.filt.2.normalized, ordination = PCoA, color = "SampleType", title = ""
) + scale_fill_manual (values = colors) + geom_point(aes(color = SampleType), alpha = 0.7, size =4) + theme_bw()
print(PCoA_Plot)

#Used this if you want to plot with ellipse of 0.95 confidence intervals.
PCoA_Plot + stat_ellipse(aes(lty=factor(SampleType)), type = "t", level = 0.95) + scale_linetype_manual(values=c(1,2,3,4,5,6)) + theme_bw()

#You could also use shape = "" to plot. Let's use shape = "Stream".
PCoA_Plot <- plot_ordination(physeq = ps.bac.filt.2.normalized, ordination = PCoA, color = "SampleType", shape = "Stream", title = ""
) + scale_fill_manual (values = colors) + geom_point(aes(color = SampleType), alpha = 0.7, size =4) + theme_bw()
print(PCoA_Plot)

set.seed(1) #to get repeatable p-values

#Calculate bray curtis distance matrix (i.e. - betadiversity)
Bray_Curtis <- phyloseq::distance(ps.bac.filt.2.normalized, method = "bray")

#Make a data frame from the sample data
Sample_DF <- data.frame(sample_data(ps.bac.filt.2.normalized))

#Perform a multivariate analysis, specifically a PERMANOVA, using the function adonis from the vegan package. 
adonis(Bray_Curtis ~ SampleType, data = Sample_DF) #p-value = 0.001, so we reject the null hypothesis that our samples/treatments have similar bacterial community composition

#Beta-dispersion measures the compositional variation of the microbiome among a group of samples
Beta_Dispersion <- betadisper(Bray_Curtis, Sample_DF$SampleType)
permutest(Beta_Dispersion) #p-value = 0.001, meaning that we can reject the null hypothesis that our groups have the same dispersion. 
                           #Basically is telling us that our results from the PERMANOVA might be due to differences in dispersion.

####Constrained Ordination and plot ordination

#Ordinate
CAP <- ordinate(
  physeq = ps.bac.filt.2.normalized, #you can use CCA,RDA, or others constrained ordination analysis
  method = "CAP",
  distance = Bray_Curtis,
  formula = ~ Time + Depth) #our response vairiable is ASV

#To select the environmental factors or predictors that contribute to a greater extent to the variance observed on the bacterial community. 
#Here we will plot both regardless, but if you have multiple environmental covariates and you only want to model significant ones use this.
ordistep(CAP)

#Plot constrained ordination
CAP_Plot <- plot_ordination(
  physeq = ps.bac.filt.2.normalized,
  ordination = CAP, color = "SampleID", title = "",
  axes = c(1,2)
) + scale_fill_manual(values = phylum_colors) #default is type = "samples", but if you include type = "split" you will create a split plot and you could include taxa (color = "Genus") in one side and samples (shape = "SampleID") in the other
                                              #if you include type = "biplot" it will plot taxa and samples in the same plot
#Add the environmental variables as arrows
Arrowenv <- vegan::scores(CAP, display = 'bp') 

#Add labels, make a data.frame
Arrowdf <- data.frame(labels = rownames(Arrowenv), Arrowenv)

#Define the arrow aesthetic mapping
Arrowmap <- aes(xend = CAP1,
                yend = CAP2,
                x = 0,
                y = 0,
                shape = NULL,
                color = NULL,
                label = labels)

#To make labels the same angle of the arrow. 
Arrowdf = transform(Arrowdf, radians = atan(CAP2/CAP1))
Arrowdf = transform(Arrowdf, angle = 360 * radians/(2 * pi))
#Quadrants II, III, IV
Arrowdf$quad234 <- apply(Arrowdf[, c("CAP1", "CAP2")], 1, function(x) {
  any(x < 0)
})
Arrowdf$quad4 <- apply(Arrowdf[, c("CAP1", "CAP2")], 1, function(x) {
  all(x < 0)
})
#If quadrant II, III, IV, add 180 to angle
if (any(Arrowdf$quad234)) {
  Arrowdf[Arrowdf$quad234, "angle"] <- Arrowdf[Arrowdf$quad234, "angle"] + 
    180
}
#If quadrant IV, add additional 180
if (any(Arrowdf$quad4)) {
  Arrowdf[Arrowdf$quad4, "angle"] <- Arrowdf[Arrowdf$quad4, "angle"] + 180
}
#For printing text, we want to flip if its greater than 90
if (any(Arrowdf$angle > 90)) {
  Arrowdf[Arrowdf$angle > 90, "angle"] <- Arrowdf[Arrowdf$angle > 90, "angle"] - 
    180
}

Labelmap <- aes(x = 1.3 * CAP1,
                y = 1.3 * CAP2,
                shape = NULL,
                color = NULL,
                label = labels,
                angle = angle)

Arrowhead = arrow(length = unit(0.01, "npc"))

#Make a new graphic
CAP_Plot + 
  geom_segment(
    mapping = Arrowmap,
    size = .5,
    data = Arrowdf,
    color = "black",
    arrow = Arrowhead
  ) + 
  geom_text(
    mapping = Labelmap,
    size = 2,
    data = Arrowdf,
    show.legend = FALSE
  )

#######Differential Abundance Analysis

#Negative binomial in microbiome differential abundance testing (https://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html)
BiocManager::install("DESeq2")
library("DESeq2"); packageVersion("DESeq2")

#Convert from phyloseq-format to DeSeqDataSet. For this, use counts not proportions.
ps.bac.filt.agglo.2 <- tax_glom(ps.bac.filt.2, taxrank = "Genus")

DeSeq2 <- phyloseq_to_deseq2(ps.bac.filt.agglo.2, ~ SampleType)
prune_samples(sample_sums(DeSeq2) > 500, DeSeq2) ##optional and parameters will change based on your data.

#Calculate  geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans = apply(counts(DeSeq2), 1, gm_mean)
DeSeq2 = estimateSizeFactors(DeSeq2, geoMeans = geoMeans)
DeSeq2 = DESeq(DeSeq2, fitType = "local")

##Create a table of the results of the tests.
results_DeSeq2 <- results(DeSeq2)
summary(results_DeSeq2)

results_DeSeq2 <- results_DeSeq2[order(results_DeSeq2$padj, na.last = NA), ]
alpha <- 0.05
sigtab <- results_DeSeq2[(results_DeSeq2$padj < alpha), ]
sigtab <- cbind(as(sigtab, "data.frame"), as(tax_table(ps.bac.filt.2)[rownames(sigtab), 
], "matrix"))
head(sigtab)
print(sigtab)

#Let’s look at just the OTUs that were significantly enriched. First, cleaning up the table a little for legibility.
posigtab <- sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab <- posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
head(posigtab)
print(posigtab)

##To make a plot
ggplot(sigtab, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

####Metacoder
install.packages("devtools")
devtools::install_github("grunwaldlab/metacoder")
library(metacoder)

#Extract abundance matrix from the phyloseq object
ASV_Table <- as(otu_table(ps.bac.filt.2), "matrix")
Tax_Table <- as(tax_table(ps.bac.filt.2), "matrix")
Metadata <- as(sample_data(ps.bac.filt.2), "matrix")

#transpose if necessary
ASV_Table_t <- t(ASV_Table)

#Coerce to data.frame
ASV_df = as.data.frame(ASV_Table_t)
Tax_df = as.data.frame(Tax_Table)
Metadata = as.data.frame(Metadata)
Metadata$SampleName <- rownames(Metadata)

#Merge OTU and Taxonomy tables
install.packages("data.table")
library(data.table)

ASVTAX_Table <- merge(ASV_df, Tax_df, by = 0, all = TRUE)
setnames(ASVTAX_Table, old = "Row.names", new = "ASV_ID")

###Created a taxmap object.
Tax_Map = parse_tax_data(ASVTAX_Table, class_cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), class_sep = " ",
                              class_key = c(tax_name = "taxon_name"))

####Convert to proportions
Tax_Map$data$otu_props <- calc_obs_props(Tax_Map, data = "tax_data", cols = Metadata$SampleName)

####Get per-taxon counts 
Tax_Map$data$tax_abund <- calc_taxon_abund(Tax_Map, data = "otu_props", cols = Metadata$SampleName)

####Calculate differences between groups and multiple test correction.
Tax_Map$data$diff_table <- compare_groups(Tax_Map, data = "tax_abund", cols = Metadata$SampleName, groups = Metadata$SampleType)
Tax_Map$data$diff_table$wilcox_p_value <- p.adjust(Tax_Map$data$diff_table$wilcox_p_value, method = "fdr")
range(Tax_Map$data$diff_table$wilcox_p_value, finite = TRUE)
Tax_Map$data$diff_table$log2_median_ratio[Tax_Map$data$diff_table$wilcox_p_value > 0.05] <- 0

set.seed(999)

heat_tree(Tax_Map, 
          node_label = taxon_names,
          node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
          node_color = log2_median_ratio, # A column from `obj$data$diff_table`
          node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
          node_color_range = c("cyan", "gray", "tan"), # The color palette used
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

#Let's do the same but with more than two treatment groups.

####Calculate differences between groups and multiple test correction.
Tax_Map$data$diff_table <- compare_groups(Tax_Map, data = "tax_abund", cols = Metadata$SampleName, groups = Metadata$Stream)
Tax_Map$data$diff_table$wilcox_p_value <- p.adjust(Tax_Map$data$diff_table$wilcox_p_value, method = "fdr")
range(Tax_Map$data$diff_table$wilcox_p_value, finite = TRUE)
Tax_Map$data$diff_table$log2_median_ratio[Tax_Map$data$diff_table$wilcox_p_value > 0.05] <- 0

heat_tree_matrix(Tax_Map,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = diverging_palette(),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford")

