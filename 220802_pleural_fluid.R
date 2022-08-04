# 8/2/2022; Benjamin Kwok; Segal Lab
# Microbial Signatures of Malignant Pleural Effusions
# For Github publication

# Load/Save/wd-----
setwd() # set your own wd
# save.image(file="220802_pleural_fluid/220802_pleural_fluid.RData")
load(file="220802_pleural_fluid/220802_pleural_fluid.RData")

# Load libraries
library(pacman)
pacman::p_load(ade4, vegan, phyloseq, biomformat, devtools, readr, readtext, qiime2R, ape, DESeq2, edgeR, limma, Glimma, gplots, RColorBrewer, pheatmap, ggplot2, tidyverse, ggrepel, scales, data.table, fBasics, forcats, omu, maptools, reshape2, dplyr, gridExtra, tibble, formattable, htmltools, webshot, splitstackshape, decontam, grid, cowplot, colorspace, ggpubr, ggpmisc, dunn.test, survival, survminer, lubridate, lefser, microbiomeMarker,
               install= FALSE,
               update = FALSE)

'%!in%' <- function(x,y)!('%in%'(x,y))

## To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}

#remove objects ------
rm(metadata, phy, phy_sub, phy_sub.1, phy_sub.2, phy_sub.bkg, phy_sub.skin, phy_sub.mpelung, phy_sub.mpeother, phy_sub.meso, counts.raw, metadata.raw, counts, count_table, p_bar_all, p_bar_bkg, p_bar_skin, p_bar_mpelung, p_bar_mpeother, p_bar_meso, p_ggarrange, countsum_table, OTU.table, Genus.table, OTU.rel.table, Genus.rel.table) 


# build phyloseq object------
phy<-qza_to_phyloseq(features="no-miss-table-dada2.qza", 
                     tree="rooted-tree_quality.qza", 
                     taxonomy="taxonomy.qza",
                     metadata="220802_pleural_fluid_metadata.txt")

#Give a colnames to separate different taxonomic levels
colnames(tax_table(phy))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Remove taxa with 0 abundance
phy_filter = subset_taxa(phy, rowSums(otu_table(phy)) != 0)
# otu_table()   OTU Table:         [ 11250 taxa and 244 samples ] -> otu_table()   OTU Table:         [ 8507 taxa and 244 samples ]

# Remove samples with read counts < 1000 
phy_filter = prune_samples(sample_sums(phy_filter)>=1000, phy_filter)

phy_filter_rel = transformSampleCounts(phy_filter, normalizeSample)

# Tax glom 
Genus.table = tax_glom(phy_filter, taxrank = "Genus", NArm = FALSE, bad_empty=c(NA, "", " ", "\t"))
# otu_table()   OTU Table:         [ 862 taxa and 244 samples ]


# filter and remove taxa not seen more than *5* times in at least *3%* of the samples. This protects against an OTU with small mean & trivially large C.V. ----
phy_sub_a.2 = genefilter_sample(Genus.table, filterfun_sample(function(x) x > 5), A = 0.01 * nsamples(phy_filter))
Genus.table = prune_taxa(phy_sub_a.2, Genus.table)

Genus.rel.table = transformSampleCounts(Genus.table, normalizeSample)

ps <- Genus.table
ps.rel <- Genus.rel.table

#Beta diversity for BKG vs SKin vs Pleural fluid (pooled)------
Fluid.Bray.dist = phyloseq::distance(phy_filter_rel, method="bray")
Fluid.Bray.pco = dudi.pco(d = cailliez(Fluid.Bray.dist), scannf = FALSE, nf = 3)

#plot
s.class(Fluid.Bray.pco$li, interaction(sample_data(phy_filter_rel)$Description.1), col=c('grey40','grey60','dodgerblue'))

#stats - KW
adonis <- adonis(Fluid.Bray.dist ~ Description.1, data=data.frame(sample_data(phy_filter_rel)), permutations= 10000)


### 1. bkg vs skin------
# subset to have BKG and Skin only
a <- subset_samples(phy_filter_rel, Description.1 %in% c('BKG','Skin'))
# sample_data(a)$Description.1 # Levels: BKG Skin
Fluid.Bray.dist = phyloseq::distance(a, method="bray")
adonis.1 <- adonis(Fluid.Bray.dist ~ Description.1, data=data.frame(sample_data(a)), permutations=10000)

# This is repeated for each pairwise comparison

## Alpha diversity------
Shannon_diversity_fluid = vegan::diversity(otu_table(phy_filter_rel), index = "shannon", MARGIN = 2, base = exp(1))

boxplot(Shannon_diversity_fluid
        ~sample_data(phy_filter_rel)$Description.1,
        ylab ="Shannon diversity index",
        xlab ="",
        col=c('grey40','grey80','dodgerblue'),
        names=c('Background','Skin','Pleural Fluid'),
        par(cex.axis=1),
        outline=FALSE,las=1
)

stripchart(Shannon_diversity_fluid
           ~ sample_data(phy_filter_rel)$Description.1,
           vertical=TRUE,
           add=TRUE,
           method="jitter",
           pch=20,
           col="black")

alphadiversity <- cbind(data.frame(Shannon_diversity_fluid),sample_data(phy_filter_rel)$Description.1)
colnames(alphadiversity) <- c("shannon","category")


#Stats
kw <- kruskal.test(Shannon_diversity_fluid ~ sample_data(phy_filter_rel)$Description.1)
kw$p.value

# wilcox rank sum test
rm(wilcox)
if (kw$p.value <0.05){
  wilcox <- pairwise.wilcox.test(alphadiversity$shannon,
                                 g = alphadiversity$category,
                                 p.adjust.method = "BH")}
wilcox



#Beta diversity for Benign, Paramalignant, MPE-Lung, MPE-Other, Mesothelioma------
phy.rel.nobkg <- phy_filter_rel %>%
  subset_samples(Histology %in% c("Benign.Effusion","Paramalignant","MPE.Lung","MPE.Other","MESOTHELIOMA"))

#Calculate Bray distance
Fluid.Bray.dist = phyloseq::distance(phy.rel.nobkg, method="bray")
Fluid.Bray.pco = dudi.pco(d = cailliez(Fluid.Bray.dist), scannf = FALSE, nf = 3)

# plot
s.class(Fluid.Bray.pco$li, interaction(sample_data(phy.rel.nobkg)$Histology), col=c('navyblue','khaki4','red','dodgerblue','orange'))

# stats
adonis <- adonis(Fluid.Bray.dist ~ Histology, data=data.frame(sample_data(phy.rel.nobkg)), permutations= 10000)

## doing adonis for each group------
# 1. mpelung_vs_mpeother
# 2. mpelung_vs_meso
# 3. mpeother_vs_meso
# 4. benign_vs_mpelung
# 5. benign_vs_mpeother
# 6. benign_vs_meso
# 7. cytology.negative_vs_mpelung
# 8. cytology.negative_vs_mpeother
# 9. cytology.negative_vs_mesothelioma
# 10. cytology.negative_vs_benign


# extract p-values with # adonis.2$aov.tab$'Pr(>F)'[1] or # adonis.2[["aov.tab"]][["Pr(>F)"]][1]
### 1. mpelung vs mpeother------
a <- subset_samples(phy.rel.nobkg, Histology %in% c('MPE.Lung','MPE.Other'))
Fluid.Bray.dist = phyloseq::distance(a, method="bray")
adonis.1 <- adonis(Fluid.Bray.dist ~ Histology, data=data.frame(sample_data(a)), permutations=10000)

# This is repeated for each pairwise comparison




## Alpha diversity------
Shannon_diversity_fluid = vegan::diversity(otu_table(phy.rel.nobkg), index = "shannon", MARGIN = 2, base = exp(1))

boxplot(Shannon_diversity_fluid
        ~sample_data(phy.rel.nobkg)$Histology,
        ylab ="Shannon diversity index",
        xlab ="",
        col=c('darkmagenta','khaki4','red','dodgerblue','orange'),
        names=c('Benign','Paramalignant','MPE-Lung','MPE-Other','Mesothelioma'),
        par(cex.axis=1),
        outline=FALSE,las=1
)

stripchart(Shannon_diversity_fluid
           ~ sample_data(phy.rel.nobkg)$Histology,
           vertical=TRUE,
           add=TRUE,
           method="jitter",
           pch=20,
           col="black")


alphadiversity <- cbind(data.frame(Shannon_diversity_fluid),sample_data(phy.rel.nobkg)$Histology)
colnames(alphadiversity) <- c("shannon","category")

#Stats
kw <- kruskal.test(Shannon_diversity_fluid ~ sample_data(phy.rel.nobkg)$Histology)
kw$p.value

# wilcoxon
rm(wilcox)
if (kw$p.value <0.05){
  wilcox <- pairwise.wilcox.test(alphadiversity$shannon,
                                 g = alphadiversity$category,
                                 p.adjust.method = "BH")}

wilcox$p.value

# //////////////////////////////////
# DECONTAM --------- 
# //////////////////////////////////
# Load/Save/wd
setwd() # set your own wd
getwd()
# save.image(file="220802_pleural_fluid/220802_pleural_fluid.RData")
# load(file="220802_pleural_fluid/220802_pleural_fluid.RData")

#creating factor for Histology.1 and Description.2 (surgical.tool)
sample_data(Genus.table)$Histology.1 <- factor(sample_data(Genus.table)$Histology.1, level = c("BKG","Surgical.Tool","Skin","Benign.Effusion","Paramalignant","MPE.Lung","MPE.Other","MESOTHELIOMA"))

sample_data(Genus.table)$Description.2 <- factor(sample_data(Genus.table)$Description.2, level = c("BKG","Surgical.Tool","Skin","Pleural.Fluid"))


# set threshold ------
threshold = 0.75

# set taxa name------
colnames(tax_table(Genus.table))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")


# subset into comparison groups ------
#BKG vs Fluid (MPE + Benign)
ps.bkg.fluid <- subset_samples(Genus.table, Histology.1 %in% c('BKG', 'Benign.Effusion', 'Paramalignant', 'MPE.Lung', 'MPE.Other', 'MESOTHELIOMA'))
ps.rel.bkg.fluid = transformSampleCounts(ps.bkg.fluid, normalizeSample)

#tools vs pleural fluid
ps.tools.fluid <- subset_samples(Genus.table, Histology.1 %in% c('Surgical.Tool', 'Benign.Effusion', 'Paramalignant', 'MPE.Lung', 'MPE.Other', 'MESOTHELIOMA'))
ps.rel.tools.fluid = transformSampleCounts(ps.tools.fluid, normalizeSample)


# BKG vs Fluid ------
counts.raw <- as.data.frame(otu_table(ps.bkg.fluid))
relab.raw <- as.data.frame(otu_table(ps.rel.bkg.fluid))
metadata.raw <- as.data.frame(sample_data(ps.bkg.fluid))
taxa.table <- as.data.frame(tax_table(ps.bkg.fluid))

#creates column "names" that does not have trailing ASV
taxa.table$names <- paste(ifelse(!is.na(taxa.table$OTU), paste("g_",taxa.table$Genus,"_s_",taxa.table$OTU,sep=""),
                                 ifelse(!is.na(taxa.table$Genus), paste("g_",taxa.table$Genus,sep=""),
                                        ifelse(!is.na(taxa.table$Family), paste("f_",taxa.table$Family,"_g__",sep=""),
                                               ifelse(!is.na(taxa.table$Order), paste("o_",taxa.table$Order, "_f__g__",sep=""),
                                                      ifelse(!is.na(taxa.table$Class), paste("c_",taxa.table$Class, "_o__f__g__",sep=""),
                                                             ifelse(!is.na(taxa.table$Phylum), paste("p_",taxa.table$Phylum, "_c__o__f__g__",sep=""),
                                                                    ifelse(!is.na(taxa.table$Domain), paste("k_",taxa.table$Domain, "_p__c__o__f__g__",sep=""), paste(rownames(taxa.table))))))))))

#creates column "nameASV" that includes trailing ASV
taxa.table$nameASV <- paste(ifelse(!is.na(taxa.table$OTU), paste("g_",taxa.table$Genus,"_s_",rownames(taxa.table),sep=""),
                                   ifelse(!is.na(taxa.table$Genus), paste("g_",taxa.table$Genus,"_",rownames(taxa.table),sep=""),
                                          ifelse(!is.na(taxa.table$Family), paste("f_",taxa.table$Family,"_g__",rownames(taxa.table),sep=""),
                                                 ifelse(!is.na(taxa.table$Order), paste("o_",taxa.table$Order, "_f__g__",rownames(taxa.table),sep=""),
                                                        ifelse(!is.na(taxa.table$Class), paste("c_",taxa.table$Class, "_o__f__g__",rownames(taxa.table),sep=""),
                                                               ifelse(!is.na(taxa.table$Phylum), paste("p_",taxa.table$Phylum, "_c__o__f__g__",rownames(taxa.table),sep=""),
                                                                      ifelse(!is.na(taxa.table$Domain), paste("k_",taxa.table$Domain, "_p__c__o__f__g__",rownames(taxa.table),sep=""), paste(rownames(taxa.table))))))))))

#create column ASV that is from row name
taxa.table$ASV <- rownames(taxa.table)


# names(taxa.table)

# this part of the code "removes" duplicated "names" from the "taxa.table" by adding a trailing ASV to "names." Ie, substitute "names" with "nameASV" if "names" repeats
duplicated(taxa.table$names) | duplicated(taxa.table$names, fromLast = TRUE) #this returns a logic vector the rows that are duplicated
any(duplicated(taxa.table$names) | duplicated(taxa.table$names, fromLast = TRUE))
dup <- subset(taxa.table, duplicated(taxa.table$names) | duplicated(taxa.table$names, fromLast = TRUE)) #this returns a df with rows that are duplicated

taxa.table[rownames(dup),]
rownames(taxa.table[rownames(dup),])
taxa.table[rownames(dup),]$nameASV
taxa.table[rownames(dup),]$names
# taxa.table[rownames(dup),]$names <- taxa.table[rownames(dup),]$nameASV #rename the repeated "names" to include trailing ASV values
taxa.table[rownames(dup),]$names <- paste("f_",taxa.table[rownames(dup),]$Family,"_",taxa.table[rownames(dup),]$names,sep="") # rename the repeated "names" to include one "order" above

duplicated(taxa.table$names) | duplicated(taxa.table$names, fromLast = TRUE) #check to see if any duplicates
any(duplicated(taxa.table$names)) #no longer any duplicated names


# rename sampleID with histology/category
# will use Description.2
counts <- counts.raw
relab <- relab.raw
names(counts) <- metadata.raw$Histology.1[match(names(counts.raw), metadata.raw$SampleID_1)]

names(relab) <- metadata.raw$Histology.1[match(names(relab.raw), metadata.raw$SampleID_1)]

counts <- merge(x = counts, y = taxa.table, by=0)
row.names(counts) <- counts$names 
counts <- counts[,-1]
counts <- subset(counts, select=-c(Domain, Phylum, Class, Order, Family, Genus, OTU, names, nameASV, ASV))

relab <- merge(x = relab, y = taxa.table, by=0)
row.names(relab) <- relab$names
relab <- relab [,-1]
relab <- subset(relab, select=-c(Domain, Phylum, Class, Order, Family, Genus, OTU, names, nameASV, ASV))

contam <- vector(mode = "list", length = 1)
contam[[1]]

comparison.list <- vector(mode = "list", length = 1)
comparison.list[[1]]

comparison.list[[1]] <- list() # clearing previous list
comparison.list[[1]]$all <- list(grep("Benign.Effusion|Paramalignant|MPE.Lung|MPE.Other|MESOTHELIOMA",colnames(counts)), grep("BKG",colnames(counts)))

contam <- list()
contam[[1]] <- list()
counts <- counts[,unlist(comparison.list[[1]])] # head(counts.subset) #this reorders the counts table so that BKG is at the end, pleural fluid in the front
relab <- relab[,unlist(comparison.list[[1]])] # head(relab.subset) #this reorders the counts table so that BKG is at the end, pleural fluid in the front
libsizes <- colSums(counts)


comparison.list[[1]]$all[[2]] <- ifelse(comparison.list[[1]]$all[[2]], T, F)  #change the $ to the comparison you're making!
#$Meso.BKG[[2]] = BKG, so we are labeling this as TRUE
comparison.list[[1]]$all[[1]] <- ifelse(comparison.list[[1]]$all[[1]], F, T)  #change the $ to the comparison you're making!
#$XXX.BKG[[1]] = XXX, so we are going to labbel this as FALSE (ie, FALSE, this is NOT a negative/background sample)


# extract from comparison.list to do decontam
bkg.samples <- unlist(comparison.list)

contam[[1]]$prev <- isContaminant(t(as.matrix(counts)),
                                  method="prevalence",
                                  neg=bkg.samples,
                                  threshold=threshold)

plot.df.bkg.fluid <- as.data.frame(rbind(cbind("decontam prevalence",rownames(contam[[1]]$prev)[contam[[1]]$prev$contaminant == T])))
plot.df.bkg.fluid <- plot.df.bkg.fluid[plot.df.bkg.fluid[,2] %in% taxa.table$names,]


# Tools vs Fluid ------
# The same is repeated for tools vs pleural fluid (pooled)


# merging output------
plot.df.bkg.fluid
plot.df.tools.fluid

plot.df.all <- distinct(rbind(plot.df.bkg.fluid,
                              plot.df.tools.fluid))

# plot.df.all is then merged with median relative abundances for each sample type (BKG+Tools, Skin, Pleural Fluid)

#BKG as reference for rank order------
# plot.df.all is arranged by BKG samples, then plotted

# MPE as reference for rank order-----
# repeat the above code, objects are below
# arrange and save plots below

# Skin as reference for rank order-----
# repeat the above code, objects are below
# arrange and save plots below


# write plot.df.all table -----
write.table(plot.df.all,file= paste("220802_pleural_fluid/decontam",threshold,"_plot.df.all_v13.txt",sep=""), sep="\t", col.names = NA, row.names = TRUE)


# //////////////////////////////////
# LEfSE -----
# //////////////////////////////////

# Load/Save/wd
setwd() # set your own wd
save.image(file="220802_pleural_fluid/220802_pleural_fluid.RData")
load(file="220802_pleural_fluid/220802_pleural_fluid.RData")

Genus.table
Genus.rel.table

Genus.table.noBKG <- subset_samples(Genus.table, Histology  %!in% c("BKG", "Skin","Surgical.Tool"))

Genus.rel.table.noBKG = transformSampleCounts(Genus.table.noBKG, normalizeSample)

rm(res, resCladogram, res1, resdf, merged, resc1, resc2, otuC1, otuC2, rel.otuC1, rel.otuC2, otu.to.save1, otu.to.save2, df.1.df, df.2.df, df.1.meanRA, df.2.meanRA, df.1.meanRA.save, df.2.meanRA.save, RESDF, data.table, new, test) 

## set threshold ----
pvalue.thres <- 0.05
lda.thres <- 2

ps <- Genus.table.noBKG 
colnames(tax_table(ps))=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(ps))

ps.rel <- transformSampleCounts(ps, normalizeSample)

#### Using microbiomeMarker package-----
set.seed(123) 
res <- run_lefse(ps,
                 group="Histology",
                 subgroup = NULL,
                 taxa_rank = "none", # this will give ASV values. Cannot use this with Cladogram
                 # needs to be "all" to run Cladogram (needs names listed as f__, g__, s__)
                 # needs to be "all"none" to do Bubble Plot (needs ASV) -- will somehow need to get ASV changed to names
                 transform = "identity", #log10, log10p
                 norm = "CPM", #CPM is method of choice for normalizing data for LEFSE analysis
                 kw_cutoff = 0.05,
                 wilcoxon_cutoff = 0.05,
                 lda_cutoff = 2,
                 bootstrap_n = 30, #30 is default
                 bootstrap_fraction = 2/3,
                 multigrp_strat = TRUE
)


res1 <- data.frame(marker_table(res))


# //////////////////////////////////
# Kaplan-Meier Curves-----
# //////////////////////////////////


# Load/Save/wd
setwd() # set your own wd
save.image(file="220802_pleural_fluid/220802_pleural_fluid.RData")
load(file="220802_pleural_fluid/220802_pleural_fluid.RData")


Genus.table.noBKG <- subset_samples(Genus.table, Histology %in% c("MPE.Lung","MPE.Other","MESOTHELIOMA"))
Genus.table.noBKG <- subset_samples(Genus.table.noBKG, unique_surv == 1) # unique_surv not unique
metadata1 <- data.frame(sample_data(Genus.table.noBKG))
metadata1 <- metadata1[, c("alive_or_dead", 
                           "months_from_dx_to_death",
                           "Histology")]
rownames(metadata1) <- NULL
colnames(metadata1) <- c('status','time','histo')
metadata1[metadata1$status == "n.a" & metadata1$time == "n.a",]
metadata1$time <- as.numeric(metadata1$time)
metadata1$status <- as.numeric(metadata1$status)
metadata1$histo <- factor(metadata1$histo, levels = c("MPE.Lung","MPE.Other","MESOTHELIOMA"))


metadata1 <- as_tibble(metadata1)
metadata1 <- metadata1[complete.cases(metadata1),]

metadata1 <- metadata1 %>% dplyr::group_by(histo) %>%
  mutate(count = n()) %>%
  mutate(mean = mean(time)) %>%
  mutate(median = median(time)) %>%
  ungroup() %>%
  mutate(time1 = 
           ifelse(time > 36 & status == "0", 36, 
                  ifelse(time > 36 & status == "1", 36, time))) %>%
  mutate(status1 = ifelse(time > 36 & status == "1", 0, status))


#Fit survival curve
fit <- survfit(Surv(time1+0.01, status1)~ histo, data=metadata1)
fit
surv_summary(fit)
test <- surv_summary(fit)

#Draw survival curve
p_all <- ggsurvplot(fit, 
                    data = metadata1,
                    fun = "pct",
                    ggtheme = theme_pubr(),
                    xlab= "Time (months)",
                    # title = "Overall Survival",
                    # subtitle= "Based on Kaplan-Meier Estimates",
                    palette = c("red","dodgerblue","orange"),
                    size = 1,
                    # linetype = "strata",
                    # pval = TRUE,
                    pval.coord = c(0, 10),
                    # pval.method = TRUE,
                    pval.size = 4,
                    # test.for.trend = TRUE,
                    # surv.median.line = "hv",
                    legend.title = "Groups",
                    legend.labs = c("MPE-Lung","MPE-Other","Mesothelioma"),
                    legend = "top",
                    break.x.by = 6,
                    xlim = c(0,36),
                    axes.offset = T,
                    # conf.int = T,
                    
                    ###censor ###,
                    censor = TRUE,
                    censor.shape = 124,
                    
                    
                    
                    ###RISK TABLE ###
                    risk.table = "nrisk_cumevents",
                    # cumevents = TRUE,
                    cumcensor= TRUE,
                    # tables.col = "strata", #color of numbers inside risk table
                    tables.height = 0.15,
                    fontsize = 4,
                    risk.table.y.text.col = T,
                    tables.theme = theme_cleantable(),
                    # risk.table.y.text = FALSE,
                    risk.table.title = "Number at risk (cumulative number of deaths)",
                    cumcensor.title = "Cumulative number of censored subjects"
)
p_all

pdf(file="220802_pleural_fluid/KaplanMeier.pdf", width = 8, height = 8)
  p_all
dev.off()

# pairwise comparisons for p value
surv_pvalue(fit) # overall

coxfit <- coxph(
  Surv(time1+0.01, status1)~ histo,
  data=metadata1,
  ties = 'exact'
)
summary(coxfit) # MPE.Lung as reference

metadata1 %>% mutate(histo = factor(histo, levels = c("MPE.Other", "MESOTHELIOMA", "MPE.Lung"))) %>%
  coxph(
    Surv(time1+0.01, status1)~ histo,
    data=.,
    ties = 'exact'
  ) %>%
  summary(.) # MPE.Other as reference



# //////////////////////////////////
# Survival Analysis -----
# //////////////////////////////////


# Load/Save/wd-----
setwd() # set your own wd
save.image(file="220802_pleural_fluid/220802_pleural_fluid.RData")
load(file="220802_pleural_fluid/220802_pleural_fluid.RData")


# set color of x-axis----
plot.df.all
contamlist.norank <- plot.df.all$names


#MPE-Lung -----
#split Genus.table
Genus.table.mpelung <- Genus.table %>%
  subset_samples(., Histology %in% "MPE.Lung")

#remove samples with n.a follow-up
sample_data(Genus.table.mpelung)$months_from_dx_to_death
Genus.table.mpelung <- subset_samples(Genus.table.mpelung, months_from_dx_to_death %!in% "n.a")
sample_data(Genus.table.mpelung)$months_from_dx_to_death <- as.numeric(sample_data(Genus.table.mpelung)$months_from_dx_to_death)

metadata <- sample_data(Genus.table.mpelung)
TAX <- tax_table(Genus.table.mpelung)
OTU <- otu_table(Genus.table.mpelung)

# KP median is 33.6 (dx to death)

#create new column
metadata$Xmedian.survival.mpelung <- 0
metadata[metadata$months_from_dx_to_death < 33.6 & metadata$alive_or_dead == 1, ]$Xmedian.survival.mpelung <- "Dead"
metadata[metadata$months_from_dx_to_death >= 33.6, ]$Xmedian.survival.mpelung <- "Alive"
table(metadata$Xmedian.survival.mpelung)
metadata$Xmedian.survival.mpelung <- factor(metadata$Xmedian.survival.mpelung, levels = c("Dead","Alive"))

metadata <- metadata %>% data.frame () %>% select(Histology, alive_or_dead, months_from_dx_to_death, Xmedian.survival.mpelung) %>% .[complete.cases(.),]
data.frame(table(metadata$Xmedian.survival.mpelung))

#for 33.6 months
# Var1 Freq
# 1  Dead   19
# 2 Alive   15

#build ps object again
Genus.table.mpelung <- phyloseq(sample_data(metadata), otu_table(OTU), tax_table(TAX))

Genus.rel.table.mpelung = transformSampleCounts(Genus.table.mpelung, normalizeSample)


###Beta diversity ------
#Calculate Bray distance
Fluid.Bray.dist = phyloseq::distance(Genus.rel.table.mpelung, method="bray")
Fluid.Bray.pco = dudi.pco(d = cailliez(Fluid.Bray.dist), scannf = FALSE, nf = 3)

s.class(Fluid.Bray.pco$li, interaction(sample_data(Genus.rel.table.mpelung)$Xmedian.survival.mpelung), col=c('red','darkgreen'))
# dev.off()
adonis <- adonis(Fluid.Bray.dist ~ Xmedian.survival.mpelung, data=data.frame(sample_data(Genus.rel.table.mpelung)), permutations= 10000)

adonis[["aov.tab"]][["Pr(>F)"]][1]


## Alpha diversity------
Shannon_diversity_fluid = vegan::diversity(otu_table(Genus.rel.table.mpelung), index = "shannon", MARGIN = 2, base = exp(1))

boxplot(Shannon_diversity_fluid
        ~sample_data(Genus.rel.table.mpelung)$Xmedian.survival.mpelung,
        ylab ="Shannon diversity index",
        xlab ="",
        col=c('red','darkgreen'),
        names=c('Deceased','Alive'),
        par(cex.axis=0.50),
        outline=FALSE,las=1
)

stripchart(Shannon_diversity_fluid
           ~ sample_data(Genus.rel.table.mpelung)$Xmedian.survival.mpelung,
           vertical=TRUE,
           add=TRUE,
           method="jitter",
           pch=20,
           col="black")


alphadiversity <- cbind(data.frame(Shannon_diversity_fluid),sample_data(Genus.rel.table.mpelung)$Xmedian.survival.mpelung)
colnames(alphadiversity) <- c("shannon","category")

#Stats
kw <- kruskal.test(Shannon_diversity_fluid ~ sample_data(Genus.rel.table.mpelung)$Xmedian.survival.mpelung)
kw$p.value



##LEFSE -------
## set threshold ----
pvalue.thres
lda.thres
## MPE Lung - dead vs alive
Genus.table.mpelung
Genus.rel.table.mpelung


ps <- Genus.table.mpelung 
colnames(tax_table(ps))=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(ps))
# view(tax_table(ps))
ps.rel <- transformSampleCounts(ps, normalizeSample)
# view(otu_table(ps.rel))

#Prune Data to whatever dimensions you need
# OTU.wh1 = genefilter_sample(ps, filterfun_sample(function(x) x > 5), A = 0.03 * nsamples(ps))
# ps = prune_taxa(OTU.wh1, ps)
# ps.rel = transformSampleCounts(ps, normalizeSample)
# 
# sample_data(ps)$Description.1 #MPE.Lung, MESOTHELIOMA

#### Using microbiomeMarker package-----
set.seed(123) 
res <- run_lefse(ps,
                 group="Xmedian.survival.mpelung",
                 subgroup = NULL,
                 taxa_rank = "none", # this will give ASV values. Cannot use this with Cladogram
                 # needs to be "all" to run Cladogram (needs names listed as f__, g__, s__)
                 # needs to be "all"none" to do Bubble Plot (needs ASV) -- will somehow need to get ASV changed to names
                 transform = "identity", #log10, log10p
                 norm = "CPM", #CPM is method of choice for normalizing data for LEFSE analysis
                 kw_cutoff = 0.05,
                 wilcoxon_cutoff = 0.05,
                 lda_cutoff = 2,
                 bootstrap_n = 30, #30 is default
                 bootstrap_fraction = 2/3,
                 multigrp_strat = TRUE
)


res1 <- data.frame(marker_table(res))

# repeat abive section for MPE-Other and Mesothelioma using cutoffs for deceased/alive from Kaplan-Meier survival model:
## MPE-Other 36 months
## Mesothelioma 19.6 months