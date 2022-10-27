#################################
# Analysis of 16S amplicon data #
#################################

Sys.setenv(LANG = "en") # in case system language changes to German

######              Parameters             #########


# >>> alpha diversity (univariate)

### Inverse Simpson index (also can be used to interpret eveness)

### Richness (number of ASVs)



# >>>beta diversity (multivariate)

### Bray-Curtis dissimilarity (asymmetric)

### Euclidean distance (symmetric)
# consider common 0 as similarities

### Jaccard dissimilarity (asymmetric)
# presence absence
# more conservative
# does not take into account the sequence proportions (can be good for Euks )



# what to do with these indexes
### Mantel test (multivariate with different tables)
# correlate one measure of dissimilarity of one table to another measure of another table





#              Visualization            #


### rarefaction curves 
#(to note if the seq depth must be taking in to account to interpret alpha diversity results)
#     inv simpson  
#     richness


# non metric scaling methods

### NMDS (ordination - simplify multi dimensional data)
#     BC dissimilarities

# metric scaling methods


### PCA ordination 
#     Euclidean 

### PCoA
#     NOT RESTRICTED

### beta dispersion
# calculated based in the PCoA
# heterogeneity of the community

### clusters
# hierarchical clustering
# different linkage (single, complete)





#              Comparisons              #

# start by setting colors in the visualization 

# later I can try with more statistical approach (ANOVA for example)

############



### Start analyses

### prepare environment ####

# set working directory
setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_pros")

# load packages
require(vegan)
#BiocManager::install("ALDEx2")
require(ALDEx2)
require(reshape)
require(gplots)
require(tidyverse)
#install.packages("ecodist") # for PCoA analysis
require(ecodist)

# loading and saving workspace
#load("data_exploration_all_pros.Rdata")
#save.image("data_exploration_all_pros.Rdata")

#####


### read and adjust data ####

# ASV table
ASV <- read.table(
  "otu_tab_ssu_all_pros_curated.txt", 
  h = T,
  sep = "\t",
  row.names = 1
)
# adjusting Purified_PCR_product from _ to - 
colnames(ASV) <- gsub(".", "-", colnames(ASV), fixed = TRUE)

# Metadata
META <- readxl::read_xlsx("Lists_bioinformatic_all_pros.xlsx") 
META <- META[META$Purified_PCR_product != "PR-1178" & META$Purified_PCR_product != "PR-0703",] 
# removing sample that was excluded in data_curation script (no sequences left after filtering) - PR-1178
# removing sample that was identified as a Failed PCR (detailed later in this script) - PR-0703


# removing ASVs of this sample from ASV table
META <- META %>% column_to_rownames("Purified_PCR_product") # assign PCR product as row name of META 
ASV <- ASV[, rownames(META)]

# taxonomy of microbial community
TAX <- read.table(
  "tax_tab_ssu_all_pros_curated_2.txt", 
  h = T, 
  sep = "\t", 
  row.names = 1,
  comment.char = "",
  quote = "", 
  stringsAsFactors = F
)

# sub-setting TAX according to ASV table filtered above
TAX <- TAX[rownames(ASV),]

# check data types
str(ASV)
str(TAX)
str(META)


# adjust object and data types
ASV <- as.matrix(ASV)
TAX <- as.matrix(TAX)
META$Site <- factor(META$Site, levels = c("Capo Passero", "Plemmirio", "Tel Shikmona", "Eilat", "NC"))
META$Condition <- factor(META$Condition, levels = c("Shallow_Algae", "Shallow_Rubbles", "Deep_Algae", "Deep_Rubbles", "NC"))
META$Sample_type <- factor(META$Sample_type, levels = c("filter", "sediment", "single_cell", "NC_env", "NC_single_cell"))
META$Sample_NC <- as.factor(META$Sample_NC)
levels(META$Sample_type)


# reorder data
META <- META[order(META$Site, META$Condition), ]
ASV <- ASV[, rownames(META)]

# check that table are in the correct order
all.equal(rownames(ASV), rownames(TAX)) # TRUE  
all.equal(colnames(ASV), rownames(META)) # TRUE

# setting color scheme for sites
site.color <- c("darkolivegreen4", "darkorchid1", "gold1", "firebrick", "black")
META$color_site <- META$Site
levels(META$color_site) <- site.color
META$color_site <- as.character(META$color_site)

# setting color scheme for sample type
sample.type.color <- c("dodgerblue", "blue", "red", "azure4", "dimgrey")
META$color_sample_type <- META$Sample_type
levels(META$color_sample_type) <- sample.type.color
META$color_sample_type <- as.character(META$color_sample_type)
str(META)


# calculate proportions
ASV.rel <- prop.table(ASV, 2) * 100

#######


##### define ASV.filt cut off of sequences (according to data exploration in "screen_nc_pros.R" script) #####
## Given the high dissimilarity between the technical replicates in this data set, we cannot apply the same protocol 
## as we did in the Euks data set (keeping only ASVs present in both replicates), because this caused a major modification in the
## community diversity. 
## Therefore, we should define a different threshold, but that still have a conservative approach in order to avoid false positives

## Attempt 1: Remove ASVs that don't occur with at least 0.1% in 3 replicates
# ASV.filt <- ASV[apply(ASV.rel, 1, function(x) sum(x >= 0.1) >= 3), ]
# # comparing number of ASVs kept
# dim(ASV.filt)/dim(ASV)
# #0.07589587 - only 7.5% kept

## Attempt 2: Remove ASVs that don't occur with at least 0.01% in 2 replicates
ASV.filt <- ASV[apply(ASV.rel, 1, function(x) sum(x >= 0.01) >= 2), ]
# comparing number of ASVs kept
summary(dim(ASV.filt)/dim(ASV))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4349  0.5761  0.7174  0.7174  0.8587  1.0000 

## Attempt 3: Remove ASVs that do not occur in both technical replicates OR in both replicates of another sample (of the same sample)
# presence absence approach

# create dataset where ASVs are only kept if present in both technical replicates of the same sample
# to use as one of the conditions to look at later

# let stick with Attempt 2!

#### merging ___ does not work after removing sample PR-0703 ######
# METAnoNC <- META[META$Sample_NC == "Sample",] # working only with samples 
# ASVnoNC <- ASV.filt[, rownames(METAnoNC)]
# 
# ASV.strict <- map_dfr(
#   unique(METAnoNC$Extraction_Voucher),
#   function(i) { 
#     ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i]
#     ASV.rep <- rowSums(ASV.rep[apply(ASV.rep, 1, min) > 0, ])
#     return(ASV.rep)
#   }
# ) %>% t()
# colnames(ASV.strict) <- unique(METAnoNC$Extraction_Voucher)
# ASV.strict[is.na(ASV.strict)] <- 0
# 
# ASV.merge <- map_dfr(
#   unique(METAnoNC$Extraction_Voucher),
#   function(i) { 
#     ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i]
#     ASV.rep <- ASV.rep[rowSums(ASV.rep) > 0, ]
#     tmp <- ASV.strict[, colnames(ASV.strict) != i]
#     tmp <- tmp[rowSums(tmp) > 0, ]
#     ASV.rep <- rowSums(ASV.rep[rownames(ASV.rep) %in% rownames(tmp) | apply(ASV.rep, 1, min) > 0, ])
#     return(ASV.rep)
#   }
# ) %>% t()
# colnames(ASV.merge) <- unique(METAnoNC$Extraction_Voucher)
# ASV.merge[is.na(ASV.merge)] <- 0
# 
# summary(colSums(ASV.merge)/c(by(colSums(ASVnoNC), METAnoNC$Extraction_Voucher, sum))[colnames(ASV.merge)])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3433  0.8756  0.9292  0.9118  0.9645  0.9941 

#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. # using ASV.filt to obtain ASVnoNC
# 0.4636  0.9075  0.9508  0.9373  0.9886  0.9999

# best attempt so far...
######

#### seeing changes in betadiversity on the PCR-replicate level#####
ASV.check <- map(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
    ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i]
    tmp <- ASV.strict[, colnames(ASV.strict) != i]
    tmp <- tmp[rowSums(tmp) > 0, ]
    ASV.rep[!rownames(ASV.rep) %in% rownames(tmp) & apply(ASV.rep, 1, min) == 0, ] <- 0
    return(ASV.rep)
  }
)
ASV.check_df <- do.call("cbind", ASV.check)
ASV.check_df <- ASV.check_df[, colnames(ASVnoNC)]
ASV.check_df <- ASV.check_df[rowSums(ASV.check_df) > 0, ]
ASV.check.rel <- prop.table(ASV.check_df, 2) * 100
summary(colSums(ASV.check_df)/colSums(ASVnoNC))

ASVnoNC.rel <- prop.table(ASVnoNC, 2) * 100

plot(
  vegdist(t(ASVnoNC.rel)),
  vegdist(t(ASV.check.rel))
)
abline(0, 1)

# much, much better
# there are however still some outliers,
# trying to identify them...

d1 <- melt(as.matrix(vegdist(t(ASVnoNC.rel))))
d2 <- melt(as.matrix(vegdist(t(ASV.check.rel))))

pw.combn <- combn(colnames(ASVnoNC.rel), 2)
d1sub <- map_dfr(1:ncol(pw.combn), function(x) d1[d1$X1 == pw.combn[1, x] & d1$X2 == pw.combn[2, x], ]) # does not work after removing the sample PR0703 from the data set
d2sub <- map_dfr(1:ncol(pw.combn), function(x) d2[d2$X1 == pw.combn[1, x] & d2$X2 == pw.combn[2, x], ]) # does not work after removing the sample PR0703 from the data set

dev.off()
plot(d1sub$value, d2sub$value)
abline(v = 0.95, h = 0.8)
tmp <- data.frame(d1sub, d2sub$value)
tmp.sub <- tmp[tmp$value > 0.95 & tmp$d2sub.value > 0.8, ]
tmp.sub$diff <- tmp.sub$value - tmp.sub$d2sub.value
plot(tmp.sub$value, tmp.sub$d2sub.value, col = (tmp.sub$X1 == "PR-0703" | tmp.sub$X2 == "PR-0703") + 1, pch = 16, cex = 0.5)
# those are all associated with PR-0703
METAnoNC["PR-0703", ]
METAnoNC[METAnoNC$Extraction_Voucher == "S0512", ]
colSums(ASVnoNC)[c("PR-0703", "PR-0704")]
plot(tmp$value, tmp$d2sub.value, col = (tmp$X1 == "PR-0703" | tmp$X2 == "PR-0703") + 1, pch = 16, cex = 0.5)

# remove this technical replicate (PR-0703) from dataset (in the beginning of the script) 
# and come back to analysis, do the filtering again (attempt 3), but applying drop equals false,
# otherwise code wont work because there is now one sample with only one replicate and it 
# cannot calculate a triangular matrix, with "drop = F" we allow to construct a matrix with only one column

#######

##### ASV.merge table to merge ASVs from same extraction voucher based in ASV.filt ######
METAnoNC <- META[META$Sample_NC == "Sample",] # working only with samples 
ASVnoNC <- ASV.filt[, rownames(METAnoNC)]

# create dataset where ASVs are only kept if present in both technical replicates of the same sample
# to use as one of the conditions to look at later
ASV.strict <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
    ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i, drop = F] # difference from the code before is the "drop = F"
    ASV.rep <- rowSums(ASV.rep[apply(ASV.rep, 1, min) > 0, , drop = F]) # difference from the code before is the "drop = F"
    return(ASV.rep)
  }
) %>% t()
colnames(ASV.strict) <- unique(METAnoNC$Extraction_Voucher)
ASV.strict[is.na(ASV.strict)] <- 0

ASV.merge <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
    ASV.rep <- ASVnoNC[, METAnoNC$Extraction_Voucher == i, drop = F]
    ASV.rep <- ASV.rep[rowSums(ASV.rep) > 0, , drop = F]
    tmp <- ASV.strict[, colnames(ASV.strict) != i]
    tmp <- tmp[rowSums(tmp) > 0, ]
    ASV.rep <- rowSums(ASV.rep[rownames(ASV.rep) %in% rownames(tmp) | apply(ASV.rep, 1, min) > 0, , drop = F])
    return(ASV.rep)
  }
) %>% t()
colnames(ASV.merge) <- unique(METAnoNC$Extraction_Voucher)
ASV.merge[is.na(ASV.merge)] <- 0

summary(colSums(ASV.merge)/c(by(colSums(ASVnoNC), METAnoNC$Extraction_Voucher, sum))[colnames(ASV.merge)])
# Min.   1st Qu.  Median   Mean  3rd Qu.   Max. 
# 0.6060  0.8791  0.9335  0.9170  0.9655  1.0000 

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.        # using ASV.filt to obtain ASVnoNC
# 0.6562  0.9092  0.9524  0.9422  0.9896  1.0000 

str(ASV.merge)
# num [1:22352, 1:144] 239 25 13 85 1 7 8 2 3 5 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:22352] "sq6" "sq7" "sq9" "sq14" ...
# ..$ : chr [1:144] "S0566" "S0612" "S0613" "S0514" ...


# num [1:14228, 1:144] 239 25 13 85 1 7 8 2 3 5 ...    # using ASV.filt to obtain ASVnoNC
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14228] "sq6" "sq7" "sq9" "sq14" ...
# ..$ : chr [1:144] "S0566" "S0612" "S0613" "S0514" ...

# Creating META merge table to analyze together with ASV.merge
META.merge <- META[META$Sample_NC == "Sample",] # taking META without NCs
META.merge <- META.merge[!grepl("_1", META.merge$PCR_Product),] # removing the replicated lines for the technical replicates 
#(keeping the ones with "_2" because I excluded the "_1" of extraction voucher S0512 (PR-0703) in the beginning of the script)
# this is not the best code to solve this issue, but serves for my case.

rownames(META.merge) <- META.merge$Extraction_Voucher
META.merge$DNA_yield_mean <- c(NA) # create empty column to add data of DNA concentration mean between replicates
#META.merge$DNA_yield_dif <- c(NA) # create empty column to add data of DNA concentration difference between replicates
for (i in META.merge$Extraction_Voucher){ # selecting sample by sample in a loop
  META.merge[i,"DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) # calculating the mean of DNA conc among the technical replicates
  #META.merge[i,"DNA_yield_dif"] <- max(c(dist(META$DNA_concentration[META$Extraction_Voucher == i]))) # calculating the difference between the DNA conc of the technical replicates
  #using max to make a more universal code in case one has more than two PCR replicates
}

str(META.merge)
META.merge$Sample_type <- droplevels(META.merge$Sample_type)
META.merge$Site <- droplevels(META.merge$Site)
META.merge$Condition <- droplevels(META.merge$Condition)

ASV.merge <- ASV.merge[, rownames(META.merge)] # same order than META
all.equal(colnames(ASV.merge), rownames(META.merge))
# [1] TRUE

TAX.merge <- TAX[rownames(ASV.merge),]
# confirming TAX.merge table still matches with the ASV.merge table
all.equal(rownames(ASV.merge), rownames(TAX.merge))
# [1] TRUE

ASV.rel.merge <- prop.table(ASV.merge, 2) * 100
#######

##### saving RDS files to use in next script "pattern_exploration_pros.R"#####

saveRDS(ASV.merge, "ASV_merge_pros.RDS")
saveRDS(TAX.merge, "TAX_merge_pros.RDS")
saveRDS(META.merge, "META_merge_pros.RDS")


#######


### OLD: rearrange data set to combine the technical replicates - important step before moving to the statistical analysis #####
# METAnoNC <- META[META$Sample_NC == "Sample",] # working only with samples before I decide what to do with the NCs
# 
# 
# # using the function inside map_dfr from dplyr
# require(dplyr)
# ASV.merge <- map_dfr(
#   unique(METAnoNC$Extraction_Voucher),
#   function(i) { 
#     META.rep <- METAnoNC[METAnoNC$Extraction_Voucher == i, ]
#     ASV.rep <- ASV.filt[, rownames(META.rep)] # selecting ASVs of the technical replicates from the same sample
#     temp <- ASV.rep # keeping all ASVs present in both replicates
#     ASV.rep.temp <- as.data.frame(rowSums(temp)) # summing the replicates in a new data frame
#     colnames(ASV.rep.temp) <- i  #saving column name as the extraction voucher
#     data.frame(t(ASV.rep.temp))
#     
#   })
# ASV.merge <- t(ASV.merge) # transpose again to have in the same format as the original ASV table
# ASV.merge[is.na(ASV.merge)] <- 0 # setting NAs as zero


# Observing how this cut off affected the BC dissimilarities (comparing original vs filtered data sets)

# creating a ASV.merge.original which is based in the original ASV before applying the cut off
ASV.merge.original <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) {
    META.rep <- METAnoNC[METAnoNC$Extraction_Voucher == i,]
    ASV.rep <- ASV[, rownames(META.rep), drop = F] # selecting ASVs of the technical replicates from the same sample and excluding NCs
    temp <- ASV.rep # keeping all ASVs present in both replicates
    ASV.rep.temp <- as.data.frame(rowSums(temp), drop = F) # summing the replicates in a new data frame
    colnames(ASV.rep.temp) <- i  #saving column name as the extraction voucher
    data.frame(t(ASV.rep.temp))

  })
ASV.merge.original <- t(ASV.merge.original) # transpose again to have in the same format as the original ASV table
ASV.merge.original[is.na(ASV.merge.original)] <- 0 # setting NAs as zero
ASV.rel.merge.original  <- prop.table(ASV.merge.original, 2) * 100
#
# # calculate Bray-Curtis dissimilarities
BC.merge.filtered <- vegdist(t(prop.table(ASV.merge, 2) * 100))
BC.merge.original <- vegdist(t(prop.table(ASV.merge.original, 2) * 100))
#
plot(BC.merge.original, BC.merge.filtered)
abline (0, 1)

# visualizing h clusts
META.plot <- META.merge
BC.clust <- hclust(BC.merge.filtered, method = "average")
BC.plot <- plot(BC.clust, labels = paste(META.plot$Site, META.plot$Sample_type), cex = 0.7)


#############


####### Exploring similarity between single cell and environmental samples ######

# selecting ASVs of single cells only
ASV.merge.sing <- ASV.merge[,META.merge$Sample_type == "single_cell"] # select single cell
ASV.merge.sing <- ASV.merge.sing[rowSums(ASV.merge.sing) > 0,] # excluding ASVs that are not in the subset
META.sing <- META.merge[colnames(ASV.merge.sing),]
all.equal(rownames(META.sing), colnames(ASV.merge.sing))

ASV.rel.sing <- prop.table(ASV.merge.sing, 2) * 100

# selecting ASVs of env samples only
ASV.merge.env <- ASV.merge[,META.merge$Sample_type != "single_cell"] # select env samples
ASV.merge.env <- ASV.merge.env[rowSums(ASV.merge.env) > 0,] # excluding ASVs that are not in the subset
META.env <- META.merge[colnames(ASV.merge.env),]

ASV.rel.env <- prop.table(ASV.merge.env, 2) * 100
######

### NMDS plot ####

# calculate BC dissimilarities
BC <- vegdist(t(ASV.rel.merge)) # run ASV.rel above (either ASV.rel.sing, ASV.rel.env, or ASV.rel.merge - for all together)

# Plotting h clust
META.plot <- META.merge # select META table to work with (META.sing, META.env or META.merge, for instance)

BC.clust <- hclust(BC, method = "average")
BC.plot <- plot(BC.clust, labels = META.plot$Site, cex = 0.7)


# Define groups based on equal similarity
# Those will be used in the NMDS plot
rect.hclust(BC.clust, h = 0.8) 
# 0.8 was better for sing cell and for all together
BC.groups <- cutree(BC.clust, h = 0.8)


# calculate NMDS
NMDS <- metaMDS(t(ASV.rel.merge), k = 2)
NMDS$stress
# [1] 0.1071108    all samples
# [1]     single cell
# [1]     env samples all

# [1]     sediment

# [1]     seawater
# Warning message:
#   In metaMDS(t(ASV.rel), k = 2) :
#   stress is (nearly) zero: you may have insufficient data

# build nice plot element by element
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n")

# add hulls based on equal similarity
ordihull(NMDS, META.plot$Site)

# add points for sites
# sing
points(
  NMDS$points[, 1],
  NMDS$points[, 2],
  col = META.plot$color_site,
  pch = 16,
  cex = 1
)


# add legend for sites
legend(
  "bottom",
  legend = levels(META.plot$Site),
  col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
  pch = 15,
  pt.cex = 1,
  title = "Site"
)

########

### Betadispersion ####
betadispersion <- betadisper(
  vegdist(t(ASV.rel.merge)), 
  META.plot$Site,
  sqrt.dist = T
)

# show PCoA
plot(betadispersion)

# show within group heterogeneity
boxplot(betadispersion)

#######

### PERMANOVA ######
# change for the sub set dataset
META.adonis <- META.merge

str(META.adonis)
META.adonis$Sample_type <- droplevels(META.adonis$Sample_type)
META.adonis$Site <- droplevels(META.adonis$Site)
META.adonis$Condition <- droplevels(META.adonis$Condition)

# How much of the variability in the data can be explained by the environmental /PCR variables
adonis2(t(ASV.rel.merge) ~ Site + PCR + Depth + Sample_type, data = META.adonis, sqrt.dist = T)

# Let's also have a look at pairwise comparisons for sites
site.comb <- combn(levels(META.adonis$Site), 2) # creating a triangular pairwise comparison
permanova.list <-apply(
  site.comb,
  2,
  function(x) {
    adonis2(
      t(ASV.rel.sing[, rownames(META.adonis)[META.adonis$Site %in% x]]) ~ Site, 
      data = droplevels(META.adonis[META.adonis$Site %in% x, ]), 
      sqrt.dist = T
    )
  }
)

# reformat as table
permanova.df <- data.frame(
  t(site.comb),
  do.call(
    "rbind",
    lapply(
      permanova.list,
      function(x) {
        c(x$R2[1], x$'F'[1], x$'Pr(>F)'[1])
      }
    )
  )
)
colnames(permanova.df) <- c("cond1", "cond2", "R2", "F", "p")
# adjust p-values for multiple testing
permanova.df$p.adj <- p.adjust(permanova.df$p, method = "fdr")

permanova.df
# all sites are different!

########

### Venn Diagram ####

# Plotting each site the ASVs of single cell vs sediment vs filters
# to observe if there is an intake from the environmental biota

library("ggVennDiagram")
require(ggplot2)

## to see per specific site
s <- 'Eilat' # # select site to subset in s
ASV.venn.diagram <- ASV.merge[, META.merge$Site == s] # select site s
ASV.venn.diagram <- ASV.venn.diagram[rowSums(ASV.venn.diagram) > 0, ] # remove ASVs that are not present in any sample in this sub setting
METAs <- META.merge[META.merge$Site == s,] # select site s

# create list for venn diagram setting different thresholds (if wished) according to sample type
venn.list <- vector("list", length = 3)
names(venn.list) <- c("single cell", "sediment", "filter")
tmp.asv <- ASV.venn.diagram[, METAs$Sample_type == "single_cell"]
tmp.rel <- prop.table(tmp.asv, 2) * 100
venn.list$`single cell` <- rownames(tmp.rel)[apply(tmp.rel, 1, function(x) sum(x > 0) >= 1)] # this option takes all ASVs (larger than 0 and occurs in at least 1 sample)
tmp.asv <- ASV.venn.diagram[, METAs$Sample_type == "sediment"]
tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
venn.list$sediment <- rownames(tmp.asv)
tmp.asv <- ASV.venn.diagram[, METAs$Sample_type == "filter"]
tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
venn.list$filter <- rownames(tmp.asv)

# Venn Diagram
str(venn.list)
dim(tmp.rel)

ggVennDiagram(venn.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow") + 
  ggtitle(s)


# another option to create the list
# setting same criteria to all sample types

# venn.list <- vector("list", length = length(unique(METAs$Sample_type)))
# names(venn.list) <- unique(METAs$Sample_type)
# 
# for (v in unique(METAs$Sample_type)) { # without unique also works, but take twice the time to run the code
#   temp1 <- c(rownames(METAs[METAs$Sample_type == v,]))
#   ASV.temp <- ASV.merge.site[,(colnames(ASV.merge.site) %in% temp1)] # selecting ASVs of the same sample type
#   ASV.temp <- ASV.temp[apply(ASV.temp, 1, function(x) !all(x==0)),]
#   rnames.ASV.temp <- rownames(ASV.temp)
#   venn.list[[v]] <- rnames.ASV.temp
# }



##########

### Jaccard boxplot #####
library(reshape)
require(tidyverse)

s <- 'Tel Shikmona' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s
ASV.jaccard <- ASV.merge[,rownames(METAs)] # select site s
ASV.jaccard <- ASV.jaccard[rowSums(ASV.jaccard) > 0,] # remove ASVs that are not present in any sample in this sub setting

ASV.jaccard.pa <- decostand(ASV.jaccard, "pa")

JC <- vegdist(t(ASV.jaccard.pa), method = "jaccard", binary = T)


tmp <- combn(rownames(METAs), 2) # 2 for pairwise comparison
JCmelt <- melt(as.matrix(JC))
JC.df <- map_dfr(1:ncol(tmp), function(x) JCmelt[JCmelt[, 1] ==tmp[1,x] & JCmelt[, 2] == tmp[2,x],])

JC.df$shared <- 1-JC.df$value

JC.df$groupX1 <- METAs[JC.df$X1, "Sample_type"]
JC.df$groupX2 <- METAs[JC.df$X2, "Sample_type"]


JC.df$comparison <- ifelse(
  JC.df$groupX1 == JC.df$groupX2 | !apply(JC.df[, c("groupX1", "groupX2")], 1, function(x) any(x == "single_cell")),
  "ignore",
  ifelse(
    apply(JC.df[, c("groupX1", "groupX2")], 1, function(x) any(x == "sediment")),
    "single_vs_sediment",
    "single_vs_filter"
  )
)


boxplot(shared*100 ~ comparison, data = JC.df[JC.df$comparison != "ignore",])




##########


### BC boxplot ######

#s <- 'Capo Passero' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s
ASV.BC <- ASV.merge[,rownames(METAs)] # select site s
ASV.BC <- ASV.BC[rowSums(ASV.BC) > 0,] # remove ASVs that are not present in any sample in this sub setting

BC <- vegdist(t(ASV.BC))
              
tmp <- combn(rownames(METAs), 2) # 2 for pairwise comparison
BCmelt <- melt(as.matrix(BC))
BC.df <- map_dfr(1:ncol(tmp), function(x) BCmelt[BCmelt[, 1] ==tmp[1,x] & BCmelt[, 2] == tmp[2,x],]) # removing comparisons between same sample and removing the second time the same comparison is made (x1 vs x2 is equal to x2 vs x1)

BC.df$shared <- 1-BC.df$value

BC.df$groupX1 <- METAs[BC.df$X1, "Sample_type"]
BC.df$groupX2 <- METAs[BC.df$X2, "Sample_type"]


BC.df$comparison <- ifelse(
  BC.df$groupX1 == BC.df$groupX2 | !apply(BC.df[, c("groupX1", "groupX2")], 1, function(x) any(x == "single_cell")),
  "ignore",
  ifelse(
    apply(BC.df[, c("groupX1", "groupX2")], 1, function(x) any(x == "sediment")),
    "single_vs_sediment",
    "single_vs_filter"
  )
)


boxplot(shared*100 ~ comparison, data = BC.df[BC.df$comparison != "ignore",])

##########













##########################################################
################# ANOSIMposthoc function ################# 
##########################################################

ANOSIMposthoc=function(M,E,distance="bray", padj="fdr"){
  E=droplevels(E)
  Mlist=list()
  for(i in 1:length(levels(E))){
    Mlist[[i]]=M[E==levels(E)[i],]
  }
  
  Elist=list()
  for(i in 1:length(levels(E))){
    Elist[[i]]=as.numeric(E[E==levels(E)[i]])
  }
  
  result=list(anosimR=matrix(NA,length(levels(E)),length(levels(E))),
              anosimP=matrix(NA,length(levels(E)),length(levels(E))),
              anosimPadj=matrix(NA,length(levels(E)),length(levels(E))))
  colnames(result$anosimR)=colnames(result$anosimP)=levels(E)
  rownames(result$anosimR)=rownames(result$anosimP)=levels(E)
  for(i in 1:(length(levels(E))-1)){
    for(j in (i+1):length(levels(E))){
      temp=anosim(rbind(Mlist[[i]],Mlist[[j]]),c(Elist[[i]],Elist[[j]]),distance=distance)
      result$anosimR[j,i]=temp$statistic
      result$anosimP[j,i]=temp$signif
    }
  }
  result$anosimPadj=matrix(p.adjust(as.vector(result$anosimP),method=padj,
                                    n=length(which(!is.na(as.vector(result$anosimP))))),
                           length(levels(E)),length(levels(E)))
  colnames(result$anosimPadj)=levels(E)
  rownames(result$anosimPadj)=levels(E)
  return(result)
}



