#################################
# Analysis of 18S amplicon data #
#################################

#              Parameters             #


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

######

Sys.setenv(LANG = "en") # in case system language changes to german



### prepare environment ####

# set working directory
setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_euks")

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
#load("data_exploration_all_euks.Rdata")
#save.image("data_exploration_all_euks.Rdata")

#####


### read data ####

# ASV table
ASV <- read.table(
  "otu_tab_ssu_all_euks_curated.txt", 
  h = T,
  sep = "\t",
  row.names = 1
)

META <- readxl::read_xlsx("Lists_bioinformatic_all_euks.xlsx") 

# taxonomy of microbial community
TAX <- read.table(
  "tax_tab_ssu_all_euks_curated_2.txt", 
  h = T, 
  sep = "\t", 
  row.names = 1,
  comment.char = "",
  quote = "", 
  stringsAsFactors = F
)


# check data types
str(ASV)
str(TAX)
str(META)

# adjusting Purified_PCR_product from _ to - 
colnames(ASV) <- gsub("_", "-", colnames(ASV))


# adjust object and data types
ASV <- as.matrix(ASV)
TAX <- as.matrix(TAX)
META$Site <- factor(META$Site, levels = c("Capo Passero", "Plemmirio", "Tel Shikmona", "Eilat", "NC"))
META$Condition <- factor(META$Condition, levels = c("Shallow_Algae", "Shallow_Rubbles", "Deep_Algae", "Deep_Rubbles", "NC"))
META$Sample_type <- factor(META$Sample_type, levels = c("filter", "sediment", "single_cell", "NC_env", "NC_single_cell"))
META$Sample_NC <- as.factor(META$Sample_NC)
levels(META$Sample_type)
# [1] "filter"         "sediment"       "single_cell"    "NC_env"         "NC_single_cell"

# reorder data
META <- META[order(META$Site, META$Condition), ]
META <- META %>% column_to_rownames("Purified_PCR_product") # had to do it first in my data
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
#####


#### cleaning data #####

# it was necessary to remove a PCR replicate (NC amplified by mistake - PR-1040) and an entire sample (S0533 - PR-0736 and PR-0737)
# which was a failed PCR (when we looked the rarefaction curves, it just had 400 seqs)   

samples_to_remove <-c("PR-1040", "PR-0736", "PR-0737") # creating vector with samples to remove 
META <- META[!(row.names(META) %in% samples_to_remove),] # removing row names with the name of the samples
ASV <- ASV[,rownames(META)]
# investigate if there are ASVs with zero sequences after removing these samples
sort(rowSums(ASV)) # to view sum sorted from lower to higher number
# no ASVs with zero sequences! 
# no need to change the TAX table in this case

all.equal(colnames(ASV), rownames(META))
# [1] TRUE
all.equal(rownames(ASV), rownames(TAX))
# [1] TRUE

#######

#### rearrange data set to combine the technical replicates - important step before moving to the statistical analysis #####
METAnoNC <- META[META$Sample_NC == "Sample",] # working only with samples, excluding the NCs

# creating vector list to receive loop information
ASV.rep.list <- vector("list", length = length(unique(METAnoNC$Extraction_Voucher)))
names(ASV.rep.list) <- unique(METAnoNC$Extraction_Voucher)

for (i in unique(METAnoNC$Extraction_Voucher)) { # without unique also works, but take twice the time to run the code
  META.rep <- METAnoNC[METAnoNC$Extraction_Voucher == i, ]
  replicates <- c(rownames(META.rep)) # take the technical replicates names correspondents to this extraction voucher
  ASV.rep <- ASV[,(colnames(ASV) %in% replicates)] # selecting ASVs of the technical replicates
  temp <- ASV.rep[apply(ASV.rep, 1, function(x) !any(x==0)),] # keeping only ASVs present in both replicates
  ASV.rep.temp <- as.data.frame(rowSums(temp)) # summing the replicates in a new data frame
  colnames(ASV.rep.temp) <- i  #saving column name as the extraction voucher
  ASV.rep.list[[i]] <- ASV.rep.temp
}

# converting the list into a data frame with map_dfr  
require(dplyr)
ASV.rep.df <- map_dfr(ASV.rep.list, function (x) data.frame(t(x))) # it has to be transpose so the function combine the rows and keep all columns (adding NAs when there is no value)
ASV.merge <- t(ASV.rep.df) # transpose again to have in the same format as the original ASV table
ASV.merge[is.na(ASV.merge)] <- 0 # setting NAs as zero 


# another (shorter) approach (no loop and using the entire function inside map_dfr)
require(dplyr)
ASV.merge <- map_dfr(
  unique(METAnoNC$Extraction_Voucher),
  function(i) { 
  META.rep <- METAnoNC[METAnoNC$Extraction_Voucher == i, ]
  ASV.rep <- ASV[, rownames(META.rep)] # selecting ASVs of the technical replicates from the same sample
  temp <- ASV.rep[apply(ASV.rep, 1, function(x) !any(x==0)),] # keeping only ASVs present in both replicates
  ASV.rep.temp <- as.data.frame(rowSums(temp)) # summing the replicates in a new data frame
  colnames(ASV.rep.temp) <- i  #saving column name as the extraction voucher
  data.frame(t(ASV.rep.temp))
  
})
ASV.merge <- t(ASV.rep.df) # transpose again to have in the same format as the original ASV table
ASV.merge[is.na(ASV.merge)] <- 0 # setting NAs as zero

# ASV.merge2 <- map_dfr(
#   unique(METAnoNC$Extraction_Voucher),
#   function(i) { 
#     META.rep <- METAnoNC[METAnoNC$Extraction_Voucher == i, ]
#     ASV.rep <- ASV[, rownames(META.rep)] # selecting ASVs of the technical replicates from the same sample
#     temp <- ASV.rep # keeping all ASVs present in both replicates
#     ASV.rep.temp <- as.data.frame(rowSums(temp)) # summing the replicates in a new data frame
#     colnames(ASV.rep.temp) <- i  #saving column name as the extraction voucher
#     data.frame(t(ASV.rep.temp))
#     
#   })
# ASV.merge2 <- t(ASV.rep.df) # transpose again to have in the same format as the original ASV table
# ASV.merge2[is.na(ASV.merge2)] <- 0 # setting NAs as zero

# calculating how many ASVs we lost
summary(colSums(ASV.merge)/c(by(colSums(ASV), META$Extraction_Voucher, sum))[colnames(ASV.merge)])
# function by considers the elements (1 - data, 2 - grouping factor, 3 - function to apply)
# we subset with [] just the rows in ASV that correspond to the ASVs present in the colnames of ASV.merge ( making sure we dont take from NC)
# this by function works because we make sure with all.equal that the colnames of ASV are in the same order than the rownames of META
# we use c to obtain the value as a vector (the by function gives a list) that we could use in the calculation

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9052  0.9651  0.9799  0.9754  0.9921  0.9985 

# we don't have a major loss, so we are good to go and continue the analysis

nrow(ASV.merge) #[1] 10180
nrow(ASV) # [1] 13769

# Creating META merge table to analyze together with ASV.merge
META.merge <- META[META$Sample_NC == "Sample",] # taking META without NCs
META.merge <- META.merge[!grepl("_2", META.merge$PCR_Product),] # removing the replicated lines for the technical replicates
rownames(META.merge) <- META.merge$Extraction_Voucher
META.merge$DNA_yield_mean <- c(NA) # create empty column to add data of DNA concentration mean between replicates
META.merge$DNA_yield_dif <- c(NA) # create empty column to add data of DNA concentration difference between replicates
for (i in META.merge$Extraction_Voucher){ # selecting sample by sample in a loop
  META.merge[i,"DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) # calculating the mean of DNA conc among the technical replicates
  META.merge[i,"DNA_yield_dif"] <- max(c(dist(META$DNA_concentration[META$Extraction_Voucher == i]))) # calculating the difference between the DNA conc of the technical replicates
  #using max to make a more universal code in case one has more than two PCR replicates
}

str(META.merge)
META.merge$Sample_type <- droplevels(META.merge$Sample_type)
META.merge$Site <- droplevels(META.merge$Site)
META.merge$Condition <- droplevels(META.merge$Condition)

ASV.merge <- ASV.merge[, rownames(META.merge)]
all.equal(colnames(ASV.merge), rownames(META.merge))
# [1] TRUE

# adjusting TAX table to just show the ASVs that remained in ASV.merge
TAX.merge <- TAX[rownames(ASV.merge),]
all.equal(rownames(ASV.merge), rownames(TAX.merge))
# [1] TRUE

ASV.rel.merge <- prop.table(ASV.merge, 2) * 100


# Saving names of ASVs kept, to select them on new fasta file (to do it in GitBash script "select_seqs_fasta.bash")
names.accnos <- rownames(ASV.merge)

# saving the vector file to a txt file in desktop
write(names.accnos, "names.accnos.txt")

#saving ASV.merge in desktop for lulu/swarm investigation
write.table(ASV.merge, "asv_tab_ssu_all_euks_merge.txt", quote = F, sep = "\t")

#saving TAX.merge (which is the same as META.pen)
write.table(TAX.merge, "TAX_merge_all_euks.txt", quote = F, sep = "\t")

# saving as R objects to be read in further scripts 
saveRDS(ASV.merge, "ASV_merge_euks.RDS")
saveRDS(TAX.merge, "TAX_merge_euks.RDS")
saveRDS(META.merge, "META_merge_euks.RDS")

#######


##### Removing taxa that are not possible symbionts - DO NOT DO IT FOR NOW ######

# Remove the genus Jania of the Rhodophyta division 
# Jania is one of the main macro algae in which foraminifera leaves attached to
temp <- TAX.merge[TAX.merge[, "Genus"] == "Jania",] # just looking how many seqs as there, not so many. is okay to remove.
temp <- ASV.merge[TAX.merge[, "Genus"] == "Jania",]
tmp <- META.merge[colnames(temp),] # all sites and all sample types have these seqs

# Doing the same for Metazoans (cannot possibily be an endosymbiont of foraminifera)
temp <- TAX.merge[TAX.merge[, "Division"] == "Metazoa",] 
temp <- ASV.merge[TAX.merge[, "Division"] == "Metazoa",]
tmp <- META.merge[colnames(temp),] # all sites and all sample types have these seqs

# Applying the filters to remove these TAX
TAX.filt <- TAX.merge[TAX.merge[, "Genus"] != "Jania" & TAX.merge[, "Division"] != "Metazoa",]
ASV.filt <- ASV.merge[TAX.merge[, "Genus"] != "Jania" & TAX.merge[, "Division"] != "Metazoa",]

all.equal(rownames(ASV.filt), rownames(TAX.filt))
# [1] TRUE

nrow(TAX.filt)/nrow(TAX.merge)
# 0.6908644

sum(ASV.filt)/sum(ASV.merge)
# 0.8995442


#######

##### Subsetting datasets to just investigate the diatom symbionts ######

TAX.diat <- as.data.frame(TAX.merge)
str(TAX.diat)
TAX.diat$Class <- as.factor(TAX.diat$Class)
TAX.diat <- TAX.diat[TAX.diat$Class == "Bacillariophyta",]
TAX.diat$Class <- droplevels(TAX.diat$Class)
TAX.diat <- as.matrix(TAX.diat)

# o use the code bellow - the same output
TAX.diat <- TAX.merge[TAX.merge[, "Class"] == "Bacillariophyta", ]

ASV.diat <- ASV.merge[rownames(TAX.diat),]
all.equal(rownames(ASV.diat), rownames(TAX.diat))
# [1] TRUE

META.diat <- META.merge[colnames(ASV.diat),]
all.equal(colnames(ASV.diat), rownames(META.diat))
# [1] TRUE

str(META.diat)
META.diat$Site <- droplevels(META.diat$Site)
levels(META.diat$Site)
#"Capo Passero" "Plemmirio"    "Tel Shikmona" "Eilat"  

#######

##### Subsetting diatoms data set to just observe the arapid pennate family ######
# the absolute most abundant familiy in the foraminifera microbiome
TAX.pen <- TAX.diat[TAX.diat[, "Family"] == "Araphid-pennate",]

ASV.pen <- ASV.diat[rownames(TAX.pen),]
all.equal(rownames(ASV.pen), rownames(TAX.pen))
# [1] TRUE

META.pen <- META.diat[colnames(ASV.pen),]
all.equal(colnames(ASV.pen), rownames(META.pen))
# [1] TRUE

# Saving names of ASVs kept, to select them on new fasta file (to do it in GitBash script "select_seqs_fasta.bash")
names.pennate <- rownames(ASV.pen)

# saving the vector file to a txt file in desktop
write(names.pennate, "names.pennate.txt")

#saving ASV.pen in desktop for phylogenetic investigation
write.table(ASV.pen, "asv_tab_ssu_all_penatte.txt", quote = F, sep = "\t")

#saving META.merge (which is the same as META.pen)
write.table(META.merge, "META_all_euks.txt", quote = F, sep = "\t")
#######

##### Exploring how many ASVs are in common between single cell and environmental samples ####

ASV.diat.filt.sing <- ASV.diat[,META.diat$Site == "Eilat" & META.diat$Sample_type == "single_cell"] # select site and sample type
ASV.diat.filt.sing <- ASV.diat.filt.sing[apply(ASV.diat.filt.sing, 1, function(x) !all(x==0)),] # remove ASVs that are present in this filtering setting

ASV.diat.filt.env <- ASV.diat[,META.diat$Site == "Eilat" & META.diat$Sample_type == "filter"] # select site and sample type
ASV.diat.filt.env <- ASV.diat.filt.env[apply(ASV.diat.filt.env, 1, function(x) !all(x==0)),] # remove ASVs that are present in this filtering setting

nrow(ASV.diat.filt.sing)
nrow(ASV.diat.filt.env)

# to see ASVs that are in single cell and not in env and vice-versa
rownames(ASV.diat.filt.sing)[!rownames(ASV.diat.filt.sing) %in% rownames(ASV.diat.filt.env)] # ASVs in single cell that are not in the env samples
rownames(ASV.diat.filt.env)[!rownames(ASV.diat.filt.env) %in% rownames(ASV.diat.filt.sing)] # ASVs in env sample that are not in the single cells

# ASVs in common single cell and env samples
rownames(ASV.diat.filt.sing)[rownames(ASV.diat.filt.sing) %in% rownames(ASV.diat.filt.env)]


# ASVs in common between Capo Passero single cell and Capo Passero env samples
# [1] "sq21"  "sq373" "sq38"  "sq353" "sq422" "sq654"                          # sediment
# [1] "sq21"  "sq190" "sq38"  "sq422" "sq654"                                  # sea water

# ASVs in common between Plemmirio single cell and Plemmirio env samples
# [1] "sq21"  "sq190" "sq373" "sq38"  "sq63"  "sq353" "sq947" "sq654" "sq280"
# [1] "sq21"  "sq190" "sq38"  "sq63"  "sq353" "sq654"

# ASVs in common between Tel Shikmona single cell and Tel Shikmona env samples
# [1] "sq21"   "sq190"  "sq38"   "sq358"  "sq422"  "sq1106" "sq1527" "sq280"  "sq4448"
# [1] "sq21"   "sq38"   "sq358"  "sq1527" "sq280"

# ASVs in common between Eilat single cell and Eilat env samples
# [1] "sq21" "sq38"
# [1] "sq21" "sq190" "sq38" 

# re shaping ASV table for plotting
# melt from reshape package
ASV.diat.plot <- melt(ASV.diat)
colnames(ASV.diat.plot) <- c("ASV_name", "Sample_ID", "n_seqs")
str(ASV.diat.plot)
max(ASV.diat.plot$n_seqs)
# [1] 89445
summary(ASV.diat.plot$n_seqs)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00     0.00     0.00    34.36     0.00 89445.00 

single_cell_samples <- rownames(META.diat)[META.diat$Sample_type == "single_cell"] # select only ASVs of single cell 
env_samples <- rownames(META.diat)[META.diat$Sample_type != "single_cell"] # select only ASVs of env samples

max(ASV.diat.plot[ASV.diat.plot$Sample_ID %in% single_cell_samples,]$n_seqs)
# [1] 89445

max(ASV.diat.plot[ASV.diat.plot$Sample_ID %in% env_samples,]$n_seqs)
# [1] 8214


# Plotting number of ASVs (y) per ASV type (x), one facet for the single cells (all sites) and the other for env samples (all sites)

require(ggplot2)
require(patchwork)


# env samples
p <- ggplot(ASV.diat.plot[ASV.diat.plot$Sample_ID %in% env_samples,], aes(x=ASV_name, y=n_seqs)) + 
  geom_boxplot()
p

# single cell
p1 <- ggplot(ASV.diat.plot[ASV.diat.plot$Sample_ID %in% single_cell_samples,], aes(x=ASV_name, y=n_seqs)) + 
  geom_boxplot()
p1

# using patchwork package to see plots together
p1 / p


## to see per specific site
ASV.diat.site <- ASV.diat[,META.diat$Site == "Eilat"] # select site X
ASV.diat.site <- ASV.diat.site[apply(ASV.diat.site, 1, function(x) !all(x==0)),]
ASV.diat.site.plot <- melt(ASV.diat.site)
colnames(ASV.diat.site.plot) <- c("ASV_name", "Sample_ID", "n_seqs")

METAs <- META.diat[META.diat$Site == "Eilat",] # select site X

single_cell_samples <- rownames(METAs)[METAs$Sample_type == "single_cell"] # select only ASVs of single cell 
env_samples <- rownames(METAs)[METAs$Sample_type != "single_cell"] # select only ASVs of env samples 


# env samples
p <- ggplot(ASV.diat.site.plot[ASV.diat.site.plot$Sample_ID %in% env_samples,], aes(x=ASV_name, y=n_seqs)) + 
  geom_boxplot()


# single cell
p1 <- ggplot(ASV.diat.site.plot[ASV.diat.site.plot$Sample_ID %in% single_cell_samples,], aes(x=ASV_name, y=n_seqs)) + 
  geom_boxplot()


# using patchwork package to see plots together
p1 / p


# comparing sites

 ggplot(ASV.diat.site.plot[ASV.diat.site.plot$Sample_ID %in% env_samples,], aes(x=ASV_name, y=n_seqs)) + 
  geom_boxplot()


pEil / pShi / pCap / pPle

######

####### Exploring similarity between single cell and environmental samples #####

# complete data set
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

##### only diatoms
# calculate proportions
ASV.rel.diat <- prop.table(ASV.diat, 2) * 100

# only single cell
ASV.diat.sing <- ASV.diat[,META.diat$Sample_type == "single_cell"] # select single cell
ASV.diat.sing <- ASV.diat.sing[rowSums(ASV.diat.sing) > 0,]
META.sing <- META.diat[colnames(ASV.diat.sing),]

ASV.rel.sing <- prop.table(ASV.diat.sing, 2) * 100

# only env samples
ASV.diat.env <- ASV.diat[,META.diat$Sample_type != "single_cell"] # select env samples
ASV.diat.env <- ASV.diat.env[apply(ASV.diat.env, 1, function(x) !all(x==0)),]
META.env <- META.diat[colnames(ASV.diat.env),]

ASV.rel.env <- prop.table(ASV.diat.env, 2) * 100
######


### Visualization of data 


### rarefaction curves Richness and Inverse Simpson ####

### Richness (number of ASVs)

# per site
rarecurve(t(ASV.merge), step = 100, col = META.merge$color_site, label = F, ylab = "Number of ASVs")

legend(
  'topright',
  legend = levels(META.merge$Site),
  col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick", "black"),
  pch = 15,
  pt.cex = 1.3,
  title = "Site",
  cex= 0.75 # to make legend box smaller
)

min(colSums(ASV.merge))

# per sample type
rarecurve(t(ASV.merge), step = 100, col = META.merge$color_sample_type, label = F, ylab = "Number of ASVs")

legend(
  'topright',
  legend = levels(META.merge$Sample_type),
  col = c("dodgerblue", "blue", "red", "azure4", "dimgrey"),
  pch = 15,
  pt.cex = 1.3,
  title = "Sample type",
  cex= 0.75 # to make legend box smaller
)


### Inverse Simpson index
# rarefaction curves for the inverse Simpson index, using some custom code
# for each sequencing depth, repeatedly randomly subsample the data set and calculate inverse Simpson index
knots <- 100 # pick 100 sequencing depth for rarefaction
n <- 10 # randomly subsample 10 times at each sequencing depth
# the output in a list object
invs.repeated <- vector("list", length = ncol(ASV.merge))
# repeat the process for each sample (i.e. columns in ASV)
for(i in 1:ncol(ASV.merge)) {
  print(i) # show progress (sample number)
  temp <- ASV.merge[, i] # subset to only one sample
  temp <- temp[temp > 0] # remove ASVs not occurring in that sample
  
  # build matrix to collect inverse Simpson values in
  invs.repeated[[i]] <- matrix(NA, nrow = knots, ncol = n + 1)
  
  # determine 'knots' sequencing depth for subsampling from 0 to the available sequencing depth
  temp.sub <- round(seq(0, sum(temp), length.out = knots))
  
  # save these values in the first column of the output table for this sample
  invs.repeated[[i]][, 1] <- temp.sub
  
  # repeat n times
  for(j in 1:n) {
    
    # for each sequencing depth (starting with the second - zero is not used)
    for(k in 2:knots) {
      # subsample and calculate inverse Simpson
      invs.repeated[[i]][k, j + 1] <- diversity(
        rrarefy(
          temp,
          temp.sub[k]
        ),
        "invsimpson"
      )
    }
    
    # set zero sequencing depth to zero inverse Simpson
    invs.repeated[[i]][temp.sub==0, j + 1] <- 0
  }
}

str(invs.repeated)

# get maximum for y axis
maxIndex <- max(
  sapply(
    invs.repeated,
    function(x) {
      apply(x[, -1], 1, mean)
    }
  )
)

# start plot per site
plot(
  0, 0,
  type = "n",
  xlim = c(0, max(colSums(ASV.merge))),
  ylim = c(0, max(maxIndex)),
  xlab = "Sequencing depth",
  ylab = "Inverse Simpson index",
  cex.axis = 0.7,
  las = 1,
  axes = F,
  cex.lab = 0.9,
  mgp = c(1.8, 0.4, 0)
)
axis(1, las = 1, mgp = c(1, 0.2, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
axis(2, las = 1, mgp = c(2, 0.4, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
box("plot", lwd = 0.5)
for(k in order(colSums(ASV.merge), decreasing = T)) {
  subsample <- invs.repeated[[k]][, 1]
  # show mean of the 10 time repeated random subsampling
  temp <- apply(invs.repeated[[k]][, -1], 1, mean)
  lines(
    subsample,
    temp,
    lwd = 1,
    col = META.merge$color_site[k]
  )
}
# show minimum sequencing depth in the data set
abline(v = min(colSums(ASV.merge)))



# start plot per sample type
plot(
  0, 0,
  type = "n",
  xlim = c(0, max(colSums(ASV.merge))),
  ylim = c(0, max(maxIndex)),
  xlab = "Sequencing depth",
  ylab = "Inverse Simpson index",
  cex.axis = 0.7,
  las = 1,
  axes = F,
  cex.lab = 0.9,
  mgp = c(1.8, 0.4, 0)
)
axis(1, las = 1, mgp = c(1, 0.2, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
axis(2, las = 1, mgp = c(2, 0.4, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
box("plot", lwd = 0.5)
for(k in order(colSums(ASV.merge), decreasing = T)) {
  subsample <- invs.repeated[[k]][, 1]
  # show mean of the 10 time repeated random subsampling
  temp <- apply(invs.repeated[[k]][, -1], 1, mean)
  lines(
    subsample,
    temp,
    lwd = 1,
    col = META.merge$color_sample_type[k]
  )
}
# show minimum sequencing depth in the data set
abline(v = min(colSums(ASV.merge)))


# # show only the samples with less than 1000 seq
# # start plot
# plot(
#   0, 0,
#   type = "n",
#   xlim = c(0, 1000),
#   ylim = c(0, max(maxIndex)),
#   xlab = "Sequencing depth",
#   ylab = "Inverse Simpson index",
#   cex.axis = 0.7,
#   las = 1,
#   axes = F,
#   cex.lab = 0.9,
#   mgp = c(1.8, 0.4, 0)
# )
# axis(1, las = 1, mgp = c(1, 0.2, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
# axis(2, las = 1, mgp = c(2, 0.4, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
# box("plot", lwd = 0.5)
# for(k in order(colSums(ASV.merge), decreasing = T)) {
#   subsample <- invs.repeated[[k]][, 1]
#   # show mean of the 10 time repeated random subsampling
#   temp <- apply(invs.repeated[[k]][, -1], 1, mean)
#   if (colSums(ASV) [k] < 1000) {lines(
#     subsample,
#     temp,
#     lwd = 1,
#     col = META.merge$color_site[k]
#   )}
# }
# # show minimum sequencing depth in the data set
# abline(v = min(colSums(ASV.merge)))
# 


#####

####### Plotting seq depths #####

barplot(sort(colSums(ASV.merge)))

barplot(sort(colSums(ASV.merge))[1:100], ylim = c(0, 20000))
abline(h=500)


barplot(sort(colSums(ASV.merge))[1:50], ylim = c(0, 5000))
abline(h=500)


# Filtering META data
View(META.merge[colSums(ASV.merge) < 500, ]) 
# only one sample S0533

# failed PCR? maybe remove it...



####### 

### alpha diversity (boxplot) ####

# calculate the inverse Simpson index without subsampling
invS <- diversity(t(ASV.merge), "invsimpson")

# show boxplot with points
par(mar = c(10,4,4,2) + 0.1) # adjusting the margins size to fit the x axis labels

boxplot(
  invS ~ droplevels(interaction(META.merge$Sample_type, META.merge$Site)), 
  las = 2, # 2 means that the orientation of the x axis label will be perpendicular to the x axis
  ylab = "Inverse Simpson index",
  xlab = "", 
  cex.axis = 0.7,
  col = c("dodgerblue", "blue", "red", "dodgerblue", "blue", "red", "dodgerblue", "blue", "red", "dodgerblue", "blue", "red"),
  outline = F
  
)
points(
  jitter(as.numeric(droplevels(interaction(META.merge$Sample_type, META.merge$Site)))), 
  invS,
  #bg = META.merge$Condition,
  pch = 16
)

# add legend for experimental conditions
# legend(
#   "topleft",
#   legend = levels(META.merge$Sample_type),
#   pch = 21:25,
#   pt.bg = "white",
#   pt.cex = 1.3,
#   title = "Sample type"
# )
# legend(
#   "topleft",
#   legend = levels(META.merge$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick", "black"),
#   pch = 15,
#   pt.cex = 1.3,
#   title = "Site"
# )



#####

### taxonomic composition ####

# barplot of the taxonomic composition
TAX.pooled <- vector(mode = "list", length = 8)
names(TAX.pooled) <- colnames(TAX.merge)[1:8] # kingdom to species
for (i in 1:8) {
  temp <- aggregate(
    ASV.merge,
    by = list(TAX.merge[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}
TAX.pooled.rel <- lapply(
  TAX.pooled,
  function(x) {
    prop.table(x, 2) * 100
  }
)

# PlotAbund function at the end of the script

# Here we will take a shortcut to produce a stacked barplot, only showing the most dominant taxa
# The following function can be downloaded from: https://raw.githubusercontent.com/chassenr/NGS/master/Plotting/PlotAbund.R
#source("C:/Users/chassenrueck/Documents/Repos/NGS/Plotting/PlotAbund.R")

PlotAbund(
  TAX.pooled.rel$Class, # change levels to see the different plots for each one
  abund = 2,
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T
)


# Also try other input parameters to optomize visualization
#   Change taxonomic resolution
#   Change abund parameter
#   Switch between method "percentage" and "nmost"

# Which taxa dominate the community?
# Bacillariophyta (diatoms)

# do same plot but looking only at the diatoms (Class Bacillariophyta)
TAX.pooled <- vector(mode = "list", length = 8)
names(TAX.pooled) <- colnames(TAX.diat)[1:8] # kingdom to species
for (i in 1:8) {
  temp <- aggregate(
    ASV.diat,
    by = list(TAX.diat[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}
TAX.pooled.rel <- lapply(
  TAX.pooled,
  function(x) {
    prop.table(x, 2) * 100
  }
)

PlotAbund(
  TAX.pooled.rel$Genus, # change levels to see the different plots for each one
  abund = 2,
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T
)


# Looking at the most abundant ASVs in the microbiome of the single cells

# selecting ASVs of single cells only
# all euks
ASV.merge.sing <- ASV.merge[,META.merge$Sample_type == "single_cell"] # select single cell
ASV.merge.sing <- ASV.merge.sing[rowSums(ASV.merge.sing) > 0,] # excluding ASVs that are not in the subset
META.sing <- META.merge[colnames(ASV.merge.sing),]
all.equal(rownames(META.sing), colnames(ASV.merge.sing))

ASV.rel.sing <- prop.table(ASV.merge.sing, 2) * 100

# Define ASV and TAX table based on most abundant ASVs
# Remove rare ASVs (that don't occur with at least 5% in 1 sample)
ASV.abund <- ASV.merge.sing[apply(ASV.rel.sing, 1, function(x) sum(x >= 5) >= 1), ]
nrow(ASV.abund)/nrow(ASV.merge.sing)
# 0.05275407

ASV.abund <- ASV.rel.sing[apply(ASV.rel.sing, 1, function(x) sum(x >= 5) >= 1), ]
summary(colSums(ASV.abund))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 61.65   91.70   94.55   92.97   96.13   98.03 

# removing rare ASVs of TAX table
TAX.abund <- TAX.merge[rownames(ASV.abund),]

# looking at the barplot of the taxonomic composition
TAX.pooled <- vector(mode = "list", length = 8)
names(TAX.pooled) <- colnames(TAX.abund)[1:8] # kingdom to species
for (i in 1:8) {
  temp <- aggregate(
    ASV.abund,
    by = list(TAX.abund[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}
# TAX.pooled.rel <- lapply(
#   TAX.pooled,
#   function(x) {
#     prop.table(x, 2) * 100
#   }
# )
dev.off()
barplot(TAX.pooled$Division, col = rainbow(nrow(TAX.pooled$Division)), legend.text = row.names(TAX.pooled$Division), args.legend = list(x = "right"), ylim=c(0,100))

# PlotAbund(
#   TAX.pooled.rel$Division, # change levels to see the different plots for each one
#   abund = 2,
#   method = "nmost", # or nmost
#   open.window = T,
#   plot.ratio = c(3.5, 1),
#   sort.taxa = T
# )


# looking abundant ASVs from single cell in this env dataset
ASV.abund.share <- ASV.rel.merge[rownames(ASV.abund),META.merge$Sample_type != "single_cell"] # select env samples
rowSums(ASV.abund.share)
# sq1   sq2   sq3   sq4   sq5   sq6   sq7   sq8  sq10  sq13  sq14  sq22  sq30  sq31  sq59  sq64  sq67  sq84 
# 0     0     0     0     0     0     0     0     0     0     0     0     0     0  1701     0     0     0 
# sq153 sq164   sq9  sq11  sq12  sq47  sq88 sq108 sq125 sq303  sq35  sq46  sq28  sq33  sq73  sq81  sq16  sq76 
# 300     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0 
# sq83 sq100 sq123 sq411  sq17 sq188  sq18  sq38  sq49 sq173 sq318  sq48  sq68  sq27  sq41  sq79 sq104  sq26 
# 0     0     0     0     0     0    25  6093  7623     0     0     0  2923     0    32     0     0     0 
# sq168  sq97  sq51  sq60 sq130  sq74  sq34 sq107  sq54 sq112 sq170  sq82  sq43  sq56 
# 0     0     0     0     0    96     0    78     0     0     0     0     0     0 

ASV.abund.share <- ASV.abund.share[rowSums(ASV.abund.share) > 0,] # excluding ASVs that are not in the subset
# sq59 sq153  sq18  sq38  sq49  sq68  sq41  sq74 sq107 
# 1701   300    25  6093  7623  2923    32    96    78 

# sq59       sq153        sq18        sq38        sq49        sq68        sq41        sq74       sq107 
# 3.20195702  0.56111386  0.05106530 10.41828118 12.97318894  5.39201489  0.06554073  0.16018896  0.14155697

# removing rare ASVs of TAX table
TAX.abund.share <- TAX.merge[rownames(ASV.abund.share),]

colnames(ASV.abund.share)[ASV.abund.share["sq38",] >0]
META.merge[colnames(ASV.abund.share)[ASV.abund.share["sq38",] >0],]

summary(ASV.rel.sing["sq38",])

sort(ASV.rel.sing["sq38",])

META.merge["S0494",]

colSums(ASV.merge.sing)["S0494"]
# S0494 
# 88430 


#####

### cluster diagram ####

# calculate Bray-Curtis dissimilarities
BC <- vegdist(t(ASV.rel))

# # save triangular matrix
# BC.full <- as.matrix(as.dist(BC, upper = TRUE, diag = TRUE))
# write.table(
#   BC.full,
#   "Bray_Curtis_triangular.txt",
#   sep = "\t",
#   quote = F
# )
# 
# # save as pairwise list
# BC.list <- melt(as.matrix(BC))
# # this one has still the duplicated values from the upper triangle
# # to filter out those, first extract all possible pairwise comparisons
# pw.comb <- combn(colnames(ASV), 2)
# # then only select those rows in BC.list that correspond to the pairs in pw.comb
# BC.list.clean <- do.call(
#   "rbind",
#   apply(
#     pw.comb, 
#     2, 
#     function(x) {
#       BC.list[BC.list$Var1 == x[1] & BC.list$Var2 == x[2], ]
#     }
#   )
# )
# # write to file
# write.table(
#   BC.list.clean,
#   "Bray_Curtis_pw_list.txt",
#   sep = "\t",
#   quote = F,
#   row.names = F
# )

# Hierarchical clustering

# We will use complete linkage hierarchical clustering based on Bray-Curtis dissimilarity of ASV sequence proportions
# average linkage method shown to be the same in cophenetic analysis performed in "screen_nc_euks.R" script

BC.clust <- hclust(BC, method = "average")

# change for the sub set dataset
META.plot <- META.diat 
plot(BC.clust, labels = paste(META.plot$Site, META.plot$Sample_type), cex = 0.7)



# How and why do the clusters change when choosing average linkage?

#   which linkage algorithm is better?

# Define groups based on equal similarity
# Those will be used in the NMDS plot
rect.hclust(BC.clust, h = 0.7)
BC.groups <- cutree(BC.clust, h = 0.7)



#####

### NMDS plot ####

# calculate NMDS
NMDS <- metaMDS(t(ASV.rel), k = 2)
NMDS$stress
# 0.13432   # sll samples
# 0.2577109 # only single cell

# build nice plot element by element
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n")

# add hulls based on equal similarity
ordihull(NMDS, BC.groups)

# add points for sites
points(
  NMDS$points[, 1],
  NMDS$points[, 2],
  col = META.plot$color_site,
  pch = 16,
  cex = 1
)

# add legend for sites
legend(
  "bottomleft",
  legend = levels(META.plot$Site),
  col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
  pch = 15,
  pt.cex = 1,
  title = "Site"
)

# or

# dev.off() #to reset the PAR settings 
# 
# plot(NMDS, display = "sites", type = "n")
# 
# # add hulls based on equal similarity
# ordihull(NMDS, BC.groups)
# 
# # add points for sample type
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 2],
#   col = META.merge$color_sample_type,
#   pch = 16,
#   cex = 1
# )
# 
# # add legend for sample type
# legend(
#   "bottomright",
#   legend = levels(META.merge$Sample_type),
#   col = c("dodgerblue", "blue", "red", "azure4", "dimgrey"),
#   pch = 15,
#   pt.cex = 1.3,
#   title = "Sample type"
# )

# How accurate is the NMDS ordination? (Hint: look at stress value found in NMDS$stress)
# How well represented are the cluster from the hierarchical cluster diagram?
# How well represented are the original Bray-Curtis dissimilarities?
# Describe the patterns that you see.

# IMPORTANT: Always consider the original Bray-Curtis dissimilarities in your interpretation


#####

#### Jaccard  ##########
library(reshape)
require(tidyverse)

s <- 'Capo Passero' # # select site to subset in s
METAs <- META.diat[META.diat$Site == s,] # select site s
ASV.diat.site <- ASV.diat[,rownames(METAs)] # select site s
ASV.diat.site <- ASV.diat.site[rowSums(ASV.diat.site) > 0,] # remove ASVs that are not present in any sample in this sub setting

ASV.site.pa <- decostand(ASV.diat.site, "pa")

JC <- vegdist(t(ASV.site.pa), method = "jaccard", binary = T)


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


# check for all sites
# check BC as well


######

### Betadispersion ####

betadispersion <- betadisper(
  vegdist(t(ASV.rel.merge)), 
  META.plot$Site,
  sqrt.dist = T
)

# show PCoA
plot(betadispersion)

eigenvals(betadispersion)/sum(eigenvals(betadispersion))
# PCoA1    PCoA2
#

# show within group heterogeneity
boxplot(betadispersion)


# show PCoA for every site
# abort this option cause is not the right approach to sub set the site before running the statistics
# this way the axis calculation will change each time and therefore not being comparable.

dev.off()
layout(matrix(c(1,2,3,4),2,2, byrow=T))

# Capo Passero
s <- 'Capo Passero' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

betadispersion <- betadisper(
  vegdist(t(ASV.rel.plot)), 
  METAs$Sample_type,
  sqrt.dist = T
)

eigenvals(betadispersion)/sum(eigenvals(betadispersion))
# PCoA1    PCoA2
# 4.780551 1.526000

plot(betadispersion, main=paste("Betadisper",s))
plot(betadispersion, main=paste("Betadisper",s), label = T, label.cex = 0.5)
plot(betadispersion, main=paste("Betadisper",s), label = F)

# Plemmirio
s <- 'Plemmirio' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

betadispersion <- betadisper(
  vegdist(t(ASV.rel.plot)), 
  METAs$Sample_type,
  sqrt.dist = T
)

eigenvals(betadispersion)/sum(eigenvals(betadispersion))
# PCoA1    PCoA2


plot(betadispersion, main=paste("Betadisper",s))
plot(betadispersion, main=paste("Betadisper",s), label = F)



# Tel Shikmona
s <- 'Tel Shikmona' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

betadispersion <- betadisper(
  vegdist(t(ASV.rel.plot)), 
  METAs$Sample_type,
  sqrt.dist = T
)
eigenvals(betadispersion)/sum(eigenvals(betadispersion))
# PCoA1    PCoA2


plot(betadispersion, main=paste("Betadisper",s))
plot(betadispersion, main=paste("Betadisper",s), label = F)

# Eilat
s <- 'Eilat' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

betadispersion <- betadisper(
  vegdist(t(ASV.rel.plot)), 
  METAs$Sample_type,
  sqrt.dist = T
)


eigenvals(betadispersion)/sum(eigenvals(betadispersion))
# PCoA1    PCoA2
# 0.29723432 0.12863717

plot(betadispersion, main=paste("Betadisper",s))
plot(betadispersion, main=paste("Betadisper",s), label = F)

# improve code to see all the points
#####

#### PCoA ecodist package ######
PCoA.dist <- vegdist(t(ASV.rel.merge))
PCoA.dist.pco <- pco(PCoA.dist)

# scatterplot of the first two dimensions
plot(PCoA.dist.pco$vectors[,1:2], col=as.numeric(META.merge$Sample_type),
     pch=as.numeric(META.merge$Sample_type), main="PCoA all euks", xlab="PCoA 1", ylab="PCoA 2")

# legend(
#   "bottom",
#   legend = levels(META.merge$Sample_type),
#   col = as.numeric(META.merge$Sample_type),
#   pch = as.numeric(META.merge$Sample_type),
#   pt.cex = 1,
#   title = "Sample type"
# )
# legend does not work


# do this PCoA for every site
# abort this option cause is not the right approach to sub set the site before running the statistics
# this way the axis calculation will change each time and therefore not being comparable.

dev.off()
layout(matrix(c(1,2,3,4),2,2, byrow=T))

# Capo Passero
s <- 'Capo Passero' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

PCoA.dist <- vegdist(t(ASV.rel.plot))
PCoA.dist.pco <- pco(PCoA.dist)

# scatterplot of the first two dimensions
plot(PCoA.dist.pco$vectors[,1:2], col=as.numeric(METAs$Sample_type),
     pch=as.numeric(METAs$Sample_type), main=paste("PCoA",s), xlab="PCoA 1", ylab="PCoA 2")

# Plemmirio
s <- 'Plemmirio' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

PCoA.dist <- vegdist(t(ASV.rel.plot))
PCoA.dist.pco <- pco(PCoA.dist)

# scatterplot of the first two dimensions
plot(PCoA.dist.pco$vectors[,1:2], col=as.numeric(METAs$Sample_type),
     pch=as.numeric(METAs$Sample_type), main=paste("PCoA",s), xlab="PCoA 1", ylab="PCoA 2")

# Tel Shikmona
s <- 'Tel Shikmona' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

PCoA.dist <- vegdist(t(ASV.rel.plot))
PCoA.dist.pco <- pco(PCoA.dist)

# scatterplot of the first two dimensions
plot(PCoA.dist.pco$vectors[,1:2], col=as.numeric(METAs$Sample_type),
     pch=as.numeric(METAs$Sample_type), main=paste("PCoA",s), xlab="PCoA 1", ylab="PCoA 2")

# Eilat
s <- 'Eilat' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.plot <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.plot <- ASV.plot[rowSums(ASV.plot) > 0, ] # remove ASVs that are not present in any sample in this sub setting

ASV.rel.plot <- prop.table(ASV.plot, 2) * 100

PCoA.dist <- vegdist(t(ASV.rel.plot))
PCoA.dist.pco <- pco(PCoA.dist)

# scatterplot of the first two dimensions
plot(PCoA.dist.pco$vectors[,1:2], col=as.numeric(METAs$Sample_type),
     pch=as.numeric(METAs$Sample_type), main=paste("PCoA",s), xlab="PCoA 1", ylab="PCoA 2")

#######

### ANOSIM ####


# Group 1:
compare.site <- "Capo Passero"
# Group 2:
compare.site <- "Plemmirio"
# Group 3:
compare.site <- "Tel Shikmona"
# Group 4:
compare.site <- "Eilat"

# subset data set
META.sub <- META.diat[META.diat$Site == compare.site, ]
ASV.rel.sub <- ASV.rel[, rownames(META.sub)]
ASV.rel.sub <- ASV.rel.sub[rowSums(ASV.rel.sub) > 0, ]

# calculate ANOSIM between sample type
anosim(t(ASV.rel.sub), META.sub$Sample_type)

# calculate pairwise ANOSIM statistics
# download function from: https://raw.githubusercontent.com/chassenr/Tutorials/master/R_roundtable_sequence_analysis/anosimPosthoc.R
source("C:/Users/chassenrueck/Documents/Repos/Tutorials/R_roundtable_sequence_analysis/anosimPosthoc.R") # doesn#t work
# run function anosimPosthoc in the end of the script
ANOSIMposthoc(t(ASV.rel.sub), META.sub$Sample_type)



#####

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
      t(ASV.rel.merge[, rownames(META.adonis)[META.adonis$Site %in% x]]) ~ Site, 
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

### Venn-Diagram ####

# Plotting each site the ASVs of single cell vs sediment vs filters
# to observe if there is an intake from the environmental biota

library("ggVennDiagram")
require(ggplot2)

## to see per specific site
s <- 'Eilat' # # select site to subset in s
METAs <- META.merge[META.merge$Site == s,] # select site s in META
ASV.venn.diagram <- ASV.merge[,rownames(METAs)] # select site s ASV
ASV.venn.diagram <- ASV.venn.diagram[rowSums(ASV.venn.diagram) > 0, ] # remove ASVs that are not present in any sample in this sub setting

# ceeate list for venn diagram
venn.list <- vector("list", length = length(unique(METAs$Sample_type)))
names(venn.list) <- unique(METAs$Sample_type)

for (v in unique(METAs$Sample_type)) { # without unique also works, but take twice the time to run the code
  temp1 <- c(rownames(METAs[METAs$Sample_type == v,]))
  ASV.temp <- ASV.diat.site[,(colnames(ASV.diat.site) %in% temp1)] # selecting ASVs of the same sample type
  ASV.temp <- ASV.temp[apply(ASV.temp, 1, function(x) !all(x==0)),]
  rnames.ASV.temp <- rownames(ASV.temp)
  venn.list[[v]] <- rnames.ASV.temp
}

# or, another option to write the list, setting filters before taking the ASVs names
venn.list <- vector("list", length = 3)
names(venn.list) <- c("single", "sediment", "filter")
tmp.asv <- ASV.diat.site[, METAs$Sample_type == "single_cell"]
tmp.rel <- prop.table(tmp.asv, 2) * 100
venn.list$single <- rownames(tmp.rel)[apply(tmp.rel, 1, function(x) sum(x > 0) >= 1)]
tmp.asv <- ASV.diat.site[, METAs$Sample_type == "sediment"]
tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
venn.list$sediment <- rownames(tmp.asv)
tmp.asv <- ASV.diat.site[, METAs$Sample_type == "filter"]
tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
venn.list$filter <- rownames(tmp.asv)

# Venn Diagram
str(venn.list)
dim(tmp.rel)

ggVennDiagram(venn.list, label_alpha = 0)+
  ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")




#####

### Jaccard boxplot #####
library(reshape)
require(tidyverse)

s <- 'Eilat' # # select site to subset in s
METAs <- META.diat[META.diat$Site == s,] # select site s
ASV.jaccard <- ASV.diat[,rownames(METAs)] # select site s
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
#METAs <- META.diat[META.diat$Site == s,] # select site s
ASV.BC <- ASV.diat[,rownames(METAs)] # select site s
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

### ALDEx2 (does not work) ####

BiocManager::install(c("ALDEx2"))
library("ALDEx2")
packageVersion("ALDEx2")
# '1.24.0'

# Identification of differentially enriched ASVs between treatments
# When looking at trends among individual ASVs, it is most important to account for compositionality effects
# We do this by using the R package ALDEx2

# To prepare the data for ALDEx2, it is necessary to remove rare and low-coverage ASVs,
#   i.e. those not occurring in many samples
# Statistically, those ASVs are also likely to yields false positives or constitute outliers,
#   so removing them should also improve the analysis from this aspect

# Remove ASVs that don't occur with at least 0.1% in 3 samples
ASV.sub <- ASV.merge[rownames(ASV.rel.merge), rownames(META.merge)]
ASV.sub.filt <- ASV.sub[apply(ASV.rel.merge, 1, function(x) sum(x >= 0.1) >= 3), ]
dim(ASV.sub.filt)

# check that filtering did not change beta diversity patterns
plot(
  vegdist(t(ASV.rel.merge)),
  vegdist(prop.table(t(ASV.sub.filt), 1) * 100)
)
# highly correlated

# generate CLR transformed counts
ASV.sub.filt.clr <- aldex.clr(
  ASV.sub.filt,
  META.merge$Site
)

# run ALDEx2 GLM and Kruskal Wallis non-parametric test
aldex.out <- aldex.kw(ASV.sub.filt.clr) # code takes too long to run and never finish

# look at distribution of p-values for parametric tests...
#summary(aldex.out$glm.eBH) # that took ~ 20 min
hist(aldex.out$glm.eBH, breaks = 50)
abline(v = 0.05)
sum(aldex.out$glm.eBH <= 0.05)
# How many ASVs detected as differentially enriched?

# ... and non-parametric tests
summary(aldex.out$kw.eBH)
hist(aldex.out$kw.eBH, breaks = 50)
abline(v = 0.05)
sum(aldex.out$kw.eBH <= 0.05)
# How many ASVs detected as differentially enriched?

# Due to the low sample size, it is possible that p-values are not reliable
# Especially for non-parametric tests, it may not be mathematically possible to obtain 
#   p-values low enough to not be above the significance threshold after correction

# compare parametric and non-parametric results
plot(aldex.out$kw.ep, aldex.out$glm.eBH)
abline(v = 0.05, h = 0.05)

# As a compromise, select those ASVs with...
# ... parametric corrected p-value <= 0.05
# ... non-parametric uncorrected p-value <= 0.05

# Why are both parametric and non-parametric results considered for identifying differentially enriched ASVs?

#####

### heatmap (based on Aldex2)#### 
install.packages("gplots")
library(gplots)
# is this the right package?

# Visualize outcome of differential OTU enrichment analysis as heatmap

# Subset ASV table to those ASV selected as showing strongest effects between treatments
# For better human-readability, we will plot proportions and not log ratios
select.asvs <- rownames(aldex.out)[aldex.out$kw.ep <= 0.05 & aldex.out$glm.eBH <= 0.05]
ASV.rel.diff <- ASV.rel.merge[select.asvs, ]
dim(ASV.rel.diff)

# Also filter by proportion to show only changes among dominant ASVs (here: 1%)
ASV.rel.diff.filt <- ASV.rel.diff[apply(ASV.rel.diff, 1, max) >= 1, ]
dim(ASV.rel.diff.filt)

# Also subset taxonomy table
TAX.diff.filt <- TAX[rownames(ASV.rel.diff.filt), ]

# While it is usually recommended to build a heatmap from sratch using basic plotting elements,
#   we will use a shortcut here (provided by the gplots package)
# For better visibility of small changes, short square-root transformed proportions

gplots::heatmap.2(
  sqrt(ASV.rel.diff.filt),
  labRow = paste0(TAX.diff.filt[, "genus"], " [", select.asvs, "]"),
  col = colorRampPalette(c("grey95", "darkred"))(50),
  margins = c(8, 18),
  trace = "none",
  density.info = "none",
  Colv = F,
  dendrogram = "none",
  keysize = 1,
  key.title = "",
  key.xlab = "Sequence proportion [sqrt %]"
)

# Take a look at the ecological role of the taxa which are changing the most
# How much of the total community is changing?


# testing code directly from ASV.rel.merge, not from aldex.out


gplots::heatmap.2(
  sqrt(ASV.rel.merge),
  labRow = paste0(TAX.merge[, "Family"], " [",  "]"),
  col = colorRampPalette(c("grey95", "darkred"))(50),
  margins = c(8, 18),
  trace = "none",
  density.info = "none",
  Colv = F,
  dendrogram = "none",
  keysize = 1,
  key.title = "",
  key.xlab = "Sequence proportion [sqrt %]"
)

# did not work well

# try with phyloseq package

#####

### heatmap (phyloseq) ######
library("phyloseq")
library("ggplot2")      # graphics
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("patchwork")    # to combine plots later
library("dplyr")

#create phyloseq objects
# for all euks
ASV.phy <- otu_table(ASV.merge, taxa_are_rows = TRUE)
TAX.phy <- tax_table(TAX.merge)
META.phy <- sample_data(META.merge)

# for diatoms only
ASV.phy <- otu_table(ASV.diat, taxa_are_rows = TRUE)
TAX.phy <- tax_table(TAX.diat)
META.phy <- sample_data(META.diat)

ps_euks <- phyloseq(ASV.phy, TAX.phy, META.phy)

ps_euks
# all euks
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 10180 taxa and 144 samples ]
# sample_data() Sample Data:       [ 144 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 10180 taxa by 8 taxonomic ranks ]

# selecting only diatoms
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1502 taxa and 144 samples ]
# sample_data() Sample Data:       [ 144 samples by 29 sample variables ]
# tax_table()   Taxonomy Table:    [ 1502 taxa by 8 taxonomic ranks ]

head(sample_data(ps_euks))

#Normalize number of reads in each sample using median sequencing depth.
total <- median(sample_sums(ps_euks))
standf <- function(x, t=total) round(t * (x / sum(x)))
ps_euks <- transform_sample_counts(ps_euks, standf)


# Drawing heatmap

# firs sub set the most abundant taxa

# selecting only ASVs that represent at least 5% of reads in at least one sample
ps_euks_abund <- filter_taxa(ps_euks, function(x) sum(x > total*0.05) > 0, TRUE)
ps_euks_abund


# basic heatmap plot
plot_heatmap(ps_euks, method = "NMDS", distance = "bray") 
# see documentation
# bray curtins in the proportions not counts

plot_heatmap(ps_euks_abund, method = "NMDS", distance = "bray")

# Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
#   invalid graphics state
dev.off() # to solve error above

# then run code for plot again

# showing by taxonomic group
# for all euks
plot_heatmap(ps_euks_abund, method = "NMDS", distance = "bray", 
             taxa.label = "Class", taxa.order = "Class", sample.label = "Site", sample.order = "Sample_type",
             trans=NULL, low="beige", high="red", na.value="beige")


# for diatoms
plot_heatmap(ps_euks_abund, method = "NMDS", distance = "bray", 
             taxa.label = "Family", taxa.order = "Family", sample.label = "Site", sample.order = "Sample_type",
             trans=NULL, low="beige", high="red", na.value="beige")


#####





##########################################################
################### PlotAbund function ################### 
##########################################################
PlotAbund <- function(relData, abund, 
                      margin = par()$mar,
                      method = c("nmost", "precentage"),
                      colorPalette = c("violet", "purple", "blue", "white", "darkgreen", "yellow", "red", "darkred"),
                      plot.ratio = c(3, 1),
                      open.window = T,
                      save.OTU = F,
                      sample.names = colnames(relData),
                      sort.taxa = F,
                      ...) {
  if(method == "nmost") {
    abund_names <- c()
    for (i in 1:ncol(relData)) {
      abund_names <- unique(c(abund_names, rownames(relData)[order(relData[, i], decreasing = T)][1:abund]))
    }
    abund_rel0 <- relData[abund_names, ]
  }
  if(method == "percentage") {
    abund_rel0 <- relData[apply(relData, 1, function(x) { max(x) >= abund }), ]
  }
  if (sort.taxa == T) {
    abund_rel0 <- abund_rel0[order(rownames(abund_rel0)), ]
  }
  abund_rel <- rbind(abund_rel0, 100 - colSums(abund_rel0))
  rownames(abund_rel)[nrow(abund_rel)] <- "other"
  abund_rel <- as.matrix(abund_rel)
  if (open.window == T) {
    grDevices::windows(width = 20, height = 10) # I had to call the package because my R was not getting the right function
  }
  # par(mfrow = c(1, 2), mar = margin, xpd = NA)
  layout(mat = matrix(c(1, 2), 1, 2, byrow = T), widths = plot.ratio)
  par(mar = margin)
  barplot(abund_rel,
          col = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
          ylim = c(0, 100), 
          las = 2,
          ylab = "Relative sequence abundance [%]",
          names.arg = sample.names,
          cex.names = 0.8,
          ...)
  par(mar = c(1,1,1,1))
  plot.new()
  legend("center", 
         pch = 22, 
         col = "black", 
         pt.bg = rev(c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey")),
         legend = rev(rownames(abund_rel)),
         pt.cex = 1.5, 
         cex = 0.8
  ) 
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(mfrow = c(1,1))
  
  if(save.OTU == T) {
    return(
      data.frame(
        abund_rel,
        color = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
        stringsAsFactors = F
      )
    )
  }
}







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

