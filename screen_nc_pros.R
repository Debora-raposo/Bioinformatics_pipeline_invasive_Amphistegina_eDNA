#################################
# Exploring neg controls -Prok  #
#################################

### prepare environment ####
# set working directory
setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_pros")

#load("screen_nc_pros.Rdata")
#save.image("screen_nc_pros.Rdata")

# load packages
require(vegan)
require(scales)
require(tidyverse)


### read data ####

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
META <- META[META$Purified_PCR_product != "PR-1178",] # removing sample that was excluded in data_curation script (no sequences left after filtering)

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
# [1] "filter"         "sediment"       "single_cell"    "NC_env"         "NC_single_cell"

# reorder data
META <- META[order(META$Site, META$Condition), ]
META <- META %>% column_to_rownames("Purified_PCR_product") # had to do it first in my data
ASV <- ASV[, rownames(META)]

# check that table are in the correct order
all.equal(rownames(ASV), rownames(TAX)) 
# TRUE  
all.equal(colnames(ASV), rownames(META)) 
# TRUE

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

# Cluster diagram 
BC <- vegdist(t(ASV.rel))

# complete linkage (default)
BC.clust <- hclust(BC)
cor(BC, cophenetic(BC.clust))
# [1] 0.7707682

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
#yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal)
#yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
#yes


# average linkage
BC.clust <- hclust(BC, method = "average")
cor(BC, cophenetic(BC.clust))
# [1] 0.9091584 # best method

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
#yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal) 
#yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
#yes

# ward linkage
BC.clust <- hclust(BC, method = "ward.D2")
cor(BC, cophenetic(BC.clust))
# [1] 0.6153868
# not a good representation of the original dissimilarities

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
#yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal) 
#yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
#yes




# investigate dissimilarities by PCR batch (11 samples + 1 NC)


# subset dataset to just work with same PCR
# all PCRs possible 
META$PCR <- as.factor(META$PCR)
levels(META$PCR)
# [1] "PCR-S077" "PCR-S101" "PCR-S103" "PCR-S107" "PCR-S109" "PCR-S111" "PCR-S113" "PCR-S115" "PCR-S119" "PCR-S121" "PCR-S123" "PCR-S125" "PCR-S127" "PCR-S129"
# [15] "PCR-S131" "PCR-S135" "PCR-S79"  "PCR-S83"  "PCR-S85"  "PCR-S87"  "PCR-S89"  "PCR-S91"  "PCR-S95"  "PCR-S99" 

pdf("dissimilarities_average_PCR_proks.pdf", width = 7, height = 8, onefile = T)
par(mar=c(7,4,4,1))

for (i in levels(META$PCR)) {
  META.PCR1 <- META[META$PCR == i,] # change name of PCR according to levels above and run all lines until line 145... 
  # to check if NCs are dissimilar from tthe samples (meaning, if they can be excluded) or not
  
  ASV.PCR1 <- ASV[, rownames(META.PCR1)]
  
  # calculate proportions
  ASV.rel <- prop.table(ASV.PCR1, 2) * 100
  BC <- vegdist(t(ASV.rel))
  
  # complete linkage (default)
  BC.clust <- hclust(BC, method = "average")
  plot(BC.clust, labels = paste(META.PCR1$Sample_type, META.PCR1$Extraction_Voucher), cex = 0.7, main = i) # to see if NC is dissimilar from the samples (ideal)
}
dev.off()


# Calculate pairwise distances between dissimilarities of NC (or NCs) and every sample for each PCR batch

NC_pairwise <- vector("list", length = length(levels(META$PCR))) # create a list to object to save the output. each PCR batch become a line in the list
names(NC_pairwise) <- levels(META$PCR) # save elements of list with the name of the PCR
for(i in levels(META$PCR)) { # i is each PCR
  META.PCR1 <- META[META$PCR == i, ] # create META table with just one PCR
  ASV1 <- ASV[, rownames(META.PCR1)] # subset ASV table to just show samples from the same PCR
  ASV.rel <- prop.table(as.matrix(ASV1), 2) * 100 # calculate relative proportions
  NC_pairwise[[i]] <- vector("list", length=sum(META.PCR1$Sample_NC == "NC_extraction")) # create another list to divide the lines in the list according to the number of NCs
  names(NC_pairwise[[i]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "NC_extraction"] # save this sub-lines with the name of the NC
  for(n in which(META.PCR1$Sample_NC == "NC_extraction")) { # n is the row number in which a NC is located in the META table for the PCR batch
    ASV.NC1 <- ASV.rel[, n] # subset ASV table to just show the NC
    for(s in which(META.PCR1$Sample_NC == "Sample")) { # s is the row number in which a sample is located in the META table for the PCR batch
      temp <- t(cbind(ASV.NC1, ASV.rel[, s])) # creating a transposed table to be in the shape for vegdist
      temp <- temp[, colSums(temp) > 0] # ignoring ASVs that dont occurs in the samples we are comparing (reduce the calculation requirements)
      NC_pairwise[[i]][[rownames(META.PCR1)[n]]] <- c(NC_pairwise[[i]][[rownames(META.PCR1)[n]]], vegdist(temp)) # to find the rowname (NC_extraction) according to the row number n
    }
    names(NC_pairwise[[i]][[rownames(META.PCR1)[n]]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "Sample"] # save the name of elements
  }
}

# see the minimum dissimilarities for each PCR batch
#list apply lapply
lapply(NC_pairwise, function(x){
 sapply(x, min) # simplifies de output of the inner list
})

NC_pairwise$`PCR-S125`$`PR-1134`
# PR-1137   PR-1128   PR-1138   PR-1129   PR-1135   PR-1126   PR-1136   PR-1127   PR-1130   PR-1139   PR-1140   PR-1131   PR-1132   PR-1141   PR-1133   PR-1142 
# 0.7991401 0.7978159 0.7990998 0.7965212 0.8213431 0.8218273 0.8198843 0.8193901 0.5815874 0.5733967 0.5743956 0.5842579 0.7841393 0.8399092 0.7834158 0.8385638 
## PR-1130, PR-1131, PR-1139, PR-1140 are problematic

NC_pairwise$`PCR-S127`$`PR-1175`
# PR-1171   PR-1172   PR-1164   PR-1169   PR-1165   PR-1170   PR-1173   PR-1166   PR-1167   PR-1174   PR-1179   PR-1176   PR-1177   PR-1180 
# 0.6544979 0.6509770 0.6818738 0.6916317 0.6864243 0.6906137 0.7352490 0.7932862 0.7949614 0.7370530 0.7886442 0.7886151 0.7931141 0.7897282
## PR-1164, PR-1165, PR-1171, PR-1172, PR-1169, PR-1170 are problematic    PR-1173, PR-1174 are okayish

NC_pairwise$`PCR-S129`$`PR-1208`
# PR-1202   PR-1203   PR-1211   PR-1212   PR-1198   PR-1199   PR-1209   PR-1210   PR-1204   PR-1205   PR-1213   PR-1214 
# 0.9395673 0.9426239 0.9523574 0.9526441 0.6709149 0.6687939 0.9368552 0.9377753 0.9565789 0.9575527 0.9586702 0.9597585 
## PR-1198 and PR-1199 are problematic


# the gels of these 3 PCRs look all fine. 


## next step: investigate the ASVs that are in the NC and in the samples









## identify what is causing the dissimilarites with the technical replicates
# do correlations of the dissimilarities with:
#   DNA yield
#   Sample type
# and other variables 

require(reshape)
META_per_sample <- unique(META[META$Sample_NC != "NC_extraction",c("Extraction_Voucher", "Site", "Depth", "Substrate", "Sample_type")]) 

rownames(META_per_sample) <- META_per_sample$Extraction_Voucher
ASV.rel <- prop.table(as.matrix(ASV), 2) * 100 
META_per_sample$PCR_BC <- c(NA) # create empty column to add data of BC dissimilarities between replicates
META_per_sample$DNA_yield_mean <- c(NA) # create empty column to add data of DNA concentration mean between replicates
META_per_sample$DNA_yield_dif <- c(NA) # create empty column to add data of DNA concentration difference between replicates
for (i in META_per_sample$Extraction_Voucher){ # selecting sample by sample in a loop
 ASVi <- t(ASV.rel[,META$Extraction_Voucher == i]) # sub-setting ASV.rel to just show the ones regarding the samples we want to analyse (technical replicates)
 ASVi <- ASVi[,colSums(ASVi)>0] # to ignore ASVs that are not present in the samples analyzed
 META_per_sample[i,"PCR_BC"] <- max(c(vegdist(ASVi))) # calculating the BC dissimilarities between the technical replicates and adding the result to the column "PCR_BC" in the line corresponding to the sample we are analyzing.
 #using max to have a more universal code in case one has more than two PCR replicates
 META_per_sample[i,"DNA_yield_mean"] <- mean(META$DNA_concentration[META$Extraction_Voucher == i]) # calculating the mean of DNA conc among the technical replicates
 META_per_sample[i,"DNA_yield_dif"] <- max(c(dist(META$DNA_concentration[META$Extraction_Voucher == i]))) # calculating the difference between the DNA conc of the technical replicates
 #using max to make a more universal code in case one has more than two PCR replicates
}

str(META_per_sample)
META_per_sample$Sample_type <- droplevels(META_per_sample$Sample_type)

# transform variables that are not yet as factor
META_per_sample$Site <- as.factor(META_per_sample$Site)
META_per_sample$Depth <- as.factor(META_per_sample$Depth)
META_per_sample$Substrate <- as.factor(META_per_sample$Substrate)

# visualize how factor affects distance between technical replicates
boxplot(PCR_BC ~ Sample_type, data = META_per_sample)

# visualize the effect of DNA yield in the distance between technical replicates
plot(META_per_sample$DNA_yield_mean, META_per_sample$PCR_BC, col=as.numeric(META_per_sample$Sample_type), pch=16)

plot(META_per_sample$DNA_yield_mean, META_per_sample$PCR_BC, col=as.numeric(META_per_sample$Site), pch=16)

plot(META_per_sample$DNA_yield_dif, META_per_sample$PCR_BC, col=as.numeric(META_per_sample$Sample_type), pch=16)

plot(META_per_sample$DNA_yield_dif, META_per_sample$PCR_BC, col=as.numeric(META_per_sample$Site), pch=16)


# Compute the analysis of variance
summary(aov(PCR_BC ~ DNA_yield_dif + DNA_yield_mean + Site + Depth + Sample_type, META_per_sample))
#                 Df Sum Sq Mean Sq F value   Pr(>F)    
# DNA_yield_dif    1 0.0285  0.0285   1.425   0.2347    
# DNA_yield_mean   1 0.5507  0.5507  27.556 5.88e-07 ***
# Site             3 0.2260  0.0753   3.769   0.0123 *  
# Depth            3 0.0237  0.0079   0.395   0.7569    
# Sample_type      2 2.5153  1.2577  62.934  < 2e-16 ***
# Residuals      133 2.6579  0.0200                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 

# removing variables that are not significant
summary(aov(PCR_BC ~ DNA_yield_mean + Site + Sample_type, META_per_sample))
#                  Df Sum Sq Mean Sq F value   Pr(>F)    
# DNA_yield_mean   1 0.4076  0.4076  19.976 1.63e-05 ***
# Site             3 0.2334  0.0778   3.813   0.0116 *  
# Sample_type      2 2.5656  1.2828  62.871  < 2e-16 ***
# Residuals      137 2.7954  0.0204                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# to look only single cell
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type == "single_cell",]))
#                 Df Sum Sq Mean Sq F value Pr(>F)   
# DNA_yield_dif    1 0.0698 0.06979   2.853 0.0942 . 
# DNA_yield_mean   1 0.1801 0.18014   7.365 0.0078 **
# Site             3 0.2665 0.08884   3.632 0.0154 * 
# Depth            3 0.0392 0.01307   0.534 0.6598   
# Residuals      103 2.5192 0.02446                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# to look only env samples
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type != "single_cell",]))
#                Df  Sum Sq Mean Sq F value  Pr(>F)    
# DNA_yield_dif   1 0.02503 0.02503   8.296 0.00845 ** 
# DNA_yield_mean  1 0.10608 0.10608  35.164 4.8e-06 ***
# Site            3 0.00404 0.00135   0.447 0.72206    
# Depth           3 0.00649 0.00216   0.717 0.55203    
# Residuals      23 0.06939 0.00302                    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



# is the DNA yield dependent of the sample type?
boxplot(DNA_yield_mean ~ Sample_type, data = META_per_sample)


# removing variables that are not significant
summary(aov(DNA_yield_mean ~ Sample_type, META_per_sample))
#               Df Sum Sq Mean Sq F value   Pr(>F)    
# Sample_type   2  23.95  11.973   7.629 0.000714 ***
# Residuals   141 221.27   1.569                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1





###### Another way to calculate this table #######

# extract output of vegdist as matrix
distmat <- as.matrix(BC)

# I have to combine the distance between the technical replicates (row names and column names of distmat & row names of META) 
#that come from the same sample (in META$Extraction_Voucher)

# create a table with the following columns
#Extraction_Voucher/ dist_between_replicates / DNA yield / etc

# rename columns and rows as the name of the sample (extraction voucher)
distmat.renamed <- distmat
rownames(distmat.renamed) <- META$Extraction_Voucher
colnames(distmat.renamed) <- META$Extraction_Voucher

# create data.frame with three columns (pair A, pair B and distance)
dist <- data.frame(col=colnames(distmat.renamed)[col(distmat.renamed)], 
                   row=rownames(distmat.renamed)[row(distmat.renamed)], 
                   dist=c(distmat.renamed)) 
# this gives the entire pairwise list, but I want just the pairs of same sample 

# keeping only pairs of same sample 
distpairs <- subset(dist, col==row)

# removing pairwise comparison with same technical replicate (always zero)
distpairs <- distpairs[distpairs$dist > 0, ]

# create META table without NC so it is possible to compare with the distpairs table
# necessary since the NCs are removed when we keep only samples that have techical replicates
METAnoNC <- META[META$Sample_NC == "Sample", ]
nrow(METAnoNC)#288
nrow(distpairs)#294
# six extra samples in distpairs?? which ones?

# to see samples that are in distpairs and not in METAnoNC
distpairs$col[!(distpairs$col %in% METAnoNC$Extraction_Voucher)]
#[1] "S0465" "S0633" "S0657" "S0465" "S0657" "S0633"

# those are actually 3 NCs that happened to be sequenced twice (different PCRs). We can remove for this analysis.
distpairs <- distpairs[distpairs$col != "S0465" & distpairs$col != "S0633" & distpairs$col != "S0657", ]

distpairs <- as.data.frame(distpairs)
all.equal(METAnoNC$Extraction_Voucher, distpairs$col)
# [1] TRUE
rownames(distpairs) <- rownames(METAnoNC)

distpairs$dna_yield <- METAnoNC$DNA_concentration
distpairs$site <- METAnoNC$Site
distpairs$depth <- METAnoNC$Depth
distpairs$substrate <- METAnoNC$Substrate
distpairs$sample_type <- METAnoNC$Sample_type

# transform as factor
distpairs$depth <- as.factor(distpairs$depth)
distpairs$substrate <- as.factor(distpairs$substrate)
levels(distpairs$substrate)
#[1] "Algae"   "bottom"  "Rubbles" "surface"
levels(distpairs$sample_type)
#[1] "filter"         "sediment"       "single_cell"    "NC_env"         "NC_single_cell"
# these are not the right number of levels, there shouldn't have NCs
distpairs$sample_type <- droplevels(distpairs$sample_type)
levels(distpairs$sample_type)
#[1] "filter"      "sediment"    "single_cell"
# ok, number of levels are fixed

str(distpairs)


# however, before analysing,  I should first average the technical replicates for a better statistical approach

# removing duplicated column with repeated name of sample
distpairs_aver <- subset(distpairs, select = -c(row)) 

# calculating the mean of dna yield (between technical replicates)
distpairs_aver  <- distpairs_aver %>%
  unite(col="Experiment", col, site, depth, substrate, sample_type, sep="-") %>%
  group_by(Experiment) %>%
  summarise_all(funs(n(), mean)) %>%
  separate(Experiment, c("col", "site", "depth", "substrate", "sample_type"), sep="-") %>%
  select(-"dist_n", -"dna_yield_n") 

str(distpairs_aver)

# transform as factor
distpairs_aver$site <- as.factor(distpairs_aver$site)
distpairs_aver$depth <- as.factor(distpairs_aver$depth)
distpairs_aver$substrate <- as.factor(distpairs_aver$substrate)
distpairs_aver$sample_type <- as.factor(distpairs_aver$sample_type)
#######



### Next steps 

## Remove sample S0512 - failed PCR (replicates are with very different DNA yield and NC slightly reacted)

# subsetting data sets
samples_to_remove <-c("PR-0703","PR-0704") # creating vector with samples to remove (two technical replicates)
META.clean <- META[!(row.names(META) %in% samples_to_remove),] # removing row names with the name of the samples
ASV.clean <- ASV[,rownames(META.clean)] # removing these samples also from ASV table
# investigate if there are ASVs with zero sequences after removing these samples
sort(rowSums(ASV.clean)) # to view sum sorted from lower to higher number
# there are some ASVs with zero sequences!
sum(rowSums(ASV.clean) == 0)
# [1] 163
# remove these ASVs which has zero sequences
ASV.clean2 <- ASV.clean[rowSums(ASV.clean) > 0,]
# clean TAX table by removing these ASVs too
TAX.clean <- TAX[rownames(ASV.clean2),]

all.equal(colnames(ASV.clean2), rownames(META.clean))
# [1] TRUE
all.equal(rownames(ASV.clean2), rownames(TAX.clean))
# [1] TRUE



# Remove ASVs that don't occur with at least 0.1% in 3 samples
dim(ASV)
#[1] 38566   317

ASV.filt <- ASV[apply(ASV.rel, 1, function(x) sum(x >= 0.1) >= 3), ]
dim(ASV.filt)
#[1] 2927  317

# or using  0.01%
# [1] 12777   317

# or using 0.1% in at least 2 samples
# [1] 4651  317


# making sure META and ASV.filt still match
all.equal(rownames(META), colnames(ASV.filt))

# subsetting TAX accoridng to filtered ASV table
TAX.filt <- TAX[rownames(ASV.filt),]
all.equal(rownames(TAX.filt), rownames(ASV.filt))

# calculate proportions
ASV.rel.filt <- prop.table(ASV.filt, 2) * 100









