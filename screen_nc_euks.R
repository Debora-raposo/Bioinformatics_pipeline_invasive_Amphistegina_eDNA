##################################
# Exploring neg controls - Euks  #
##################################
Sys.setenv(LANG = "en") 
### prepare environment ####
# set working directory
setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_euks")


#load("screen_nc_euks.Rdata")
#save.image("screen_nc_euks.Rdata")

# load packages
require(vegan)
require(scales)
require(tidyverse)

### read data ####

# ASV table
ASV <- read.table(
  "otu_tab_ssu_all_euks_curated.txt", 
  h = T,
  sep = "\t",
  row.names = 1
)

# Metadata
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


# it was necessray to remove some samples 
# NC amplified by mistake
samples_to_remove <-c("PR-1040") # creating vector with samples to remove (NC amplified by mistake)
META <- META[!(row.names(META) %in% samples_to_remove),] # removing row names with the name of the samples
ASV <- ASV[,rownames(META)]
# investigate if there are ASVs with zero sequences after removing these samples
sort(rowSums(ASV)) # to view sum sorted from lower to higher number
# no ASVs with zero sequences! no need to change the TAX table in this case

all.equal(colnames(ASV), rownames(META))
# [1] TRUE
all.equal(rownames(ASV), rownames(TAX))
# [1] TRUE


# calculate proportions
ASV.rel <- prop.table(ASV, 2) * 100



## Statistical analysis

# Cluster diagram 
BC <- vegdist(t(ASV.rel))

# complete linkage (default)
BC.clust <- hclust(BC)
cor(BC, cophenetic(BC.clust))
# [1] 0.95937

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
#yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal)
#yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
#yes


# average linkage
BC.clust <- hclust(BC, method = "average")
cor(BC, cophenetic(BC.clust))
# [1] 0.9798066 # best method

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
#yes
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal) 
#yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
#yes

# ward linkage
BC.clust <- hclust(BC, method = "ward.D2")
cor(BC, cophenetic(BC.clust))
# [1] 0.7050669 # not a good representation of the original dissimilarities

plot(BC.clust, labels = META$Sample_type, cex = 0.7) # to see if NCs are randomly distributed (ideal)
#no! NCs are clustered together in this linkage method
plot(BC.clust, labels = META$Extraction_Voucher, cex = 0.7) # to see if technical replicates are together (ideal) 
#yes
plot(BC.clust, labels = META$PCR, cex = 0.7) # to see if samples from same PCR are randomly distributed (ideal)
#no! PCRs seems to be clustered together

 

# investigate dissimilarities by PCR batch (11 samples + 1 NC)


# subset dataset to just work with same PCR
# all PCRs possible 
META$PCR <- as.factor(META$PCR)
levels(META$PCR)
# [1] "PCR-S074" "PCR-S100" "PCR-S102" "PCR-S106" "PCR-S108" "PCR-S110" "PCR-S112" "PCR-S118" "PCR-S120" "PCR-S122" "PCR-S124" "PCR-S126" "PCR-S128"
# [15] "PCR-S130" "PCR-S134" "PCR-S78"  "PCR-S82"  "PCR-S84"  "PCR-S86"  "PCR-S88"  "PCR-S90"  "PCR-S94"  "PCR-S98"

pdf("dissimilarities_average_PCR_euks.pdf", width = 7, height = 8, onefile = T)
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
    plot(BC.clust, labels = paste(META.PCR1$Sample_type,META.PCR1$Extraction_Voucher), cex = 0.7, main = i) # to see if NC is dissimilar from the samples (ideal)
 
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
      temp <- t(cbind(ASV.NC1, ASV.rel[, s])) # binding by column and creating a transposed table to be in the shape for vegdist
      temp <- temp[, colSums(temp) > 0] # ignoring ASVs that dont occurs in the samples we are comparing (reduce the calculation requirements)
      NC_pairwise[[i]][[rownames(META.PCR1)[n]]] <- c(NC_pairwise[[i]][[rownames(META.PCR1)[n]]], vegdist(temp)) # to find the rowname (NC_extraction) according to the row number n
    }
    names(NC_pairwise[[i]][[rownames(META.PCR1)[n]]]) <- rownames(META.PCR1)[META.PCR1$Sample_NC == "Sample"] # save the name of elements
  }
}

# see the minimum dissimilarities for each PCR batch
#list apply lapply
lapply(NC_pairwise, function(x){
  sapply(x, min) # sapply simplifies de output of the inner list
})

NC_pairwise$`PCR-S134`
# $`PR-1250` # NC
# PR-1242   PR-1243   PR-1240   PR-1241   PR-1236   PR-1237   PR-1244   PR-1246   PR-1245   PR-1247   PR-1248   PR-1238   PR-1249   PR-1239 
# 0.9697972 0.9703625 0.9572631 0.9584791 0.9538770 0.9509242 0.9161347 0.9397924 0.9134775 0.9398344 0.8633405 0.4622106 0.8652955 0.4988943 
### PR-1238 and PR-1239 are problematic

NC_pairwise$`PCR-S78`
# $`PR-0624`
# PR-0606   PR-0604   PR-0607   PR-0605   PR-0608   PR-0609   PR-0614   PR-0615   PR-0612   PR-0610   PR-0611   PR-0613   PR-0616   PR-0617   PR-0622   PR-0618 
# 0.8694478 0.8230440 0.8654556 0.8232387 0.5984669 0.6037479 0.7860039 0.7820872 0.8256205 0.8050043 0.8018250 0.8071423 0.8736349 0.8713768 0.8377466 0.9063099 
# PR-0620   PR-0619   PR-0621   PR-0623 
# 0.9054207 0.9085224 0.9067833 0.8343758 
### PR-0608 and PR-0609 are problematic

NC_pairwise$`PCR-S82`
# $`PR-0660`
# PR-0648   PR-0652   PR-0649   PR-0653   PR-0650   PR-0651   PR-0658   PR-0659   PR-0654   PR-0656   PR-0655   PR-0657 
# 0.6617718 0.6515615 0.6613881 0.6514887 0.9594049 0.9629418 0.9987503 0.9976634 0.6432839 0.6343754 0.6330457 0.6272073
### PR-0648, PR-0652, PR-0649, PR-0653, PR-0654, PR-0656, PR-0655, PR-0657 are problematic

NC_pairwise$`PCR-S84`
# $`PR-0692`
# PR-0688    PR-0690    PR-0689    PR-0691    PR-0680    PR-0681    PR-0678    PR-0679    PR-0686    PR-0687    PR-0682    PR-0683    PR-0684    PR-0685 
# 0.99542142 0.95615625 0.99467209 0.95868313 0.93646356 0.93753287 0.89813644 0.91120469 0.09199655 0.07233529 0.79505929 0.79282217 0.48517186 0.47706875 
# PR-0676    PR-0677 
# 0.98950691 0.99012796
###  PR-0687 is somewhat problematic, but close to 0.8
### PR-0684 and PR-0685 are problematic

NC_pairwise$`PCR-S88`
# $`PR-0750`
# PR-0732   PR-0733   PR-0738   PR-0739   PR-0736   PR-0735   PR-0737   PR-0734   PR-0744   PR-0746   PR-0745   PR-0747   PR-0742   PR-0743   PR-0740   PR-0741 
# 0.9944029 0.9914476 0.9930090 0.9893759 0.5789053 0.9968544 0.6120320 0.9963521 0.9933054 0.9850920 0.9884999 0.9781471 0.5549959 0.5286585 0.5843360 0.5810379 
### PR-0736, PR-0737, PR-0742, PR-0743, PR-0740, PR-0741 are problematic

NC_pairwise$`PCR-S90`
# $`PR-0790`
# PR-0780   PR-0781   PR-0784   PR-0782   PR-0785   PR-0783   PR-0778   PR-0779   PR-0776   PR-0777   PR-0788   PR-0786   PR-0789   PR-0787   PR-0770   PR-0771 
# 0.4192687 0.4141188 0.9830532 0.5583104 0.9737703 0.5600334 0.9776692 0.9858800 0.5265999 0.5263060 0.5566847 0.4656797 0.5575210 0.4663329 0.8806385 0.8796387 
# PR-0774   PR-0775 
# 0.9552567 0.9174148 
### PR-0780, PR-0781, PR-0782, PR-0783, PR-0776, PR-0777, PR-0788, PR-0786, PR-0789, PR-0787 are problematic



### to do
### look gel - if it doesnt look good (NC reacting, different bands for technical replicates), remove samples
# gel looked good. so there is no need to remove samples in this step

### look ASVs present in the NCs that are also in the samples






## identify what the factors that have major influence in the dissimilarites between the technical replicates
# do correlations of the dissimilarities with:
#   DNA yield
#  site
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
#                 Df  Sum Sq Mean Sq F value   Pr(>F)
# DNA_yield_dif    1 0.00146 0.00146   0.681  0.41067    
# DNA_yield_mean   1 0.00000 0.00000   0.000  0.99605    
# Site             3 0.03484 0.01161   5.432  0.00148 ** 
# Depth            3 0.01444 0.00481   2.251  0.08532 .  
# Sample_type      2 0.10354 0.05177  24.216 1.08e-09 ***
# Residuals      133 0.28433 0.00214                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# removing variables that are not significant
summary(aov(PCR_BC ~ Site + Sample_type, META_per_sample))
#               Df  Sum Sq Mean Sq F value   Pr(>F)    
# Site          3 0.03484 0.01161   5.296  0.00174 ** 
# Sample_type   2 0.10116 0.05058  23.066 2.28e-09 ***
# Residuals   138 0.30261 0.00219                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# to look only single cell
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type == "single_cell",]))
#                 Df  Sum Sq  Mean Sq F value Pr(>F)  
# DNA_yield_dif    1 0.00037 0.000370   0.153 0.6966  
# DNA_yield_mean   1 0.00271 0.002712   1.121 0.2922  
# Site             3 0.02582 0.008608   3.558 0.0169 *
# Depth            3 0.01843 0.006144   2.540 0.0606 .
# Residuals      103 0.24917 0.002419                 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# to look only env samples
summary(aov(PCR_BC ~ DNA_yield_dif+ DNA_yield_mean + Site + Depth, META_per_sample[META_per_sample$Sample_type != "single_cell",]))
#                Df   Sum Sq   Mean Sq F value Pr(>F)  
# DNA_yield_dif   1 0.001159 0.0011595   1.456 0.2398  
# DNA_yield_mean  1 0.000841 0.0008405   1.056 0.3149  
# Site            3 0.004206 0.0014021   1.761 0.1827  
# Depth           3 0.006363 0.0021210   2.664 0.0718 .
# Residuals      23 0.018314 0.0007963                 
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



# is the DNA yield dependent of the sample type?
boxplot(DNA_yield_mean ~ Sample_type, data = META_per_sample)


# removing variables that are not significant
summary(aov(DNA_yield_mean ~ Sample_type, META_per_sample))
#              Df Sum Sq Mean Sq F value   Pr(>F)    
# Sample_type   2  17.03   8.514   7.935 0.000543 ***
# Residuals   141 151.30   1.073                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 


####### Another way to calculate this table ######


# extract output of vegdist as matrix
distmat <- as.matrix(BC)

# Combine the distance between the technical replicates (row names and column names of distmat & row names of META) 
#that come from the same sample (in META$Extraction_Voucher)

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
#[1] "S0465" "S0465" "S0633" "S0657" "S0633" "S0657"

# those are actually 3 NCs that happened to be sequenced twice (different PCRs). We can remove for this analysis.
distpairs <- distpairs[distpairs$col != "S0465" & distpairs$col != "S0633" & distpairs$col != "S0657", ]

distpairs <- as.data.frame(distpairs)
all.equal(METAnoNC$Extraction_Voucher, distpairs$col)
#[1] TRUE
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
# we need to fix this number of levels
distpairs$sample_type <- droplevels(distpairs$sample_type)
levels(distpairs$sample_type)
#[1] "filter"      "sediment"    "single_cell"
# now is correct

str(distpairs)


# before continuing, I should first average (or sum) the technical replicates

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



##### Remove sample S0512? #######
# S0512 seems like a failed PCR in the proks and must be removed there. Shall I also remove here in the euks datasets? 
#(replicates are with very different DNA yield in Proks and NC slightly reacted)
samples_to_remove <-c("PR-0686","PR-0687") # creating vector with samples to remove (two technical replicates of S0512)
META.clean <- META[!(row.names(META) %in% samples_to_remove),] # removing row names with the name of the samples
ASV.clean <- ASV[,rownames(META.clean)] # removing these samples also from ASV table
# investigate if there are ASVs with zero sequences after removing these samples
sort(rowSums(ASV.clean)) # to view sum sorted from lower to higher number
# there are two ASVs with zero sequences!
sum(rowSums(ASV.clean) == 0)
# [1] 2
# remove these ASVs which has zero sequences
ASV.clean2 <- ASV.clean[rowSums(ASV.clean) > 0,]
# clean TAX table by removing these ASVs too
TAX.clean <- TAX[rownames(ASV.clean2),]

all.equal(colnames(ASV.clean2), rownames(META.clean))
# [1] TRUE
all.equal(rownames(ASV.clean2), rownames(TAX.clean))
# [1] TRUE
#######


