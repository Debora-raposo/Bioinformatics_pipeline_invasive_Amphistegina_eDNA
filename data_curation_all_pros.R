#################################################
#                                               #
#  Data curation    SSU     all      pros       #
#                                               #
#################################################



############# Working on local drive ######

setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_pros")
#save.image("data_curation_all_pros.RData")
#load("data_curation_all_pros.RData")
require(vegan)
require(readxl)

# load OTU and TAX tables
OTU <- read.table("otu_tab_ssu_all_pros.txt", h = T, sep = "\t")
TAX <- read.table("tax_tab_ssu_all_pros.txt", h = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
TAX.boot <- read.table("tax_bootstrap_ssu_all_pros.txt", h = T, sep = "\t")
META <- data.frame(read_xlsx("Lists_bioinformatic_all_pros.xlsx"))
#NSEQ <- read.table("nSeq_dada2_ssu_all_pros.txt", h = T, sep = "\t")

# setting Purified_PCR_product as rownames to META, using gsub to change from _ to . so it matches the OTU table
rownames(META) <- gsub("_", ".", META$Purified_PCR_product)

# to set in same order as occurs as columns in OTU
META <- META[colnames(OTU), ]

all.equal(rownames(META), colnames(OTU)) # [1] TRUE



# Analyzing bootstrap ranges

apply(TAX.boot, 2, summary)
# Kingdom    Phylum    Class     Order    Family     Genus
# Min.     10.00000   5.00000   5.0000   2.00000   2.00000   1.00000
# 1st Qu. 100.00000 100.00000 100.0000  99.00000  98.00000  71.00000
# Median  100.00000 100.00000 100.0000 100.00000 100.00000  97.00000
# Mean     99.99101  98.80495  98.3601  95.05006  93.20681  83.08141
# 3rd Qu. 100.00000 100.00000 100.0000 100.00000 100.00000 100.00000
# Max.    100.00000 100.00000 100.0000 100.00000 100.00000 100.00000

apply(TAX.boot,2, function(x) quantile(x, seq(0, 1, 0.05)))
#       Kingdom Phylum Class Order Family Genus
# 0%        10      5     5     2      2     1
# 5%       100     96    92    62     52    31
# 10%      100    100    99    84     74    43
# 15%      100    100   100    94     87    52
# 20%      100    100   100    98     94    61
# 25%      100    100   100    99     98    71
# 30%      100    100   100   100     99    79
# 35%      100    100   100   100    100    86
# 40%      100    100   100   100    100    91
# 45%      100    100   100   100    100    95
# 50%      100    100   100   100    100    97
# 55%      100    100   100   100    100    99
# 60%      100    100   100   100    100   100
# 65%      100    100   100   100    100   100
# 70%      100    100   100   100    100   100
# 75%      100    100   100   100    100   100
# 80%      100    100   100   100    100   100
# 85%      100    100   100   100    100   100
# 90%      100    100   100   100    100   100
# 95%      100    100   100   100    100   100
# 100%     100    100   100   100    100   100

# 70 seems to be a okay cutoff. we will loose around 25% of ASVs, but it is okay.

# apllying bootstrap
TAX.clean <- as.matrix(TAX[, -7])
TAX.clean[is.na(TAX.clean)] <- "NA"
TAX.clean[TAX.boot < 70] <- "NA" # cutoff observed above in the tables
TAX.clean <- TAX.clean[rowSums(OTU) >= 2 & TAX.clean[, 1] == "Bacteria" & TAX.clean[, 4] != "Chloroplast" & TAX.clean[, 5] != "Mitochondria", ]
#the >=2 means that we only keep ASVs that occurs at least twice, this way we exclude singletons
OTU.clean <- OTU[rownames(TAX.clean), ] # select only rownames that were kept in Tax.clean and apply it as OTU.clean
dim(OTU.clean)# [1] 38566   318
nrow(OTU.clean)/nrow(OTU) # 94.2%
summary(colSums(OTU.clean)/colSums(OTU))
# Min.     1st Qu.  Median    Mean   3rd Qu.    Max. 
# 0.00000 0.07888   0.16770 0.35027 0.56834    0.99925 

# check if there are samples with no sequences left. sorting in ascending order
sort(colSums(OTU.clean)/colSums(OTU))
# sample  PR.1178 has no sequences

# View only the META from samples with less than 5% of sequences kept
# this code only works because META was adjusted in the same order of OTU previously (line 27)
View(META[colSums(OTU.clean)/colSums(OTU) < 0.05, ])

# Order in ascending order and view the first 50 (percentage of sequences kept)
View(META[order(colSums(OTU.clean)/colSums(OTU))[1:50], ])

# Order in ascending order and view the first 50 (absolute number of sequences kept)
View(META[order(colSums(OTU.clean))[1:50], ])

# subsetting to remove samples with zero sequences (PR.1178)
OTU.clean2 <- OTU.clean[, colSums(OTU.clean) > 0 ]
META.clean <- META[colnames(OTU.clean2), ]



# calculate rarefaction curves
plot.in <- OTU.clean2
maxIndex <- colSums(decostand(plot.in, "pa"))
plot(
  0, 0,
  type = "n",
  xlim = c(0, max(colSums(plot.in))),
  ylim = c(0, max(maxIndex)),
  ylab = "Number of OTUs (Proks)",
  xlab = "Sequencing depth",
  cex.axis = 0.7,
  las = 1,
  axes = F,
  cex.lab = 0.9,
  mgp = c(1.8, 0.4, 0)
)
axis(1, las = 1, mgp = c(1, 0.2, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
axis(2, las = 1, mgp = c(2, 0.4, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
box("plot", lwd = 0.5)
for(k in order(colSums(plot.in), decreasing = T)) {
  temp <- rarefy(
    plot.in[, k],
    sample = round(seq(0, sum(plot.in[, k]), length.out = 100))
  )
  lines(
    round(seq(0, sum(plot.in[, k]), length.out = 100)),
    temp,
    lwd = 1, col = META.clean$Color[k]
  )
}


# Looks okay


#Get nSeqs summary
track <- colSums(OTU.clean)
track <- data.frame(track)
colnames(track) <- c("tabled")






# Curating taxonomic assignments
# saving as "unclassified" instead of "NA"

# set the NA as characters to NA 
TAX.temp <- TAX.clean
TAX.temp[TAX.temp == "NA"] <- NA


#################### function curate_taxpath ######################

#################### function 
curate_taxpath <- function(tax) {
  out <- tax
  out[out == "uncultured"] <- NA
  k <- ncol(out) - 1
  for (i in 2:k) {
    if (sum(is.na(out[, i])) > 1) {
      test <- out[is.na(out[, i]), ]
      for (j in 1:nrow(test)) {
        if (sum(is.na(test[j, i:(k + 1)])) == length(test[j, i:(k + 1)])) {
          test[j, i] <- paste(test[j, (i - 1)], "_unclassified", sep = "")
          test[j, (i + 1):(k + 1)] <- test[j, i]
        }
      }
      out[is.na(out[, i]), ] <- test
    }
    if (sum(is.na(out[, i])) == 1) {
      test <- out[is.na(out[, i]), ]
      if (sum(is.na(test[i:(k + 1)])) == length(test[i:(k + 1)])) {
        test[i] <- paste(test[(i - 1)], "_unclassified", sep = "")
        test[(i + 1):(k + 1)] <- test[i]
      }
      out[is.na(out[, i]),] <- test
    }
  }
  out[is.na(out[, (k + 1)]), (k + 1)] <- paste(out[is.na(out[, (k + 1)]), k], "_unclassified", sep = "")
  return(out)
}


###

TAX.curated <- curate_taxpath(TAX.temp) # to apply the function
# the TAX.curated has the names saved as "unclassified" instead of "NA"


# Exporting tables

write.table(track, "nSeq_dada2_ssu_all_pros_curated.txt", quote = F, sep = "\t")
write.table(OTU.clean2, "otu_tab_ssu_all_pros_curated.txt", quote = F, sep = "\t")
write.table(TAX.curated, "tax_tab_ssu_all_pros_curated.txt", quote = F, sep = "\t")

