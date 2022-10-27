#################################################
#                                               #
#  Data curation    SSU     all      euks       #
#                                               #
#################################################



############# Working on local drive ######

setwd("C:/Users/draposo/Documents/PhD/Chapter2/R_out_euks")
#save.image("data_curation_all_euks.RData")
#load("data_curation_all_euks.RData")
require(vegan)
require(readxl)

# load OTU and TAX tables
OTU <- data.frame(read.table("otu_tab_ssu_all_euks.txt", h = T, sep = "\t"))
TAX <- read.table("tax_tab_ssu_all_euks.txt", h = T, sep = "\t", stringsAsFactors = F, comment.char = "", quote = "")
TAX.boot <- read.table("tax_bootstrap_ssu_all_euks.txt", h = T, sep = "\t")
META <- data.frame(read_xlsx("Lists_bioinformatic_all_euks.xlsx"))
#NSEQ <- read.table("nSeq_dada2_ssu_all_euks.txt", h = T, sep = "\t")


# setting Purified_PCR_product as rownames to META, to match with the OTU table later
rownames(META) <- META$Purified_PCR_product

# adjusting names of OTU tables (all with underscore)
colnames(OTU) <- gsub(".", "_", colnames(OTU), fixed = T)


# to set in same order as occurs as columns in OTU
META <- META[colnames(OTU), ]

all.equal(rownames(META), colnames(OTU)) # [1] TRUE


# Analyzing bootstrap ranges

# with entire dataset
apply(TAX.boot, 2, summary)
#           Kingdom Supergroup  Division     Class     Order    Family     Genus   Species
# Min.      1.00000    1.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000
# 1st Qu. 100.00000   98.00000  95.00000  79.00000  61.00000  50.00000  25.00000  21.00000
# Median  100.00000  100.00000 100.00000 100.00000 100.00000  99.00000  63.00000  53.00000
# Mean     95.93787   86.01946  84.65126  81.31967  77.85298  75.09076  59.15145  53.91155
# 3rd Qu. 100.00000  100.00000 100.00000 100.00000 100.00000 100.00000  98.00000  92.00000
# Max.    100.00000  100.00000 100.00000 100.00000 100.00000 100.00000 100.00000 100.00000

# applying cut off observed in boostrap plot 

# remove seqs shorter than 353bp and longer than 410bp
TAX.cut <- as.matrix(TAX.boot[, -9])
TAX.cut <- TAX.cut[TAX$seqlen >= 353 & TAX$seqlen <= 410, ]
apply(TAX.cut, 2, summary)
#           Kingdom Supergroup  Division     Class     Order    Family     Genus  Species
# Min.      1.00000    1.00000   1.00000   0.00000   0.00000   0.00000   0.00000   0.0000
# 1st Qu. 100.00000  100.00000 100.00000  98.00000  90.00000  79.00000  36.00000  31.0000
# Median  100.00000  100.00000 100.00000 100.00000 100.00000 100.00000  73.00000  60.0000
# Mean     98.75751   93.08513  91.90282  89.10536  85.75213  82.85584  65.18317  59.4063
# 3rd Qu. 100.00000  100.00000 100.00000 100.00000 100.00000 100.00000  99.00000  95.0000
# Max.    100.00000  100.00000 100.00000 100.00000 100.00000 100.00000 100.00000 100.0000

# by this range we can define that the bootstrap of 70 is a good option. it was the same used 
#for the proks dataset and it means that we will remove 25% of the with <70 Genus bootstrap



# Applying cut off, removing singletons and setting bootstrap
TAX.clean <- as.matrix(TAX[, -9])
TAX.clean[TAX.boot < 70] <- "NA" # 70 bootstrap
TAX.clean <- TAX.clean[TAX$seqlen >= 353 & TAX$seqlen <= 410 & rowSums(OTU) >= 2 & TAX.clean[, 1] == "Eukaryota", ]
OTU.clean <- OTU[rownames(TAX.clean), ]
dim(OTU.clean)# [1] 13769   318
nrow(OTU.clean)/nrow(OTU) # 85.9%
summary(colSums(OTU.clean)/colSums(OTU))
# Min.     1st Qu.  Median    Mean   3rd Qu.    Max. 
# 0.8159  0.9832   0.9967    0.9886  0.9995  1.0000 




# NSEQ$classified <- colSums(OTU.clean)
# NSEQ.perc <- data.frame(round(apply(NSEQ, 2, function(x) x/NSEQ$Demux) * 100, 2))
# summary(NSEQ.perc$classified)
# # 
# summary(NSEQ$classified)
# #


# calculate rarefaction curves
plot.in <- OTU.clean
maxIndex <- colSums(decostand(plot.in, "pa"))
plot(
  0, 0,
  type = "n",
  xlim = c(0, max(colSums(plot.in))),
  ylim = c(0, max(maxIndex)),
  ylab = "Number of OTUs (Euks)",
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
    lwd = 1, col = META$Color[k]
  )
}

# looks okay

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

##### not doing this for now #########
# set detection threshold to exclude transient/food signal
# OTU.rel <- prop.table(as.matrix(OTU.clean), 2) * 100
# hist(OTU.rel[OTU.rel > 0], breaks = 200000, xlim = c(0, 0.1))
# abline(v = 0.01)
# hist(OTU.clean[OTU.clean > 0], breaks = 200000, xlim = c(0, 100))
# 
# # which proportion are 3 sequences translating to
# quantile(OTU.rel[OTU.clean == 3], seq(0, 1, 0.1))
# # 0.01% seems to be a good filter...
# # check again after excluding negative controls from the assessment
# OTU.filt <- OTU.clean
# OTU.filt[OTU.rel < 0.01] <- 0
# OTU.filt <- OTU.filt[rowSums(OTU.filt) >= 2, ]
# TAX.filt <- TAX.curated[rownames(OTU.filt), ]
# 
# # redraw rarefaction curves
# plot.in <- OTU.filt
# maxIndex <- colSums(decostand(plot.in, "pa"))
# plot(
#   0, 0,
#   type = "n",
#   xlim = c(0, max(colSums(plot.in))),
#   ylim = c(0, max(maxIndex)),
#   ylab = "Number of OTUs (Euks - filtered)",
#   xlab = "Sequencing depth",
#   cex.axis = 0.7,
#   las = 1,
#   axes = F,
#   cex.lab = 0.9,
#   mgp = c(1.8, 0.4, 0)
# )
# axis(1, las = 1, mgp = c(1, 0.2, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
# axis(2, las = 1, mgp = c(2, 0.4, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
# box("plot", lwd = 0.5)
# for(k in order(colSums(plot.in), decreasing = T)) {
#   temp <- rarefy(
#     plot.in[, k],
#     sample = round(seq(0, sum(plot.in[, k]), length.out = 100))
#   )
#   lines(
#     round(seq(0, sum(plot.in[, k]), length.out = 100)),
#     temp,
#     lwd = 1, col = META$Color[k]
#   )
# }
# all saturated (with the exception of neg controls)

# NSEQ$final <- colSums(OTU.filt)
# NSEQ.perc$final <- round(NSEQ$final / NSEQ$Demux * 100, 2)
# summary(NSEQ.perc$final)
# quantile(NSEQ.perc$final, seq(0, 1, 0.1))


write.table(track, "nSeq_dada2_ssu_all_euks_curated.txt", quote = F, sep = "\t")
write.table(OTU.clean, "otu_tab_ssu_all_euks_curated.txt", quote = F, sep = "\t")
write.table(TAX.curated, "tax_tab_ssu_all_euks_curated.txt", quote = F, sep = "\t")
