#################################################
#                                               #
#  NSEQ tables    SSU     all      euks         #
#                                               #
#################################################



############# Working on local drive ######

setwd("C:/Users/draposo/Documents/PhD/Chapter2/NSEQ")
#save.image("NSEQ_all_pros.RData")
#load("NSEQ_all_pros.RData")



# load OTU and TAX tables
NSEQ_2nd_run <- read.table("nSeq_dada2_ssu_2nd_run_pros_curated.txt", h = T, sep ="\t")
clip_den_3rd_run <- read.table("nSeq_dada2_ssu_3rd_run_pros.txt", h = T, sep = "\t")
merged_nochim <- read.table("nSeq_dada2_ssu_all_pros.txt", h = T, sep = "\t")
tabled <- read.table("nSeq_dada2_ssu_all_pros_curated.txt", h = T, sep = "\t")


# adjusting names of samples (all with dash)
rownames(tabled) <- gsub(".", "-", rownames(tabled), fixed = T)

# selecting only clipped, filtered and denoised columns of NSEQ table from 2nd run
clip_den_2nd_run <- NSEQ_2nd_run[,2:5]

# combining 1st run and 3rd run tables
clip_den <- merge(clip_den_2nd_run, clip_den_3rd_run, by=0, all=TRUE) # "by=0" stands for by row.names and "all=TRUE" to include all samples, even the ones not duplicated
clip_den[is.na(clip_den)] <- 0  # replace NA values by zero
clip_den$Clipped <- rowSums(clip_den[ , c("Clipped.x","Clipped.y")], na.rm=TRUE)
clip_den$Filtered <- rowSums(clip_den[ , c("Filtered.x","Filtered.y")], na.rm=TRUE)
clip_den$Denoised_fwd <- rowSums(clip_den[ , c("Denoised_fwd.x","Denoised_fwd.y")], na.rm=TRUE)
clip_den$Denoised_rev <- rowSums(clip_den[ , c("Denoised_rev.x","Denoised_rev.y")], na.rm=TRUE)
rownames(clip_den) <- clip_den$Row.names # saving rownames again
clip_den_final <- clip_den[,10:13] # selecting only the final summed columns

# combining in the same table
NSEQ <- cbind(clip_den_final[rownames(tabled),], merged_nochim[rownames(tabled),], tabled) # making sure to have sample names in same order

# NSEQ perc
NSEQ.perc <- data.frame(round(apply(NSEQ, 2, function(x) x/NSEQ$Clipped) * 100, 2))
summary(NSEQ.perc$tabled)

# identify how many more sequences were added by resequencing

# load metadata
META <- readxl::read_xlsx("C:/Users/draposo/Documents/PhD/Chapter2/R_out_pros/Lists_bioinformatic_all_pros.xlsx") 
META$Sample_NC <- as.factor(META$Sample_NC)
str(META)

colSums(NSEQ) # all samples

# setting NSEQ and META in same order
META <- as.data.frame(META) # setting first as df
rownames(META) <- META$Purified_PCR_product
META <- META[rownames(NSEQ),]

table(META$Sample_type)
# filter         NC_env NC_single_cell       sediment    single_cell 
# 32             11             19             32            224 

by(NSEQ$tabled, META$Sample_type, summary)
# META$Sample_type: filter
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 14826   17603   23992   24420   29783   39444 
# ------------------------------------------------------------------------------------------------ 
#   META$Sample_type: NC_env
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   338.5   584.0   719.5   837.0  1921.0 
# ------------------------------------------------------------------------------------------------ 
#   META$Sample_type: NC_single_cell
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 47.0   177.5   341.0  2150.3  1013.5 27561.0 
# ------------------------------------------------------------------------------------------------ 
#   META$Sample_type: sediment
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 12460   20284   22855   23556   26657   34723 
# ------------------------------------------------------------------------------------------------ 
#   META$Sample_type: single_cell
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 197    1973    4244    5999    7635   54345 

by(NSEQ$tabled, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))
by(NSEQ$nochim, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))

all.equal(rownames(NSEQ.perc), rownames(META)) # TRUE


by(NSEQ.perc$Filtered, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))

by(NSEQ.perc$Denoised_rev, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))

by(NSEQ.perc$merged, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))

by(NSEQ.perc$nochim, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))

by(NSEQ.perc$tabled, META$Sample_type, function(x) quantile (x, seq(0,1,0.05)))




write.table(NSEQ, "nSeq_dada2_ssu_all_pros_curated_complete.txt", quote = F, sep = "\t")
