#################################################
#                                               #
#  NSEQ tables    SSU     all      euks         #
#                                               #
#################################################



############# Working on local drive ######

setwd("C:/Users/draposo/Documents/PhD/Chapter2/NSEQ")
#save.image("NSEQ_all_euks.RData")
#load("NSEQ_all_euks.RData")



# load OTU and TAX tables

NSEQ_1st_run <- read.table("nSeq_dada2_ssu_1st_run_euks_curated.txt", h = T, sep ="\t")
clip_den_3rd_run <- read.table("nSeq_dada2_ssu_3rd_run_euks.txt", h = T, sep = "\t")
merged_nochim <- read.table("nSeq_dada2_ssu_all_euks.txt", h = T, sep = "\t")
tabled <- read.table("nSeq_dada2_ssu_all_euks_curated.txt", h = T, sep = "\t")


# adjusting names of samples (all with dash)
rownames(NSEQ_1st_run) <- gsub("_", "-", rownames(NSEQ_1st_run))
rownames(merged_nochim) <- gsub("_", "-", rownames(merged_nochim))
rownames(tabled) <- gsub("_", "-", rownames(tabled))


# selecting only clipped, filtered and denoised columns of NSEQ table from 1st run
clip_den_1st_run <- NSEQ_1st_run[,2:5]

# combining 1st run and 3rd run tables
clip_den <- rbind(clip_den_1st_run, clip_den_3rd_run)

# combining in the same table
NSEQ <- cbind(clip_den[rownames(tabled),], merged_nochim[rownames(tabled),], tabled) # making sure to have sample names in same order

# NSEQ perc
NSEQ.perc <- data.frame(round(apply(NSEQ, 2, function(x) x/NSEQ$Clipped) * 100, 2))
summary(NSEQ.perc$tabled)


write.table(NSEQ, "nSeq_dada2_ssu_all_euks_curated_complete.txt", quote = F, sep = "\t")
