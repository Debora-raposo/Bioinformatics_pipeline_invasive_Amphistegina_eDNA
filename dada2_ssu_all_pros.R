# GitBash, so copy output files in my local computer
#scp chh@10.223.24.140:/storage/hdd1/chh/Debora_amplicons/Single_pros/seqtab_2ndRun_pros* ./
  
#/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_pros

# This script is an example of the basic dada2 workflow.
# See: https://benjjneb.github.io/dada2/tutorial.html
# adapted for highly variable fragment lengths

# load packages
require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)
require(Biostrings)
require(scales)

packageVersion("dada2")
# 1.16.0

# save and load workspace

setwd("/storage/hdd1/chh/Debora_amplicons/all_pros")
# save.image("dada2_ssu_all_pros.Rdata")



# Load objects from 3rdRun_pros and single_cell_pros


# 3rdRUn
seqtab.rf.rc_3rdRun <- readRDS("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_pros/seqtab.rf.rc_3rdRun_pros.RDS")
seqtab.fr_3rdRun <- readRDS("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_pros/seqtab.fr_3rdRun_pros.RDS")

# single cell
seqtab.rf.rc_single <- readRDS("/storage/hdd1/chh/Debora_amplicons/Single_pros/seqtab.rf.rc_single_pros.RDS")
seqtab.fr_single <- readRDS("/storage/hdd1/chh/Debora_amplicons/Single_pros/seqtab.fr_single_pros.RDS")



# Merge sequence tables
seqtab <- mergeSequenceTables( 
  seqtab.fr_3rdRun,
  seqtab.fr_single,
  seqtab.rf.rc_3rdRun,
  seqtab.rf.rc_single,
  repeats = "sum"
)
dim(seqtab) #  318 156115
# 


# Remove chimeras
# This may remove quite a bit of ASVs, but only a small fraction of your total sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 120, verbose = TRUE, minFoldParentOverAbundance = 2)
ncol(seqtab.nochim)/ncol(seqtab)
# about 54.82241% of ASVs classified as non-chimeric
summary(rowSums(seqtab.nochim)/rowSums(seqtab))
# Min.   1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7178  0.9724  0.9927  0.9587  0.9954  1.0000

# 50%(median) of the samples stll presents 99% of the sequences after chimera removal


# Inspect ASV length distribution
table(nchar(colnames(seqtab.nochim))) # to see how many counts of unique sequences for each length
table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))) # to see how many counts of total sequences for each length

# Check unusual sequence lengths
# uniquesToFasta(
#   seqtab.nochim[,which((nchar(colnames(seqtab.nochim)) < 398 & nchar(colnames(seqtab.nochim)) != 260) | nchar(colnames(seqtab.nochim)) > 429)], 
#   "check_ssu_all_pros.fasta"
# )

#we can skip since it was the same picks than tsingle cell


# Remove potential junk sequences and singletons (sequence that only occurs once)
# dada does not generate singletons, any singletons are introduced in the merging step
# Adjust range of sequence lengths based on expected length of marker gene fragment
seqtab.nochim2 <- seqtab.nochim[, colSums(seqtab.nochim) > 1 & nchar(colnames(seqtab.nochim)) >= 398 & nchar(colnames(seqtab.nochim)) <= 429]
dim(seqtab.nochim2) # 318 40927 ASVs
ncol(seqtab.nochim2)/ncol(seqtab) # 26%
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.8761  0.9918  0.9986  0.9861  0.9994  1.0000


#Get nSeqs summary
track <- cbind(
  rowSums(seqtab),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("merged", "nochim")
track <- data.frame(track)

write.table(track, "nSeq_dada2_ssu_all_pros.txt", quote = F, sep = "\t")

# Taxonomic classification
# Available options:
#   DECIPHER (IdTaxa), 
#   RDP
#   Blast
#   silvangs
#   etc.

# For now, we will use RDP
# I am disabling the bootstrap filtering, but saving the bootstrap values
# so that we can manually filter by bootstrap later
tax <- assignTaxonomy(
  seqtab.nochim2, 
  "/storage/hdd6/DB/Dada2/Silva/v138.1/silva_nr99_v138.1_train_set.fa.gz", # name of the training fasta provided by dada2 for RDP???
  multithread = 120,
  minBoot = 0, #not applying treshold now
  outputBootstraps = T,
  taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
)
# remember to save your workspace regularly

# continue from here



# don't do any filtering yet
# save output as is
write.table(track, "nSeq_dada2_ssu_all_pros.txt", quote = F, sep = "\t")
otu.print <- t(seqtab.nochim2) # t of tab separated?
rownames(otu.print) <- paste("sq", 1:ncol(seqtab.nochim2), sep = "")
write.table(otu.print, "otu_tab_ssu_all_pros.txt", quote = F, sep = "\t")
write.table(seqtab.nochim2, "otu_tab_with_seqs_ssu_all_pros.txt", quote = F, sep = "\t")
uniquesToFasta(seqtab.nochim2, "dada2_unique_ssu_all_pros.fasta")
tax.print <- tax$tax
rownames(tax.print) <- paste("sq", 1:nrow(tax$tax), sep = "")
all.equal(rownames(tax.print), rownames(otu.print)) 
#[1] TRUE
write.table(
  data.frame(
    tax.print, 
    seqlen = nchar(colnames(seqtab.nochim2))
  ),
  "tax_tab_ssu_all_pros.txt", 
  quote = F, 
  sep = "\t"
)
boot.print <- tax$boot
rownames(boot.print) <- paste("sq", 1:nrow(tax$boot), sep = "")
all.equal(rownames(boot.print), rownames(otu.print))
#[1] TRUE
write.table(boot.print, "tax_bootstrap_ssu_all_pros.txt", quote = F, sep = "\t")



# save files to my local drive ("C:/Users/draposo/.../R_out_pros")


# otu_tab_ssu_all_pros.txt
# otu_tab_with_seqs_ssu_all_pros.txt
# dada2_unique_ssu_all_pros.fasta
# tax_tab_ssu_all_pros.txt
# tax_bootstrap_ssu_all_pros.txt


# go to next script for data analysis
# next script: "data_curation_all_pros.RData"



