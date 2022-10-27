# GitBash, to copy output files in my local computer
#scp chh@10.223.24.140:/storage/hdd1/chh/Debora_amplicons/all_euks/ErrorProfiles_ssu_env_euks2* ./
  


# This script is an example of the basic dada2 workflow.
# See: https://benjjneb.github.io/dada2/tutorial.html
# adapted for highly variable fragment lengths


# save and load workspace
setwd("/storage/hdd1/chh/Debora_amplicons/all_euks")
# save.image("dada2_ssu_all_euks.Rdata")

# load packages
require(dada2)
require(ShortRead)
require(ggplot2)
require(gridExtra)
require(Biostrings)
require(scales)


# Load objects from 3rdRun_pros and single_cell_pros


# 3rdRUn
seqtab.rf.rc_3rdRun <- readRDS("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_euks/seqtab.rf.rc_3rdRun_euks.RDS")
seqtab.fr_3rdRun <- readRDS("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_euks/seqtab.fr_3rdRun_euks.RDS")

# single cell
seqtab.rf.rc_single <- readRDS("/storage/hdd1/chh/Debora_amplicons/Single_euks/seqtab.rf.rc_single_euks.RDS")
seqtab.fr_single <- readRDS("/storage/hdd1/chh/Debora_amplicons/Single_euks/seqtab.fr_single_euks.RDS")



# Merge sequence tables
seqtab <- mergeSequenceTables( 
  seqtab.fr_3rdRun,
  seqtab.fr_single,
  seqtab.rf.rc_3rdRun,
  seqtab.rf.rc_single,
  repeats = "sum"
)
dim(seqtab) #  318 83803




# Remove chimeras
# This may remove quite a bit of ASVs, but only a small fraction of your total sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = 170, verbose = TRUE, minFoldParentOverAbundance = 2)
ncol(seqtab.nochim)/ncol(seqtab)
# about 27.5% of ASVs classified as non-chimeric
summary(rowSums(seqtab.nochim)/rowSums(seqtab))
# Min.   1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6799  0.8963  0.9878  0.9414  0.9929  0.9981
# Looks good!

# Inspect ASV length distribution
table(nchar(colnames(seqtab.nochim))) # to see how many counts of unique sequences for each length
table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))) # to see how many counts of total sequences for each length

# range seems to be between 368 and 410
# when analyzing only single cell was between 369 and 396

# # Check unusual sequence lengths
uniquesToFasta(
  seqtab.nochim[, sample(which(nchar(colnames(seqtab.nochim)) < 368 | nchar(colnames(seqtab.nochim)) > 410), 100)],
  "check_ssu_all_euks.fasta"
)

# not conclusive by looking BLAST (NCBI)

# Remove potential junk sequences and singletons
# dada does not generate singletons, any singletons are introduced in the merging step
# Adjust range of sequence lengths based on expected length of marker gene fragment
seqtab.nochim2 <- seqtab.nochim[, colSums(seqtab.nochim) > 1 & nchar(colnames(seqtab.nochim)) >= 260 ]
dim(seqtab.nochim2) # 318 16032 ASVs
ncol(seqtab.nochim2)/ncol(seqtab) # 19.1%
summary(rowSums(seqtab.nochim2)/rowSums(seqtab.nochim))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.9703  0.9987  0.9999  0.9984  1.0000  1.0000
# ok!


#Get nSeqs summary
track <- cbind(
  rowSums(seqtab),
  rowSums(seqtab.nochim)
)
colnames(track) <- c("merged", "nochim")
track <- data.frame(track)

write.table(track, "nSeq_dada2_ssu_all_euks.txt", quote = F, sep = "\t")

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
  "/storage/hdd1/chh/Debora_amplicons/Single_euks/pr2_version_4.14.0_SSU_dada2.fasta.gz",
  multithread = 120,
  minBoot = 0,
  outputBootstraps = T,
  taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species")
)
# remember to save your workspace regularly



# don't do any filtering yet
# save output as is
write.table(track, "nSeq_dada2_ssu_all_euks.txt", quote = F, sep = "\t")
otu.print <- t(seqtab.nochim2)
rownames(otu.print) <- paste("sq", 1:ncol(seqtab.nochim2), sep = "")
write.table(otu.print, "otu_tab_ssu_all_euks.txt", quote = F, sep = "\t")
write.table(seqtab.nochim2, "otu_tab_with_seqs_ssu_all_euks.txt", quote = F, sep = "\t")
uniquesToFasta(seqtab.nochim2, "dada2_unique_ssu_all_euks.fasta")
tax.print <- tax$tax
rownames(tax.print) <- paste("sq", 1:nrow(tax$tax), sep = "")
all.equal(rownames(tax.print), rownames(otu.print))
#[1] TRUE
write.table(
  data.frame(
    tax.print, 
    seqlen = nchar(colnames(seqtab.nochim2))
  ),
  "tax_tab_ssu_all_euks.txt", 
  quote = F, 
  sep = "\t"
)
boot.print <- tax$boot
rownames(boot.print) <- paste("sq", 1:nrow(tax$boot), sep = "")
all.equal(rownames(boot.print), rownames(otu.print))
#[1] TRUE
write.table(boot.print, "tax_bootstrap_ssu_all_euks.txt", quote = F, sep = "\t")


# further data curation
# remove seqs of unexpected lengths based also on bootstrap
pdf("length_vs_genus_bootstrap_all_euks_4.pdf")
plot(
  tax$boot[, 6], 
  nchar(colnames(seqtab.nochim2)),
  pch = 16,
  col = rgb(red = 0, green = 0, blue = 0, alpha = rescale(sqrt(colSums(seqtab.nochim2)), to = c(0.1, 1))),
  cex = rescale(sqrt(colSums(seqtab.nochim2)), to = c(0.3, 1))
)
abline(h = c(353, 410), v = 70)
dev.off()

# by this plot we can observe that the best cut off is between 353-410.




# save files from server to  local drive ("C:/Users/draposo/.../R_out_euks")

# otu_tab_ssu_all_euks.txt
# otu_tab_with_seqs_ssu_all_euks.txt
# dada2_unique_ssu_all_euks.fasta
# tax_tab_ssu_all_euks.txt
# tax_bootstrap_ssu_all_euks.txt
# length_vs_genus_bootstrap_all_euks_4.pdf

# go to next script for data analysis
# next script: "data_curation_all_euks.RData"
