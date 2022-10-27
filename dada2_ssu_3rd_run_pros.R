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

setwd("/storage/hdd1/chh/Debora_amplicons/Environmental/R_out_pros")
# save.image("dada2_ssu_3rdRun_pros.Rdata")

# Libraries analyzed in the 3rd run
#Lib 7  -     10%       -      proks      - single cell
#Lib 8  -     10%       -      proks      - single cell
#Lib 13 -     20%       -      euks       - env (sediment & filter)
#Lib 14 -     20%       -      euks       - env (sediment & filter)
#Lib 15 -     20%       -      proks      - env (sediment & filter)
#Lib 16 -     20%       -      proks      - env (sediment & filter)


# specify path to input fastq files
path <- "../Clipped" # to take from level before
fnFs.fr <- sort(list.files(path, pattern="clip_fr_R1.fastq", full.names = TRUE))
fnRs.fr <- sort(list.files(path, pattern="clip_fr_R2.fastq", full.names = TRUE))
fnFs.rf <- sort(list.files(path, pattern="clip_rf_R1.fastq", full.names = TRUE))
fnRs.rf <- sort(list.files(path, pattern="clip_rf_R2.fastq", full.names = TRUE))
# example of clipped name files:  PR-1265_clip_fr_R2.fastq

# Extract sample names
sample.names <-  sapply(strsplit(basename(fnFs.fr), "_"), function(x) x[1])
#or the code bellow ?
#sample.names <- sapply(strsplit(basename(fnFs.fr), "_"), function(x) paste(x[1:2], collapse = "_"))

# set names
names(fnFs.fr) <- sample.names
names(fnRs.fr) <- sample.names
names(fnFs.rf) <- sample.names
names(fnRs.rf) <- sample.names


# separate here
pros <- scan("../sample_names_7_8_15_16.txt", what="character", sep="\n")

# to only select the ones with same names than pros
fnFs.fr <- fnFs.fr[pros] 
fnRs.fr <- fnRs.fr[pros] 
fnFs.rf <- fnFs.rf[pros] 
fnRs.rf <- fnRs.rf[pros] 

# run only after using sample names before
sample.names <- pros

# to check how many products
length(fnFs.fr)

# pros

# quality check
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_quality_check.R")
quality_check(
  c(fnFs.fr, fnFs.rf),
  c(fnRs.fr, fnRs.rf),
  file_base = "QualityProfile_ssu_3rdRun_pros"
)

# Place filtered files in Filtered/ subdirectory
filtFs.fr <- file.path("Filtered", paste0(sample.names, "_filt_fr_R1.fastq"))
filtRs.fr <- file.path("Filtered", paste0(sample.names, "_filt_fr_R2.fastq"))
filtFs.rf <- file.path("Filtered", paste0(sample.names, "_filt_rf_R1.fastq"))
filtRs.rf <- file.path("Filtered", paste0(sample.names, "_filt_rf_R2.fastq"))
names(filtFs.fr) <- sample.names
names(filtRs.fr) <- sample.names
names(filtFs.rf) <- sample.names
names(filtRs.rf) <- sample.names



# Run trimming with optimal parameters
# same as the ones used for single cell pros

filt.out.fr <- filterAndTrim(
  fwd = fnFs.fr, 
  filt = filtFs.fr, 
  rev = fnRs.fr, 
  filt.rev = filtRs.fr,
  truncLen = c(260, 200),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 150 # number of cores in the server
)
filt.out.rf <- filterAndTrim(
  fwd = fnFs.rf, 
  filt = filtFs.rf, 
  rev = fnRs.rf, 
  filt.rev = filtRs.rf,
  truncLen = c(260, 200),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 150
)

# Repeat quality check after trimming
quality_check(
  c(filtFs.fr, filtFs.rf),
  c(filtRs.fr, filtRs.rf),
  file_base = "QualityProfileFiltered_ssu_3rdRun_pros"
)

#  

# Learn error rates
# It is generally not necessary to increase the number of nbases used for the error estimation
# It is possible that with 10 rounds (MAX_CONSIST), the algorithm for learning the errors won't converge
# Increasing MAX_CONSIST will lead to longer run times, and may only marginally improve error estimation
# I would not recommend setting MAX_CONSIST higher than 15
errF.fr <- learnErrors(filtFs.fr, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.fr <- learnErrors(filtRs.fr, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errF.rf <- learnErrors(filtFs.rf, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.rf <- learnErrors(filtRs.rf, multithread = 150, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
# it is a good idea to save your workspace here

# check convergence of error estimation and plot error profiles
pdf("ErrorProfiles_ssu_3rdRun_pros.pdf")
barplot(log10(dada2:::checkConvergence(errF.fr) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.fr) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(errF.rf) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.rf) + 1), main = "Convergence_rev")
plotErrors(errF.fr, nominalQ = TRUE)
plotErrors(errR.fr, nominalQ = TRUE)
plotErrors(errF.rf, nominalQ = TRUE)
plotErrors(errR.rf, nominalQ = TRUE)
dev.off()
#  

# save.image("dada2_ssu_3rdRun_pros.Rdata")


# Dereplicate and denoise samples
# This step takes a while...
# For large data set (e.g. full HiSeq lane), I strongly recommend pool = "pseudo"
# I would not use pool = FALSE as this will strongly impact (i.e. lower) your alpha diversity,
# which seems to be rather an artifact of the change in parameters than any true signal
dadaFs.fr <- dada(filtFs.fr, err = errF.fr, multithread = 150, pool = TRUE)
# 
dadaRs.fr <- dada(filtRs.fr, err = errR.fr, multithread = 150, pool = TRUE)
# 
dadaFs.rf <- dada(filtFs.rf, err = errF.rf, multithread = 150, pool = TRUE)
# 
dadaRs.rf <- dada(filtRs.rf, err = errR.rf, multithread = 150, pool = TRUE)
# 
# it is a good idea to save your workspace here

# work until here today (18.10)


# Merge reads
mergers.fr <- mergePairs(
  dadaFs.fr,
  filtFs.fr, 
  dadaRs.fr, 
  filtRs.fr, 
  minOverlap = 10,
  verbose = TRUE
)
mergers.rf <- mergePairs(
  dadaFs.rf,
  filtFs.rf, 
  dadaRs.rf, 
  filtRs.rf, 
  minOverlap = 10,
  verbose = TRUE
)


# Create sequence table
seqtab.fr <- makeSequenceTable(mergers.fr)
seqtab.rf <- makeSequenceTable(mergers.rf)
dim(seqtab.fr) # [1]   157 52080
dim(seqtab.rf) # [1]   157 109952


# This is the step at which separate denoising runs should be combined
# (e.g. if data comes from different sequencer runs or lanes, 
# or if fwd-rev and rev-fwd orientation were processed separately) 

# Generate reverse complement of rf
seqtab.rf.rc <- seqtab.rf
colnames(seqtab.rf.rc) <- rc(colnames(seqtab.rf))


# Get nSeqs summary
getN <- function(x) sum(getUniques(x))
track <- cbind(
  filt.out.fr + filt.out.rf,
  sapply(dadaFs.fr, getN) + sapply(dadaFs.rf, getN), 
  sapply(dadaRs.fr, getN) + sapply(dadaRs.rf, getN) 
)
colnames(track) <- c("Clipped", "Filtered", "Denoised_fwd", "Denoised_rev")
rownames(track) <- c(sample.names)
track <- data.frame(track)

write.table(track, "nSeq_dada2_ssu_3rd_run_pros.txt", quote = F, sep = "\t")


# save objects
saveRDS(seqtab.rf.rc, "seqtab.rf.rc_3rdRun_pros.RDS")
saveRDS(seqtab.fr,  "seqtab.fr_3rdRun_pros.RDS")

# Go to next script to merge with single cell sequences
# next script: dada2_ssu_all_pros.R

