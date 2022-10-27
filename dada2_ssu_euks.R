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
setwd("/storage/hdd1/chh/Debora_amplicons/Single_euks")
# save.image("dada2_ssu_euks.Rdata")

# specify path to input fastq files
path <- "Clipped"
fnFs.fr <- sort(list.files(path, pattern="clip_fr_R1.fastq", full.names = TRUE))
fnRs.fr <- sort(list.files(path, pattern="clip_fr_R2.fastq", full.names = TRUE))
fnFs.rf <- sort(list.files(path, pattern="clip_rf_R1.fastq", full.names = TRUE))
fnRs.rf <- sort(list.files(path, pattern="clip_rf_R2.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs.fr), "_"), function(x) paste(x[1:2], collapse = "_"))
# optional: use mixedsort of package gtools to get full alphanumeric sort

# quality check
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_quality_check.R")
quality_check(
  c(fnFs.fr, fnFs.rf),
  c(fnRs.fr, fnRs.rf),
  file_base = "QualityProfile_ssu_euks"
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

# Considerations for trimming:
# this is a 2x300bp run
# expected max length: 400bp (?)
# min overlap: 30bp
# It is recommended to trim to just enough for the required length for sufficient overlap
# Caution: don't remove too much

# Define ranges for truncLen
range_truncLen <- matrix(
  c(280, 160,
    270, 170,
    260, 180,
    250, 190),
  nrow = 4,
  ncol = 2,
  byrow = T
)
# Based on the expected maximum fragment length, the trimming could be even stricter
# Feel free to adjust the parameter ranges further

# Define ranges for maxEE
range_maxEE <- matrix(
  c(2, 2,
    2, 3,
    3, 3),
  nrow = 3,
  ncol = 2,
  byrow = T
)

# Run parameter optimization
# This is quite time consuming and should only be attempted on a server with as many cores as you have samples (or at least 20)
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/dada2_screen_settings.R")
screen_filt_settings_fr <- screen_settings(
  sample.names, 
  fnFs.fr, 
  fnRs.fr, 
  range_maxEE, 
  range_truncLen, 
  cpus = 200
)
screen_filt_settings_rf <- screen_settings(
  sample.names, 
  fnFs.rf, 
  fnRs.rf, 
  range_maxEE, 
  range_truncLen, 
  cpus = 200
)

# This is just a gut feeling, but I would optimize for the following criteria:
#   small difference between 10 and 90 percentile of retained reads
#   high total proportion of retained reads
#   most stringent maxEE that does not result in severe loss of reads
pdf("Parameter_screening_ssu_euks.pdf", width = 7, height = 7)
plot(
  screen_filt_settings_fr[, "prop.total"],
  screen_filt_settings_fr[, "q90"] - screen_filt_settings_fr[, "q10"],
  col = rep(rainbow(nrow(range_maxEE)), nrow(range_truncLen)),
  pch = 16
)
text(
  screen_filt_settings_fr[, "prop.total"],
  screen_filt_settings_fr[, "q90"] - screen_filt_settings_fr[, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
plot(
  screen_filt_settings_rf[, "prop.total"],
  screen_filt_settings_rf[, "q90"] - screen_filt_settings_rf[, "q10"],
  col = rep(rainbow(nrow(range_maxEE)), nrow(range_truncLen)),
  pch = 16
)
text(
  screen_filt_settings_rf[, "prop.total"],
  screen_filt_settings_rf[, "q90"] - screen_filt_settings_rf[, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
dev.off()

# Run trimming with optimal parameters
filt.out.fr <- filterAndTrim(
  fwd = fnFs.fr, 
  filt = filtFs.fr, 
  rev = fnRs.fr, 
  filt.rev = filtRs.fr,
  truncLen = c(260, 180),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 200
)
filt.out.rf <- filterAndTrim(
  fwd = fnFs.rf, 
  filt = filtFs.rf, 
  rev = fnRs.rf, 
  filt.rev = filtRs.rf,
  truncLen = c(260, 180),
  maxN = 0,
  minQ = 2,
  maxEE = c(3, 3), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = 200
)

# Repeat quality check after trimming
quality_check(
  c(filtFs.fr, filtFs.rf),
  c(filtRs.fr, filtRs.rf),
  file_base = "QualityProfileFiltered_ssu_euks"
)

#  Looks ok...

# Learn error rates
# It is generally not necessary to increase the number of nbases used for the error estimation
# It is possible that with 10 rounds (MAX_CONSIST), the algorithm for learning the errors won't converge
# Increasing MAX_CONSIST will lead to longer run times, and may only marginally improve error estimation
# I would not recommend setting MAX_CONSIST higher than 15
errF.fr <- learnErrors(filtFs.fr, multithread = 70, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.fr <- learnErrors(filtRs.fr, multithread = 70, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errF.rf <- learnErrors(filtFs.rf, multithread = 70, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
errR.rf <- learnErrors(filtRs.rf, multithread = 70, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
# it is a good idea to save your workspace here

# check convergence of error estimation and plot error profiles
pdf("ErrorProfiles_ssu_euks.pdf")
barplot(log10(dada2:::checkConvergence(errF.fr) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.fr) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(errF.rf) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(errR.rf) + 1), main = "Convergence_rev")
plotErrors(errF.fr, nominalQ = TRUE)
plotErrors(errR.fr, nominalQ = TRUE)
plotErrors(errF.rf, nominalQ = TRUE)
plotErrors(errR.rf, nominalQ = TRUE)
dev.off()

# There are several instances, where the predicted error frequencies are lower than the observed ones
# Manually set error frequencies to value of quality score 40, if predicted lower than that
# errF.fr_mod <- errF.fr
# errR.fr_mod <- errR.fr
# errF.rf_mod <- errF.rf
# errR.rf_mod <- errR.rf
# errF.fr_mod$err_out <- t(apply(getErrors(errF.fr), 1, function(x) ifelse(x < x[41], x[41], x)))
# errR.fr_mod$err_out <- t(apply(getErrors(errR.fr), 1, function(x) ifelse(x < x[41], x[41], x)))
# errF.rf_mod$err_out <- t(apply(getErrors(errF.rf), 1, function(x) ifelse(x < x[41], x[41], x)))
# errR.rf_mod$err_out <- t(apply(getErrors(errR.rf), 1, function(x) ifelse(x < x[41], x[41], x)))
# pdf("ErrorProfiles_ssu_euks_mod.pdf")
# plotErrors(errF.fr_mod, nominalQ = TRUE)
# plotErrors(errR.fr_mod, nominalQ = TRUE)
# plotErrors(errF.rf_mod, nominalQ = TRUE)
# plotErrors(errR.rf_mod, nominalQ = TRUE)
# dev.off()

# Still not the best
# Try alternative Loess function
source("/storage/hdd1/chh/Repos/Tutorials/Dada2_workshop_UniHB/loessErrfun2.R")
err2F.fr <- learnErrors(filtFs.fr, errorEstimationFunction = loessErrfun2, multithread = 32, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
err2R.fr <- learnErrors(filtRs.fr, errorEstimationFunction = loessErrfun2, multithread = 32, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
err2F.rf <- learnErrors(filtFs.rf, errorEstimationFunction = loessErrfun2, multithread = 32, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
err2R.rf <- learnErrors(filtRs.rf, errorEstimationFunction = loessErrfun2, multithread = 32, randomize = TRUE, verbose = 1, MAX_CONSIST = 10)
pdf("ErrorProfiles_ssu_euks_2.pdf")
plotErrors(err2F.fr, nominalQ = TRUE)
plotErrors(err2R.fr, nominalQ = TRUE)
plotErrors(err2F.rf, nominalQ = TRUE)
plotErrors(err2R.rf, nominalQ = TRUE)
dev.off()

# This is much better :)

# Dereplicate and denoise samples
# This step takes a while...
# For large data set (e.g. full HiSeq lane), I strongly recommend pool = "pseudo"
# I would not use pool = FALSE as this will strongly impact (i.e. lower) your alpha diversity,
# which seems to be rather an artifact of the change in parameters than any true signal
dadaFs.fr <- dada(filtFs.fr, err = err2F.fr, multithread = 70, pool = TRUE)
# 243 samples were pooled: 3943346 reads in 245978 unique sequences.
dadaRs.fr <- dada(filtRs.fr, err = err2R.fr, multithread = 70, pool = TRUE)
# 243 samples were pooled: 3943346 reads in 189017 unique sequences.
dadaFs.rf <- dada(filtFs.rf, err = err2F.rf, multithread = 70, pool = TRUE)
# 243 samples were pooled: 4035401 reads in 220285 unique sequences.
dadaRs.rf <- dada(filtRs.rf, err = err2R.rf, multithread = 70, pool = TRUE)
# 243 samples were pooled: 4035401 reads in 308435 unique sequences.
# it is a good idea to save your workspace here

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
dim(seqtab.fr) # 243 6517
dim(seqtab.rf) # 243 5662

# This is the step at which separate denoising runs should be combined
# (e.g. if data comes from different sequencer runs or lanes, 
# or if fwd-rev and rev-fwd orientation were processed separately) 

# Generate reverse complement of rf
seqtab.rf.rc <- seqtab.rf
colnames(seqtab.rf.rc) <- rc(colnames(seqtab.rf))


# saving to merge later with 3rd run sequences
saveRDS(seqtab.rf.rc, "seqtab.rf.rc_single_euks.RDS")
saveRDS(seqtab.fr,  "seqtab.fr_single_euks.RDS")

# doing the next steps in another script now...("dada2_ssu_all_euks.Rdata")


