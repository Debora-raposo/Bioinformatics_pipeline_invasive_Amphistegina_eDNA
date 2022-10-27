###################
#                 #
# SSU prokaryotes #
#                 #
###################


# set working directory on server
WDIR="/storage/.../Single_pros"
cd $WDIR

# activate conda environment
conda activate trimming
# cutadapt 3.1

# create directory for logfiles
mkdir Logfiles

# CAUTION: It is not recommended to combine the output from several sequencer runs. Some argue, lanes should also be processed independently.
# If different sequencing batches of the same project are processed independently, they can be combined into one sample-by-asv table after merging (Before Chimera detection).
# This also applies to the separate processing of reads in the wrong orientation.
# CAUTION: If in a split run, the second run/lane only has a small number of sequences, it may be better to not denoise separately, as a certain sequencing depth is required for denoising.

# As each library only contains a limited number of samples (with a sequencing depth <100K per sample), we will process all libraries in a project together regardless of sequencing run.
# Otherwise, the number of sequences per library may not be sufficient for good denoising

# prepare index files
# split index information by library
# and prepare index fasta for demultiplexing
# make sure that there are no spaces in the table with the index and primer information
for library in $(sed '1d' Prok_single_cell_Debora.txt | cut -f6 | sort | uniq)
do
  sed '1d' Prok_single_cell_Debora.txt | awk -v FS="\t" -v OFS="\t" -v library=${library} '$6 == library' > "Library_"${library}"_index_primer.txt"
done
# double-check index combinations in R 


### Demultiplex using cutadapt version 3.0
# demux utility of cutadapt does not work for library prep strategy

# create directory for demultiplexed files
mkdir Demux

# as illumina adapters were ligated regardless of orientation, demultiplex both fwd-rev and rev-fwd
# treat them separately until after denoising and merging

# run demultiplexing with 2 threads per sample
# fwd-rev and rev-fwd orientation
for i in $(ls -1 Library_*)
do
  LIB=$(echo "${i}" | sed -e 's/^Library_//' -e 's/_index_primer.txt$//')
  LANE=$(expr ${LIB} - 6)
  while read line
  do
    SID=$(echo "${line}" | cut -f1)
    FWD=$(echo "${line}" | cut -f2 | sed 's/^/\^/')
    REV=$(echo "${line}" | cut -f4 | sed 's/^/\^/')
    OLP=7 # index length minus 1
    
    cutadapt -j 2 -O ${OLP} --no-indels -e 0 -g ${FWD} -G ${REV} --discard-untrimmed -o Demux/${SID}"_fr_R1.fastq.gz" -p Demux/${SID}"_fr_R2.fastq.gz" Libraries/"Library"${LIB}"_S"${LANE}"_L001_R1_001.fastq.gz" Libraries/"Library"${LIB}"_S"${LANE}"_L001_R2_001.fastq.gz" > Logfiles/${SID}".demux_fr.log" 2>&1
    
    cutadapt -j 2 -O ${OLP} --no-indels -e 0 -g ${REV} -G ${FWD} --discard-untrimmed -o Demux/${SID}"_rf_R1.fastq.gz" -p Demux/${SID}"_rf_R2.fastq.gz" Libraries/"Library"${LIB}"_S"${LANE}"_L001_R1_001.fastq.gz" Libraries/"Library"${LIB}"_S"${LANE}"_L001_R2_001.fastq.gz" > Logfiles/${SID}".demux_rf.log" 2>&1 
  done < ${i}
done

# further processing steps are easier and faster without gz
gunzip Demux/*.gz


### Primer clipping using cutadapt version 3.0
# get list of demux sample names
cd Demux
ls -1v *_fr_R1.fastq | sed 's/_fr_R1\.fastq//' > ../sample_names.txt
cd ..
wc -l sample_names.txt
# 243 samples

# create directory for primer-clipped files
mkdir Clipped

# set clipping parameters
FWD="^CCTACGGGNGGCWGCAG"
REV="^GACTACHVGGGTATCTAATCC"
OFWD=$(expr ${#FWD} - 2)
OREV=$(expr ${#REV} - 2)
ERROR=0.16

# remove primer from each sample
while read line
do
  # process fwd-rev orientation
  cutadapt -j 2 --no-indels -e ${ERROR} -g "${FWD};o=${OFWD}" -G "${REV};o=${OREV}" -m 50 --discard-untrimmed -o Clipped/${line}"_clip_fr_R1.fastq" -p Clipped/${line}"_clip_fr_R2.fastq" Demux/${line}"_fr_R1.fastq" Demux/${line}"_fr_R2.fastq" > Logfiles/${line}".clip_fr.log" 2>&1
  # process rev-fwd orientation
  cutadapt -j 2 --no-indels -e ${ERROR} -g "${REV};o=${OREV}" -G "${FWD};o=${OFWD}" -m 50 --discard-untrimmed -o Clipped/${line}"_clip_rf_R1.fastq" -p Clipped/${line}"_clip_rf_R2.fastq" Demux/${line}"_rf_R1.fastq" Demux/${line}"_rf_R2.fastq" > Logfiles/${line}".clip_rf.log" 2>&1
done < sample_names.txt


### Count sequences for each step and switch to dada2 in R
ls -1v Demux/*_fr_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste <(ls -1v Demux/*_fr_R1.fastq | xargs -n1 basename | sed 's/_fr_R1\.fastq//') - > tmp1
ls -1v Clipped/*_clip_fr_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste tmp1 - > tmp2
echo -e 'SID\tDemux\tClipped' | cat - tmp2 > nSeqs_fr.txt
rm tmp*
ls -1v Demux/*_rf_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste <(ls -1v Demux/*_rf_R1.fastq | xargs -n1 basename | sed 's/_rf_R1\.fastq//') - > tmp1
ls -1v Clipped/*_clip_rf_R1.fastq | xargs wc -l | grep -v "total" | awk '{print $1/4}' | paste tmp1 - > tmp2
echo -e 'SID\tDemux\tClipped' | cat - tmp2 > nSeqs_rf.txt
rm tmp*

