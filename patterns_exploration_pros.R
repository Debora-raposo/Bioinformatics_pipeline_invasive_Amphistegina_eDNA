#################################
# Analysis of 16S amplicon data #
#################################

Sys.setenv(LANG = "en") # in case system language changes to German

##### prepare environment ####

# set working directory
setwd("G:/My Drive/Research/PhD/Data/Data_analysis/Metabarcoding_analyses/R_out_pros") # Windows
setwd("~/Library/CloudStorage/GoogleDrive-deboraposo@gmail.com/My Drive/Research/PhD/Data/Data_analysis/Metabarcoding_analyses/R_out_pros") # Mac

# load packages
require(vegan)
#BiocManager::install("ALDEx2")
require(ALDEx2)
require(reshape)
require(gplots)
require(tidyverse)
#install.packages("treemap")
require(treemap)
#install.packages("treemapify")
require(treemapify)
require(ggplot2)
require(reshape2) # for melt function and more
require(dplyr) # for dcast and melt functions
require(patchwork) # to combine ggplot plots
require(ecodist) # for PCoA
# install.packages("rcartocolor")
require(rcartocolor) # color blind pallete 
require("UpSetR") # upset plots

# loading and saving workspace
load("patterns_exploration_pros.Rdata")
#save.image("patterns_exploration_pros.Rdata")

#####
##### read and sort data #####

# R objects generated after data cleaning/merging in "data_exploration_all_pros.R"

ASV <- readRDS("ASV_merge_pros.RDS")
TAX <- readRDS("TAX_merge_pros.RDS")
META <- readRDS("META_merge_pros.RDS")

# sort rownames according to sample type in META (to separate single cell from env samples)
META <- META[order(META$Sample_type, META$Site),]

# re-setting color scheme for Sample Type
color.sample.type <- c("dodgerblue", "forestgreen", "red3")
META$color_sample_type <- META$Sample_type
levels(META$color_sample_type) <- color.sample.type
META$color_sample_type <- as.character(META$color_sample_type)

shape.sample.type <- c(6,10,8)
META$shp_sample_type <- META$Sample_type
levels(META$shp_sample_type) <- shape.sample.type
META$shp_sample_type <- as.numeric(as.character(META$shp_sample_type))

#change levels order of site 
META$Site <- factor(META$Site, levels = c("Eilat","Tel Shikmona", "Plemmirio","Capo Passero" ))
color.site <- c("gold1","salmon","darkorchid1","blue")
META$color_site_2 <- META$Site
levels(META$color_site_2) <- color.site
META$color_site_2 <- as.character(META$color_site_2)


# applying new order to ASV and TAX tables
ASV <- ASV[,rownames(META)]
TAX <- TAX[rownames(ASV),]

all.equal(colnames(ASV), rownames(META))
# [1] TRUE
all.equal(rownames(ASV), rownames(TAX))
# [1] TRUE

ASV.rel <- prop.table(ASV, 2) * 100

nrow(ASV)

# selecting only foraminifera microbiome
META.sing <- META[META$Sample_type == "single_cell",]
ASV.sing <- ASV[,rownames(META.sing)]
ASV.sing <- ASV.sing[rowSums(ASV.sing) > 0,] # remove ASVs with zero seq
TAX.sing <- TAX[rownames(ASV.sing),]

#####

##### rarefaction curve per sample type ######
META.plot <- META
META.plot$Sample_type <- factor(META.plot$Sample_type, levels=c("filter", "sediment", "single_cell"), 
                                labels = c("Seawater", "Sediment", "Foraminifera"))

levels(META.plot$Sample_type)

dev.off()
rarecurve(t(ASV), step = 100, col = META.plot$color_sample_type, label = F, ylab = "Number of ASVs")

legend(
  'topright',
  legend = levels(META.plot$Sample_type),
  bty = "n",
  col = color.sample.type,
  pch = 15,
  pt.cex = 2,
  title = "Sample type",
  cex= 1.5 # to adjust size of legend box 
)


##### 
##### rarefaction curve per site ######
rarecurve(t(ASV), step = 100, col = META.plot$color_site_2, label = F, ylab = "Number of ASVs")

# add legend for sites
legend(
  "topright",
  legend = levels(META.plot$Site),
  col = color.site,
  bty = "n",
  pch = 15,
  pt.cex = 2,
  title = "Site",
  cex= 1.5 # to adjust size of legend box 
)

#####
##### rarefaction curve - Inverse Simpson index  #####
# rarefaction curves for the inverse Simpson index, using some custom code
# for each sequencing depth, repeatedly randomly subsample the data set and calculate inverse Simpson index
knots <- 100 # pick 100 sequencing depth for rarefaction
n <- 10 # randomly subsample 10 times at each sequencing depth
# the output in a list object
invs.repeated <- vector("list", length = ncol(ASV.clean))
# repeat the process for each sample (i.e. columns in ASV)
for(i in 1:ncol(ASV.clean)) {
  print(i) # show progress (sample number)
  temp <- ASV.clean[, i] # subset to only one sample
  temp <- temp[temp > 0] # remove ASVs not occurring in that sample
  
  # build matrix to collect inverse Simpson values in
  invs.repeated[[i]] <- matrix(NA, nrow = knots, ncol = n + 1)
  
  # determine 'knots' sequencing depth for subsampling from 0 to the available sequencing depth
  temp.sub <- round(seq(0, sum(temp), length.out = knots))
  
  # save these values in the first column of the output table for this sample
  invs.repeated[[i]][, 1] <- temp.sub
  
  # repeat n times
  for(j in 1:n) {
    
    # for each sequencing depth (starting with the second - zero is not used)
    for(k in 2:knots) {
      # subsample and calculate inverse Simpson
      invs.repeated[[i]][k, j + 1] <- diversity(
        rrarefy(
          temp,
          temp.sub[k]
        ),
        "invsimpson"
      )
    }
    
    # set zero sequencing depth to zero inverse Simpson
    invs.repeated[[i]][temp.sub==0, j + 1] <- 0
  }
}

str(invs.repeated)

# get maximum for y axis
maxIndex <- max(
  sapply(
    invs.repeated,
    function(x) {
      apply(x[, -1], 1, mean)
    }
  )
)


# start plot per sample type
META.plot <- META
META.plot$Sample_type <- factor(META.plot$Sample_type, levels=c("filter", "sediment", "single_cell"), 
                                labels = c("Seawater", "Sediment", "Foraminifera"))

plot(
  0, 0,
  type = "n",
  xlim = c(0, max(colSums(ASV.clean))),
  ylim = c(0, max(maxIndex)),
  xlab = "Sequencing depth",
  ylab = "Inverse Simpson index",
  cex.axis = 0.7,
  las = 1,
  axes = F,
  cex.lab = 0.9,
  mgp = c(1.8, 0.4, 0)
)
axis(1, las = 1, mgp = c(1, 0.2, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
axis(2, las = 1, mgp = c(2, 0.4, 0), tcl = -0.3, cex.axis = 0.7, lwd = 0.5)
box("plot", lwd = 0.5)
for(k in order(colSums(ASV.clean), decreasing = T)) {
  subsample <- invs.repeated[[k]][, 1]
  # show mean of the 10 time repeated random subsampling
  temp <- apply(invs.repeated[[k]][, -1], 1, mean)
  lines(
    subsample,
    temp,
    lwd = 1,
    col = META.plot$color_sample_type[k]
  )
}
# show minimum sequencing depth in the data set
abline(v = min(colSums(ASV.clean)))

#####

##### Cleaning taxa with 'unclassified' or similar from the data set ######
View(unique(TAX.sing[, 1:4]))
# check Phylums present
table(TAX.sing[,"Phylum"])
# Acidobacteriota              Actinobacteriota         Bacteria_unclassified                  Bacteroidota 
# 53                           273                            88                           798 
# Bdellovibrionota                 Calditrichota              Campylobacterota                   Chloroflexi 
# 46                             7                            56                            47 
# Cyanobacteria                  Dadabacteria                  Deinococcota                  Dependentiae 
# 116                            30                            15                             1 
# Desulfobacterota             Entotheonellaeota                Fibrobacterota                    Firmicutes 
# 133                             3                             1                           186 
# Fusobacteriota               Gemmatimonadota              Latescibacterota              Margulisbacteria 
# 64                            16                            14                             1 
# Marinimicrobia (SAR406 clade)                        MBNT15                   Myxococcota                         NB1-j 
# 2                             1                            78                            65 
# Nitrospinota                  Nitrospirota               Patescibacteria               Planctomycetota 
# 5                             7                            66                           439 
# Proteobacteria  SAR324 clade(Marine group B)                 Spirochaetota                   Sumerlaeota 
# 3031                             7                            30                             3 
# Thermotogota             Verrucomicrobiota                         WPS-2                           WS2 
# 1                           168                             3                             1 

# no "unclassified" phylum, like in the euks dataset. no need to clean here.

table(TAX[,"Phylum"]) 
# one phylum to remove "Bacteria_unclassified"

table(TAX.sing[,"Class"]) # many classes to remove

table(TAX.sing[,"Order"]) # many orders to remove 

phylum_to_remove <- c("Bacteria_unclassified")

classes_to_remove <- c("Actinobacteriota_unclassified",  "Latescibacterota_unclassified", "Margulisbacteria_unclassified",
                       "Proteobacteria_unclassified", "SAR324 clade(Marine group B)_unclassified", "WPS-2_unclassified",
                       "Acidobacteriota_unclassified", "Marinimicrobia (SAR406 clade)_unclassified", "NB1-j_unclassified",
                       "Planctomycetota_unclassified", "WS2_unclassified", "Bacteria_unclassified", "Desulfobacterota_unclassified",
                       "Firmicutes_unclassified", "MBNT15_unclassified", "Myxococcota_unclassified", "Patescibacteria_unclassified", 
                       "Sva0485_unclassified", "Zixibacteria_unclassified")


# removing these taxa
TAX.clean <- TAX[!TAX[, "Phylum"] %in% phylum_to_remove, ] 
TAX.clean <- TAX.clean[!TAX.clean[, "Class"] %in% classes_to_remove, ] 

table(TAX.clean[,"Class"])
grep(pattern="_unclassified",x=names(table(TAX.clean[,"Class"])),value=TRUE)
# ok

# now using this table to remove the orders
# orders_to_remove <- c("028H05-P-BN-P5_unclassified", "Alphaproteobacteria Incertae Sedis","Bacilli_unclassified", 
#                       "Desulfobacterota_unclassified", "Gammaproteobacteria Incertae Sedis", "Gracilibacteria_unclassified",
#                       "Margulisbacteria_unclassified",  "OM190_unclassified", "Parcubacteria_unclassified", 
#                       "PAUC43f marine benthic group_unclassified", "Pla3 lineage_unclassified", "Planctomycetota_unclassified",
#                       "SAR202 clade",  "Sericytochromatia_unclassified", "Subgroup 21_unclassified", "Alphaproteobacteria_unclassified",
#                       "bacteriap25_unclassified", "Bacteroidia_unclassified", "BD7-11_unclassified", "Clostridia_unclassified",
#                       "Cyanobacteriia_unclassified", "Gammaproteobacteria_unclassified", "Latescibacterota_unclassified",
#                       "Marinimicrobia (SAR406 clade)_unclassified", "NB1-j_unclassified", "Pla4 lineage_unclassified",
#                       "Polyangia_unclassified", "Proteobacteria_unclassified", "SAR324 clade(Marine group B)_unclassified",
#                       "Subgroup 22_unclassified", "Sumerlaeia_unclassified", "vadinHA49_unclassified", "WPS-2_unclassified",
#                       "Acidobacteriota_unclassified", "Bacilli_unclassified", "Firmicutes_unclassified", "KD4-96_unclassified",
#                       "MBNT15_unclassified", "P9X2b3D02_unclassified", "Subgroup 26_unclassified", "Verrucomicrobiae_unclassified", 
#                       "BD2-11 terrestrial group_unclassified", "WS2_unclassified", "4-29-1_unclassified", "ABY1_unclassified",
#                       "AT-s3-28_unclassified", "Desulfuromonadia_unclassified", "Dojkabacteria_unclassified", "JG30-KF-CM66_unclassified",
#                       "Kiritimatiellae_unclassified", "Phycisphaerae_unclassified", "Planctomycetes_unclassified", "Thermodesulfovibrionia_unclassified",
#                       "V2072-189E03_unclassified"
# )

# removing these orders
#TAX.clean <- TAX.clean[!TAX.clean[, "Order"] %in% orders_to_remove, ] 
# not removing unclassified Orders  anymore, decided to clean only on the Phylum and Class level

ASV.clean <- ASV[rownames(TAX.clean),]
nrow(ASV)
nrow(ASV.clean)

ASV.sing.clean <- ASV.clean[,rownames(META.sing)]
ASV.sing.clean <- ASV.sing.clean[rowSums(ASV.sing.clean) > 0,] # remove ASVs with zero seq
TAX.sing.clean <- TAX.clean[rownames(ASV.sing.clean),]

# check classes remaining after cleaning
table(TAX.sing.clean[,"Class"])
#table(TAX.sing.clean[,"Order"])
#table(TAX.clean[, "Order"])


#####
##### number of sequences retained ######

# foraminifera
nrow(ASV.sing.clean)   
# [1] 5649


# sediment
ASV.sed <- ASV.clean[,META$Sample_type == "sediment"]
ASV.sed <- ASV.sed[rowSums(ASV.sed) >0, ] 
nrow(ASV.sed)
# [1] 11750
AX.sed <- TAX.clean[rownames(ASV.sed),]


# seawater
ASV.filters <- ASV.clean[,META$Sample_type == "filter"]
ASV.filters <- ASV.filters[rowSums(ASV.filters) >0, ] # excluding lines with no ASV
nrow(ASV.filters)
# [1] 9863
TAX.filters <- TAX.clean[rownames(ASV.filters),]

#####

##### Fig donut plotting all taxa (donut and tree map plots) #####

# creating grouping vector based on Phylum

tax.group <- TAX.clean[ ,"Phylum"]
tax.group[TAX.clean[,"Phylum"] == "Proteobacteria"] <- ifelse(
  TAX.clean[TAX.clean[,"Phylum"] == "Proteobacteria", "Class"] == "Alphaproteobacteria", 
  "Proteobacteria (Alpha)",
  ifelse(
    TAX.clean[TAX.clean[,"Phylum"] == "Proteobacteria", "Class"] == "Gammaproteobacteria", 
    "Proteobacteria (Gamma)",
    "Proteobacteria (other)"
  )
)


ASV.group <- aggregate(
  ASV.clean,
  by = list(tax.group),
  FUN = sum
)
rownames(ASV.group) <- ASV.group$Group.1
ASV.group <- ASV.group[, -1]

all.equal(rownames(META), colnames(ASV.clean))
# [1] TRUE

ASV.group.st <- aggregate(
  t(ASV.group),# transpose, aggregate always work on the rows
  by = list(META$Sample_type),
  FUN = sum 
)
rownames(ASV.group.st) <- ASV.group.st$Group.1
ASV.group.st <- t(ASV.group.st[, -1])

ASV.group.st.rel <- prop.table(ASV.group.st, 2)

ASV.group.st.rel.abund <- ASV.group.st.rel[apply(ASV.group.st.rel, 1, max) > 0.01,]

ASV.group.st.rel.plot <- rbind(ASV.group.st.rel.abund, 1 - colSums(ASV.group.st.rel.abund))
rownames(ASV.group.st.rel.plot)[nrow(ASV.group.st.rel.plot)] <- "Other phyla" # rownames is a vector, so we use nrow to het the final row


carto_pal(n= 14, name = "Safe") 
# this palette has only 12 colors, have to add two more
safe_plus2colors <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#F0E442", "#E69F00", "#888888")

# plotting
par(mfrow=c(1,4), mar = rep(0,4))  # to have three pie charts (one per sample type) and the legend on the side
for (i in 1:3){  # to make one per sample type
  pie(
    ASV.group.st.rel.plot[,i], 
    col = safe_plus2colors, 
    border = NA, 
    labels = NA, 
    clockwise = T
  )
}
plot.new()
legend(
  "center",
  legend = rownames(ASV.group.st.rel.plot),
  col = safe_plus2colors,
  pch = 15,
  cex = 1,
  pt.cex = 2,
  ncol = 1,
  bty = "n",
  title = "Phylum",
  title.adj = 0,
  title.cex = 1.5 # works on 4.2.1 R version on Mac
)
dev.off()


#####

##### Fig donut plot for each individual site ######
# filtering one site per time
s <- c("Eilat")
s <- c("Tel Shikmona")
s <- c("Plemmirio")
s <- c("Capo Passero")
METAs <- META[META$Site == s,]
ASV.all.s <- ASV.clean[,rownames(METAs)]
ASV.all.s <- ASV.all.s[rowSums(ASV.all.s) > 0,] # excluding ASVs not present in this subset

ASV.sw.s <- ASV.all.s[,rownames(METAs[METAs$Sample_type == "filter",])] 
ASV.sw.s <- ASV.sw.s[rowSums(ASV.sw.s) > 0,]

ASV.sed.s <- ASV.all.s[,rownames(METAs[METAs$Sample_type == "sediment",])] 
ASV.sed.s <- ASV.sed.s[rowSums(ASV.sed.s) > 0,] 

ASV.sing.s <- ASV.all.s[,rownames(METAs[METAs$Sample_type == "single_cell",])] 
ASV.sing.s <- ASV.sing.s[rowSums(ASV.sing.s) > 0,] 

TAX.clean.s <- TAX.clean[rownames(ASV.all.s),]
  
# creating grouping vector based on Phylum

tax.group <- TAX.clean.s[ ,"Phylum"]
tax.group[TAX.clean.s[,"Phylum"] == "Proteobacteria"] <- ifelse(
  TAX.clean.s[TAX.clean.s[,"Phylum"] == "Proteobacteria", "Class"] == "Alphaproteobacteria", 
  "Proteobacteria (Alpha)",
  ifelse(
    TAX.clean.s[TAX.clean.s[,"Phylum"] == "Proteobacteria", "Class"] == "Gammaproteobacteria", 
    "Proteobacteria (Gamma)",
    "Proteobacteria (other)"
  )
)

ASV.group <- aggregate(
  ASV.all.s,
  by = list(tax.group),
  FUN = sum
)
rownames(ASV.group) <- ASV.group$Group.1
ASV.group <- ASV.group[, -1]

all.equal(rownames(METAs), colnames(ASV.all.s))
# [1] TRUE

ASV.group.st <- aggregate(
  t(ASV.group),# transpose, aggregate always work on the rows
  by = list(METAs$Sample_type),
  FUN = sum 
)
rownames(ASV.group.st) <- ASV.group.st$Group.1
ASV.group.st <- t(ASV.group.st[, -1])

ASV.group.st.rel <- prop.table(ASV.group.st, 2)

# keep only taxa that appeared in general plot (to have corresponding colors for each in the different plots)
ASV.group.st.rel.taxa.general <-  ASV.group.st.rel[rownames(ASV.group.st.rel) %in% c(rownames(ASV.group.st.rel.plot)),]

# calculating porportion of "Other Phyla" based on taxa that was kicked out in the subsettiing above
ASV.group.st.rel.other.phyla <- t(as.data.frame(colSums(ASV.group.st.rel) - colSums(ASV.group.st.rel.taxa.general))) 

ASV.group.st.rel.plot.s <- rbind(ASV.group.st.rel.taxa.general, ASV.group.st.rel.other.phyla)
rownames(ASV.group.st.rel.plot.s)[nrow(ASV.group.st.rel.plot.s)] <- "Other phyla" # rownames is a vector, so we use nrow to het the final row


#  palette with 14 colors, based on Safe palette from carto_pal pkg
safe_plus2colors <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#F0E442", "#E69F00", "#888888")

# plotting
par(mfrow=c(1,4), mar = rep(0,4))  # to have three pie charts (one per sample type) and the legend on the side
for (i in 1:3){  # to make one per sample type
  pie(
    ASV.group.st.rel.plot.s[,i], 
    col = safe_plus2colors, 
    border = NA, 
    labels = NA, 
    clockwise = T,
  )
  
}
plot.new()
legend(
  "center",
  legend = rownames(ASV.group.st.rel.plot.s),
  col = safe_plus2colors,
  pch = 15,
  cex = 1,
  pt.cex = 2,
  ncol = 1,
  bty = "n",
  title = paste("Phylum - ", s),
  title.adj = 0,
  title.cex = 1.5 # works on 4.2.1 R version on Mac
)


#####

##### NMDS all data - to plot sample type ######
ASV.clean.rel <- prop.table(ASV.clean, 2) * 100

# calculate NMDS
NMDS <- metaMDS(t(ASV.clean.rel), k = 2)
NMDS$stress
# 0.1088548 # all taxa

# set meta table
META.plot <- META

META.plot$Sample_type <- factor(META.plot$Sample_type, levels = c("filter", "sediment", "single_cell"), 
                                labels = c("Seawater", "Sediment", "Foraminifera"))

# build nice plot element by element
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n")

# add hulls based on sites

#ordihull(NMDS, META.plot$Sample_type)

# add points for sample type
points(
  NMDS$points[, 1],
  NMDS$points[, 2],
  col = META.plot$color_sample_type, #or per site
  #col = META.plot$color_site_2, 
  pch = 8, # or show sample types as shapes and sites as colors
  # pch = META.plot$shp_sample_type, 
  cex = 1
)
# add legend of stress value
legend("topright", "Stress = 0.11", bty = "n", cex = 1) 

# # add legend for sample type
legend(
  "bottom",
  legend = levels(META.plot$Sample_type),
  col = color.sample.type,
  pch = 8,
  bty = "n",
  pt.cex = 2,
  cex = 1.2
  # title = "Sample type"
)

# or

# add legend for sites and sample types
legend("bottom", legend = c(levels(META.plot$Site), levels(META.plot$Sample_type)), 
       pch = c(15, 15, 15, 15, 6, 10, 8), 
       col = c("gold1", "salmon","darkorchid1","blue","black","black","black"), 
       bty = "n", 
       pt.cex = c(2,2,2,2,1.5,1.5,1.5)) 

#####
##### inv simpson alpha diversity (boxplot) ####

META.plot <- META
META.plot$Sample_type <- factor(META.plot$Sample_type, levels=c("filter", "sediment", "single_cell"), 
                                labels = c("Seawater", "Sediment", "Foraminifera"))

levels(META.plot$Sample_type)
color.sample.type

# try with ggplot using log of axis
invS <- diversity(t(ASV.clean), "invsimpson")
# preparing data for ggplot
invSdf <- as.data.frame(invS)

all.equal(rownames(META.plot), rownames(invSdf))
# [1] TRUE
META.plot$invS <- invSdf$invS

scaleFUN <- function(x) sprintf("%.0f", x)

ggplot(META.plot, aes(y=invS, x=Sample_type, fill = Sample_type)) + 
  geom_boxplot()+
  geom_jitter()+
  scale_y_continuous(trans='log', labels = scaleFUN) +
  scale_fill_manual(values = c("dodgerblue", "forestgreen", "red3"))+
  facet_grid(~Site) +
  theme_minimal()+
  theme(strip.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10)
  )


#####
##### ANOSIM all proks in all sample types ####

# calculate ANOSIM between sample type
anosim(t(ASV.clean.rel), META$Sample_type)

#anosim(t(ASV.clean.rel), META$Site)

# calculate pairwise ANOSIM statistics
# download function from: https://raw.githubusercontent.com/chassenr/Tutorials/master/R_roundtable_sequence_analysis/anosimPosthoc.R
source("C:/Users/chassenrueck/Documents/Repos/Tutorials/R_roundtable_sequence_analysis/anosimPosthoc.R") # doesn't work, run the one in the script (end)
# run function anosimPosthoc in the end of the script
ANOSIMposthoc(t(ASV.clean.rel), META$Sample_type)

# ANOSIMposthoc(t(ASV.clean.rel), META$Site)

# to constrain each site per time:
# Group 1:
compare.site <- "Capo Passero"
# Group 2:
compare.site <- "Plemmirio"
# Group 3:
compare.site <- "Tel Shikmona"
# Group 4:
compare.site <- "Eilat"

# subset data set
META.sub <- META[META$Site == compare.site, ]
ASV.rel.sub <- ASV.clean.rel[, rownames(META.sub)]
ASV.rel.sub <- ASV.rel.sub[rowSums(ASV.rel.sub) > 0, ]

# calculate ANOSIM between sample type
anosim(t(ASV.rel.sub), META.sub$Sample_type)

# run function anosimPosthoc in the end of the script
ANOSIMposthoc(t(ASV.rel.sub), META.sub$Sample_type)


#####

##### upset plot - bacterial intersections for all sites together ####
# to see intersection of bacterial between sample types
?upset()

# all sites together
nrow(ASV.clean)
3317/13670
# create ASV.table per sample type
st <- c("filter")
st <- c("sediment")
st <- c("single_cell")

METAst <- META[META$Sample_type == st,]
ASV.all.st <- ASV.clean[,rownames(METAst)]
ASV.all.st <- ASV.all.st[rowSums(ASV.all.st) > 0,] # excluding ASVs not present in this subset

ASV.sw <- ASV.all.st[,rownames(METAst[METAst$Sample_type == "filter",])] 
ASV.sw <- ASV.sw[rowSums(ASV.sw) > 0,]

ASV.sed <- ASV.all.st[,rownames(METAst[METAst$Sample_type == "sediment",])] 
ASV.sed <- ASV.sed[rowSums(ASV.sed) > 0,] 

ASV.sing <- ASV.all.st[,rownames(METAst[METAst$Sample_type == "single_cell",])] 
ASV.sing <- ASV.sing[rowSums(ASV.sing) > 0,] 


# create list with ASVs name per sample type
upsetlist <- list(rownames(ASV.sw), rownames(ASV.sed), rownames(ASV.sing))
names(upsetlist) <- c("Seawater", "Sediment", "Foraminifera")

# plot upset
upset.plot <- upset(fromList(upsetlist), sets = c("Foraminifera", "Sediment", "Seawater"), 
                    keep.order = T,  mainbar.y.label = "Sample type Intersections", 
                    main.bar.color = c("forestgreen", "red3","dodgerblue", "black", "black","black", "black"),
                    # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
                    sets.bar.color = c("red3","forestgreen", "dodgerblue"),
                    sets.x.label = "Bacterial ASVs", text.scale = c(2,2,2,1.5,2,2)
)

upset.plot

#####
##### upset plot - filtering one site per time for Sup Mat ####
# filtering one site per time
s <- c("Eilat")
s <- c("Tel Shikmona")
s <- c("Plemmirio")
s <- c("Capo Passero")
METAs <- META[META$Site == s,]
ASV.all.s <- ASV.clean[,rownames(METAs)]
ASV.all.s <- ASV.all.s[rowSums(ASV.all.s) > 0,] # excluding ASVs not present in this subset

ASV.sw.s <- ASV.all.s[,rownames(METAs[METAs$Sample_type == "filter",])] 
ASV.sw.s <- ASV.sw.s[rowSums(ASV.sw.s) > 0,]

ASV.sed.s <- ASV.all.s[,rownames(METAs[METAs$Sample_type == "sediment",])] 
ASV.sed.s <- ASV.sed.s[rowSums(ASV.sed.s) > 0,] 

ASV.sing.s <- ASV.all.s[,rownames(METAs[METAs$Sample_type == "single_cell",])] 
ASV.sing.s <- ASV.sing.s[rowSums(ASV.sing.s) > 0,] 


# make upset list
upsetlist.s <- list(rownames(ASV.sw.s), rownames(ASV.sed.s), rownames(ASV.sing.s))
names(upsetlist.s) <- c("Seawater", "Sediment", "Foraminifera")

# plot upset
p1 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 3500,
            main.bar.color = c("forestgreen", "dodgerblue", "red3", "black", "black","black", "black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Bacterial ASVs", text.scale = c(2,2,2,1.5,2,2)
)
p1 # eilat

p2 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 3500,
            main.bar.color = c("forestgreen", "red3","dodgerblue", "black", "black","black","black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Bacterial ASVs", text.scale = c(2,2,2,1.5,2,2)
)
p2 # tel shikmona

p3 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 3500,
            main.bar.color = c("forestgreen", "dodgerblue", "red3","black", "black","black","black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Bacterial ASVs", text.scale = c(2,2,2,1.5,2,2)
)
p3 # plemmirio

p4 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 3500,
            main.bar.color = c("forestgreen", "dodgerblue", "red3","black", "black","black", "black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Bacterial ASVs", text.scale = c(2,2,2,1.5,2,2)
)
p4 # capo passero

#####

##### NMDS single cell only - to plot sites  ######

# calculate NMDS
ASV.sing.clean.rel <- prop.table(ASV.sing.clean, 2) * 100
NMDS <- metaMDS(t(ASV.sing.clean.rel), k = 3)
NMDS$stress
#  0.1758467 # single cell


# set meta table
META.plot <- META.sing

META.plot$Depth_cat <- factor(META.plot$Depth_cat, levels=c("Shallow", "Deep"))
META.plot$Substrate <- as.factor(META.plot$Substrate)

str(META.plot)

############################## plotting site as colors, substrate as shape and depth as size

# 1 and 2
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n", choices = c(1,2))

# add hulls based on sites
#ordihull(NMDS, META.plot$Site, choices = c(1,2))

# add points for sites
points(
  NMDS$points[, 1],
  NMDS$points[, 2],
  bg = META.plot$color_site_2,
  pch = as.numeric(META.plot$Substrate)+20,
  cex = as.numeric(META.plot$Depth_cat)
)
# add legend of stress value
legend("bottomright", "stress = 0.18", bty = "n", cex = 1.3) 

# 1 and 3
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n", choices = c(1,3))

# add hulls based on sites
#ordihull(NMDS, META.plot$Site, choices = c(1,3))

# add points for sites
points(
  NMDS$points[, 1],
  NMDS$points[, 3],
  bg = META.plot$color_site_2,
  pch = as.numeric(META.plot$Substrate)+20,
  cex = as.numeric(META.plot$Depth_cat)
)
# add legend of stress value
legend("bottomright", "stress = 0.18", bty = "n", cex = 1.3) 

# 2 and 3
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n", choices = c(2,3))

# add hulls based on sites
#ordihull(NMDS, META.plot$Site, choices = c(2,3))

# add points for sites
points(
  NMDS$points[, 2],
  NMDS$points[, 3],
  bg = META.plot$color_site_2,
  pch = as.numeric(META.plot$Substrate)+20,
  cex = as.numeric(META.plot$Depth_cat)
)

# add legend of stress value
legend("bottomright", "stress = 0.18", bty = "n", cex = 1.3) 

# add legend for sites, substrates and depth

legend("bottomright", legend = c(levels(META.plot$Site), levels(META.plot$Substrate), levels(META.plot$Depth_cat)), 
       pch = c(15, 15, 15, 15, 21, 22, 16, 16), 
       col = c("gold1", "salmon","darkorchid1","blue","black","black","black","black"), 
       bty = "n", 
       pt.cex = c(2,2,2,2,2,2,1,2)) 



##### 
##### PERMANOVA and posthoc for bacterial taxa in foraminifera ######

META.adonis <- META.sing # change for the sub set dataset
str(META.adonis)

META.adonis$Depth_cat <- as.factor(META.adonis$Depth_cat)
META.adonis$Substrate <- as.factor(META.adonis$Substrate)

# PERMANOVA per site 
adonis2(t(ASV.sing.clean.rel) ~ Site * Depth_cat * Substrate, data = META.adonis, sqrt.dist = T)
# 
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = t(ASV.sing.clean.rel) ~ Site * Depth_cat * Substrate, data = META.adonis, sqrt.dist = T)
#                           Df SumOfSqs R2      F       Pr(>F)    
# Site                       3    5.025 0.10677 4.5816  0.001 ***
# Depth_cat                  1    0.502 0.01068 1.3743  0.014 *  
# Substrate                  1    0.961 0.02042 2.6291  0.001 ***
# Site:Depth_cat             3    1.387 0.02947 1.2645  0.002 ** 
# Site:Substrate             2    1.763 0.03746 2.4112  0.001 ***
# Depth_cat:Substrate        1    0.558 0.01186 1.5267  0.004 ** 
# Site:Depth_cat:Substrate   2    1.038 0.02205 1.4195  0.001 ***
# Residual                  98   35.827 0.76128                  
# Total                    111   47.061 1.00000                  
# ---
#  

# all factor are significant, although site play a slighty higher effect


# Pairwise comparison of site
site.comb <- combn(levels(META.adonis$Site), 2) # creating a triangular pairwise comparison
permanova.list <-apply(
  site.comb,
  2,
  function(x) {
    adonis2(
      t(ASV.sing.clean.rel[, rownames(META.adonis)[META.adonis$Site %in% x]]) ~ Site, 
      data = droplevels(META.adonis[META.adonis$Site %in% x, ]), 
      sqrt.dist = T
    )
  }
)

# reformat as table
permanova.df <- data.frame(
  t(site.comb),
  do.call(
    "rbind",
    lapply(
      permanova.list,
      function(x) {
        c(x$R2[1], x$'F'[1], x$'Pr(>F)'[1])
      }
    )
  )
)
colnames(permanova.df) <- c("cond1", "cond2", "R2", "F", "p")
# adjust p-values for multiple testing
permanova.df$p.adj <- p.adjust(permanova.df$p, method = "fdr")

permanova.df
# all sites are different!



# Pairwise comparison of depth_cat
depth.comb <- combn(levels(META.adonis$Depth_cat), 2) # creating a triangular pairwise comparison
permanova.list <-apply(
  depth.comb,
  2,
  function(x) {
    adonis2(
      t(ASV.sing.clean.rel[, rownames(META.adonis)[META.adonis$Depth_cat %in% x]]) ~ Depth_cat, 
      data = droplevels(META.adonis[META.adonis$Depth_cat %in% x, ]), 
      sqrt.dist = T
    )
  }
)

# reformat as table
permanova.df <- data.frame(
  t(depth.comb),
  do.call(
    "rbind",
    lapply(
      permanova.list,
      function(x) {
        c(x$R2[1], x$'F'[1], x$'Pr(>F)'[1])
      }
    )
  )
)
colnames(permanova.df) <- c("cond1", "cond2", "R2", "F", "p")
# adjust p-values for multiple testing
permanova.df$p.adj <- p.adjust(permanova.df$p, method = "fdr")

permanova.df
# the two depths categoroies are different!


# Pairwise comparison of substrate
subst.comb <- combn(levels(META.adonis$Substrate), 2) # creating a triangular pairwise comparison
permanova.list <-apply(
  subst.comb,
  2,
  function(x) {
    adonis2(
      t(ASV.sing.clean.rel[, rownames(META.adonis)[META.adonis$Substrate %in% x]]) ~ Substrate, 
      data = droplevels(META.adonis[META.adonis$Substrate %in% x, ]), 
      sqrt.dist = T
    )
  }
)

# reformat as table
permanova.df <- data.frame(
  t(subst.comb),
  do.call(
    "rbind",
    lapply(
      permanova.list,
      function(x) {
        c(x$R2[1], x$'F'[1], x$'Pr(>F)'[1])
      }
    )
  )
)
colnames(permanova.df) <- c("cond1", "cond2", "R2", "F", "p")
# adjust p-values for multiple testing
permanova.df$p.adj <- p.adjust(permanova.df$p, method = "fdr")

permanova.df
# the two substrate types are different!
levels()



# Pairwise comparison of site*depth*substrate

META.adonis$interaction <- droplevels(interaction(META.adonis$Site,META.adonis$Depth_cat, META.adonis$Substrate))
levels(META.adonis$interaction)

posthoc.comb <- combn(levels(META.adonis$interaction), 2) # creating a triangular pairwise comparison
permanova.list <-apply(
  posthoc.comb,
  2,
  function(x) {
    adonis2(
      t(ASV.rel.sing[, rownames(META.adonis)[META.adonis$interaction %in% x]]) ~ interaction, 
      data = droplevels(META.adonis[META.adonis$interaction %in% x, ]), 
      sqrt.dist = T
    )
  }
)

# reformat as table
permanova.df <- data.frame(
  t(posthoc.comb),
  do.call(
    "rbind",
    lapply(
      permanova.list,
      function(x) {
        c(x$R2[1], x$'F'[1], x$'Pr(>F)'[1])
      }
    )
  )
)
colnames(permanova.df) <- c("cond1", "cond2", "R2", "F", "p")
# adjust p-values for multiple testing
permanova.df$p.adj <- p.adjust(permanova.df$p, method = "fdr")

View(permanova.df)

# the two substrate types are different!
##### 







##### Old plots and ideas
##### betadispersion within foraminifera, to check sites ####
betadispersion <- betadisper(
  vegdist(t(ASV.sing.clean.rel)), 
  META.sing$Site,sqrt.dist = T
)

par(mar= c(7,4,1,1))
boxplot(betadispersion, ylim = c(0,1), las = 2)

# plotting on ggplot2
str(betadispersion)
betadispersion.df <- as.data.frame(betadispersion$distances)

META.plot <- META.sing 
all.equal(rownames(META.plot), rownames(betadispersion.df))
# [1] TRUE
META.plot$betadisp <- betadispersion.df$`betadispersion$distances`

ggplot(META.plot, aes(y=betadisp, x=Site, fill = Site)) + 
  geom_boxplot()+
  geom_jitter(shape=8, position=position_jitter(0.2), alpha=.4)+
  coord_cartesian(ylim = c(0, 1))+
  scale_fill_manual(values = c("gold1", "salmon", "darkorchid1", "blue"))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10)
  )

#####








############################################################################################################################################

##### looking the dominant ASVs within the foraminifera bacterial microbiome - didnt work well, no important dominant ASVs #####


# selecting ASVs of single cells only
ASV.sing <- ASV[,META$Sample_type == "single_cell"] # select single cell
ASV.sing <- ASV.sing[rowSums(ASV.sing) > 0,] # excluding ASVs that are not in the subset
META.sing <- META[colnames(ASV.sing),] # selecting only rows that represents the single cell in META table
all.equal(rownames(META.sing), colnames(ASV.sing))
ASV.rel.sing <- prop.table(ASV.sing, 2) * 100

TAX.sing <- TAX[rownames(ASV.sing),] # selecting TAX just of single cells


# plotting barplot with PlotAbund

p <- PlotAbund(
  ASV.rel.sing, # change levels to see the different plots for each one
  abund = 10, # selecting the 3 most dominant taxa in each sample and showing their proportions in the other samples
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T,
  save.OTU = T
)

dim(p)
# 376 113
# 375(376 - 1 row "other") most abundant ASVs 

summary(t(p[,-ncol(p)])[,"other"]) # percentage for others
# Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
# 2.337  22.393  30.598  30.566  38.612  73.766 

abundant.ASVs<- rownames(p)[-nrow(p)]

# saving the vector file to a txt file in desktop (to use to select these seqs from fasta file)
#write(abundant.ASVs, "abundant.diatoms.foraminifera.txt")

# obtaining ASV table just for these abundant ASVs
ASV.abundant <- ASV.sing[abundant.ASVs,]

# same for the tax table
TAX.abundant <- TAX.sing[abundant.ASVs,]

# How much percent of the diversity these ASVs represent?
sum(ASV.abundant)/sum(ASV.sing) # 0.65
nrow(TAX.abundant)/nrow(TAX.sing) # 0.06

# 65% of the number of ASVs from foraminifera but only 6% of the diversity
##### 
##### plotting foraminiferal abundant microbiome data as tree map ##########
TAX.pooled.abundant <- vector(mode = "list", length = 6)
names(TAX.pooled.abundant) <- colnames(TAX.abundant)[1:6] 
for (i in 1:6) {
  temp <- aggregate(
    ASV.abundant,
    by = list(TAX.abundant[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.pooled.abundant[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# take table with TAX.pooled.abundant values within the Class levels 
data.abundant <- melt(TAX.pooled.abundant$Class) # change level as wished
colnames(data.abundant) <- c("Class", "Extraction_Voucher", "Value")

# join metadata
require(dplyr)
data.abundant.meta <- merge(x=data.abundant, y=META, by="Extraction_Voucher") 

# create table with sum of value per group
data.abundant.site <- dcast(data.abundant.meta, Site + Class ~ ., sum, value.var = "Value") # this way is better for plotting
colnames(data.abundant.site) <- c("Site", "Class", "Value")

# plotting tree map
ggplot(data.abundant.site, aes(area = Value, fill = Class, label = Class)) +
  geom_treemap() +
  facet_grid(. ~ Site)

# --> the most freq Class in the sicilian sites: Gammaproteobacteria
# --> the most freq Class in the israeli sites: Alphaproteobacteria


##### 
##### plotting all data from the foraminifera, not only the abundandt ########

TAX.sing.pooled <- vector(mode = "list", length = 6)
names(TAX.sing.pooled) <- colnames(TAX.sing)[1:6] 
for (i in 1:6) {
  temp <- aggregate(
    ASV.sing,
    by = list(TAX.sing[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.sing.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# take table with TAX.sing.pooled values within the Class levels 
data.sing <- melt(TAX.sing.pooled$Class) # change level as wished
colnames(data.sing) <- c("Class", "Extraction_Voucher", "Value")

# join metadata
data.sing.meta <- merge(x=data.sing, y=META, by="Extraction_Voucher") 

# create table with sum of value per group
data.sing.meta.site <- dcast(data.sing.meta, Site + Class ~ ., sum, value.var = "Value") # this way is better for plotting
colnames(data.sing.meta.site) <- c("Site", "Class", "Value")

# plotting tree map
ggplot(data.sing.meta.site, aes(area = Value, fill = Class, label = Class)) +
  geom_treemap() +
  facet_grid(. ~ Site)

# the most freq Class in the sicilian sites: Gammaproteobacteria
# the most freq Class in the israeli sites: Alphaproteobacteria
##### 
##### tree map for environmental samples to see if they also differ - do not clean or select abundant data ######


# selecting ASVs of single cells only
ASV.env <- ASV[,META$Sample_type != "single_cell"] # select environmental samples
ASV.env <- ASV.env[rowSums(ASV.env) > 0,] # excluding ASVs that are not in the subset
META.env <- META[colnames(ASV.env),] # selecting only rows that represents the single cell in META table
all.equal(rownames(META.env), colnames(ASV.env))

TAX.env <- TAX[rownames(ASV.env),] 

TAX.env.pooled <- vector(mode = "list", length = 6)
names(TAX.env.pooled) <- colnames(TAX.env)[1:6] 
for (i in 1:6) {
  temp <- aggregate(
    ASV.env,
    by = list(TAX.env[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.env.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# take table with TAX.env.pooled values within the Class levels 
data.env <- melt(TAX.env.pooled$Class) # change level as wished
colnames(data.env) <- c("Class", "Extraction_Voucher", "Value")

# join metadata
data.env.meta <- merge(x=data.env, y=META, by="Extraction_Voucher") 

# create table with sum of value per group
data.env.meta.site <- dcast(data.env.meta, Site + Class ~ ., sum, value.var = "Value") # this way is better for plotting
colnames(data.env.meta.site) <- c("Site", "Class", "Value")

# plotting tree map
ggplot(data.env.meta.site, aes(area = Value, fill = Class, label = Class)) +
  geom_treemap() +
  facet_grid(. ~ Site)

# the most freq Class in the sicilian sites: Gammaproteobacteria, Alphaproteobacteria and Bacteroidia
# the most freq Class in the israeli sites: Alphaproteobacteria, Bacteroidia and Cyanobacteriia



# # observing from which taxa these ASVs are and their frequency per site
# 
# # aggregating per site
# agg.env <- aggregate(t(ASV.env),
#                  by = list(META.env$Site),
#                  FUN = sum)
# rownames(agg.env) <- agg.env$Group.1
# agg.env <- agg.env[,-1] # removing first column "Group.1"
# agg.env.t <- t(agg.env)
# 
# all.equal(rownames(TAX.env), rownames(agg.env.t))
# # [1] TRUE
# 
# # combining with TAX table
# TAX.env.freq <-  merge(agg.env.t, TAX.env, by=0, all=TRUE) # by=0 is the same as by="row.names"
# 
# 
# 
##### 
##### NMDS all taxa - to plot sample types  ######
# calculate NMDS

NMDS <- metaMDS(t(ASV.rel), k = 2)
NMDS$stress
# 0.1072931 # all taxa

# set meta table
META.plot <- META

# build nice plot element by element
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n")

# add hulls based on sites

ordihull(NMDS, META.plot$color_sample_type)

# add points for sample type
points(
  NMDS$points[, 1],
  NMDS$points[, 2],
  col = META.plot$color_sample_type,
  pch = 16,
  cex = 1
)

# add legend for sample type
legend(
  "bottom",
  legend = levels(META.plot$Sample_type),
  col = c("dodgerblue", "blue", "red"),
  pch = 15,
  pt.cex = 1,
  title = "Sample type"
)

##### 
##### NMDS all taxa - to plot site #####
# calculate NMDS
NMDS <- metaMDS(t(ASV.rel), k = 2)
NMDS$stress
#  0.1072932 # all taxa

# set meta table
META$Substrate <- as.factor(META$Substrate)
META$Depth_cat <- as.factor(META$Depth_cat)
META.plot <- META

# build nice plot element by element
dev.off() #to reset the PAR settings I configured earlier

plot(NMDS, display = "sites", type = "n")

# add hulls based on sample type

ordihull(NMDS, META.plot$Sample_type)
# add points for sites
points(
  NMDS$points[, 1],
  NMDS$points[, 2],
  bg = META.plot$color_site,
  pch = as.numeric(META.plot$Substrate)+20,
  cex = as.numeric(META.plot$Depth_cat)
)

# add legend for sites
legend(
  "bottom",
  legend = levels(META.plot$Site),
  col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
  pch = 15,
  pt.cex = 1,
  title = "Site"
)
# add legend for substrates
legend(
  "top",
  legend = levels(META.plot$Substrate),
  col = "black",
  pch = c(21, 22),
  pt.cex = 1,
  title = "Substrate"
)

# add legend for depth
legend(
  "center",
  legend = levels(META.plot$Depth_cat),
  col = "black",
  pch = 16,
  pt.cex = c(1,2),
  title = "Depth"
)

##### inv simpson all sample together and plotting each sample type per time ####

META.plot <- META
META.plot$Sample_type <- factor(META.plot$Sample_type, levels=c("filter", "sediment", "single_cell"), 
                                labels = c("Seawater", "Sediment", "Foraminifera"))

levels(META.plot$Sample_type)
color.sample.type

# calculate the inverse Simpson index without subsampling
invS <- diversity(t(ASV.clean), "invsimpson")

# show boxplot with points
par(mar = c(10,4,4,2) + 0.1) # adjusting the margins size to fit the x axis labels

boxplot(
  invS ~ droplevels(interaction(META.plot$Site, META.plot$Sample_type)), 
  las = 2, # 2 means that the orientation of the x axis label will be perpendicular to the x axis
  ylab = "Inverse Simpson index",
  xlab = "", 
  ylim = c(0,600),
  cex.axis = 0.7,
  col = c("dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue", "forestgreen", "forestgreen", "forestgreen", "forestgreen", "red", "red", "red", "red"),
  outline = F
  
)
points(
  jitter(as.numeric(droplevels(interaction(META.plot$Site, META.plot$Sample_type)))), 
  invS,
  #bg = META.merge$Condition,
  pch = 1,
)
legend(
  "topright",
  legend = levels(META.plot$Sample_type),
  col = color.sample.type,
  pch = 15,
  pt.cex = 2,
  bty = "n",
  title = "Sample type"
)
# doesnt look so good diversity is too different between sample tyopes and when comparing with the same plot for the euks

# plot one sample type per time, to compare euks vs proks

# create ASV.table per sample type
# plotting with ggplot for better adjustment of Y axis

# foraminifera
st <- c("single_cell")
METAst <- META[META$Sample_type == st,]
ASV.all.st <- ASV.clean[,rownames(METAst)]
ASV.all.st <- ASV.all.st[rowSums(ASV.all.st) > 0,] # excluding ASVs not present in this subset

ASV.sing <- ASV.all.st[,rownames(METAst[METAst$Sample_type == "single_cell",])] 
ASV.sing <- ASV.sing[rowSums(ASV.sing) > 0,] 

# calculate the inverse Simpson index without subsampling
invS <- diversity(t(ASV.sing), "invsimpson")

# preparing data for ggplot
invSdf <- as.data.frame(invS)
all.equal(rownames(METAst), rownames(invSdf))
# [1] TRUE
METAst$invS <- invSdf$invS

plot_invS_sing <- ggplot(METAst, aes(y=invS, x=Site)) + 
  geom_boxplot(fill = "red3")+
  geom_jitter()+
  scale_y_continuous(trans='sqrt', limits = c(0,100)) +
  ylab(paste("Inverse Simpson index\n\n", "Foraminifera\n"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 15)
  )

plot_invS_sing 

# seawater
st <- c("filter")
METAst <- META[META$Sample_type == st,]
ASV.all.st <- ASV.clean[,rownames(METAst)]
ASV.all.st <- ASV.all.st[rowSums(ASV.all.st) > 0,] # excluding ASVs not present in this subset

ASV.sw <- ASV.all.st[,rownames(METAst[METAst$Sample_type == "filter",])] 
ASV.sw <- ASV.sw[rowSums(ASV.sw) > 0,]

# calculate the inverse Simpson index without subsampling
invS <- diversity(t(ASV.sw), "invsimpson")

# preparing data for ggplot
invSdf <- as.data.frame(invS)
all.equal(rownames(METAst), rownames(invSdf))
# [1] TRUE
METAst$invS <- invSdf$invS

plot_invS_sw <- ggplot(METAst, aes(y=invS, x=Site)) + 
  geom_boxplot(fill = "dodgerblue")+
  geom_jitter()+
  scale_y_continuous(trans='sqrt', limits = c(0,150)) +
  ylab(paste("Inverse Simpson index\n\n", "Seawater\n"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 15)
  )

plot_invS_sw 

# # sediment
st <- c("sediment")
METAst <- META[META$Sample_type == st,]
ASV.all.st <- ASV.clean[,rownames(METAst)]
ASV.all.st <- ASV.all.st[rowSums(ASV.all.st) > 0,] # excluding ASVs not present in this subset

ASV.sed <- ASV.all.st[,rownames(METAst[METAst$Sample_type == "sediment",])] 
ASV.sed <- ASV.sed[rowSums(ASV.sed) > 0,] 

# calculate the inverse Simpson index without subsampling
invS <- diversity(t(ASV.sed), "invsimpson")

# preparing data for ggplot
invSdf <- as.data.frame(invS)
all.equal(rownames(METAst), rownames(invSdf))
# [1] TRUE
METAst$invS <- invSdf$invS

plot_invS_sed <- ggplot(METAst, aes(y=invS, x=Site)) + 
  geom_boxplot(fill = "forestgreen")+
  geom_jitter()+
  scale_y_continuous(trans='sqrt', limits = c(0,550)) +
  ylab(paste("Inverse Simpson index\n\n", "Sediment\n"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 15, hjust = 0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size = 15)
  )

plot_invS_sed 


#####
##### bar plot - taxonomic composition all taxa ####

# creating list to plot data
TAX.pooled <- vector(mode = "list", length = 6)
names(TAX.pooled) <- colnames(TAX.clean)[1:6] # kingdom to Genus
for (i in 1:6) {
  temp <- aggregate(
    ASV.clean,
    by = list(TAX.clean[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# calculating relative proportions based on list
TAX.pooled.rel <- lapply(
  TAX.pooled,
  function(x) {
    prop.table(x, 2) * 100
  }
)

barplot(TAX.pooled.rel$Phylum, col = rainbow(nrow(TAX.pooled.rel$Phylum)), 
        legend.text = row.names(TAX.pooled.rel$Phylum), args.legend = list(x = "right"), 
        ylim=c(0,100))

barplot(TAX.pooled.rel$Order, col = rainbow(nrow(TAX.pooled.rel$Order)), 
        legend.text = row.names(TAX.pooled.rel$Order), args.legend = list(x = "right"), 
        ylim=c(0,100))


# Showing the most dominant taxa with PlotAbund function
# PlotAbund function at the end of the script

PlotAbund(
  TAX.pooled.rel$Phylum, # change levels to see the different plots for each one
  abund = 3, # selecting the 3 most dominant taxa in each sample and showing their proportions in the other samples
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T
)


PlotAbund(
  TAX.pooled.rel$Order, # change levels to see the different plots for each one
  abund = 5, # selecting the 5 most dominant taxa in each sample and showing their proportions in the other samples
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T
)


#####
##### calculating the most abundant taxa out of the clean dataset of foraminifera only ######

ASV.sing.clean.rel <- prop.table(ASV.sing.clean, 2) * 100

b <- PlotAbund(
  ASV.sing.clean.rel, # change levels to see the different plots for each one
  abund = 3, # selecting the 3 most dominant taxa in each sample and showing their proportions in the other samples
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T,
  save.OTU = T
)

dim(b)
# [1]  93 113 # 92 ASVs, using the 3 most abundant

nrow(b)/nrow(ASV.sing.clean)*100 # 1.78% of the total number of ASVS

summary(t(b[,-ncol(b)])[,"other"]) # too see percentage for others

# 3 most abundant taxa
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.359  39.089  48.277  48.402  60.547  87.044 

# These ASVs are 52% of the total number of ASVs in the foraminifera samples, and represents only 1.78%  of the community diversity
# this means that the percentage for others is too high (48%), consider more than 3 abundant ASVs

# using the 20 more abundant 
b <- PlotAbund(
  ASV.sing.clean.rel, # change levels to see the different plots for each one
  abund = 20, # selecting the 3 most dominant taxa in each sample and showing their proportions in the other samples
  method = "nmost", # or nmost
  open.window = T,
  plot.ratio = c(3.5, 1),
  sort.taxa = T,
  save.OTU = T
)

dim(b)
# [1]  766 113 # 765 ASVs, using the 20 most abundant

nrow(b)/nrow(ASV.sing.clean)*100 # 14.69% of the total number of ASVS

summary(t(b[,-ncol(b)])[,"other"]) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.411  11.523  17.757  18.262  23.716  59.221 

# the proportion of "others" is still high (18%) even using 20 most abundant ASVs, this reveals that the ASVs a high variable i nthe prokaryotic microbiome

# obtaining ASV table just for these abundant ASVs

abundant.taxa <- rownames(b)[-nrow(b)]
ASV.abundant.taxa <- ASV.sing.clean[abundant.taxa,]

# same for the tax table
TAX.abundant.taxa <- TAX.sing.clean[abundant.taxa,]

# How much percent of the diversity these ASVs represent?
sum(ASV.abundant.taxa)/sum(ASV.sing.clean) # 0.78
nrow(TAX.abundant.taxa)/nrow(TAX.sing.clean) # 0.14
#####
#### old donut plot code (for ggplot plot) ####
# pooling tax data as list
TAX.pooled <- vector(mode = "list", length = 6)
names(TAX.pooled) <- colnames(TAX.clean)[1:6] # kingdom to genus
for (i in 1:6) {
  temp <- aggregate(
    ASV.clean,
    by = list(TAX.clean[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# preparing data for plotting

# take table with TAX.pooled values within the Phylum levels 
data <- melt(TAX.pooled$Phylum) # change level as wished
colnames(data) <- c("Phylum", "Extraction_Voucher", "Value")

# join metadata
data.meta <- merge(x=data, y=META, by="Extraction_Voucher") 

# create table with sum of value per group
data.meta.st <- dcast(data.meta, Sample_type + Phylum ~ ., sum, value.var = "Value") # this way is better for plotting
colnames(data.meta.st) <- c("Sample_type", "Phylum", "Value")

### donuts plot 
data.meta.st.don <- data.meta.st # make copy of data

## Compute percentages per sample type

# select sample type (st) to plot
st <- c("filter")
st <- c("sediment")
st <- c("single_cell")
data.meta.st.don.i <- data.meta.st.don[data.meta.st.don$Sample_type == st,] 

# include row for Class Bacillariophyta
META.st <- META[META$Sample_type == st,]
ASV.st <- ASV.clean[,rownames(META.st)]
ASV.st <- ASV.st[rowSums(ASV.st) > 0,] # remove ASVs with zero seq
TAX.st <- TAX.clean[rownames(ASV.st),]

# calculating how many ASVs from class alphaproteobacteria and gammaproteobacteria within this sample type
TAX.st.alpha <- TAX.st[TAX.st[, "Class"] == "Alphaproteobacteria", ] 
ASV.st.alpha <- ASV.st[rownames(TAX.st.alpha),]

class_alpha <- sum(ASV.st.alpha)
class_alpha_df <- data.frame(col1=st,
                             col2="Proteobacteria (Alpha)",
                             col3=class_alpha)

TAX.st.gamma <- TAX.st[TAX.st[, "Class"] == "Gammaproteobacteria", ] 
ASV.st.gamma <- ASV.st[rownames(TAX.st.gamma),]

class_gamma <- sum(ASV.st.gamma)
class_gamma_df <- data.frame(col1=st,
                             col2="Proteobacteria (Gamma)",
                             col3=class_gamma)


# subtracting to obtain how many Proteobacteria are not from alpha and gamma classes
sum_alpha_gamma <- class_alpha_df$col3 + class_gamma_df$col3

proteobacteria_other <- data.meta.st.don.i[data.meta.st.don.i$Phylum == "Proteobacteria",]$Value - sum_alpha_gamma
proteobacteria_other_df <- data.frame(col1=st, 
                                      col2="Proteobacteria (other classes)", 
                                      col3=proteobacteria_other)
# summarize in a table
proteobacteria_df <- rbind(class_alpha_df, class_gamma_df, proteobacteria_other_df)
colnames(proteobacteria_df) <- c("Sample_type", "Phylum", "Value")

# remove roteobacteria row from original table (to then add the dat divided in the three categories we just calculated)
data.meta.st.don.i <- data.meta.st.don.i[data.meta.st.don.i$Phylum != "Proteobacteria",]

# combining now the ochrophyta values splitted as we want
data.meta.st.don.i <- rbind(data.meta.st.don.i, proteobacteria_df)

# sorting in ascending order
str(data.meta.st.don.i)
data.meta.st.don.i$Phylum <- as.character(data.meta.st.don.i$Phylum) # has to convert as character before using order()
data.meta.st.don.i <- data.meta.st.don.i[order(data.meta.st.don.i$Phylum),] 


# calculate percentages of taxa within sample type
data.meta.st.don.i$fraction = data.meta.st.don.i$Value / sum(data.meta.st.don.i$Value) # calculate percentages of taxa within sample type

# summarize taxa with abundance bellow 0.01 as "Other divisions"
data.meta.st.don.i.rare <- data.meta.st.don.i[data.meta.st.don.i$fraction < 0.01,]
sumValue.other.divisions <-  sum(data.meta.st.don.i.rare$Value)
data.other.divisions <-  sum(data.meta.st.don.i.rare$fraction)
data.other.divisions_df <- data.frame(col1=st, 
                                      col2="Other phyla", 
                                      col3=sumValue.other.divisions,
                                      col4=data.other.divisions)
colnames(data.other.divisions_df) <- c("Sample_type", "Phylum", "Value", "fraction")

# removing rare taxa of original table, for plotting 
data.meta.st.don.i.abund <- data.meta.st.don.i[data.meta.st.don.i$fraction >= 0.01,]

# combining with value of "other divisions"
data.meta.st.don.i <- rbind(data.meta.st.don.i.abund, data.other.divisions_df)

sum(data.meta.st.don.i$fraction)
#[1] 1 - it worked!

# Compute the cumulative percentages (top of each rectangle)
data.meta.st.don.i$ymax = cumsum(data.meta.st.don.i$fraction)

# Compute the bottom of each rectangle
data.meta.st.don.i$ymin = c(0, head(data.meta.st.don.i$ymax, n=-1))

# Make the plot
d1 <- ggplot(data.meta.st.don.i, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Phylum)) +
  geom_rect() +
  #scale_fill_brewer(palette = 4)+
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15))
#theme(legend.position = "none")
#ggtitle("Seawater") +
# theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

d2 <- ggplot(data.meta.st.don.i, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Phylum)) +
  geom_rect() +
  #scale_fill_brewer(palette = 4)+
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15))
#theme(legend.position = "none")
#ggtitle("Sediment")+
# theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

d3 <- ggplot(data.meta.st.don.i, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Phylum)) +
  geom_rect() +
  #scale_fill_brewer(palette = 4)+
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15))
#ggtitle("Foraminifera")+
# theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))+ 
# guides(shape = guide_legend(override.aes = list(size = 0.5))) 


donut.plot <- d1 + d2 + d3

#####





##########################################################
################### PlotAbund function ################### 
##########################################################
PlotAbund <- function(relData, abund, 
                      margin = par()$mar,
                      method = c("nmost", "precentage"),
                      colorPalette = c("violet", "purple", "blue", "white", "darkgreen", "yellow", "red", "darkred"),
                      plot.ratio = c(3, 1),
                      open.window = T,
                      save.OTU = F,
                      sample.names = colnames(relData),
                      sort.taxa = F,
                      ...) {
  if(method == "nmost") {
    abund_names <- c()
    for (i in 1:ncol(relData)) {
      abund_names <- unique(c(abund_names, rownames(relData)[order(relData[, i], decreasing = T)][1:abund]))
    }
    abund_rel0 <- relData[abund_names, ]
  }
  if(method == "percentage") {
    abund_rel0 <- relData[apply(relData, 1, function(x) { max(x) >= abund }), ]
  }
  if (sort.taxa == T) {
    abund_rel0 <- abund_rel0[order(rownames(abund_rel0)), ]
  }
  abund_rel <- rbind(abund_rel0, 100 - colSums(abund_rel0))
  rownames(abund_rel)[nrow(abund_rel)] <- "other"
  abund_rel <- as.matrix(abund_rel)
  if (open.window == T) {
    grDevices::windows(width = 20, height = 10) # I had to call the package because my R was not getting the right function
  }
  # par(mfrow = c(1, 2), mar = margin, xpd = NA)
  layout(mat = matrix(c(1, 2), 1, 2, byrow = T), widths = plot.ratio)
  par(mar = margin)
  barplot(abund_rel,
          col = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
          ylim = c(0, 100), 
          las = 2,
          ylab = "Relative sequence abundance [%]",
          names.arg = sample.names,
          cex.names = 0.8,
          ...)
  par(mar = c(1,1,1,1))
  plot.new()
  legend("center", 
         pch = 22, 
         col = "black", 
         pt.bg = rev(c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey")),
         legend = rev(rownames(abund_rel)),
         pt.cex = 1.5, 
         cex = 0.8
  ) 
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(mfrow = c(1,1))
  
  if(save.OTU == T) {
    return(
      data.frame(
        abund_rel,
        color = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
        stringsAsFactors = F
      )
    )
  }
}





##########################################################
################# ANOSIMposthoc function ################# 
##########################################################

ANOSIMposthoc=function(M,E,distance="bray", padj="fdr"){
  E=droplevels(E)
  Mlist=list()
  for(i in 1:length(levels(E))){
    Mlist[[i]]=M[E==levels(E)[i],]
  }
  
  Elist=list()
  for(i in 1:length(levels(E))){
    Elist[[i]]=as.numeric(E[E==levels(E)[i]])
  }
  
  result=list(anosimR=matrix(NA,length(levels(E)),length(levels(E))),
              anosimP=matrix(NA,length(levels(E)),length(levels(E))),
              anosimPadj=matrix(NA,length(levels(E)),length(levels(E))))
  colnames(result$anosimR)=colnames(result$anosimP)=levels(E)
  rownames(result$anosimR)=rownames(result$anosimP)=levels(E)
  for(i in 1:(length(levels(E))-1)){
    for(j in (i+1):length(levels(E))){
      temp=anosim(rbind(Mlist[[i]],Mlist[[j]]),c(Elist[[i]],Elist[[j]]),distance=distance)
      result$anosimR[j,i]=temp$statistic
      result$anosimP[j,i]=temp$signif
    }
  }
  result$anosimPadj=matrix(p.adjust(as.vector(result$anosimP),method=padj,
                                    n=length(which(!is.na(as.vector(result$anosimP))))),
                           length(levels(E)),length(levels(E)))
  colnames(result$anosimPadj)=levels(E)
  rownames(result$anosimPadj)=levels(E)
  return(result)
}






