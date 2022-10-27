#################################
# Analysis of 18S amplicon data #
#################################

Sys.setenv(LANG = "en") # in case system language changes to german

##### prepare environment ####

# set working directory
setwd("G:/My Drive/Research/PhD/Data/Data_analysis/Metabarcoding_analyses/R_out_euks") # Windows
setwd("~/Library/CloudStorage/GoogleDrive-deboraposo@gmail.com/My Drive/Research/PhD/Data/Data_analysis/Metabarcoding_analyses/R_out_euks") # Mac


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
require(rcartocolor) # color blind pallete "Safe"
require("UpSetR") # upset plots
citation("vegan")

# loading and saving workspace
load("patterns_exploration_clean_euks.Rdata")
#save.image("patterns_exploration_clean_euks.Rdata")

#####
##### read and sort data #####
# R objects generated after data cleaning/merging in "data_exploration_all_euks.R"

ASV <- readRDS("ASV_merge_euks.RDS")
TAX <- readRDS("TAX_merge_euks.RDS")
META <- readRDS("META_merge_euks.RDS")

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

# saving META.sing RDS file
#saveRDS(META.sing, "METAsing.RDS")

##### 

##### rarefaction curve - to plot sample type #####

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
##### rarefaction curve - to plot site ######
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

##### Cleaning taxa from single cell and environmental data sets ######

# check divisions present in foraminifera that could not be symbionts be substrate 
table(TAX.sing[,"Division"])

# Alveolata_unclassified                Apicomplexa                   Cercozoa                Chlorophyta 
# 2                          4                         11                         50 
# Ciliophora             Dinoflagellata     Eukaryota_unclassified                      Fungi 
# 8                         70                         64                          6 
# Lobosa              Mesomycetozoa                    Metazoa                 Ochrophyta 
# 3                          5                         75                        871 
# Opisthokonta_unclassified          Prasinodermophyta                Pseudofungi                 Rhodophyta 
# 1                          2                         10                         27 
# Sagenista Stramenopiles_unclassified               Streptophyta 
# 47                         15                         18 


# there are some taxa that should be remove!

# unclassified Divisions - does not add in the discussion since we cannot even identify the division

# Rhodophyta Division - all species identified are from pluricelular organisms / macro algae
#therefore they represent the algae in which the foraminifera lives attached to, unlikely symbionts

# Embryophycea Class, Streptophyta Division - large group of green plant, the taxa found in our dataset are from land plants or unclassified.
#therefore, they are not likely symbionts / does not add in the discussion

# Metazoa Division - too big to be symbionts


# Fungi Division - the taxa we found are normally found on animal and human skin . we can remove it

# Division Chlorophyta - all species identified are from pluricelular organisms / macro algae
#therefore they represent the algae in which the foraminifera lives attached to, unlikely symbionts

# Pseudo Fungi Division - only taxa that are considered parasites of macro algae. since this taxa did not occur in eilat
# the only site we dont see macro algae, we believe this indicates that the source of this ASVs are indeed from the environment and
# therefore, they should be removed.
# Ciliophora Division - the taxa found  represents benthic ciliates that are also too big to ve symbionts (~250 ?m). therefore, remove it


# check divisions to remove for entire dataset
table(TAX[,"Division"])

# make list of Division to remove
divisions_to_remove <- c("Alveolata_unclassified", "Archaeplastida_unclassified", "Ciliophora", "Chlorophyta", "Eukaryota_unclassified",
                             "Metazoa", "Fungi", "Hacrobia_unclassified", "Pseudofungi", "Opisthokonta_unclassified", "Opisthokonta_X", "Rhodophyta",
                             "Stramenopiles_unclassified", "Stramenopiles_X", "Streptophyta"
)


# removing these division
TAX.clean <- TAX[!TAX[, "Division"] %in% divisions_to_remove_all, ] 

# check classes to remove after already cleaning divisions
grep(pattern="_unclassified",x=names(table(TAX.clean[,"Class"])),value=TRUE)

grep(pattern="_X",x=names(table(TAX.clean[,"Class"])),value=TRUE)

# make list of Classes to remove
classes_to_remove <- c("Ochrophyta_unclassified", "Apicomplexa_unclassified", "Cercozoa_unclassified", "Choanoflagellida_unclassified",
                       "Dinoflagellata_unclassified", "Haptophyta_unclassified", "Opalozoa_unclassified", "Sagenista_unclassified",  
                       "Ancoracystida_X", "Centroheliozoa_X", "Cercozoa_X", "Choanoflagellida_X", "Dinophyta_X", "Fornicata_X",
                       "Lobosa_X", "Picozoa_X", "Telonemia_X"
)

# removing these classes
TAX.clean <- TAX.clean[!TAX.clean[, "Class"] %in% classes_to_remove, ] 
#TAX.clean <- TAX.clean[!TAX.clean[, "Family"] %in% family_to_remove, ] 

View(unique(TAX.clean[, 1:4])) # ok, no unclassified o _X Division and Classes

# checking families
# grep(pattern="_unclassified",x=names(table(TAX.clean[,"Family"])),value=TRUE)
# [1] "Acantharea_2_unclassified"           "Acantharea_E_unclassified"           "Acantharea_F_unclassified"           "Acantharea_unclassified"            
# [5] "Acanthoecida_unclassified"           "Apicomplexa_unclassified"            "Apusomonadidae_Group-2_unclassified" "Bicoecea_unclassified"              
# [9] "Centroheliozoa_X_unclassified"       "Cercozoa_unclassified"               "Chlorarachnida_unclassified"         "Choanoflagellatea_unclassified"     
# [13] "Choanoflagellida_unclassified"       "Chrysophyceae_X_unclassified"        "Coccidiomorphea_unclassified"        "Coccolithales_unclassified"         
# [17] "Cryptofilida_unclassified"           "Cryptophyceae_unclassified"          "Dictyochophyceae_X_unclassified"     "Dino-Group-II_unclassified"         
# [21] "Dinoflagellata_unclassified"         "Dinophyceae_unclassified"            "Ebriida_unclassified"                "Eimeriida_unclassified"             
# [25] "Endomyxa-Phytomyxea_unclassified"    "Euglyphida_unclassified"             "Eugregarinorida_unclassified"        "Filosa-Granofilosea_unclassified"   
# [29] "Filosa-Imbricatea_unclassified"      "Filosa-Thecofilosea_unclassified"    "Gonyaulacales_unclassified"          "Gymnodiniales_unclassified"         
# [33] "Haptophyta_unclassified"             "Ichthyosphonida_unclassified"        "Labyrinthulomycetes_unclassified"    "Labyrinthulomycetes_X_unclassified" 
# [37] "MAST-12_unclassified"                "MAST-3_unclassified"                 "MAST-4_unclassified"                 "MAST-7_unclassified"                
# [41] "Opalozoa_unclassified"               "Pelagophyceae_unclassified"          "Peridiniales_unclassified"           "Prymnesiales_unclassified"          
# [45] "Prymnesiophyceae_unclassified"       "Pterocystida_unclassified"           "Sagenista_unclassified"              "Sarcinochrysidales_unclassified"    
# [49] "Spumellaria_unclassified"            "Suessiales_unclassified"             "Syndiniales_unclassified"            "Syracosphaerales_unclassified"      
# [53] "Tubulinea_unclassified"              "Vampyrellida_unclassified"  
# 
# grep(pattern="_X",x=names(table(TAX.clean[,"Family"])),value=TRUE)
# [1] "Acantharea_A_X"                     "Acanthoecida_X"                     "Ancoracystida_XXX"                  "Anoecales_X"                       
# [5] "Apusomonadidae_Group-1_XX"          "Apusomonadidae_Group-2B_X"          "Calcihaptophycidae_X"               "Centroheliozoa_X_unclassified"     
# [9] "Centroheliozoa_XXX"                 "Cercozoa_XXX"                       "Chlorarachnida_X"                   "Choanoflagellida_XX_Clade-3"       
# [13] "Chrysophyceae_X_unclassified"       "Cryomonadida_X"                     "Cryptomonadales_X"                  "Dictyochophyceae_X_unclassified"   
# [17] "Dictyochophyceae_XX"                "Dino-Group-I_X"                     "Dino-Group-II_X"                    "Dino-Group-III_X"                  
# [21] "Dino-Group-V_X"                     "Dinophyceae_XX"                     "Dinophyta_XXX"                      "Ebriida_X"                         
# [25] "Euglyphida_X"                       "Filosa_XX"                          "Filosa-Imbricatea_XX"               "Filosa-Thecofilosea_XX"            
# [29] "Goniomonadales_X"                   "Gregarinomorphea_X_GRE3_X"          "Gregarinomorphea_XX"                "Haptophyta_Clade_HAP2_XX"          
# [33] "Haptophyta_Clade_HAP3_XX"           "Haptophyta_Clade_HAP4_XX"           "Haptophyta_Clade_HAP5_XX"           "Ichthyosphonida_X"                 
# [37] "Katablepharidales_X"                "Labyrinthulomycetes_X_LAB1/6/8"     "Labyrinthulomycetes_X_LAB14"        "Labyrinthulomycetes_X_LAB7"        
# [41] "Labyrinthulomycetes_X_unclassified" "Mantamonadida_XX"                   "Marimonadida_X"                     "MAST-10_XX"                        
# [45] "MAST-11_XX"                         "MAST-12_XX"                         "MAST-12A_X"                         "MAST-12B_X"                        
# [49] "MAST-12D_X"                         "MAST-12E_X"                         "MAST-20_XX"                         "MAST-3_XX"                         
# [53] "MAST-3A_X"                          "MAST-3B_X"                          "MAST-3C_X"                          "MAST-3D_X"                         
# [57] "MAST-3E_X"                          "MAST-3F_X"                          "MAST-3I_X"                          "MAST-3J_X"                         
# [61] "MAST-3K_X"                          "MAST-3L_X"                          "MAST-4A_X"                          "MAST-4B_X"                         
# [65] "MAST-4C_X"                          "MAST-4D_X"                          "MAST-4E_X"                          "MAST-6_XX"                         
# [69] "MAST-7B_X"                          "MAST-7C_X"                          "MAST-7D_X"                          "MAST-7E_X"                         
# [73] "MAST-8A_X"                          "MAST-8B_X"                          "MAST-8C_X"                          "MAST-9A_X"                         
# [77] "MAST-9C_X"                          "MAST-9D_X"                          "Microhelida_X"                      "MOCH-1_XX"                         
# [81] "MOCH-2_XX"                          "MOCH-3_XX"                          "MOCH-4_XX"                          "MOCH-5_XX"                         
# [85] "Novel-clade-9_X"                    "Parmales_X"                         "Peridiniales_X"                     "Perkinsida_XX"                     
# [89] "Phaeophyceae_XX"                    "Picozoa_XXX"                        "Prymnesiales_X"                     "Prymnesiophyceae_Clade_D_X"        
# [93] "Prymnesiophyceae_Clade_E_X"         "Prymnesiophyceae_XX"                "Pterocystida_X"                     "RAD-A_Xa"                          
# [97] "RAD-A_Xf"                           "RAD-B_X_Group-II"                   "RAD-B_X_Group-IVe"                  "Raphidophyceae_XX"                 
# [101] "Syndiniales_XX"                     "Syracosphaerales_X"                 "Vampyrellida_X"                     "Ventricleftida_X"     
# 

#family_to_remove <- c("Bacillariophyta_X_unclassified") 
# not removing unclassified families anymore, decided to clean only on the Division and Class level

# too many families to remove manualy
# not a good idea, might change too much the data set


# applying cleaning to ASV
ASV.clean <- ASV[rownames(TAX.clean),]

# create clean foraminifera table
ASV.sing.clean <- ASV.clean[,rownames(META.sing)]
ASV.sing.clean <- ASV.sing.clean[rowSums(ASV.sing.clean) > 0,] # remove ASVs with zero seq
TAX.sing.clean <- TAX.clean[rownames(ASV.sing.clean),]


# save RDS files to upload in next scripts
saveRDS(ASV.clean, "ASV_merge_clean_euks.RDS")
saveRDS(TAX.clean, "TAX_merge_clean_euks.RDS")


#####
##### number of sequences retained ######

# foraminifera
nrow(ASV.sing.clean)   
# [1] 1005

nrow(ASV)
nrow(ASV.clean)
# sediment
ASV.sed <- ASV.clean[,META$Sample_type == "sediment"]
ASV.sed <- ASV.sed[rowSums(ASV.sed) >0, ] 
nrow(ASV.sed)
# [1] 1539
TAX.sed <- TAX.clean[rownames(ASV.sed),]



# seawater
ASV.filters <- ASV.clean[,META$Sample_type == "filter"]
ASV.filters <- ASV.filters[rowSums(ASV.filters) >0, ] # excluding lines with no ASV
nrow(ASV.filters)
# [1] 3085
TAX.filters <- TAX.clean[rownames(ASV.filters),]


#####
##### proportion of Ochrophyta and Bacillariophyta in foraminifera samples ####
# frist checking how much they represent in the Division Ochrophyta
TAX.alg <- TAX.sing.clean[TAX.sing.clean[, "Division"] == "Ochrophyta", ] 
TAX.diat <- TAX.sing.clean[TAX.sing.clean[, "Class"] == "Bacillariophyta", ]

nrow(TAX.alg)/nrow(TAX.sing.clean)
# [1] 0.8571 - Ochrophyta is 86% of the foraminiferal microbiome diversity
ASV.alg <- ASV.sing.clean[rownames(TAX.alg),]
ASV.alg <- ASV.alg[rowSums(ASV.alg) > 0,]
sum(ASV.alg)/sum(ASV.sing.clean)
# 0.9936947 # Ochrophyta represents 99.4% of the total number of sequences

nrow(TAX.diat)/nrow(TAX.alg)
# [1] 0.997 - diatoms are 99.7% of the Ochrophyta division we can work only with them!

# checking araphid pennates (most abundant family)
TAX.araphid <- TAX.diat[TAX.diat[, "Family"] == "Araphid-pennate", ] 

ASV.araphid <- ASV.alg[rownames(TAX.araphid),]
ASV.araphid <- ASV.araphid[rowSums(ASV.araphid) > 0,]
sum(ASV.araphid)/sum(ASV.alg) 
# 0.9881067 # Araphid family represents 98.8% of the total number of sequences

#####

##### Fig donut plotting all taxa (donut and tree map plots) #####

# creating grouping vector based on division

tax.group <- TAX.clean[ ,"Division"]
tax.group[TAX.clean[,"Division"] == "Ochrophyta"] <- ifelse(
  TAX.clean[TAX.clean[,"Division"] == "Ochrophyta", "Class"] == "Bacillariophyta", 
  "Ochrophyta (Bacillariophyta)",
  "Ochrophyta (other)"
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
rownames(ASV.group.st.rel.plot)[nrow(ASV.group.st.rel.plot)] <- "Other divisions" # rownames is a vector, so we use nrow to het the final row

# plotting
par(mfrow=c(1,4), mar = rep(0,4)) # to have three pie charts (one per sample type) and the legend on the side

for (i in 1:3){ # to make one per sample type
  pie(
    ASV.group.st.rel.plot[,i], 
    col = carto_pal(nrow(ASV.group.st.rel.plot), "Safe"), 
    border = NA, 
    labels = NA, 
    clockwise = T
  )
}
plot.new() # so legend enter as a new plot
legend(
  "center",
  legend = rownames(ASV.group.st.rel.plot),
  col = carto_pal(nrow(ASV.group.st.rel.plot), "Safe"),
  pch = 15,
  cex = 1,
  pt.cex = 2,
  ncol = 1,
  bty = "n",
  title = "Division",
  title.adj = 0,
  title.cex = 1.5
)
dev.off()
?legend()
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

tax.group <- TAX.clean.s[ ,"Division"]
tax.group[TAX.clean.s[,"Division"] == "Ochrophyta"] <- ifelse(
  TAX.clean.s[TAX.clean.s[,"Division"] == "Ochrophyta", "Class"] == "Bacillariophyta", 
  "Ochrophyta (Bacillariophyta)",
  "Ochrophyta (other)"
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

# calculating proportion of "Other Phyla" based on taxa that was kicked out in the subsetting above
ASV.group.st.rel.other.phyla <- t(as.data.frame(colSums(ASV.group.st.rel) - colSums(ASV.group.st.rel.taxa.general))) 

ASV.group.st.rel.plot.s <- rbind(ASV.group.st.rel.taxa.general, ASV.group.st.rel.other.phyla)
rownames(ASV.group.st.rel.plot.s)[nrow(ASV.group.st.rel.plot.s)] <- "Other phyla" # rownames is a vector, so we use nrow to het the final row

# plotting
par(mfrow=c(1,4), mar = rep(0,4)) # to have three pie charts (one per sample type) and the legend on the side

for (i in 1:3){ # to make one per sample type
  pie(
    ASV.group.st.rel.plot.s[,i], 
    col = carto_pal(nrow(ASV.group.st.rel.plot.s), "Safe"), 
    border = NA, 
    labels = NA, 
    clockwise = T
  )
}
plot.new() # so legend enter as a new plot
legend(
  "center",
  legend = rownames(ASV.group.st.rel.plot.s),
  col = carto_pal(nrow(ASV.group.st.rel.plot.s), "Safe"),
  pch = 15,
  cex = 1,
  pt.cex = 2,
  ncol = 1,
  bty = "n",
  title = paste("Division - ", s),
  title.adj = 0,
  title.cex = 1.5
)


#####
##### boxplots of rare taxa for forams ####

Ochrophyta_to_remove <- c("Ochrophyta")

# removing these taxa
tax.rare <- TAX.sing.clean[!TAX.sing.clean[ ,"Division"] %in% Ochrophyta_to_remove,]
ASV.rare.rel <- ASV.sing.clean.rel[rownames(tax.rare),]
tax.rare <- tax.rare[, "Division"] # convert to vector with Divisions to use on aggregate function bellow

ASV.rare <- aggregate(
  ASV.rare.rel,
  by = list(tax.rare),
  FUN = sum
)
rownames(ASV.rare) <- ASV.rare$Group.1
ASV.rare <- t(ASV.rare[, -1])

# I need a df with 4 columns: Sample, Site, Division and Proportion

data.barplot <- as.data.frame(ASV.rare)
all.equal(rownames(data.barplot), rownames(META.sing))
# [1] TRUE
data.barplot$Site <- META.sing$Site

# melt to four columns (Sample, Site, Division and Proportion)
data.barplot <- rownames_to_column(data.barplot)

data.barplot <- melt(data.barplot,  id.vars = c("rowname", "Site"), # try id.vars = c("family_id", "age_mother"), before: id.vars = "rowname"
                                   measure.vars = c("Apicomplexa", "Cercozoa", "Dinoflagellata", "Lobosa", 
                                                    "Mesomycetozoa", "Prasinodermophyta", "Sagenista"))
colnames(data.barplot) <- c("Sample", "Site", "Division", "Proportion")
str(data.barplot)

# taking note of colors in the same order that they appeared in donut plot
carto_pal(name = "Safe")
"#88CCEE" # Apicomplexa # from Safe
"#CC6677" # Cercozoa # from Safe
"#332288" # Dinoflag # from Safe
"#D3D3D3" # light gray Lobosa # new colors for rare taxa
"#808080" # gray Mesomyc # new colors for rare taxa
"#000000" # Black Prasino # new colors for rare taxa
"#6699CC" # Sagenista # from Safe

# box plot for divisions and facets as sites
str(data.barplot)
scaleFUN <- function(x) sprintf("%.2f", x)

plot5 <- ggplot(data.barplot, aes(x = Division, y = Proportion, fill = Division)
) +
  geom_boxplot() +
  labs(y = "Rare taxa in foraminifera (%)\n") +
  facet_grid(~Site) +
  scale_y_continuous(trans='log', labels = scaleFUN)+
  scale_fill_manual(values= c("#88CCEE", "#CC6677", "#332288", "#D3D3D3", "#808080", "#000000", "#6699CC"))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        axis.ticks.x = element_line(size = 0.5)
  )
plot5

#####

##### Fig NMDS all data - to plot sample type ######
ASV.clean.rel <- prop.table(ASV.clean, 2) * 100

# calculate NMDS
NMDS <- metaMDS(t(ASV.clean.rel), k = 2)
NMDS$stress
# 0.1229098 # all taxa

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
legend("topright", "Stress = 0.12", bty = "n", cex = 1) 

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

# calculate the inverse Simpson index without subsampling
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
##### ANOSIM all euks in all sample types ####

# calculate ANOSIM between sample type
anosim(t(ASV.clean.rel), META$Sample_type)

#anosim(t(ASV.clean.rel), META$Site)
?anosim()
# calculate pairwise ANOSIM statistics
# download function from: https://raw.githubusercontent.com/chassenr/Tutorials/master/R_roundtable_sequence_analysis/anosimPosthoc.R
source("C:/Users/chassenrueck/Documents/Repos/Tutorials/R_roundtable_sequence_analysis/anosimPosthoc.R") # doesn't work, run the one in the script (end)
# run function anosimPosthoc in the end of the script
ANOSIMposthoc(t(ASV.clean.rel), META$Sample_type)

ANOSIMposthoc(t(ASV.clean.rel), META$Site)


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

##### upset plot - diatoms intersections for all sites together  ####
?upset()

# constraining the diatoms only
TAX.all.diat <- TAX.clean[TAX.clean[, "Class"] == "Bacillariophyta", ]
ASV.all.diat <- ASV.clean[rownames(TAX.all.diat),]
nrow(ASV.all.diat)
# [1] 1502

# create list with ASVs name per sample type
ASV.diat.sed <- ASV.all.diat[,rownames(META[META$Sample_type == "sediment",])] 
ASV.diat.sed <- ASV.diat.sed[rowSums(ASV.diat.sed) > 0,] 

ASV.diat.sw <- ASV.all.diat[,rownames(META[META$Sample_type == "filter",])] 
ASV.diat.sw <- ASV.diat.sw[rowSums(ASV.diat.sw) > 0,] 

ASV.diat.sing <- ASV.all.diat[,rownames(META[META$Sample_type == "single_cell",])] 
ASV.diat.sing <- ASV.diat.sing[rowSums(ASV.diat.sing) > 0,] 

upsetlist <- list(rownames(ASV.diat.sw), rownames(ASV.diat.sed), rownames(ASV.diat.sing))
names(upsetlist) <- c("Seawater", "Sediment", "Foraminifera")

# plot upset
upset.plot <- upset(fromList(upsetlist), sets = c("Foraminifera", "Sediment", "Seawater"), 
                    keep.order = T,  mainbar.y.label = "Sample type Intersections", 
                    main.bar.color = c("red3","forestgreen", "dodgerblue", "black", "black","black"),
                    # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
                    sets.bar.color = c("red3","forestgreen", "dodgerblue"),
                    sets.x.label = "Diatoms' ASVs", text.scale = c(2,2,2,1.5,2,2)
                    )

upset.plot
#####
##### upset plot - filtering one site per time for Sup Mat  ####

s <- c("Eilat")
s <- c("Tel Shikmona")
s <- c("Plemmirio")
s <- c("Capo Passero")
METAs <- META[META$Site == s,]
ASV.all.diat.s <- ASV.all.diat[,rownames(METAs)]
ASV.all.diat.s <- ASV.all.diat.s[rowSums(ASV.all.diat.s) > 0,] # excluding ASVs not present in this subset

# the rest is the same
ASV.diat.sed.s <- ASV.all.diat.s[,rownames(METAs[METAs$Sample_type == "sediment",])] 
ASV.diat.sed.s <- ASV.diat.sed.s[rowSums(ASV.diat.sed.s) > 0,] 

ASV.diat.sw.s <- ASV.all.diat.s[,rownames(METAs[METAs$Sample_type == "filter",])] 
ASV.diat.sw.s <- ASV.diat.sw.s[rowSums(ASV.diat.sw.s) > 0,] 

ASV.diat.sing.s <- ASV.all.diat.s[,rownames(METAs[METAs$Sample_type == "single_cell",])] 
ASV.diat.sing.s <- ASV.diat.sing.s[rowSums(ASV.diat.sing.s) > 0,] 

# make upset list
upsetlist.s <- list(rownames(ASV.diat.sw.s), rownames(ASV.diat.sed.s), rownames(ASV.diat.sing.s))
names(upsetlist.s) <- c("Seawater", "Sediment", "Foraminifera")

# plot upset
p1 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
                    keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 600,
                    main.bar.color = c("red3","dodgerblue", "forestgreen", "black", "black","black"),
                    # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
                    sets.bar.color = c("red3","forestgreen", "dodgerblue"),
                    sets.x.label = "Diatoms' ASVs", text.scale = c(2,2,2,1.5,2,2)
            )
p1 # eilat

p2 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 600,
            main.bar.color = c("red3","dodgerblue", "forestgreen", "black", "black","black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Diatoms' ASVs", text.scale = c(2,2,2,1.5,2,2)
            )
p2 # tel shikmona

p3 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 600,
            main.bar.color = c("red3","forestgreen", "dodgerblue", "black", "black","black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Diatoms' ASVs", text.scale = c(2,2,2,1.5,2,2)
            )
p3 # plemmirio

p4 <- upset(fromList(upsetlist.s), sets = c("Foraminifera", "Sediment", "Seawater"), 
            keep.order = T,  mainbar.y.label = paste(s, "\nSample Type Intersections"), mainbar.y.max = 600,
            main.bar.color = c("red3","forestgreen", "dodgerblue", "black", "black","black", "black"),
            # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
            sets.bar.color = c("red3","forestgreen", "dodgerblue"),
            sets.x.label = "Diatoms' ASVs", text.scale = c(2,2,2,1.5,2,2)
)
p4 # capo passero


#####
##### upset plot - all euks intersections for all sites together  ####
?upset()
33/4842
11/1502
# constraining the diatoms only
nrow(ASV.all.diat)
# [1] 1502
nrow(ASV.clean)
# [1] 4842 # much more ASVs then when constraining only the diatoms

# create list with ASVs name per sample type
ASV.clean.sed <- ASV.clean[,rownames(META[META$Sample_type == "sediment",])] 
ASV.clean.sed <- ASV.clean.sed[rowSums(ASV.clean.sed) > 0,] 

ASV.clean.sw <- ASV.clean[,rownames(META[META$Sample_type == "filter",])] 
ASV.clean.sw <- ASV.clean.sw[rowSums(ASV.clean.sw) > 0,] 

ASV.clean.sing <- ASV.clean[,rownames(META[META$Sample_type == "single_cell",])] 
ASV.clean.sing <- ASV.clean.sing[rowSums(ASV.clean.sing) > 0,] 

upsetlist <- list(rownames(ASV.clean.sw), rownames(ASV.clean.sed), rownames(ASV.clean.sing))
names(upsetlist) <- c("Seawater", "Sediment", "Foraminifera")

# plot upset
upset.plot <- upset(fromList(upsetlist), sets = c("Foraminifera", "Sediment", "Seawater"), 
                    keep.order = T,  mainbar.y.label = "Sample type Intersections", 
                    main.bar.color = c("red3","forestgreen", "dodgerblue", "black", "black","black", "black"),
                    # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
                    sets.bar.color = c("red3","forestgreen", "dodgerblue"),
                    sets.x.label = "Eukaryotes ASVs", text.scale = c(2,2,2,1.5,2,2)
)

upset.plot
#####
##### NMDS diatoms single-cell only - to plot sites  #######
ASV.diat.rel <- prop.table(ASV.sing.diat, 2) * 100
# calculate NMDS
NMDS <- metaMDS(t(ASV.diat.rel), k = 3) # k number of dimensions
NMDS$stress
#  0.1732803 # single cell diatoms all
# stress is lower than with 2 dimensions (0.2581434)

# set meta table
# META.sing$Depth_cat <- factor(META.sing$Depth_cat, levels=c("Shallow", "Deep"))
# META.sing$Substrate <- as.factor(META.sing$Substrate)

META.plot <- META.sing
str(META.plot)


# plotting site as colors, substrate as shape and depth as size

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
legend("bottomright", "stress = 0.17", bty = "n", cex = 1.3)

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
legend("bottomright", "stress = 0.17", bty = "n", cex = 1.3)

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
legend("bottomright", "stress = 0.17", bty = "n", cex = 1.3)

# add legend for sites, substrates and depth

legend("bottomright", legend = c(levels(META.plot$Site), levels(META.plot$Substrate), levels(META.plot$Depth_cat)),
       pch = c(15, 15, 15, 15, 21, 22, 16, 16),
       col = c("gold1", "salmon","darkorchid1","blue","black","black","black","black"),
       bty = "n",
       pt.cex = c(2,2,2,2,2,2,1,2))

##### 
##### PERMANOVA and posthoc for diatoms in foraminifera ######
# change for the sub set dataset
META.adonis <- META.sing
str(META.adonis)
META.adonis$Site <- droplevels(META.adonis$Site)
META.adonis$Condition <- droplevels(META.adonis$Condition)
?adonis2()
# with all factors

adonis2(t(ASV.diat.rel) ~ Site * Substrate * Depth_cat, data = META.adonis, sqrt.dist = T)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = t(ASV.diat.rel) ~ Site * Substrate * Depth_cat, data = META.adonis, sqrt.dist = T)
# Df SumOfSqs      R2       F Pr(>F)    
# Site                       3   16.627 0.34900 18.8811  0.001 ***
# Substrate                  1    0.194 0.00407  0.6610  0.805    
# Depth_cat                  1    0.248 0.00520  0.8445  0.571    
# Site:Substrate             2    0.449 0.00942  0.7648  0.782    
# Site:Depth_cat             3    0.644 0.01352  0.7314  0.868    
# Substrate:Depth_cat        1    0.395 0.00829  1.3450  0.188    
# Site:Substrate:Depth_cat   2    0.611 0.01283  1.0413  0.376    
# Residual                  97   28.472 0.59766                   
# Total                    110   47.640 1.00000                   
# ---
#   

# depth and substrate not significant

# PERMANOVA per site
adonis2(t(ASV.diat.rel) ~ Site, data = META.adonis, sqrt.dist = T)


# Pairwise comparison of site
# Let's also have a look at pairwise comparisons for sites
site.comb <- combn(levels(META.adonis$Site), 2) # creating a triangular pairwise comparison
permanova.list <-apply(
  site.comb,
  2,
  function(x) {
    adonis2(
      t(ASV.diat.rel[, rownames(META.adonis)[META.adonis$Site %in% x]]) ~ Site, 
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
# cond1        cond2        R2        F     p p.adj
# 1 Capo Passero    Plemmirio 0.1815072 13.52723 0.001 0.001
# 2 Capo Passero Tel Shikmona 0.2565076 21.04522 0.001 0.001
# 3 Capo Passero        Eilat 0.3179143 20.97412 0.001 0.001
# 4    Plemmirio Tel Shikmona 0.2415915 19.75014 0.001 0.001
# 5    Plemmirio        Eilat 0.2957898 19.32140 0.001 0.001
# 6 Tel Shikmona        Eilat 0.3263612 22.28586 0.001 0.001

# all sites are different!


##### 

##### betadispersion within foraminifera, to check sites ####
betadispersion <- betadisper(
  vegdist(t(ASV.diat.rel)), 
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
        axis.text.y = element_blank(), # removing the text because it will be combined with the proks plot
  )



#####
##### frequency of the dominant diatoms per site, for haplotype network  #############

#  calculating the most abundant taxa out of the clean dataset

nrow(ASV.diat.rel)

b <- PlotAbund(
  ASV.diat.rel, # change levels to see the different plots for each one
  abund = 3, # selecting the 3 most dominant taxa in each sample and showing their proportions in the other samples
  method = "nmost", # or nmost
  open.window = F,
  plot.ratio = c(3.5, 1),
  sort.taxa = T,
  save.OTU = T
)

dim(b)
# 66 112    # 3 most abundant ASVs per sample


summary(t(b[,-ncol(b)])[,"other"]) # percentage for others

# 3 most abundant taxa
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9905  2.1524  3.2786  5.3463  6.4209 38.6009 

abundant.taxa <- rownames(b)[-nrow(b)]

# saving the vector file to a txt file in desktop (to use to select these seqs from fasta file)
#write(abundant.taxa, "abundant.taxa.foraminifera.txt")

# obtaining ASV table just for these abundant ASVs
ASV.abundant.taxa <- ASV.sing.clean[abundant.taxa,]

# same for the tax table
TAX.abundant.taxa <- TAX.sing.clean[abundant.taxa,]

# How much percent of the diversity these ASVs represent?
sum(ASV.abundant.taxa)/sum(ASV.sing.clean) # 0.9471
nrow(TAX.abundant.taxa)/nrow(TAX.sing.clean) # 0.072


# first rarefy to avoid interpretation based on the sequencing depth 

# selecting just the same family for haplotype network (main family Araphid pennate)
TAX.abundant.araphid <- TAX.abundant.taxa[TAX.abundant.taxa[, "Family"] == "Araphid-pennate", ] 
abundant.diatoms.araphid <- rownames(TAX.abundant.araphid)
write(abundant.diatoms.araphid, "abundant.diatoms.araphid.foraminifera.txt")# saving the vector file to a txt file in desktop (to use to select these seqs from fasta file)

ASV.abundant.araphid <- ASV.abundant.taxa[abundant.diatoms.araphid,]

min(rowSums(t(ASV.abundant.taxa))) # [1] 25565
max(rowSums(t(ASV.abundant.taxa))) # [1] 125494

min(rowSums(t(ASV.abundant.araphid))) # 1433
max(rowSums(t(ASV.abundant.araphid))) # 125367

dim(ASV.abundant.taxa)


# rarefying to the minimum sequence depth
?rarefy()

ASV.abundant.raref <- rrarefy(t(ASV.abundant.taxa), min(rowSums(t(ASV.abundant.taxa))))
ASV.abundant.araphid.raref <- rrarefy(t(ASV.abundant.araphid), min(rowSums(t(ASV.abundant.araphid))))

min(colSums(t(ASV.abundant.raref))) # [1] 25565
max(colSums(t(ASV.abundant.raref))) # [1] 25565

min(colSums(t(ASV.abundant.araphid.raref))) # [1] 1433
max(colSums(t(ASV.abundant.araphid.raref))) # [1] 1433

min(colSums(ASV.abundant.raref)) # minimum is not zero, no ASV lost
min(colSums(ASV.abundant.araphid.raref)) # minimum is not zero, no ASV lost

dim(ASV.abundant.raref)

min(ASV.abundant.taxa[ASV.abundant.taxa > 0]) # to see the minimum seq depth that is not zero

all.equal(rownames(META.sing), rownames(ASV.abundant.raref))
# [1] TRUE

agg <- aggregate(ASV.abundant.raref,
                 by = list(META.sing$Site),
                 FUN = sum)
saveRDS(agg, "ASV.abundant.raref.agg.site.RDS")



all.equal(rownames(META.sing), rownames(ASV.abundant.araphid.raref))
# [1] TRUE

agg2 <- aggregate(ASV.abundant.araphid.raref,
                  by = list(META.sing$Site),
                  FUN = sum)
saveRDS(agg2, "ASV.abundant.araphid.raref.agg.site.RDS")

sum(ASV.abundant.araphid.raref)
# [1] 159063
sum(ASV.abundant.araphid)
# [1] 6866113

##### 





##### NMDS eukaryotes single-cell only - to plot sites ####

ASV.clean.sing.rel <- prop.table(ASV.clean.sing, 2) * 100

# calculate NMDS
NMDS <- metaMDS(t(ASV.clean.sing.rel), k = 3) # k number of dimensions
NMDS$stress
#  0.1732298 # single cell eukaryotes all
# stress is lower than with 2 dimensions (0.2577527)

META.plot <- META.sing
str(META.plot)


# plotting site as colors, substrate as shape and depth as size

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
legend("bottomright", "stress = 0.17", bty = "n", cex = 1.3)

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
legend("bottomright", "stress = 0.17", bty = "n", cex = 1.3)

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
legend("bottomright", "stress = 0.17", bty = "n", cex = 1.3)

# add legend for sites, substrates and depth

legend("bottomright", legend = c(levels(META.plot$Site), levels(META.plot$Substrate), levels(META.plot$Depth_cat)),
       pch = c(15, 15, 15, 15, 21, 22, 16, 16),
       col = c("gold1", "salmon","darkorchid1","blue","black","black","black","black"),
       bty = "n",
       pt.cex = c(2,2,2,2,2,2,1,2))


#####



#### To to later ###########################################################################################################################
##### investigating samples that are not clustered together in the unifrac distances S0494 and S0688 ####
samples_to_investigate <- c("S0494", "S0688") # S0494 on the 29.3 already identified as an odd sample

ASV.filt <- as.matrix(ASV.diat.sing[, samples_to_investigate])
ASV.filt <- as.matrix(ASV.filt[rowSums(ASV.filt) >0,])
TAX.filt <- TAX.diat[rownames(ASV.filt),]

unique(TAX.filt[,"Family"])
# NEXT
# observe how are these samples in the barplot with taxonomic composition per sample
# use similar code then for Fig 1B

#look only the single cell Diatoms table
all.equal(rownames(TAX.diat), rownames(ASV.diat.sing))
# [1] TRUE
all.equal(rownames(META.sing), colnames(ASV.diat.sing))
# [1] TRUE

TAX.sing.diat.pooled <- vector(mode = "list", length = 8)
names(TAX.sing.diat.pooled) <- colnames(TAX.diat)[1:8] # kingdom to species
for (i in 1:8) {
  temp <- aggregate(
    ASV.diat.sing,
    by = list(TAX.diat[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.sing.diat.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# take table with TAX.sing.diat.pooled values within the Family levels 
data.diat.family <- melt(TAX.sing.clean.pooled$Family) #  <--- changed to table with all taxa to test
colnames(data.diat.family) <- c("Family", "Sample", "Value")

# transform to relative proportions per sample
# dcast per Sample ~ Division
data.diat.family.sample <- dcast(data.diat.family, Sample ~ Family)

# Sample as rownames
data.diat.family.sample <- as.data.frame(data.diat.family.sample)
data.diat.family.sample <- column_to_rownames(data.diat.family.sample, 'Sample')

str(data.diat.family.sample)

# transform in matrix to calculate proportion table
data.diat.family.sample <- as.matrix(data.diat.family.sample)

# calculate proportions
data.diat.family.sample.rel <- prop.table(data.diat.family.sample, 1) * 100

# include Site information from META
all.equal(rownames(data.diat.family.sample.rel), rownames(META.sing))
# [1] TRUE


# add Site information
data.diat.family.sample.rel <- as.data.frame(data.diat.family.sample.rel)
data.diat.family.sample.rel$Site <- META.sing$Site

# melt to four columns (Sample, Site, Family and Proportion)
data.diat.family.sample.rel <- as.data.frame(data.diat.family.sample.rel)
data.diat.family.sample.rel <- rownames_to_column(data.diat.family.sample.rel)

data.diat.family.sample.rel <- melt(data.diat.family.sample.rel,  id.vars = c("rowname", "Site"), 
                                    measure.vars = c("Araphid-pennate", "Polar-centric-Mediophyceae", "Raphid-pennate"))

data.diat.family.sample.rel <- melt(data.diat.family.sample.rel,  id.vars = c("rowname", "Site"), 
                                    measure.vars = c("Abeoformidae_Group_MAIP_1","Acanthamoebidae","Araphid-pennate","Cercozoa_unclassified",
                                                     "Dermamoebidae","Dictyochophyceae_XX","Dinophyceae_unclassified",
                                                     "Gymnodiniaceae","Ichthyosphonida_unclassified","Ichthyosphonida_X","Labyrinthulaceae",
                                                     "Labyrinthulomycetes_X_LAB1/6/8","Lecudinidae","op14-lineage","Phaeophyceae_XX",
                                                     "Polar-centric-Mediophyceae","Prasinococcaceae","Protaspa-lineage","Pseudoperkinsidae",
                                                     "Raphid-pennate","Suessiaceae","Suessiales_unclassified","Symbiodiniaceae",
                                                     "Thraustochytriaceae","Vampyrellida_X","Vermistella-lineage"))
colnames(data.diat.family.sample.rel) <- c("Sample", "Site", "Family", "Proportion")
str(data.diat.family.sample.rel)




# plotting for Shikmona only 

ggplot(data.diat.family.sample.rel[data.diat.family.sample.rel$Site == "Tel Shikmona",], aes(fill=Family, y=Proportion, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Tel Shikmona")+
  labs(y = "Proportion of all families in Foraminifera (%)")+
  # ylim(0,100)+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90)
  )




#####
##### Mantel test - example from R workshop ####

# Comparing the trends in beta diversity, i.e. the dissimilarity matrices from 2 multivariate objects.
# As example, we will use the prokaryotic and eukaryotic community tables.
# For this, we will have to ensure that both tables contain data from the same samples.

otu_pro_sub_euk <- otu_pro_rel[, colnames(otu_euk_rel)]
otu_pro_sub_euk <- otu_pro_sub_euk[rowSums(otu_pro_sub_euk) > 0, ]

# Run Mantel test
mantel(
  vegdist(t(otu_euk_rel)),
  vegdist(t(otu_pro_sub_euk))
)

# Show correlation
plot(
  vegdist(t(otu_euk_rel)),
  vegdist(t(otu_pro_sub_euk))
)


#####
##### cross domain interactions  ####
#####




##### Old plots and ideas






############################################################################################################################################
##### proportion of diatoms in environmental samples #######

# foraminifera
TAX.sing.diat <- TAX.sing.clean[TAX.sing.clean[, "Class"] == "Bacillariophyta", ]
ASV.sing.diat <- ASV.sing.clean[rownames(TAX.sing.diat),]

nrow(TAX.sing.diat)/nrow(TAX.sing.clean)
# 0.8699612
# diatoms represents 86% of the diversity in the foraminifera eukaryotic microbiome 


# sea water
temp <- ASV[,rownames(META[META$Sample_type == "filter",])]
temp <- temp[rowSums(temp) > 0,]
TAX.temp <- TAX[rownames(temp),]
TAX.diat <- TAX.temp[TAX.temp[, "Class"] == "Bacillariophyta", ]
ASV.diat <- temp[rownames(TAX.diat),]

nrow(TAX.diat)/nrow(TAX.temp)
# 0.06807679
# diatoms represents 6.8% of the diversity in the sea water eukaryotic community



# sediment
temp <- ASV[,rownames(META[META$Sample_type == "sediment",])]
temp <- temp[rowSums(temp) > 0,]
TAX.temp <- TAX[rownames(temp),]
TAX.diat <- TAX.temp[TAX.temp[, "Class"] == "Bacillariophyta", ]
ASV.diat <- temp[rownames(TAX.diat),]

nrow(TAX.diat)/nrow(TAX.temp)
# 0.09319555
# diatoms represents 9% of the diversity in the sediment eukaryotic community

##### 
##### another (much longer way to code for the donut and bar plot #####

# ) pooling tax data as list
TAX.pooled <- vector(mode = "list", length = 8)
names(TAX.pooled) <- colnames(TAX.clean)[1:8] # kingdom to species
for (i in 1:8) {
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

# take table with TAX.pooled values within the Division levels 
data <- melt(TAX.pooled$Division) # change level as wished
colnames(data) <- c("Division", "Extraction_Voucher", "Value")

# join metadata
data.meta <- merge(x=data, y=META, by="Extraction_Voucher") 

# create table with sum of value per group
data.meta.st <- dcast(data.meta, Sample_type + Division ~ ., sum, value.var = "Value") # this way is better for plotting
colnames(data.meta.st) <- c("Sample_type", "Division", "Value")

### tree map ? - not good
# 
# # tree map per Division and Sample taype ##
# 
# ggplot(data.meta.st, aes(area = Value, fill = Division, label = Division)) +
#   geom_treemap() +
#   facet_grid(. ~ Sample_type)
# # too many categories, not good to read


### donuts plot is better
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

# calculating how many ASVs from class bacillariophyta within this sample type
TAX.st.diat <- TAX.st[TAX.st[, "Class"] == "Bacillariophyta", ] 
ASV.st.diat <- ASV.st[rownames(TAX.st.diat),]

Ochrophyta_class_bacilariophyta <- sum(ASV.st.diat)
Ochrophyta_class_bacilariophyta_df <- data.frame(col1=st, 
                                                 col2="Ochrophyta (Bacillariophyta)", 
                                                 col3=Ochrophyta_class_bacilariophyta)

# subtracting to obtain how many ochrophytas are not from bacillariophyta
Ochrophyta_not_bacilariophyta <- data.meta.st.don.i[data.meta.st.don.i$Division == "Ochrophyta",]$Value - Ochrophyta_class_bacilariophyta_df$col3
Ochrophyta_not_bacilariophyta_df <- data.frame(col1=st, 
                                                 col2="Ochrophyta (other classes)", 
                                                 col3=Ochrophyta_not_bacilariophyta)
# summarize in a table
Ochrophyta_df <- rbind(Ochrophyta_class_bacilariophyta_df, Ochrophyta_not_bacilariophyta_df)
colnames(Ochrophyta_df) <- c("Sample_type", "Division", "Value")

# remove ochrophyta row from original table (to then add the dat divided in the two categories we just calculated, 
# bacillariophyta and not bacillariophyta)
data.meta.st.don.i <- data.meta.st.don.i[data.meta.st.don.i$Division != "Ochrophyta",]

# combining now the ochrophyta values splitted as we want
data.meta.st.don.i <- rbind(data.meta.st.don.i, Ochrophyta_df)

# sorting in ascending order
str(data.meta.st.don.i)
data.meta.st.don.i$Division <- as.character(data.meta.st.don.i$Division) # has to convert as character before using order()
data.meta.st.don.i <- data.meta.st.don.i[order(data.meta.st.don.i$Division),] 


# calculate percentages of taxa within sample type
data.meta.st.don.i$fraction = data.meta.st.don.i$Value / sum(data.meta.st.don.i$Value) # calculate percentages of taxa within sample type

# summarize taxa with abundance bellow 0.01 as "Other divisions"
data.meta.st.don.i.rare <- data.meta.st.don.i[data.meta.st.don.i$fraction < 0.01,]
sumValue.other.divisions <-  sum(data.meta.st.don.i.rare$Value)
data.other.divisions <-  sum(data.meta.st.don.i.rare$fraction)
data.other.divisions_df <- data.frame(col1=st, 
                                      col2="Other divisions", 
                                      col3=sumValue.other.divisions,
                                      col4=data.other.divisions)
colnames(data.other.divisions_df) <- c("Sample_type", "Division", "Value", "fraction")

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
d1 <- ggplot(data.meta.st.don.i, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Division)) +
  geom_rect() +
  #scale_fill_brewer(palette = 4)+
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15))
  #theme(legend.position = "none") #+
  #ggtitle("Seawater")+
  #theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))


d2 <- ggplot(data.meta.st.don.i, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Division)) +
  geom_rect() +
  #scale_fill_brewer(palette = 4)+
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void()+
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15))
  #theme(legend.position = "none")#+
  #ggtitle("Sediment")+
  #theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))

d3 <- ggplot(data.meta.st.don.i, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Division)) +
  geom_rect() +
  #scale_fill_brewer(palette = 4)+
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.text=element_text(size=13),
        legend.title=element_text(size=15))
  #ggtitle("Foraminifera")+
  #theme(plot.title = element_text(size = 30, face = "bold", hjust = 0.5))+ 
  # guides(shape = guide_legend(override.aes = list(size = 0.5))) 


donut.plot <- d1 + d2 + d3
donut.plot




# OLD and long code -  barplots with foraminifera rare taxa 
# bars per individual samples, clustered by site

# tax pooled
TAX.sing.clean.pooled <- vector(mode = "list", length = 8)
names(TAX.sing.clean.pooled) <- colnames(TAX.sing.clean)[1:8] # kingdom to species
for (i in 1:8) {
  temp <- aggregate(
    ASV.sing.clean,
    by = list(TAX.sing.clean[, i]),
    FUN = sum
  )
  rownames(temp) <- temp$Group.1
  TAX.sing.clean.pooled[[i]] <- as.matrix(temp[, -1])
  rm(temp)
}

# take table with TAX.sing.clean.pooled values within the Division levels 
data.sing.clean <- melt(TAX.sing.clean.pooled$Division) # change level as wished
colnames(data.sing.clean) <- c("Division", "Sample", "Value")

# transform to relative proportions per sample
# dcast per Sample ~ Division
data.sing.clean.sample <- dcast(data.sing.clean, Sample ~ Division)

# Sample as rownames
data.sing.clean.sample <- as.data.frame(data.sing.clean.sample)
data.sing.clean.sample <- column_to_rownames(data.sing.clean.sample, 'Sample')

str(data.sing.clean.sample)

# transform in matrix to calculate proportion table
data.sing.clean.sample <- as.matrix(data.sing.clean.sample)

# calculate proportions
data.sing.clean.sample.rel <- prop.table(data.sing.clean.sample, 1) * 100

# include Site information from META
all.equal(rownames(data.sing.clean.sample.rel), rownames(META.sing))
# [1] TRUE

# add Site information
data.sing.clean.sample.rel <- as.data.frame(data.sing.clean.sample.rel)
data.sing.clean.sample.rel$Site <- META.sing$Site

# melt to four columns (Sample, Site, Division and Proportion)
data.sing.clean.sample.rel <- as.data.frame(data.sing.clean.sample.rel)
data.sing.clean.sample.rel <- rownames_to_column(data.sing.clean.sample.rel)

data.sing.clean.sample.rel <- melt(data.sing.clean.sample.rel,  id.vars = c("rowname", "Site"), # try id.vars = c("family_id", "age_mother"), before: id.vars = "rowname"
                                   measure.vars = c("Apicomplexa", "Cercozoa", "Dinoflagellata", "Lobosa", 
                                                    "Mesomycetozoa", "Ochrophyta", "Prasinodermophyta", "Sagenista"))
colnames(data.sing.clean.sample.rel) <- c("Sample", "Site", "Division", "Proportion")
str(data.sing.clean.sample.rel)


#remove the Ochrophyta division
data.sing.clean.sample.rel.rare <- data.sing.clean.sample.rel[data.sing.clean.sample.rel$Division != "Ochrophyta", ]

# plotting for each site separated


plot1 <- ggplot(data.sing.clean.sample.rel.rare[data.sing.clean.sample.rel.rare$Site == "Eilat",], aes(fill=Division, y=Proportion, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Eilat")+
  labs(y = "Rare taxa in foraminifera samples (%)\n")+
  scale_y_continuous(trans='sqrt', limits = c(0,16)) +
  scale_fill_manual(values= carto_pal(levels(data.sing.clean.sample.rel.rare$Division), "Safe"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
  )

plot1 

plot2 <- ggplot(data.sing.clean.sample.rel.rare[data.sing.clean.sample.rel.rare$Site == "Tel Shikmona",], aes(fill=Division, y=Proportion, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Tel Shikmona")+
  scale_y_continuous(trans='sqrt', limits = c(0,16)) +
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
  )

plot3 <- ggplot(data.sing.clean.sample.rel.rare[data.sing.clean.sample.rel.rare$Site == "Capo Passero",], aes(fill=Division, y=Proportion, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Capo Passero")+
  scale_y_continuous(trans='sqrt', limits = c(0,16)) +
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
  )

plot4 <- ggplot(data.sing.clean.sample.rel.rare[data.sing.clean.sample.rel.rare$Site == "Plemmirio",], aes(fill=Division, y=Proportion, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  ggtitle("Plemmirio") + 
  scale_y_continuous(trans='sqrt', limits = c(0,16)) +
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15)
  )

plot1+plot2+plot3+plot4 + plot_layout(widths = c(1, 2, 2, 2), guides = "collect") 

#####
##### barplot of all sample types - diatoms ASVs only - all sites #####

# filtering only diatoms
TAX.all.diat <- TAX.clean[TAX.clean[, "Class"] == "Bacillariophyta", ]
ASV.all.diat <- ASV.clean[rownames(TAX.all.diat),]
nrow(ASV.all.diat)
# 1502 ASVs in total

ASV.diat.sed <- ASV.all.diat[,rownames(META[META$Sample_type == "sediment",])] 
ASV.diat.sed <- ASV.diat.sed[rowSums(ASV.diat.sed) > 0,] 
a <- nrow(ASV.diat.sed)
# [1] 478

ASV.diat.sw <- ASV.all.diat[,rownames(META[META$Sample_type == "filter",])] 
ASV.diat.sw <- ASV.diat.sw[rowSums(ASV.diat.sw) > 0,] 
b <- nrow(ASV.diat.sw)
# [1] 383

ASV.diat.sing <- ASV.all.diat[,rownames(META[META$Sample_type == "single_cell",])] 
ASV.diat.sing <- ASV.diat.sing[rowSums(ASV.diat.sing) > 0,] 
c <- nrow(ASV.diat.sing)
# [1] 866

##### Exploring how many ASVs are in common between sed and filter samples
ASV.diat.env <- as.data.frame(rownames(ASV.diat.sw)[rownames(ASV.diat.sw) %in% rownames(ASV.diat.sed)])
d <- nrow(ASV.diat.env)
# [1] 201

# saving the ASVs in common between env samples
ASVsinCommon <- c(ASV.diat.env$`rownames(ASV.diat.sw)[rownames(ASV.diat.sw) %in% rownames(ASV.diat.sed)]`)


##### Exploring how many ASVs are in common between single cell and environmental samples
# ASVs in common single cell and env samples
rownames(ASV.diat.sing)[rownames(ASV.diat.sing) %in% ASVsinCommon]
# [1] "sq21"   "sq190"  "sq38"   "sq63"   "sq353"  "sq358"  "sq422"  "sq947"  "sq1527" "sq654"  "sq280" 
e <- nrow(as.data.frame(rownames(ASV.diat.sing)[rownames(ASV.diat.sing) %in% ASVsinCommon]))
# [1] 11


# create table with these counts
tab <- matrix(c(a, b, c, d, e, "Sed", "SW", "Foram", "Sed&SW", "Foram&Sed&SW"), nrow = 5, ncol = 2, byrow=F)
colnames(tab) <- c('DiatomASVs','Sample')
tab <- as.table(tab)
tab <- as.data.frame.matrix(tab)
tab$DiatomASVs <- as.numeric(tab$DiatomASVs)
tab$Sample <- factor(tab$Sample, levels = c("Sed", "SW", "Foram", "Sed&SW", "Foram&Sed&SW"))

# plot barplots
ggplot(tab, aes(x=Sample, y=DiatomASVs, fill=Sample)) + 
  geom_bar(stat = "identity")+
  geom_text(aes(label=DiatomASVs), vjust=-0.25)+
  scale_fill_manual(values = c("forestgreen", "dodgerblue", "red3", "darkslategray4", "red3") ) +
  theme_minimal()+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_blank())+
  ggtitle("All sites")


#####
##### barplot of all sample types - diatoms ASVs only - per site #####

# filtering one site per time
# s <- c("Plemmirio")
# METAs <- META[META$Site == s,]
# ASV.all.diat.s <- ASV.all.diat[,rownames(METAs)]
# ASV.all.diat.s <- ASV.all.diat.s[rowSums(ASV.all.diat.s) > 0,] # excluding ASVs not present in this subset
# 
# # the rest is the same
# ASV.diat.sed <- ASV.all.diat.s[,rownames(METAs[METAs$Sample_type == "sediment",])] 
# ASV.diat.sed <- ASV.diat.sed[rowSums(ASV.diat.sed) > 0,] 
# a <- nrow(ASV.diat.sed)
# 
# ASV.diat.sw <- ASV.all.diat.s[,rownames(METAs[METAs$Sample_type == "filter",])] 
# ASV.diat.sw <- ASV.diat.sw[rowSums(ASV.diat.sw) > 0,] 
# b <- nrow(ASV.diat.sw)
# 
# ASV.diat.sing <- ASV.all.diat.s[,rownames(METAs[METAs$Sample_type == "single_cell",])] 
# ASV.diat.sing <- ASV.diat.sing[rowSums(ASV.diat.sing) > 0,] 
# c <- nrow(ASV.diat.sing)
# 
# # ASVs  in common between sed and filter samples
# ASV.diat.env <- as.data.frame(rownames(ASV.diat.sw)[rownames(ASV.diat.sw) %in% rownames(ASV.diat.sed)])
# d <- nrow(ASV.diat.env)
# 
# # saving the ASVs in common between env samples
# ASVsinCommon <- c(ASV.diat.env$`rownames(ASV.diat.sw)[rownames(ASV.diat.sw) %in% rownames(ASV.diat.sed)]`)
# 
# # ASVs in common single cell and env samples
# rownames(ASV.diat.sing)[rownames(ASV.diat.sing) %in% ASVsinCommon]
# e <- nrow(as.data.frame(rownames(ASV.diat.sing)[rownames(ASV.diat.sing) %in% ASVsinCommon]))
# 
# 
# # make table
# tab <- matrix(c(a, b, c, d, e, "Sed", "SW", "Foram", "Sed&SW", "Foram&Sed&SW"), nrow = 5, ncol = 2, byrow=F)
# colnames(tab) <- c('DiatomASVs','Sample')
# tab <- as.table(tab)
# tab <- as.data.frame.matrix(tab)
# tab$DiatomASVs <- as.numeric(tab$DiatomASVs)
# tab$Sample <- factor(tab$Sample, levels = c("Sed", "SW", "Foram", "Sed&SW", "Foram&Sed&SW"))
# 
# # plot barplots
# p4 <- ggplot(tab, aes(x=Sample, y=DiatomASVs, fill=Sample)) + 
#   geom_bar(stat = "identity")+
#   geom_text(aes(label=DiatomASVs), vjust=-0.25)+
#   scale_fill_manual(values = c("forestgreen", "dodgerblue", "red3", "darkslategray4", "red3") ) +
#   theme_minimal()+
#   theme(legend.position="none",
#         panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.title.x = element_blank())+
#   ggtitle(s)
# 
# p1 + p2 + p3 + p4
#   
#####
##### barplot of rare taxa for each foraminifera individual ####

Ochrophyta_to_remove <- c("Ochrophyta")

# removing these taxa
tax.rare <- TAX.sing.clean[!TAX.sing.clean[ ,"Division"] %in% Ochrophyta_to_remove,]
ASV.rare.rel <- ASV.sing.clean.rel[rownames(tax.rare),]
tax.rare <- tax.rare[, "Division"] # convert to vector with Divisions to use on aggregate function bellow

ASV.rare <- aggregate(
  ASV.rare.rel,
  by = list(tax.rare),
  FUN = sum
)
rownames(ASV.rare) <- ASV.rare$Group.1
ASV.rare <- t(ASV.rare[, -1])

# I need a df with 4 columns: Sample, Site, Division and Proportion

data.barplot <- as.data.frame(ASV.rare)
all.equal(rownames(data.barplot), rownames(META.sing))
# [1] TRUE
data.barplot$Site <- META.sing$Site

# melt to four columns (Sample, Site, Division and Proportion)
data.barplot <- rownames_to_column(data.barplot)

data.barplot <- melt(data.barplot,  id.vars = c("rowname", "Site"), # try id.vars = c("family_id", "age_mother"), before: id.vars = "rowname"
                     measure.vars = c("Apicomplexa", "Cercozoa", "Dinoflagellata", "Lobosa", 
                                      "Mesomycetozoa", "Prasinodermophyta", "Sagenista"))
colnames(data.barplot) <- c("Sample", "Site", "Division", "Proportion")
str(data.barplot)

# taking note of colors in the same order that they appeared in donut plot
carto_pal(name = "Safe")
"#88CCEE" # Apicomplexa # from Safe
"#CC6677" # Cercozoa # from Safe
"#332288" # Dinoflag # from Safe
"#D3D3D3" # light gray Lobosa # new colors for rare taxa
"#808080" # gray Mesomyc # new colors for rare taxa
"#000000" # Black Prasino # new colors for rare taxa
"#6699CC" # Sagenista # from Safe


# plotting one site per time, to combine later
i <- "Eilat"
plot1 <- ggplot(data.barplot[data.barplot$Site == i,], 
                aes(x = Sample, 
                    y = Proportion, 
                    fill = Division)
) +
  geom_col() +
  ggtitle(i) + 
  labs(y = "Rare taxa in foraminifera (%)\n") +
  scale_y_continuous(trans='sqrt', limits = c(0,16))+
  scale_fill_manual(values= c("#88CCEE", "#CC6677", "#332288", "#D3D3D3", "#808080", "#000000", "#6699CC"))+
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 4, angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 11, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.4, 'cm')
  )

i <- "Tel Shikmona"
plot2 <- ggplot(data.barplot[data.barplot$Site == i,], 
                aes(x = Sample, 
                    y = Proportion, 
                    fill = Division)
) +
  geom_col() +
  ggtitle(i) + 
  scale_y_continuous(trans='sqrt', limits = c(0,16))+
  scale_fill_manual(values= c("#88CCEE", "#CC6677", "#332288", "#D3D3D3", "#808080", "#000000", "#6699CC"))+
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 4, angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 11, hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, 'cm')
  )

i <- "Plemmirio"
plot3 <- ggplot(data.barplot[data.barplot$Site == i,], 
                aes(x = Sample, 
                    y = Proportion, 
                    fill = Division)
) +
  geom_col() +
  ggtitle(i) + 
  scale_y_continuous(trans='sqrt', limits = c(0,16))+
  scale_fill_manual(values= c("#88CCEE", "#CC6677", "#332288", "#D3D3D3", "#808080", "#000000", "#6699CC"))+
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 4, angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 11, hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, 'cm')
  )

i <- "Capo Passero"
plot4 <- ggplot(data.barplot[data.barplot$Site == i,], 
                aes(x = Sample, 
                    y = Proportion, 
                    fill = Division)
) +
  geom_col() +
  ggtitle(i) + 
  scale_y_continuous(trans='sqrt', limits = c(0,16))+
  scale_fill_manual(values= c("#88CCEE", "#CC6677", "#332288", "#D3D3D3", "#808080", "#000000", "#6699CC"))+
  theme_minimal()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 4, angle = 90),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 11, hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, 'cm')
  )

# Combining with patchwork

# svglite::svglite(
#   filename="Fig_barplot_rare_euks.svg",
#   width = 970,
#   height = 315,
#   ) # for mac, problem saving as .svg with mouse

plot1+plot2+plot3+plot4 + 
  plot_layout(
    widths = c(1, 2, 2, 2), 
    guides = "collect"
  ) 
# dev.off()

#####
##### inv simpson (boxplot) all sample together and plotting for each sample type per time ####

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

# try plot one sample type per time, to compare euks vs proks

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
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
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
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
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
  theme_minimal()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 15, hjust = 0.95),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15)
  )

plot_invS_sed 
#####
##### Intersection plot - Upset plot for foraminifera diatoms - to show sites in the plot ####

ASV.diat.eilat <- ASV.diat[,rownames(META.sing[META.sing$Site == "Eilat",])] 
ASV.diat.eilat <- ASV.diat.eilat[rowSums(ASV.diat.eilat) > 0,] 

ASV.diat.shik <- ASV.diat[,rownames(META.sing[META.sing$Site == "Tel Shikmona",])] 
ASV.diat.shik <- ASV.diat.shik[rowSums(ASV.diat.shik) > 0,] 

ASV.diat.plem <- ASV.diat[,rownames(META.sing[META.sing$Site == "Plemmirio",])] 
ASV.diat.plem <- ASV.diat.plem[rowSums(ASV.diat.plem) > 0,] 

ASV.diat.capo <- ASV.diat[,rownames(META.sing[META.sing$Site == "Capo Passero",])] 
ASV.diat.capo <- ASV.diat.capo[rowSums(ASV.diat.capo) > 0,] 

# create list with ASVs name per sample type
upsetlist.sing <- list(rownames(ASV.diat.eilat), rownames(ASV.diat.shik), rownames(ASV.diat.plem), rownames(ASV.diat.capo))
names(upsetlist.sing) <- c("Eilat", "Tel Shikmona", "Plemmirio", "Capo Passero")

# plot upset
upset.plot.sing <- upset(fromList(upsetlist.sing), order.by = "freq",  mainbar.y.label = "Site Intersections", 
                         main.bar.color = c("salmon", "black","darkorchid1", "black", "blue","black","gold1",
                                            "black","black","black","black","black","black","black","black"),
                         # matrix.color = c("red3", "forestgreen", "darkgrey", "dodgerblue", "orange", "grey"),
                         sets.bar.color = c("salmon", "darkorchid1", "blue", "gold1"),
                         sets.x.label = "Diatoms' ASVs", text.scale = c(2,1.7,2,1.5,2,2)
)
upset.plot.sing


#####
##### NMDS all sample types  - diatoms ASVs only  ######
ASV.all.diat.rel <- prop.table(ASV.all.diat, 2) * 100

# calculate NMDS
NMDS <- metaMDS(t(ASV.all.diat.rel), k = 2)
NMDS$stress
#  0.1335469 # all diatoms

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
  pch = 8,
  cex = 1
)

# add legend of stress value
legend("bottomleft", "stress = 0.13", bty = "n", cex = 1) 

# add legend for sample type
legend(
  "bottom",
  legend = levels(META.plot$Sample_type),
  col = color.sample.type,
  pch = 15,
  pt.cex = 1,
  title = "Sample type"
)

##### 
##### constraining foraminifera microbiome to abundant diatoms only #######

# unique(TAX.sing.diat[,"Family"])
# 
# # "Araphid-pennate"                "Raphid-pennate"                 "Polar-centric-Mediophyceae"     "Bacillariophyta_X_unclassified"
# 
# sum(TAX.sing.diat[, "Family"] == "Bacillariophyta_X_unclassified")
# # [1] 17 ASVs 
# 
# nrow(TAX.sing.diat[TAX.sing.diat[, "Family"] == "Bacillariophyta_X_unclassified",])/nrow(TAX.sing.diat)
# # [1] 0.01963048 --> less then 0.02%
# 
# 
# 
# ASV.sing.diat
# #### pooling tax data as list
# TAX.pooled.diat <- vector(mode = "list", length = 8)
# names(TAX.pooled.diat) <- colnames(TAX.sing.diat)[1:8] # kingdom to species
# for (i in 1:8) {
#   temp <- aggregate(
#     ASV.diat,
#     by = list(TAX.sing.diat[, i]),
#     FUN = sum
#   )
#   rownames(temp) <- temp$Group.1
#   TAX.pooled.diat[[i]] <- as.matrix(temp[, -1])
#   rm(temp)
# }
# 
# 
# # show table with summary of ASVs or proportions per family of Class Bacillariophyta
# 
# 
# # take table with TAX.pooled.abundant values within the Family levels 
# data.diat <- melt(TAX.pooled.diat$Family) # change level as wished
# colnames(data.diat) <- c("Family", "Extraction_Voucher", "Value")
# 
# # join metadata
# data.diat.meta <- merge(x=data.diat, y=META.sing, by="Extraction_Voucher") 
# 
# # create table with sum of value per Site
# data.diat.site <- dcast(data.diat.meta, Site + Family ~ ., sum, value.var = "Value") # this way is better for plotting
# colnames(data.diat.site) <- c("Site", "Family", "Value")
# 
# 
# # transform in relative proportions
# # dcast per Sample ~ Division
# data.diat.site.dcast <- dcast(data.diat.site, Site ~ Family)
# 
# # Sample as rownames
# data.diat.site.dcast <- as.data.frame(data.diat.site.dcast)
# data.diat.site.dcast <- column_to_rownames(data.diat.site.dcast, 'Site')
# 
# # transform in matrix to calculate proportion table
# data.diat.site.dcast <- as.matrix(data.diat.site.dcast)
# 
# # calculate proportions
# data.diat.site.rel <- prop.table(data.diat.site.dcast, 1) * 100
# 
# # back to data frame
# data.diat.site.rel <- as.data.frame(data.diat.site.rel)
# data.diat.site.rel <- rownames_to_column(data.diat.site.rel)
# 
# # export as excel file
# require("writexl")
# write_xlsx(data.diat.site.rel,"diatoms_families_prop_foraminifera_microbiome_cleaned.xlsx")
# # Araphid pennate is the absolutely dominant family





all.equal(rownames(ASV.diat), rownames(ASV.sing.diat))

#####
##### Calculating frequency of foraminifera diatom ASVs per site and saving RDS files ####

# all.equal(rownames(META.sing), colnames(ASV.diat.sing))
# # [1] TRUE
# 
# min(colSums(ASV.diat.sing)) # [1] 27220
# max(colSums(ASV.diat.sing)) # [1] 140554
# 
# # rarefy
# 
# ASV.diat.raref <- rrarefy(t(ASV.diat.sing), min(rowSums(t(ASV.diat.sing))))
# 
# min(colSums(t(ASV.diat.raref))) # [1] 27220
# max(colSums(t(ASV.diat.raref))) # [1] 27220
# 
# min(colSums(ASV.diat.raref)) # if minimum is not zero, no ASV was lost- minimum is zero, so we lost ASVs
# 
# dim(ASV.diat.raref)
# min(ASV.diat.raref[ASV.diat.raref > 0]) # to see the minimum seq depth that is not zero
# 
# diat.freq.site <- aggregate(ASV.diat.raref,
#                  by = list(META.sing$Site),
#                  FUN = sum)

# saving in the directory
#saveRDS(diat.freq.site, "ASV.diat.raref.agg.site.RDS")


#####
##### Bray curtis vs geographical distances - pairwise comparison ####

# 
# require(spaa)
# bc.dist <- vegdist(t(ASV.diat.rel))
# str(bc.dist)
# bc <- dist2list(bc.dist)
# 
# nrow(bc)
# 
# require(geosphere) # to calculate distance between two points
# 
# geo_dist <- read.table("Sites_and_coord.txt", 
#                        h = T,
#                        sep = "\t",
#                        row.names = 1)
# # same order than META
# geo_dist <- geo_dist[rownames(META.sing),]
# all.equal(rownames(geo_dist), rownames(META.sing))
# # [1] TRUE
# 
# sample.comb <- combn(rownames(geo_dist), 2)
# combn(1:6, 2) # just to understand how t works. it already remove the duplicate combinations and doesnt make comb among the same
# 
# # pairwise comparisons
# geo_dist_out <- data.frame(
#   t(sample.comb),
#   geo_dist = map_dbl(
#     1:ncol(sample.comb),
#     function(x) {
#       distm(c(t(geo_dist[sample.comb[1, x], c("Long", "Lat")])), c(t(geo_dist[sample.comb[2, x], c("Long", "Lat")])), fun = distHaversine)
#     }
#   ),
#   bc_dis = map_dbl(
#     1:ncol(sample.comb),
#     function(x) {
#       bc[bc$col == sample.comb[1, x] & bc$row == sample.comb[2, x], 3]
#     }
#   ),
#   site_comparison = map_chr(
#     1:ncol(sample.comb),
#     function(x) {
#       paste(geo_dist[sample.comb[1, x], "Site"], geo_dist[sample.comb[2, x], "Site"])
#     }
#   )
# )
# 
# # convert to km
# geo_dist_out$geo_dist_km <- geo_dist_out$geo_dist/1000
# 
# 
# levels(as.factor(geo_dist_out$site_comparison))
# # [1] "Capo Passero Capo Passero" "Capo Passero Eilat"        "Capo Passero Plemmirio"    "Capo Passero Tel Shikmona" "Eilat Eilat"              
# # [6] "Plemmirio Eilat"           "Plemmirio Plemmirio"       "Plemmirio Tel Shikmona"    "Tel Shikmona Eilat"        "Tel Shikmona Tel Shikmona"
# 
# 
# 
# # plotting controlling settings
# geo_dist_plot <- geo_dist_out
# 
# geo_dist_plot$site_comparison <- as.factor(geo_dist_plot$site_comparison)
# 
# 
# geo_dist_plot$site_comparison <- factor(geo_dist_plot$site_comparison, 
#                                         levels = c("Eilat Eilat", "Tel Shikmona Tel Shikmona", "Plemmirio Plemmirio", "Capo Passero Capo Passero",  
#                                                    "Tel Shikmona Eilat", "Plemmirio Eilat", "Capo Passero Eilat","Plemmirio Tel Shikmona",
#                                                    "Capo Passero Tel Shikmona", "Capo Passero Plemmirio"),
#                                         labels = c("Eilat vs Eilat", "Tel Shikmona vs Tel Shikmona", "Plemmirio vs Plemmirio", "Capo Passero vs Capo Passero",  
#                                                    "Eilat vs Tel Shikmona", "Eilat vs Plemmirio", "Eilat vs Capo Passero",  "Tel Shikmona vs Plemmirio",
#                                                    "Tel Shikmona vs Capo Passero", "Capo Passero vs Plemmirio")
# )
# str(geo_dist_plot)
# 
# levels(geo_dist_plot$site_comparison)
# # [1] "Eilat vs Eilat"               "Tel Shikmona vs Tel Shikmona" "Plemmirio vs Plemmirio"       "Capo Passero vs Capo Passero"
# # [5] "Eilat vs Tel Shikmona"        "Eilat vs Plemmirio"           "Eilat vs Capo Passero"        "Tel Shikmona vs Plemmirio"   
# # [9] "Tel Shikmona vs Capo Passero" "Capo Passero vs Plemmirio" 
# 
# # setting color scheme and shape for type of comparison
# color.type <- c("gold1", "salmon", "darkorchid1", "blue", "yellow3","turquoise", "turquoise4", "lightsteelblue3", "lightsteelblue4", "lightsteelblue1")
# geo_dist_plot$color_comparison_type <- geo_dist_plot$site_comparison
# levels(geo_dist_plot$color_comparison_type) <- color.type
# geo_dist_plot$color_comparison_type <- as.character(geo_dist_plot$color_comparison_type)
# 
# dev.off()
# 
# plot(
#   jitter(geo_dist_plot$geo_dist_km, 20), 
#   geo_dist_plot$bc_dis,
#   pch = 16, 
#   col = adjustcolor(geo_dist_plot$color_comparison_type, alpha.f = 0.1), 
#   cex = 0.5,
#   xlab="Geographical distance (km)", ylab="Bray Curtis dissimilarity"
# )
# 
# 
# # add legend for comparison type
# legend(
#   "bottom",
#   legend = levels(geo_dist_plot$site_comparison),
#   col = color.type,
#   pch = 16,
#   bty = "n", 
#   pt.cex = 2,
#   cex = 1.2
#   # title = "Sample type"
# )

#####
##### NMDS within foraminifera abundant taxa ######

# ASV.abundant.diat # ASV table just for abundant ASVs
# ASV.abundant.diat.rel <- prop.table(ASV.abundant.diat, 2) * 100
# 
# ASV.abundant.araphid # ASV table just for abundant ASVs from Araphid pennate family only
# ASV.abundant.araphid.rel <- prop.table(ASV.abundant.araphid, 2) * 100
# 
# 
# # calculate NMDS for abundant taxa
# NMDS <- metaMDS(t(ASV.abundant.diat.rel), k = 2) # k number of dimensions
# NMDS$stress
# #  0.232288 # single cell abundant diatoms
# # stress is high
# 
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(1,2))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(1,2))
# 
# # add points for sites
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 2],
#   col = META.plot$color_site,
#   pch = 16,
#   cex = 1
# )
# 
# # add legend for sites
# legend(
#   "right",
#   legend = levels(META.plot$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
#   pch = 15,
#   pt.cex = 1,
#   title = "Site"
# )
# 
# # calculate NMDS for abundant araphid pennates
# NMDS <- metaMDS(t(ASV.abundant.araphid.rel), k = 2) # k number of dimensions
# NMDS$stress
# #  0.2271361# single cell abundant diatoms
# # stress is high
# 
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(1,2))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(1,2))
# 
# # add points for sites
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 2],
#   col = META.plot$color_site,
#   pch = 16,
#   cex = 1
# )
# 
# # add legend for sites
# legend(
#   "right",
#   legend = levels(META.plot$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
#   pch = 15,
#   pt.cex = 1,
#   title = "Site"
# )
##### 
##### testing NMDS with 3 dimensions to reduce stress #######
# ASV.diat.rel <- prop.table(ASV.diat, 2) * 100
# # calculate NMDS
# NMDS <- metaMDS(t(ASV.diat.rel), k = 3) # k number of dimensions
# NMDS$stress
# #  0.1506323 # single cell abundant diatoms
# # stress is lower than with 2 dimensions, but still not bellow 0.1
# 
# #  0.1732888 # single cell diatoms all
# 
# # set meta table
# META.sing$Depth_cat <- factor(META.sing$Depth_cat, levels=c("Shallow", "Deep"))
# META.sing$Substrate <- as.factor(META.sing$Substrate)
# 
# META.plot <- META.sing
# str(META.plot)
# 
# 
# ############################## plotting sites as colors
# 
# # 1 and 2
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(1,2))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(1,2))
# 
# # add points for sites
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 2],
#   col = META.plot$color_site,
#   pch = 16,
#   cex = 1
# )
# 
# # add legend for sites
# legend(
#   "right",
#   legend = levels(META.plot$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
#   pch = 15,
#   pt.cex = 1,
#   title = "Site"
# )
# 
# # 1 and 3
# 
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(1,3))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(1,3))
# 
# # add points for sites
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 3],
#   col = META.plot$color_site,
#   pch = 16,
#   cex = 1
# )
# 
# # add legend for sites
# legend(
#   "right",
#   legend = levels(META.plot$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
#   pch = 15,
#   pt.cex = 1,
#   title = "Site"
# )
# 
# 
# # 2 and 3
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(2,3))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(2,3))
# 
# # add points for sites
# points(
#   NMDS$points[, 2],
#   NMDS$points[, 3],
#   col = META.plot$color_site,
#   pch = 16,
#   cex = 1
# )
# 
# # add legend for sites
# legend(
#   "right",
#   legend = levels(META.plot$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
#   pch = 15,
#   pt.cex = 1,
#   title = "Site"
# )
# 
# 
# 
# 
# 
# ############################## plotting site as colors, substrate as shape and depth as size
# 
# # 1 and 2
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(1,2))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(1,2))
# 
# # add points for sites
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 2],
#   bg = META.plot$color_site,
#   pch = as.numeric(META.plot$Substrate)+20,
#   cex = as.numeric(META.plot$Depth_cat)
# )
# 
# # add legend for sites
# legend(
#   "right",
#   legend = levels(META.plot$Site),
#   col = c("darkolivegreen4", "darkorchid1", "gold1", "firebrick"),
#   pch = 15,
#   pt.cex = 1,
#   title = "Site"
# )
# # add legend for substrates
# legend(
#   "topright",
#   legend = levels(META.plot$Substrate),
#   col = "black",
#   pch = c(21, 22),
#   pt.cex = 1,
#   title = "Substrate"
# )
# 
# # add legend for depth
# legend(
#   "bottomright",
#   legend = levels(META.plot$Depth_cat),
#   col = "black",
#   pch = 16,
#   pt.cex = c(1,2),
#   title = "Depth"
# )
# 
# # 1 and 3
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(1,3))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(1,3))
# 
# # add points for sites
# points(
#   NMDS$points[, 1],
#   NMDS$points[, 3],
#   bg = META.plot$color_site,
#   pch = as.numeric(META.plot$Substrate)+20,
#   cex = as.numeric(META.plot$Depth_cat)
# )
# 
# # 2 and 3
# dev.off() #to reset the PAR settings I configured earlier
# 
# plot(NMDS, display = "sites", type = "n", choices = c(2,3))
# 
# # add hulls based on sites
# ordihull(NMDS, META.plot$Site, choices = c(2,3))
# 
# # add points for sites
# points(
#   NMDS$points[, 2],
#   NMDS$points[, 3],
#   bg = META.plot$color_site,
#   pch = as.numeric(META.plot$Substrate)+20,
#   cex = as.numeric(META.plot$Depth_cat)
# )

##### 
##### PERMANOVA for abundant taxa in foraminifera ######
# # change for the sub set dataset
# META.adonis <- META.sing
# str(META.adonis)
# META.adonis$Site <- droplevels(META.adonis$Site)
# META.adonis$Condition <- droplevels(META.adonis$Condition)
# 
# # with all factors
# adonis2(t(ASV.abundant.diat.rel) ~ Site + Substrate + Depth_cat, data = META.adonis, sqrt.dist = T)
# # Permutation test for adonis under reduced model
# # Terms added sequentially (first to last)
# # Permutation: free
# # Number of permutations: 999
# # 
# # adonis2(formula = t(ASV.abundant.diat.rel) ~ Site + Substrate + Depth_cat, data = META.adonis, sqrt.dist = T)
# #            Df SumOfSqs      R2       F Pr(>F)    
# # Site        3   17.184 0.36188 20.1272  0.001 ***
# # Substrate   1    0.177 0.00372  0.6208  0.807    
# # Depth_cat   1    0.242 0.00510  0.8512  0.595    
# # Residual  105   29.882 0.62929                   
# # Total     110   47.485 1.00000                   
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# # depth and substrate not significant
# 
# # PERMANOVA per site
# adonis2(t(ASV.abundant.diat.rel) ~ Site, data = META.adonis, sqrt.dist = T)
# # Permutation test for adonis under reduced model
# # Terms added sequentially (first to last)
# # Permutation: free
# # Number of permutations: 999
# # 
# # adonis2(formula = t(ASV.abundant.diat.rel) ~ Site, data = META.adonis, sqrt.dist = T)
# #           Df SumOfSqs      R2      F Pr(>F)    
# # Site       3   17.184 0.36188 20.227  0.001 ***
# # Residual 107   30.301 0.63812                  
# # Total    110   47.485 1.00000                  
# # ---
# #   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# 
# # Pairwise comparison of site
# # Let's also have a look at pairwise comparisons for sites
# site.comb <- combn(levels(META.adonis$Site), 2) # creating a triangular pairwise comparison
# permanova.list <-apply(
#   site.comb,
#   2,
#   function(x) {
#     adonis2(
#       t(ASV.diat.rel[, rownames(META.adonis)[META.adonis$Site %in% x]]) ~ Site, 
#       data = droplevels(META.adonis[META.adonis$Site %in% x, ]), 
#       sqrt.dist = T
#     )
#   }
# )
# 
# # reformat as table
# permanova.df <- data.frame(
#   t(site.comb),
#   do.call(
#     "rbind",
#     lapply(
#       permanova.list,
#       function(x) {
#         c(x$R2[1], x$'F'[1], x$'Pr(>F)'[1])
#       }
#     )
#   )
# )
# colnames(permanova.df) <- c("cond1", "cond2", "R2", "F", "p")
# # adjust p-values for multiple testing
# permanova.df$p.adj <- p.adjust(permanova.df$p, method = "fdr")
# 
# permanova.df
# # cond1        cond2        R2        F     p p.adj
# # 1 Capo Passero    Plemmirio 0.1885339 14.17258 0.001 0.001
# # 2 Capo Passero Tel Shikmona 0.2656133 22.06251 0.001 0.001
# # 3 Capo Passero        Eilat 0.3302634 22.19060 0.001 0.001
# # 4    Plemmirio Tel Shikmona 0.2525999 20.95423 0.001 0.001
# # 5    Plemmirio        Eilat 0.3101513 20.68129 0.001 0.001
# # 6 Tel Shikmona        Eilat 0.3406128 23.76175 0.001 0.001
# 
# # all sites are different!
# 




##### 
##### Plotting abundant taxa #####

# TAX.pooled.abundant <- vector(mode = "list", length = 8)
# names(TAX.pooled.abundant) <- colnames(TAX.abundant.taxa)[1:8] # kingdom to species
# for (i in 1:8) {
#   temp <- aggregate(
#     ASV.abundant.taxa,
#     by = list(TAX.abundant.taxa[, i]),
#     FUN = sum
#   )
#   rownames(temp) <- temp$Group.1
#   TAX.pooled.abundant[[i]] <- as.matrix(temp[, -1])
#   rm(temp)
# }
# 
# # take table with TAX.pooled.abundant values within the Division levels 
# data.abund <- melt(TAX.pooled.abundant$Division) # change level as wished
# colnames(data.abund) <- c("Division", "Sample", "Value")
# 
# # # join metadata
# # data.abund.meta <- merge(x=data.abund, y=META, by="Extraction_Voucher") 
# # 
# # # dcast Sample 
# # data.abund.meta.sample <- dcast(data.abund.meta, Extraction_Voucher + Division ~ ., sum, value.var = "Value") # this way is better for plotting
# # colnames(data.abund.meta.sample) <- c("Sample", "Division", "Value")
# 
# # transform to relative proportions 
# 
# # dcast per Sample ~ Division
# data.abund.sample <- dcast(data.abund, Sample ~ Division)
# 
# # Sample as rownames
# data.abund.sample <- as.data.frame(data.abund.sample)
# data.abund.sample <- column_to_rownames(data.abund.sample, 'Sample')
# 
# str(data.abund.sample)
# 
# # transform in matrix to calculate propportion table
# data.abund.sample <- as.matrix(data.abund.sample)
# 
# # calculate proportions
# data.abund.sample.rel <- prop.table(data.abund.sample, 1) * 100
# 
# # include Site information from META
# all.equal(rownames(data.abund.sample.rel), rownames(META.sing))
# 
# data.abund.sample.rel <- as.data.frame(data.abund.sample.rel)
# data.abund.sample.rel$Site <- META.sing$Site
# 
# # melt to four columns (Sample, Site, Division and Proportion)
# data.abund.sample.rel <- as.data.frame(data.abund.sample.rel)
# data.abund.sample.rel <- rownames_to_column(data.abund.sample.rel)
# 
# data.abund.sample.rel <- melt(data.abund.sample.rel,  id.vars = c("rowname", "Site"), 
#                               measure.vars = c("Apicomplexa", "Cercozoa", "Dinoflagellata", 
#                                                "Ochrophyta", "Sagenista"))
# colnames(data.abund.sample.rel) <- c("Sample", "Site", "Division", "Proportion")
# str(data.abund.sample.rel)
# 
# 
# #remove the Ochrophyta division
# data.abund.sample.rel.rare <- data.abund.sample.rel[data.abund.sample.rel$Division != "Ochrophyta", ]
# 
# # plotting for each site separated
# 
# p1 <- ggplot(data.abund.sample.rel.rare[data.abund.sample.rel.rare$Site == "Capo Passero",], aes(fill=Division, y=Proportion, x=Sample)) + 
#   geom_bar(position="stack", stat="identity")+
#   ggtitle("Capo Passero")+
#   labs(y = "Proportion of the rare taxa (%)")+
#   ylim(0,8)+
#   theme_minimal()+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.text.x = element_text(angle = 90)
#   )
# 
# p2 <- ggplot(data.abund.sample.rel.rare[data.abund.sample.rel.rare$Site == "Plemmirio",], aes(fill=Division, y=Proportion, x=Sample)) + 
#   geom_bar(position="stack", stat="identity")+
#   ggtitle("Plemmirio")+
#   ylim(0,8)+
#   theme_minimal()+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 90)
#   )
# 
# p3 <- ggplot(data.abund.sample.rel.rare[data.abund.sample.rel.rare$Site == "Tel Shikmona",], aes(fill=Division, y=Proportion, x=Sample)) + 
#   geom_bar(position="stack", stat="identity")+
#   ggtitle("Tel Shikmona")+
#   ylim(0,8)+
#   theme_minimal()+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 90)
#   )
# 
# p4 <- ggplot(data.abund.sample.rel.rare[data.abund.sample.rel.rare$Site == "Eilat",], aes(fill=Division, y=Proportion, x=Sample)) + 
#   geom_bar(position="stack", stat="identity")+
#   ggtitle("Eilat") + 
#   ylim(0,8)+
#   theme_minimal()+
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x = element_text(angle = 90)
#   )
# 
# p1+p2+p3+p4 + plot_layout(widths = c(2, 2, 2, 1), guides = "collect") 


# does it makes sense to analyse these taxa? likely is better just filter the class bacillariophyta

##### 
##### Venn-Diagram ####

# Plotting each site the ASVs of single cell vs sediment vs filters
# # to observe if there is an intake from the environmental biota
# 
# library("ggVennDiagram")
# require(ggplot2)
# 
# 
# # preparing data
# venn.list <- vector("list", length = 3)
# names(venn.list) <- c("foraminifera", "sediment", "seawater")
# tmp.asv <- ASV.clean[, META$Sample_type == "single_cell"]
# tmp.rel <- prop.table(tmp.asv, 2) * 100
# venn.list$foraminifera <- rownames(tmp.rel)[apply(tmp.rel, 1, function(x) sum(x > 0) >= 1)]# remove ASVs that are not present in any sample in this sub setting
# tmp.asv <- ASV.clean[, META$Sample_type == "sediment"]
# tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
# venn.list$sediment <- rownames(tmp.asv)
# tmp.asv <- ASV.clean[, META$Sample_type == "filter"]
# tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
# venn.list$seawater <- rownames(tmp.asv)
# 
# str(venn.list)
# dim(tmp.rel)
# 
# # Plot Venn Diagram
# ggVennDiagram(venn.list, label_alpha = 0)+
#   ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")
# 
# # 33 ASVs in common in all sample types. 
# 
# # but how many from the diatoms?
# 
# # preparing data for diatoms - looking class bacillaryophita between environmental and foraminiferal samples
# 
# TAX.sed.diat <- TAX.sed[TAX.sed[, "Class"] == "Bacillariophyta", ]
# ASV.sed.diat <- ASV.sed[rownames(TAX.sed.diat),]
# 
# TAX.filters.diat <- TAX.filters[TAX.filters[, "Class"] == "Bacillariophyta", ]
# ASV.filters.diat <- ASV.filters[rownames(TAX.filters.diat),]
# 
# 
# venn.list <- vector("list", length = 3)
# names(venn.list) <- c("foraminifera", "sediment", "seawater")
# tmp.asv <- ASV.diat
# tmp.rel <- prop.table(tmp.asv, 2) * 100
# venn.list$foraminifera <- rownames(tmp.rel)[apply(tmp.rel, 1, function(x) sum(x > 0) >= 1)]# remove ASVs that are not present in any sample in this sub setting
# tmp.asv <- ASV.sed.diat
# tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
# venn.list$sediment <- rownames(tmp.asv)
# tmp.asv <- ASV.filters.diat
# tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
# venn.list$seawater <- rownames(tmp.asv)
# 
# # Venn Diagram
# str(venn.list)
# dim(tmp.rel)
# 
# ggVennDiagram(venn.list, label_alpha = 0)+
#   ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")
# 
# 
# ### to see per specific site
# s <- 'Plemmirio' # # select site to subset in s
# METAs <- META[META$Site == s,] # select site s in META
# 
# # or, another option to write the list, setting filters before taking the ASVs names
# venn.list <- vector("list", length = 3)
# names(venn.list) <- c("foraminifera", "sediment", "seawater")
# tmp.asv <- ASV.diat[, rownames(METAs[METAs$Sample_type == "single_cell",])]
# tmp.rel <- prop.table(tmp.asv, 2) * 100
# venn.list$foraminifera <- rownames(tmp.rel)[apply(tmp.rel, 1, function(x) sum(x > 0) >= 1)]
# tmp.asv <- ASV.sed.diat[, rownames(METAs[METAs$Sample_type == "sediment",])]
# tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
# venn.list$sediment <- rownames(tmp.asv)
# tmp.asv <- ASV.filters.diat[, rownames(METAs[METAs$Sample_type == "filter",])]
# tmp.asv <- tmp.asv[rowSums(tmp.asv) > 0, ]
# venn.list$seawater <- rownames(tmp.asv)
# 
# # plot 
# ggVennDiagram(venn.list, label_alpha = 0)+
#   ggplot2::scale_fill_gradient(low="lightblue",high = "yellow")
# 

#####
##### plotting this data as tree map - options that were not good ######
# 
# ##### tree map for foraminifeta microbiome all taxa (cleaned)
# ASV.sing.clean <- ASV.sing[rownames(TAX.sing.clean),]
# 
# TAX.sing.clean.pooled <- vector(mode = "list", length = 8)
# names(TAX.sing.clean.pooled) <- colnames(TAX.sing.clean)[1:8] # kingdom to species
# for (i in 1:8) {
#   temp <- aggregate(
#     ASV.sing.clean,
#     by = list(TAX.sing.clean[, i]),
#     FUN = sum
#   )
#   rownames(temp) <- temp$Group.1
#   TAX.sing.clean.pooled[[i]] <- as.matrix(temp[, -1])
#   rm(temp)
# }
# 
# # take table with TAX.sing.clean.pooled values within the Division levels 
# data.sing.clean <- melt(TAX.sing.clean.pooled$Division) # change level as wished
# colnames(data.sing.clean) <- c("Division", "Extraction_Voucher", "Value")
# 
# # join metadata
# require(dplyr)
# data.sing.clean.meta <- merge(x=data.sing.clean, y=META, by="Extraction_Voucher") 
# 
# # create table with sum of value per group
# data.sing.clean.site <- dcast(data.sing.clean.meta, Site + Division ~ ., sum, value.var = "Value") # this way is better for plotting
# colnames(data.sing.clean.site) <- c("Site", "Division", "Value")
# 
# ggplot(data.sing.clean.site, aes(area = Value, fill = Division, label = Division)) +
#   geom_treemap() +
#   facet_grid(. ~ Site)
# 
# # Basically the foraminifera microbiome is only composed by the Division Ochrophyta, try another visualization
# 
# 
# # donuts?
# 
# ## Compute percentages
# 
# data.sing.clean.site.don <- data.sing.clean.site
# 
# data.sing.clean.site.don.CP <- data.sing.clean.site.don[data.sing.clean.site.don$Site == "Eilat",]
# 
# data.sing.clean.site.don.CP$fraction = data.sing.clean.site.don.CP$Value / sum(data.sing.clean.site.don.CP$Value)
# 
# 
# # Compute the cumulative percentages (top of each rectangle)
# data.sing.clean.site.don.CP$ymax = cumsum(data.sing.clean.site.don.CP$fraction)
# 
# # Compute the bottom of each rectangle
# data.sing.clean.site.don.CP$ymin = c(0, head(data.sing.clean.site.don.CP$ymax, n=-1))
# 
# # Make the plot
# ggplot(data.sing.clean.site.don.CP, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Division)) +
#   geom_rect() +
#   coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
#   xlim(c(2, 4)) # Try to remove that to see how to make a pie chart
# 
# ## not good yet
# 
# # pie chart?
# 
# # Basic piechart
# ggplot(data.sing.clean.site.don, aes(x="", y=Value, fill=Division)) +
#   geom_bar(stat="identity", width=1) +
#   coord_polar("y", start=0)+
#   facet_grid(. ~ Site)
# 
# ## also not good
# 
# # bar plot ?
# 
# # stacked bar plot
# ggplot(data.sing.clean.site, aes(fill=Division, y=Value, x=Site)) + 
#   geom_bar(position="fill", stat="identity")
# 
# # individuals bar plots
# ggplot(data.sing.clean.site, aes(fill=Division, y=Value, x=Site)) + 
#   geom_bar(position="dodge", stat="identity")
# 
# # transform to relative proportions before
# 
# data.sing.clean.meta.temp <- data.sing.clean.site
# 
# # dcast per Site ~ Division
# data.sing.clean.meta.temp <- dcast(data.sing.clean.meta.temp, Site ~ Division)
# 
# # site as rownames
# data.sing.clean.meta.temp <- as.data.frame(data.sing.clean.meta.temp)
# data.sing.clean.meta.temp <- column_to_rownames(data.sing.clean.meta.temp, 'Site')
# 
# str(data.sing.clean.meta.temp)
# 
# # transform in matrix to calculate propportion table
# data.sing.clean.meta.temp <- as.matrix(data.sing.clean.meta.temp)
# 
# # calculate proportions
# data.sing.clean.rel <-prop.table(data.sing.clean.meta.temp, 1) * 100
# 
# # melt to two columns
# data.sing.clean.rel <- as.data.frame(data.sing.clean.rel)
# data.sing.clean.rel <- rownames_to_column(data.sing.clean.rel)
# 
# data.sing.clean.rel <- melt(data.sing.clean.rel,  id.vars = "rowname", 
#                             measure.vars = c("Apicomplexa", "Cercozoa", "Dinoflagellata", "Lobosa", 
#                                              "Mesomycetozoa", "Ochrophyta", "Prasinodermophyta", "Sagenista"))
# colnames(data.sing.clean.rel) <- c("Site", "Division", "Value")
# str(data.sing.clean.rel)
# 
# data.sing.clean.rel$Site <- factor(data.sing.clean.rel$Site, levels = c("Capo Passero", "Plemmirio", "Tel Shikmona", "Eilat"))
# 
# 
# # individuals bar plots
# p1 <- ggplot(data.sing.clean.rel, aes(fill=Division, y=Value, x=Site)) + 
#   geom_bar(position="dodge", stat="identity")+
#   scale_y_continuous(name="Relative proportion of taxa (%)", limits=c(0, 100))+
#   theme_minimal()+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         panel.grid.minor = element_blank())
# p1 
# 
# p2 <- ggplot(data.sing.clean.rel, aes(fill=Division, y=Value, x=Site)) + 
#   geom_bar(position="dodge", stat="identity") +
#   ylim(0, 0.5)+
#   scale_y_continuous(name="< 0.5%", limits=c(0, 0.5))+
#   theme_minimal()+
#   theme(panel.grid.minor = element_blank())
# 
# p2
# 
# require(patchwork)
# p1/p2 + plot_layout(widths = c(1, 1), heights = c(3, 1), guides = "collect") 
# 
# # not good, first plot is meaningless... 
#####
##### try only the less abundant data and plot stacked bar plots ######
# 
# #remove the Ochrophyta division
# data.sing.clean.rel.rare <- data.sing.clean.rel[data.sing.clean.rel$Division != "Ochrophyta", ]
# 
# # stacked bar plot
# ggplot(data.sing.clean.rel.rare, aes(fill=Division, y=Value, x=Site)) + 
#   geom_bar(position="stack", stat="identity")
# 
# # this one is easier to read!
##### how are these abundant ASVs representend in the environmental samples? ########
# 
# ASV.filters <- ASV[,META$Sample_type == "filter"]
# ASV.sed <- ASV[,META$Sample_type == "sediment"]
# 
# # selecting only the abundant ASVs in the foraminifera from the environmental samples
# ASV.abundant.filters <- ASV.filters[abundant.diatoms,]
# ASV.abundant.sed <- ASV.sed[abundant.diatoms,]
# 
# # excluding lines with no ASV
# ASV.filters <- ASV.filters[rowSums(ASV.filters) >0, ] 
# ASV.sed <- ASV.sed[rowSums(ASV.sed) >0, ] 
# ASV.abundant.filters <- ASV.abundant.filters[rowSums(ASV.abundant.filters) >0, ] 
# # only two ASVs left
# ASV.abundant.sed <- ASV.abundant.sed[rowSums(ASV.abundant.sed) >0, ] 
# # only two ASVs left
# 
# # sq28 and sq31
# 
# 
# # What is the proportion of these ASVs in the eDNA?
# 
# # seawater
# sum(ASV.abundant.filters)/sum(ASV.filters) * 100 # 0.18% 
# # they are rare in the seawater eDNA
# nrow(ASV.abundant.filters)/nrow(ASV.filters) * 100 # 0.04%
# 
# # sediment
# sum(ASV.abundant.sed)/sum(ASV.sed) * 100 # 0.98%
# # they are rare in the sedimentary eDNA
# nrow(ASV.abundant.sed)/nrow(ASV.sed) * 100 # 0.04%
# 
# 
# # which taxa are these???
# TAX.sing.and.env <- TAX[rownames(ASV.sing.and.env),]
# 
# # Family Raphid-pennate
# # Genus Psammodictyon (sq21) and Navicula (sq38)
# 
# 
# # and in the foraminifera? are they also rare? 
# ASV.sing.and.env <- ASV.abundant[rownames(ASV.abundant.sed),]
# 
# sum(ASV.sing.and.env)/sum(ASV.abundant) * 100 # 0.49%
# nrow(ASV.sing.and.env)/nrow(ASV.abundant) * 100 # 3.17%
# # yes, they are rare also in the foraminifera
# # like we saw in the tree map from before in the topic "constraining dominant ASVs within foraminifera microbiome to diatoms only"
# 

##### 






############################################################################################################################################



##########################################################
################### PlotAbund function ################### 
##########################################################
PlotAbund <- function(relData, abund, 
                      margin = par()$mar,
                      method = c("nmost", "percentage"),
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




