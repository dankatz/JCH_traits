#Do plant traits correlate with DD reported in Comita et al. 2014, Comita et al. 2010, Johnson et al. 2012, Katz 2016?
#this script assembles the species and traits required for TRY downloads

#set up work environment
library(ggplot2)
library(dplyr)
library(metafor)
library(tidyr)
library(data.table)

rm(list = ls())

########################################################
#traits to request from TRY
########################################################

############which traits to include:
#selected traits, based on what JG thought was important and what traits had the largest amount
#of publicly available data in TRY
selectedtraits <- c( #trait description, number of observations on the species level 
  11,  #Leaf area per leaf dry mass (specific leaf area, SLA), 39
  1,  #Leaf area, 39 
  12,  #Leaf lifespan (longevity), 13
  146,  #Leaf carbon/nitrogen (C/N) ratio, 12 
  14,  #Leaf nitrogen (N) content per leaf dry mass, 40
  50, #Leaf nitrogen (N) content per leaf area, ?
  4,  #Stem dry mass per stem fresh volume (stem specific d, 55
  138,  #Seed number per reproducton unit, 5
  59,  #Plant lifespan (longevity), 13
  
  #other potentially useful traits,
  18,  #Plant height, 33 
  15,  #Leaf phosphorus (P) content per leaf dry mass, 28
  51,  #Leaf phosphorus (P) per leaf area 15
  2,   #Leaf texture (sclerophylly, physical strength), ?
  46,	#Leaf thickness, ?
  53, #	Leaf photosynthesis rate per leaf area
  7,	#Plant mycorrhizal type
  1030	#Plant mycorrhizal infection intensity
  #not enough data for:   #growth rate,   #leaf water content,   #leaf phenolics,
)
writeClipboard(toString(selectedtraits)) #paste this into the TRY web interface

# Here are all the traits that were downloaded in the first round of data exploration
# relevant traits that were measured for more than ~5000 species: 
# relevanttraits <- c(28, 42, 18, 231, 26, 1263, 11, 1, 21, 47, 38, 37, 159, 55, 343, 153, 14, 95, 33, 1140, 46, 43, 1111, 59, 604, 
#                     197, 825, 8, 2, 154, 4, 31, 7, 15, 163, 50, 520, 218, 145, 1229, 318, 13, 677, 609, 403, 413, 48, 78, 53, 125,
#                     138, 66, 981, 1258, 823, 892, 56, 982, 603, 98, 51, 199, 155, 587, 27, 232, 334, 40, 679, 45, 131, 341, 357, 
#                     1146, 1194, 1008, 358, 1137, 146, 24, 239, 719, 41, 350, 89, 230, 344, 324, 325, 827, 596, 238, 1132, 819, 
#                     1131, 193, 1135, 30, 54, 570, 233, 128, 345, 12, 1138)






########################################################
#genera/species to request from TRY
########################################################



#################
#getting genera that are in Comita et al. 2014
#################

#load in and prepare data from Comita et al.
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
cdata <- read.csv("jec12232-sup-0001-TableS1.csv")
ES <- escalc(measure = "OR", ai = cdata$Nosurvivorsnear_high, bi = cdata$No_died_near_high,
         ci = cdata$No_survivors_far_low, di = cdata$No_died_far_low, var.names = c("OR_mean", "OR_var"))

#list of species included in meta-analysis
cdata2 <- select(cdata, Species) %>% distinct() %>% arrange(Species)
names(cdata2) <- "AccSpeciesName"

#load in data from the TRY species list
tryspecieslist <- read.table("TryAccSpecies.txt", sep = "\t", header = TRUE)
names(tryspecieslist)

#combine the two
overlap_csp_trysp <- left_join(cdata2, tryspecieslist, by = "AccSpeciesName")

#conduct checks on the species that aren't turning up any matches
#make sure that the species wasn't misspelled

#alternate strategy for reducing missingness and dealing with it: get all species within a focal genus:
genuslist <- data.frame(gsub( " .*$", "", cdata2$AccSpeciesName))
names(genuslist) <- "genus"
genuslist <- distinct(genuslist, genus)

tryspecieslist$genus <- gsub( " .*$", "", tryspecieslist$AccSpeciesName)

overlap_cg_tryg <- left_join(genuslist, tryspecieslist, by = "genus")

#output the species numbers for the TRY data query
writeClipboard(toString(overlap_cg_tryg$AccSpeciesID))




#################
#getting genera that are in Comita et al. 2010 and Johnson et al. 2012
#################

#load in the TRY species list
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
tryspecieslist <- read.table("TryAccSpecies.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
tryspecieslist$genus <- gsub( " .*$", "", tryspecieslist$AccSpeciesName)
names(tryspecieslist)

#load in species that were included in Comita et al. 2010
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data/Comita et al. 2010")
comita_sptable <- read.csv("Comita_TableS1.csv", stringsAsFactors = FALSE)
comita_sp <- select(comita_sptable, genus, species)

#load in species that were included in Johnson et al.
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data/Johnson et al. 2012 FIA")
johnson_sptable <- read.csv("Johnson_TableS1.csv", stringsAsFactors = FALSE)
johnson_sp <- select(johnson_sptable, genus, species)

#combine the two species lists
sptodownload <- rbind(comita_sp, johnson_sp)
sptodownload <- mutate(sptodownload, AccSpeciesName = paste(sptodownload$genus, sptodownload$species)) %>%
  select(AccSpeciesName) %>%
  distinct()
overlap_csp_trysp <- left_join(sptodownload, tryspecieslist, by = "AccSpeciesName")

#get a list at the genus level (an alternate strategy for reducing missingness and dealing with it: get all species within a focal genus)
gentodownload <- rbind(comita_sp, johnson_sp)
gentodownload <- select(gentodownload, genus) %>%
  distinct()

overlap_cg_tryg <- left_join(gentodownload, tryspecieslist, by = "genus") 

#getting it down to 9,999 sp (otherwise I'd get all spp, which would be ~10x the size and less easy to work with)
overlap_cg_tryg_2 <- filter(overlap_cg_tryg, !is.na(ObsNum) & ObsNum > 1) #not getting sp that only have one trait measured in the database

#output the species numbers for the TRY data query
writeClipboard(toString(overlap_cg_tryg_2$AccSpeciesID))






#################
#getting genera for plant pop. growth rates paper (Katz 2016)
#################

rm(list = ls())

#load in data from the TRY species list
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
tryspecieslist <- read.table("TryAccSpecies.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
names(tryspecieslist)

#load in data from the plant population growth metaanalysis
setwd("Q:/Ibanez Lab/Dan Katz/DISSERTATION/chapter 1 -  herbivores and plant populations meta analysis")
sp <- read.csv("species.csv", stringsAsFactors = FALSE)
names(sp) <- "AccSpeciesName"
sp$AccSpeciesName[sp$AccSpeciesName == "Jacobaea vulgaris"] <- "Senecio jacobaea" #there was a name change
sp$AccSpeciesName[sp$AccSpeciesName == "Opuntia imbricata"] <- "Cylindropuntia imbricata" #there was a name change

#combine the two
overlap_csp_trysp <- left_join(sp, tryspecieslist, by = "AccSpeciesName")

#alternate strategy for reducing missingness and dealing with it: get all species within a focal genus:
genuslist <- data.frame(gsub( " .*$", "", sp$AccSpeciesName))
names(genuslist) <- "genus"
genuslist <- distinct(genuslist, genus)

tryspecieslist$genus <- gsub( " .*$", "", tryspecieslist$AccSpeciesName)
overlap_cg_tryg <- left_join(genuslist, tryspecieslist, by = "genus")

#output the species numbers for the TRY data query
writeClipboard(toString(overlap_cg_tryg$AccSpeciesID))





