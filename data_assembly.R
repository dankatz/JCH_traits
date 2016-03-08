#assembling the relevant datasets
library(dplyr)


##########
#getting species codes for relevant species from TRY (i.e., those in Comita et al. 2014) 
##########
rm(list = ls())

#load in and prepare data from Comita et al.
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
cdata <- read.csv("jec12232-sup-0001-TableS1.csv")

cdata2 <- select(cdata, Species) %>% distinct() %>% arrange(Species)
names(cdata2) <- "AccSpeciesName"

#load in data from the TRY species list
tryspecieslist <- read.table("TryAccSpecies.txt", sep = "\t", header = TRUE)
names(tryspecieslist)

#combine the two
test <- left_join(cdata2, tryspecieslist, by = "AccSpeciesName")

#conduct checks on the species that aren't turning up any matches
#make sure that the species wasn't misspelled

#alternate strategy for reducing missingness and dealing with it: get all species within a focal genus:
genuslist <- data.frame(gsub( " .*$", "", cdata2$AccSpeciesName))
names(genuslist) <- "genus"
genuslist <- distinct(genuslist, genus)

tryspecieslist$genus <- gsub( " .*$", "", tryspecieslist$AccSpeciesName)

test2 <- left_join(genuslist, tryspecieslist, by = "genus")
#output the species numbers for the TRY data query
spnameoutput <- writeClipboard(toString(test2$AccSpeciesID))


#what are the available traits
trytraitlist <- read.table("try_traitlist.txt", sep = "\t", header = TRUE)
head(trytraitlist)

#relevant traits that were measured for more than ~5000 species: 
#trait ids:
# 28, 42, 18, 231, 26, 1263, 11, 1, 21, 47, 38, 37, 159, 55, 343, 153, 14, 95, 33, 1140, 46, 43, 1111, 59, 604, 
# 197, 825, 8, 2, 154, 4, 31, 7, 15, 163, 50, 520, 218, 145, 1229, 318, 13, 677, 609, 403, 413, 48, 78, 53, 125,
# 138, 66, 981, 1258, 823, 892, 56, 982, 603, 98, 51, 199, 155, 587, 27, 232, 334, 40, 679, 45, 131, 341, 357, 
# 1146, 1194, 1008, 358, 1137, 146, 24, 239, 719, 41, 350, 89, 230, 344, 324, 325,, 827, 596, 238, 1132, 819, 
# 1131, 193, 1135, 30, 54, 570, 233, 128, 345, 12, 1138
