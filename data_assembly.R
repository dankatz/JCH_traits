#assembling the relevant datasets
library(dplyr)

##########
#getting species codes for relevant species from TRY (i.e., those in Comita et al. 2014) 
##########

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


#output the species numbers for the TRY data query


