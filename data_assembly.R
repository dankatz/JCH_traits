#assembling the relevant datasets
library(dplyr)
library(data.table)
library(metafor)
library(ggplot2)

##########
#getting species codes for relevant species from TRY (i.e., those in Comita et al. 2014) 
##########
rm(list = ls())

#load in and prepare data from Comita et al.
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
cdata <- read.csv("jec12232-sup-0001-TableS1.csv")
ES <- escalc(measure = "OR", ai = cdata$Nosurvivorsnear_high, bi = cdata$No_died_near_high,
         ci = cdata$No_survivors_far_low, di = cdata$No_died_far_low, var.names = c("OR_mean", "OR_var"))

cdata$OR_mean <- as.numeric(ES$OR_mean)
cdata$OR_var <- as.numeric(ES$OR_var)

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


#what are the available traits
trytraitlist <- read.table("try_traitlist.txt", sep = "\t", header = TRUE)
head(trytraitlist)

#relevant traits that were measured for more than ~5000 species: 
relevanttraits <- c(28, 42, 18, 231, 26, 1263, 11, 1, 21, 47, 38, 37, 159, 55, 343, 153, 14, 95, 33, 1140, 46, 43, 1111, 59, 604, 
197, 825, 8, 2, 154, 4, 31, 7, 15, 163, 50, 520, 218, 145, 1229, 318, 13, 677, 609, 403, 413, 48, 78, 53, 125,
138, 66, 981, 1258, 823, 892, 56, 982, 603, 98, 51, 199, 155, 587, 27, 232, 334, 40, 679, 45, 131, 341, 357, 
1146, 1194, 1008, 358, 1137, 146, 24, 239, 719, 41, 350, 89, 230, 344, 324, 325, 827, 596, 238, 1132, 819, 
1131, 193, 1135, 30, 54, 570, 233, 128, 345, 12, 1138)

relevanttraitlist <- trytraitlist[trytraitlist$TraitID %in% relevanttraits,]

#adding in downloaded public trait data
traits <- fread("TRYtraitpublicdownload.txt") #takes ~4 min
head(traits)
names(traits)


#graphing the various traits; note that this is for data exploration NOT for model selection or data dredging

smalltraits2 <- filter(traits, TraitID == 1229)
smalltraits2_sp <- group_by(smalltraits2, SpeciesName) %>% summarize(loopedvar = mean(StdValue))
smalltraits2_sp$Species <- smalltraits2_sp$SpeciesName

cdata_loopedvar <- left_join(cdata, smalltraits2_sp, by = "Species")

print(nrow(subset(cdata_loopedvar, !is.na(loopedvar))))

ggplot(cdata_loopedvar, aes(x = loopedvar, y = OR_mean, ymax = OR_mean + OR_var, ymin = OR_mean - OR_var)) + 
  geom_pointrange() + theme_bw() + geom_smooth(method = "lm") + xlab(unique(smalltraits2$TraitName))



#parred down list of traits
relevanttraits <- c(18, #plant height
                    1, #sla
                    231, 26, 1263, 11, 1, 21, 47, 38, 37, 159, 55, 343, 153, 14, 95, 33, 1140, 46, 43, 1111, 59, 604, 
                    197, 825, 8, 2, 154, 4, 31, 7, 15, 163, 50, 520, 218, 145, 1229, 318, 13, 677, 609, 403, 413, 48, 78, 53, 125,
                    138, 66, 981, 1258, 823, 892, 56, 982, 603, 98, 51, 199, 155, 587, 27, 232, 334, 40, 679, 45, 131, 341, 357, 
                    1146, 1194, 1008, 358, 1137, 146, 239, 719, 41, 350, 89, 230, 344, 324, 325, 827, 596, 238, 1132, 819, 
                    1131, 193, 1135, 30, 54, 570, 233, 128, 345, 12, 1138)




