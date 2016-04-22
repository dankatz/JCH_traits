#Do plant traits correlate with DD reported in Johnson et al. 2012
#this script has data exploration

#set up work environment
library(ggplot2)
library(dplyr)
library(metafor)
library(tidyr)
library(data.table)

rm(list = ls())


#################
#load in data
#################

#load in the publicly available trait data from TRY that were downloaded on 4/21/16 (see "TRY_alldatarequests.R")
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
traits <- fread("TRYtrait_Comita2010_Johnson2012.txt") #takes ~2 min

traits_selected <- filter(traits, ErrorRisk < 4) #some TRY data are flagged as unreliable; removing those records here
traits_selected$genus <- gsub( " .*$", "", traits_selected$SpeciesName)
traits_selected$genus <- tolower(traits_selected$genus)

#creating short hand names for the traits
traits_selected$t <- NA 
traits_selected$t[traits_selected$TraitID == 11] <- "leaf_SLA"
traits_selected$t[traits_selected$TraitID == 1] <- "leaf_area"
traits_selected$t[traits_selected$TraitID == 12] <- "leaf_lifespan"
traits_selected$t[traits_selected$TraitID == 146] <- "leaf_CNratio"
traits_selected$t[traits_selected$TraitID == 14] <- "leaf_N"
traits_selected$t[traits_selected$TraitID == 50] <- "leafN_area" 
traits_selected$t[traits_selected$TraitID == 4] <- "stem_density"
traits_selected$t[traits_selected$TraitID == 138] <- "seed_n"
traits_selected$t[traits_selected$TraitID == 59] <- "plant_lifespan"
traits_selected$t[traits_selected$TraitID == 18] <- "plant_height"
traits_selected$t[traits_selected$TraitID == 15] <- "leaf_P"
traits_selected$t[traits_selected$TraitID == 51] <- "leafP_area" 
traits_selected$t[traits_selected$TraitID == 2] <- "leaf_texture" 
traits_selected$t[traits_selected$TraitID == 46] <- "leaf_thickness" 
traits_selected$t[traits_selected$TraitID == 7] <- "myco_type" 
traits_selected$t[traits_selected$TraitID == 1030] <- "myco_infection_intensity" 

#load in data from Johnson et al. 2012
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data/Johnson et al. 2012 FIA")
ddj <- read.csv("Johnson_TableS1.csv")
ddj$AccSpeciesName <- paste(ddj$genus, ddj$species, sep = " ")
ddj$genus <- tolower(ddj$genus)

#################
#linking databases 
#################

#join the trait data to the C data, based on species names
traits_by_sp <- filter(traits_selected, AccSpeciesName %in% as.character(ddj$AccSpeciesName)
                       & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(sp = AccSpeciesName, trait = t) %>%
  summarise(value = median(StdValue)) %>%
  spread(trait, value) %>% #get trait into wide format before joining
  mutate(AccSpeciesName = as.character(sp))
ddj_traits_sp <- left_join(ddj, traits_by_sp, by = "AccSpeciesName")


#join the trait data to the C data, based on genus
traits_by_genus <- filter(traits_selected, genus %in% as.character(ddj$genus)
                          & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(genus, trait = t) %>%
  summarise(value = median(StdValue)) %>%
  spread(trait, value)#get trait into wide format before joining
ddj_traits_genus <- left_join(ddj, traits_by_genus, by = "genus")




#################
#data exploration 
#################

names(ddj_traits_genus)
str(ddj_traits_genus)

ggplot(ddj_traits_genus, aes(x = log(leafP_area), y = DD, ymin = DD - DD_SE, ymax = DD + DD_SE)) + 
  geom_point() + geom_errorbar() + theme_bw() + geom_smooth()
# [9] "leaf_area"                "leaf_CNratio"             "leaf_lifespan"            "leaf_N"                  
# [13] "leaf_P"                   "leaf_SLA"                 "leaf_thickness"           "leafN_area"              
# [17] "leafP_area"               "myco_infection_intensity" "plant_height"             "plant_lifespan"          
# [21] "seed_n"                   "stem_density"             NA      

ggplot(ddj_traits_sp, aes(x = n_cells, y = DD, ymin = DD - DD_SE, ymax = DD + DD_SE)) + 
  geom_point() + geom_errorbar() + theme_bw() + geom_smooth()



