#Do plant traits correlate with DD reported in Comita et al. 2010?
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
    
  #load in data from Comita et al. 2010
    setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data/Comita et al. 2010")
    comita_table1 <- read.csv("Comita_TableS1.csv", stringsAsFactors = FALSE)
    comita_table2 <- read.csv("Comita_TableS2.csv", stringsAsFactors = FALSE)
    ddc <- left_join(comita_table1, comita_table2, by = "sp_code")
    ddc$AccSpeciesName <- paste(ddc$genus, ddc$species, sep = " ")
    ddc$genus <- tolower(ddc$genus)
    
    
#################
#linking databases 
#################

#join the trait data to the C data, based on species names
traits_by_sp <- filter(traits_selected, AccSpeciesName %in% as.character(ddc$AccSpeciesName)
                       & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(sp = AccSpeciesName, trait = t) %>%
  summarise(value = median(StdValue)) %>%
  spread(trait, value) %>% #get trait into wide format before joining
  mutate(AccSpeciesName = as.character(sp))
ddc_traits_sp <- left_join(ddc, traits_by_sp, by = "AccSpeciesName")

#join the trait data to the C data, based on genus
traits_by_genus <- filter(traits_selected, genus %in% as.character(ddc$genus)
                          & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(genus, trait = t) %>%
  summarise(value = median(StdValue)) %>%
  spread(trait, value)#get trait into wide format before joining
ddc_traits_genus <- left_join(ddc, traits_by_genus, by = "genus")


#################
#data exploration 
#################

names(ddc_traits_genus)
str(ddc_traits_genus)

ggplot(ddc_traits_genus, aes(x = stem_density, y = conspecificadult, 
                             ymin = conspecificadult - conspecificadult_sd, 
                             ymax = conspecificadult + conspecificadult_sd)) + 
  geom_point() + geom_errorbar(alpha = 0.5) + theme_bw() + geom_smooth()

# [9] "leaf_area"                "leaf_CNratio"             "leaf_lifespan"            "leaf_N"                  
# [13] "leaf_P"                   "leaf_SLA"                 "leaf_thickness"           "leafN_area"              
# [17] "leafP_area"               "myco_infection_intensity" "plant_height"             "plant_lifespan"          
# [21] "seed_n"                   "stem_density"             N     

ggplot(ddc_traits_sp, aes(x = leaf_N, 
                          y = conspecificadult, 
                             ymin = conspecificadult - conspecificadult_sd, 
                             ymax = conspecificadult + conspecificadult_sd)) + 
  geom_point() + geom_errorbar(alpha = 0.5) + theme_bw() + geom_smooth()


ggplot(ddc_traits_sp, aes(x = log(leaf_P * leaf_SLA), y = conspecificadult, 
                             ymin = conspecificadult - conspecificadult_sd, 
                             ymax = conspecificadult + conspecificadult_sd)) + 
  geom_point() + geom_errorbar(alpha = 0.5) + theme_bw() + geom_smooth()






