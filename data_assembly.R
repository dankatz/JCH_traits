#assembling the privately available trait data, pulling missing data from the genus, 
#analyzing using only the relevant traits

#set up work environment
library(ggplot2)
library(dplyr)
library(metafor)
library(tidyr)
library(data.table)

rm(list = ls())


#load in and prepare data from Comita et al.
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
cdata <- read.csv("jec12232-sup-0001-TableS1.csv")
ES <- escalc(measure = "OR", ai = cdata$Nosurvivorsnear_high, bi = cdata$No_died_near_high,
             ci = cdata$No_survivors_far_low, di = cdata$No_died_far_low, var.names = c("OR_mean", "OR_var"))

cdata$OR_mean <- as.numeric(ES$OR_mean)
cdata$OR_var <- as.numeric(ES$OR_var)
cdata$genus <- gsub( " .*$", "", cdata$Species)
cdata$genus <- tolower(cdata$genus)


#load in the publicly available trait data from TRY
traits <- fread("TRYtraitprivatedownload.txt") #takes ~2 min

#selected traits, based on what JG thought was important and what traits had the largest amount
#of publicly available data in TRY
selectedtraits <- c( #trait description, number of observations on the species level 
  11,  #Leaf area per leaf dry mass (specific leaf area, SLA), 39
  1,  #Leaf area, 39 
  12,  #Leaf lifespan (longevity), 13
  146,  #Leaf carbon/nitrogen (C/N) ratio, 12 
  14,  #Leaf nitrogen (N) content per leaf dry mass, 40
  4,  #Stem dry mass per stem fresh volume (stem specific d, 55
  138,  #Seed number per reproducton unit, 5
  59,  #Plant lifespan (longevity), 13
  
  #other potentially useful traits,
  18,  #Plant height, 33 
  15,  #Leaf phosphorus (P) content per leaf dry mass, 28
  51  #Leaf phosphorus (P) per leaf area 15
  #not enough data for:   #growth rate,   #leaf water content,   #leaf phenolics,
)

traits_selected <- filter(traits, TraitID %in% selectedtraits & ErrorRisk < 4)
traits_selected$genus <- gsub( " .*$", "", traits_selected$SpeciesName)
traits_selected$genus <- tolower(traits_selected$genus)

#join the trait data to the C data, based on species names
traits_by_sp <- filter(traits_selected, AccSpeciesName %in% as.character(cdata$Species) & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(sp = AccSpeciesName, trait = TraitName) %>%
  summarise(value = median(StdValue)) %>%
  spread(trait, value) %>% #get trait into wide format before joining
  mutate(Species = as.character(sp))
cdata_traits <- left_join(cdata, traits_by_sp, by = "Species")
 
####### 
#summary of traits at the genus level and comparing to the species level traits
#######

#graphs of a particular trait at the genus level
filter(traits_selected, TraitID == 51 ) %>%
  ggplot(aes(x = genus, y = StdValue)) + geom_boxplot() + geom_jitter(alpha = 0.1) + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#graphs of correlation between species level data and genus level data
#join the trait data to the C data, based on species names
traits_by_sp_genus <- filter(traits_selected, AccSpeciesName %in% as.character(cdata$Species) & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(sp = AccSpeciesName, trait = TraitName, genus) %>%
  summarise(species_median = median(StdValue),
            species_mean = mean(StdValue),
            species_sd = sd(StdValue),
            species_n = n())

genus_level_traits <- group_by(traits_selected, trait = TraitName, genus) %>%
  summarise(genus_median = median(StdValue, na.rm = TRUE),
            genus_mean = mean(StdValue, na.rm = TRUE),
            genus_sd = sd(StdValue, na.rm = TRUE),
            genus_n = n())

traits_by_sp_genus <- left_join(traits_by_sp_genus, genus_level_traits, by = c("trait", "genus"))

ggplot(traits_by_sp_genus, aes(x = genus_mean, y = species_mean)) +  theme_bw() + 
  geom_errorbar(aes(ymin = species_mean - species_sd, ymax = species_mean + species_sd), color = "light gray") + 
  geom_errorbarh(aes(xmin = genus_mean - genus_sd, xmax = genus_mean + genus_sd), color = "light gray") + 
  geom_point() + ylab("species mean + SD") + xlab("genus mean + SD") + 
  facet_wrap(~trait, scales = "free")


#########
#creating the data file for the meta-analysis
#########

#getting short hand names for selected traits
  traits_selected$t <- NA
  traits_selected$t[traits_selected$TraitID == 11] <- "leaf_SLA"
  traits_selected$t[traits_selected$TraitID == 1] <- "leaf_area"
  traits_selected$t[traits_selected$TraitID == 12] <- "leaf_lifespan"
  traits_selected$t[traits_selected$TraitID == 146] <- "leaf_CNratio"
  traits_selected$t[traits_selected$TraitID == 14] <- "leaf_N"
  traits_selected$t[traits_selected$TraitID == 4] <- "stem_density"
  traits_selected$t[traits_selected$TraitID == 138] <- "seed_n"
  traits_selected$t[traits_selected$TraitID == 59] <- "plant_lifespan"
  traits_selected$t[traits_selected$TraitID == 18] <- "plant_height"
  traits_selected$t[traits_selected$TraitID == 15] <- "leaf_P"
  traits_selected$t[traits_selected$TraitID == 51] <- "leafP_area"  #missing _ prevents overlapping name w/ leaf_P

###########adding traits at the species level to mdata
  mdata <- cdata
  
  #trait mean
  traits_by_sp_mean <- filter(traits_selected, AccSpeciesName %in% as.character(cdata$Species) & !is.na(StdValue) & !is.na(TraitID)) %>%
    group_by(AccSpeciesName, t) %>%
    summarise(sp_trait_mean = mean(StdValue)) %>%
    spread(t, sp_trait_mean) %>% #get trait into wide format before joining
    mutate(Species = as.character(AccSpeciesName))
  
    traits_by_sp_mean$AccSpeciesName <- NULL
    name_rename <- names(traits_by_sp_mean)
    name_rename <- c(paste(name_rename[1:11], "_mean","_s",sep =""),"Species")
    setnames(traits_by_sp_mean, names(traits_by_sp_mean), name_rename)
    
    mdata <- left_join(mdata, traits_by_sp_mean, by = "Species")
    
    #trait sd
    traits_by_sp_sd <- filter(traits_selected, AccSpeciesName %in% as.character(cdata$Species) & !is.na(StdValue) & !is.na(TraitID)) %>%
      group_by(AccSpeciesName, trait = t) %>%
      summarise(sp_trait_sd = sd(StdValue)) %>%
      spread(trait, sp_trait_sd) %>% #get trait into wide format before joining
      mutate(Species = as.character(AccSpeciesName))
    
    traits_by_sp_sd$AccSpeciesName <- NULL
    name_rename <- names(traits_by_sp_sd)
    name_rename <- c(paste(name_rename[1:11], "_sd","_s",sep =""),"Species")
    setnames(traits_by_sp_sd, names(traits_by_sp_sd), name_rename)
    
    mdata <- left_join(mdata, traits_by_sp_sd, by = "Species")


###########adding traits at the genus level to mdata
    #trait mean
    genus_level_traits_mean <- group_by(traits_selected, trait = t, genus) %>%
      summarise(genus_mean = mean(StdValue, na.rm = TRUE)) %>%
      spread(trait, genus_mean) #get trait into wide format before joining
    
    name_rename <- names(genus_level_traits_mean)
    name_rename <- c("genus",paste(name_rename[2:12], "_mean","_g",sep =""))
    setnames(genus_level_traits_mean, names(genus_level_traits_mean), name_rename)
    
    mdata <- left_join(mdata, genus_level_traits_mean, by = "genus")
    
    #trait sd
    genus_level_traits_sd <- group_by(traits_selected, trait = t, genus) %>%
      summarise(genus_sd = sd(StdValue, na.rm = TRUE)) %>%
      spread(trait, genus_sd) #get trait into wide format before joining
    
    name_rename <- names(genus_level_traits_sd)
    name_rename <- c("genus",paste(name_rename[2:12], "_sd","_g",sep =""))
    setnames(genus_level_traits_sd, names(genus_level_traits_sd), name_rename)
    
    mdata <- left_join(mdata, genus_level_traits_sd, by = "genus")
    
    #trait n
    genus_level_traits_n <- group_by(traits_selected, trait = t, genus) %>%
      summarise(genus_n = n()) %>%
      spread(trait, genus_n) #get trait into wide format before joining
    
    name_rename <- names(genus_level_traits_n)
    name_rename <- c("genus",paste(name_rename[2:12], "_n","_g",sep =""))
    setnames(genus_level_traits_n, names(genus_level_traits_n), name_rename)
    
    mdata <- left_join(mdata, genus_level_traits_n, by = "genus")
    
    
    #converting latitude to decimal degrees
    library(stringr)
    library(sp)
    deg <- substr(as.character(mdata$Latitude), 1,2)
    min <- substr(as.character(mdata$Latitude), 4,5)
    sec <- substr(as.character(mdata$Latitude), 7,11)
    direction <- substr(as.character(mdata$Latitude), 13,13)
    lat_dms <- char2dms(from = paste(deg, "d", min, "m", sec, "s", direction, sep = ""), chd = "d", chm = "m", chs = "s")
    lat_dd <- as.numeric(lat_dms)
    mdata$lat_dd <- lat_dd
    
    #converting longitude to decimal degrees
    longstr <- as.character(mdata$Longitude)
    longstr[4] <- "88:59:09.59:W" #manuallly fixing a format mistake
    longstr[79] <- "107:23:06.42:W" #manuallly fixing a format mistake
    longstr <- str_pad(longstr, width = max(nchar(longstr)), side = "left", pad = " ") #the strings are different lengths
    deg <- as.numeric(substr(longstr, 1,3))
    min <- substr(longstr, 5,6)
    sec <- substr(longstr, 8,12)
    direction <- substr(longstr, 14,14)
    long_dms <- char2dms(from = paste(deg, "d", min, "m", sec, "s", direction, sep = ""), chd = "d", chm = "m", chs = "s")
    long_dd <- as.numeric(long_dms)
    mdata$long_dd <- long_dd
    
    ##########saving file so it doesn't have to be re-constructed each time
    setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/meta analysis")
    write.csv(mdata, "mdata.csv")
