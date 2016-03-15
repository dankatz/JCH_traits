#preliminary meta-analysis using the publicly available data, pulling missing data from the genus
#analyzing using only the relevant traits

#set up work environment
library(ggplot2)
library(rjags)
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


#load in the publicly available trait data from TRY
traits <- fread("TRYtraitpublicdownload.txt") #takes ~4 min

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
  51  #Leaf photosynthesis rate per leaf area 15
  #not enough data for:   #growth rate,   #leaf water content,   #leaf phenolics,
)

traits_selected <- filter(traits, TraitID %in% selectedtraits & ErrorRisk < 4)
  
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

traits_selected$genus <- gsub( " .*$", "", traits_selected$SpeciesName)
traits_selected$genus <- tolower(traits_selected$genus)

#graphs of a particular trait at the genus level
filter(traits_selected, TraitID == 1 & ErrorRisk < 4) %>%
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



#####################
#meta-analysis
#####################

setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/meta analysis")
##########getting data into the right format
mdata <- cdata_traits

mdata$study <-  as.numeric(as.factor(mdata$Author))
##########meta-analysis


sink("model.txt")
cat("   
    model{
   
    ###### data simulation 
   
    ###### model
    for(i in 1:n_obs){
      OR_mean[i] ~ dnorm(OR[i], 1/OR_var[i])

      OR[i] <- study_fe[study[i]] #+ lifestage_fe[lifestage[i]]  
#                alpha_1 * sla[species[i]] +
#                alpha_2 * 
    }

    ###### priors
    for(s in 1:nstudies){study_fe[s] ~ dnorm(0, 0.0001)}

    }
    ",fill=TRUE)
sink() 


jags <- jags.model('model.txt',
                   data = list(
                     
                  #indices
                     'n_obs' = nrow(mdata),
                     'study' = mdata$study,
                     'nstudies' = max(mdata$study),
                     
                  #Odds ratio
                    'OR_mean' = mdata$OR_mean,
                    'OR_var' = mdata$OR_var
                    
                     
                   ),
                   n.chains = 3,
                   n.adapt = 100)


update(jags,n.iter=100) #update(jags,n.iter=1000) 
mcmc_samples <- coda.samples(jags, variable.names=c("study_fe", "OR"),  n.iter=5000)

plot(mcmc_samples)  
#dic.jags <- (dic.samples(jags, 50000, "pD")); dic.jags #linear model: pD = 179.8, quadratic pD = 174.7, exponential = 172.7

resultsall<-summary(mcmc_samples)    
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param2<-substr(results$parameter,1,2)
results$param7<-substr(results$parameter,1,7)
results$param12<-substr(results$parameter,1,12)

#predicted versus observed
mdata$predicted <-  results$Mean[results$param2=="OR"]
ggplot(mdata, aes(x = predicted, y = OR_mean)) + geom_point() + geom_smooth() + theme_bw()
summary(lm((mdata$predicted ) ~ mdata$OR_mean))
##########results


