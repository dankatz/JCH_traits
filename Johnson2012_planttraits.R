#Do plant traits correlate with DD reported in Johnson et al. 2012
#this script has data exploration

#set up work environment
library(ggplot2)
library(dplyr)
library(metafor)
library(tidyr)
library(data.table)
library(stringr)

rm(list = ls())


#################
#load in data
#################

#load in the privately available trait data from TRY that were downloaded on 5/09/16 (see "TRY_alldatarequests.R")
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")
traits <- fread("TRYtrait_Comita2010_Johnson2012_private.txt") #takes ~2 min

traits_selected <- filter(traits, ErrorRisk < 4) #some TRY data are flagged as unreliable; removing those records here
traits_selected$genus <- gsub( " .*$", "", traits_selected$SpeciesName)
traits_selected$genus <- tolower(traits_selected$genus)

#NOTE: NEED TO WRITE A FUNCTION TO REMOVE THE AUTHORITY FROM SPECIES NAMES (IT IS OCCASIONALLY INCLUDED)
  
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
traits_selected$t[traits_selected$TraitID == 53] <- "leaf_photosyn" 


#load in data from Johnson et al. 2012
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data/Johnson et al. 2012 FIA")
ddj <- read.csv("Johnson_TableS1_plus.csv")
ddj$AccSpeciesName <- paste(ddj$genus, ddj$species, sep = " ")
ddj$genus <- tolower(ddj$genus)

#a manual check revealed a typo in a species name between Johnson 2012 and TRY; fixing it here
#ddj$AccSpeciesName[ddj$AccSpeciesName == "Quercus margarettae"] <- "Quercus margaretta"



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
#formatting data for analysis
#################

#select species that overlap in both Johnson and traits and get mean, sd, and n for each
  traits_by_sp <- filter(traits_selected, AccSpeciesName %in% as.character(ddj$AccSpeciesName)
                       & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(AccSpeciesName, trait = t, genus) %>%
  summarise(trait_median_s = median(StdValue),
            trait_mean_s = mean(StdValue),
            trait_sd_s = sd(StdValue),
            trait_n_s = n()) 

#select genera that overlap in both Johnson and traits and get mean, sd, and n for each
traits_by_genus <- filter(traits_selected, genus %in% as.character(ddj$genus)
                       & !is.na(StdValue) & !is.na(TraitID)) %>%
  group_by(genus, trait = t) %>%
  summarise(trait_median_g = median(StdValue),
            trait_mean_g = mean(StdValue),
            trait_sd_g = sd(StdValue),
            trait_n_g = n())  

#combine genera and species traits into one dataframe (long)
ddj_traits <- left_join(ddj, traits_by_genus, by = c("genus"))
ddj_traits <- left_join(ddj_traits, traits_by_sp, by = c("AccSpeciesName", "genus", "trait"))


#checking how well traits are conserved at the genus vs. species level  
  ggplot(ddj_traits, aes(x = log(trait_mean_g), xmin = log(trait_mean_g - trait_sd_g), xmax = log(trait_mean_g + trait_sd_g),
                         y = log(trait_mean_s), ymin = log(trait_mean_s - trait_sd_s), ymax = log(trait_mean_s + trait_sd_s))) + 
  geom_point() + geom_errorbar(alpha = 0.2, color = "gray") + geom_errorbarh(alpha = 0.2, color = "gray") + 
  facet_wrap(~trait) + theme_bw()
  
#the same, but for a particular trait
  filter(ddj_traits, trait == "leaf_N") %>%
  ggplot(aes(x = trait_mean_g, xmin = trait_mean_g - trait_sd_g, xmax = trait_mean_g + trait_sd_g,
             y = trait_mean_s, ymin = trait_mean_s - trait_sd_s, ymax = trait_mean_s + trait_sd_s)) + 
    geom_point() + geom_errorbar(alpha = 0.7, color = "gray") + geom_errorbarh(alpha = 0.7, color = "gray") + theme_bw()
  

#using species trait data when present, genera data otherwise
  ddj_traits$trait_median <- ddj_traits$trait_median_g
  ddj_traits$trait_mean <- ddj_traits$trait_mean_g
  ddj_traits$trait_sd <- ddj_traits$trait_sd_g
  ddj_traits$trait_n <- ddj_traits$trait_n_g
  ddj_traits$trait_sp_or_g <- "genus"
  
  for(i in 1:nrow(ddj_traits)){
    if(!is.na(ddj_traits$trait_median_s[i])){ddj_traits$trait_median[i] <- ddj_traits$trait_median_s[i] }
    if(!is.na(ddj_traits$trait_mean_s[i])){ddj_traits$trait_mean[i] <- ddj_traits$trait_mean_s[i] }
    if(!is.na(ddj_traits$trait_sd_s[i])){ddj_traits$trait_sd[i] <- ddj_traits$trait_sd_s[i] }
    if(!is.na(ddj_traits$trait_n_s[i])){ddj_traits$trait_n[i] <- ddj_traits$trait_n_s[i] }
    if(!is.na(ddj_traits$trait_mean_s[i])){ddj_traits$trait_sp_or_g[i] <- "species" }
   }
  
  ddj_traits$DD_SD <- ddj_traits$DD_SE * sqrt(ddj_traits$n_cells)
  
  

#saving traits so don't need to reload all the TRY data if workspace is cleared
  setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")  
  write.csv(ddj_traits, "ddj_traits.csv")  
  
  
  
  

#################
#data exploration 
#################

names(ddj_traits_genus)
str(ddj_traits_genus)

ggplot(ddj_traits_genus, aes(x = leaf_N, y = DD, ymin = DD - DD_SE, ymax = DD + DD_SE)) + 
  geom_point() + geom_errorbar() + theme_bw() + geom_smooth(method = "lm")
# "leaf_area"                "leaf_CNratio"            
# "leaf_lifespan"            "leaf_N"                   "leaf_P"                   "leaf_photosyn"            "leaf_SLA"                
# "leaf_thickness"           "leafN_area"               "leafP_area"               "myco_infection_intensity" "plant_height"            
# "plant_lifespan"           "seed_n"                   "stem_density"    

ggplot(ddj_traits_genus, aes(x = leaf_N, y = DD, ymin = DD - DD_SE, ymax = DD + DD_SE, 
     col = leaf_N)) + 
  scale_color_gradient2(low = "red", mid = "blue", high = "green", midpoint = .5) +
  geom_point() + geom_errorbar() + theme_bw() + geom_smooth(method = "lm")


ggplot(ddj_traits_sp, aes(x = leaf_N, y = DD, ymin = DD - DD_SE, ymax = DD + DD_SE, 
  col = leaf_N)) + 
  scale_color_gradient2(low = "red", mid = "blue", high = "green", midpoint = .5) +
  geom_point() + geom_errorbar() + theme_bw() + geom_smooth(method = "lm")






filter(ddj_traits, trait == "leaf_N" ) %>%
ggplot(aes(x = trait_mean_s, xmin = trait_mean_s - trait_sd_s, xmax = trait_mean_s + trait_sd_s,
           y = DD, ymin = DD - DD_SD, ymax = DD + DD_SD, color = exclude)) +  
  geom_errorbar(alpha = 0.1) + geom_errorbarh(alpha = 0.1) +
  geom_point() + theme_bw() + geom_smooth(method = "lm", se= FALSE) + 
  xlab("leaf Nitrogen (mean + SD)") + ylab("CNDD (mean + SD)") +
  geom_abline(intercept = 0, slope = 0, lty = 2) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~genus)



fit <- lm(ddj_traits_sp$DD ~ (ddj_traits_sp$leaf_N))

summary(fit)

ddj_traits_subset <- subset(ddj_traits, trait == "leaf_N" & trait_sp_or_g == "species")

fit <- lm(ddj_traits_subset$DD ~ ddj_traits_subset$trait_mean_s)
summary(fit)
plot(ddj_traits_subset$trait_mean_s, ddj_traits_subset$DD, abline(fit))

plot(ddj_traits_subset$n_cells, ddj_traits_subset$DD)

names(ddj_traits_subset)
  
  
  
  
  
  
###################################################################################
# model for a single variable at a time
###################################################################################
  setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/data")  
  ddj_traits <- read.csv("ddj_traits.csv")
  setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/meta analysis/fia")  
  library(rjags)
  
  
#####################
#formatting data
#####################


  varlist <- unique(ddj_traits$trait)
  
  #loop through all variables to see which ones look most important
  for(v in 1:15){  
  varsubset <- varlist[v]  #varsubset <- "leaf_N"
  
  ddj_traits_subset <- subset(ddj_traits, trait == varsubset & trait_sp_or_g == "species")
  #traits_selected_subset <- subset(traits_selected, trait == varsubset)
  
  
#####################
#meta-analysis
#####################
  
  sink("model.txt")
  cat(" 
      data{

      for(i in 1:n_obs){ 

      #missing data
      variable1_mean[i] ~ dunif(variable1_mean_min, variable1_mean_max)
      variable1_sd[i] ~ dunif(variable1_sd_min, variable1_sd_max)

      #obsv'd data
        variable1[i] ~ dnorm(variable1_mean[i], (1/variable1_sd[i]^2))
        variable1_p[i] <- max(variable1[i], 0.001) #to prevent it from going below zero
      }

      }

      model{
      


      ###### missing response variable
      for(i in 1:n_obs){
      #dd[i] ~ dnorm(dd_mean[i], (1/dd_sd[i])^2 )   
      dd_sd[i] ~ dunif(dd_sd_min, dd_sd_max)  
      }
      
      ###### model
      for(i in 1:n_obs){
      dd_mean[i] ~ dnorm(dd[i], (1/dd_sd[i])^2 )
      
      dd[i] <-  alpha_1 * variable1_p[i] + intercept + residual[i]
      
      dd_pred[i] <- alpha_1 * variable1_p[i] + intercept
      }
      
      ###### priors
      alpha_1 ~ dnorm(0, 0.0001)
      intercept ~ dnorm(0, 0.0001)
      for(i in 1:n_obs){residual[i] ~ dnorm(0, 1/(sigma * sigma))}
      sigma ~ dunif(0, 10)
      
      ########## simulation for figures
      
      #######predictions for alpha_1
        for(p in 1:npredlevels){
           sim_alpha_1[p] <-  alpha_1 * predlevels[p] + intercept 
        }

      }
      ",fill=TRUE)
  sink() 
  
  jags <- jags.model('model.txt',
                     data = list(
                       
                     #indices
                       'n_obs' = nrow(ddj_traits_subset),
                       
                     #response variable
                       'dd_mean' = ddj_traits_subset$DD,
                       'dd_sd' = ddj_traits_subset$DD_SD,
                       'dd_sd_min' = min(ddj_traits_subset$DD_SD, na.rm = TRUE),
                       'dd_sd_max' = max(ddj_traits_subset$DD_SD, na.rm = TRUE),
                       
                     #covariates 
                       'variable1_mean' = ddj_traits_subset$trait_mean,
                       'variable1_sd' = ddj_traits_subset$trait_sd + 0.0001,
                       'variable1_mean_min' = min(ddj_traits_subset$trait_mean, na.rm = TRUE),
                       'variable1_mean_max' = max(ddj_traits_subset$trait_mean, na.rm = TRUE),
                       'variable1_sd_min' = min(ddj_traits_subset$trait_sd + 0.0001, na.rm = TRUE),
                       'variable1_sd_max' = max(ddj_traits_subset$trait_sd + 0.0001, na.rm = TRUE),
#                       'variable1_mean_overall' = mean(mdata_subset$leafP_area_mean, na.rm = TRUE),

                       'npredlevels' = 20,
                       'predlevels' = seq(from = min(ddj_traits_subset$trait_mean), to = max(ddj_traits_subset$trait_mean), length.out = 20)
                       
#                        'variable2_mean' = mdata_subset$plant_height_mean,
#                        'variable2_sd' = mdata_subset$plant_height_sd,
#                        'variable2_mean_min' = min(mdata_subset$plant_height_mean, na.rm = TRUE),
#                        'variable2_mean_max' = max(mdata_subset$plant_height_mean, na.rm = TRUE),
#                        'variable2_sd_min' = min(mdata_subset$plant_height_sd, na.rm = TRUE),
#                        'variable2_sd_max' = max(mdata_subset$plant_height_sd, na.rm = TRUE),
#                        'variable2_mean_overall' = mean(mdata_subset$plant_height_mean, na.rm = TRUE)
                       
                     ),
                     n.chains = 3,
                     n.adapt = 100)
  
  #parameter of interest results
  update(jags,n.iter=100) #update(jags,n.iter=20000) 
  mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1", "intercept"),  n.iter=5000)
  plot(mcmc_samples)  
  resultsall<-summary(mcmc_samples)    
  results_param <- data.frame(resultsall$statistics,resultsall$quantiles) 
  results_param$parameter<-row.names(results_param)
  results_param$param2<-substr(results_param$parameter,1,2)
  

  #modelfit
  mcmc_samples_modelfit <- coda.samples(jags, variable.names=c("dd_pred"),  n.iter=15000)
  resultsall<-summary(mcmc_samples_modelfit)    
  results<-data.frame(resultsall$statistics,resultsall$quantiles) 
  results$parameter<-row.names(results)
  
 
  modelfit <- lm(ddj_traits_subset$DD ~ results$Mean )
  modelfitsummary <- summary(modelfit)
  modelR2 <- round(modelfitsummary$r.squared, 3)
  plot(ddj_traits_subset$DD, results$Mean, abline(modelfit))
#   
#   
#   #test
#   mcmc_samples <- coda.samples(jags, variable.names=c("variable1_p"),  n.iter=5000)
#   resultsall<-summary(mcmc_samples)    
#   results<-data.frame(resultsall$statistics,resultsall$quantiles) 
#   results$parameter<-row.names(results)
#   
#   plot(results$Mean, ddj_traits_subset$trait_mean)  
#   
#   test <- data.frame(cbind(results$Mean, results$SD, ddj_traits_subset$trait_mean, ddj_traits_subset$trait_sd))
#   names(test) <- c("pred_mean", "pred_sd","obsvd_mean","obsvd_sd")
#   ggplot(test, aes(x = obsvd_mean, xmin = obsvd_mean - obsvd_sd, xmax = obsvd_mean + obsvd_sd,
#                    y = pred_mean, ymin = pred_mean - pred_sd, ymax = pred_mean + pred_sd)) + 
#   geom_point() + geom_errorbar() + geom_errorbarh() + theme_bw() + geom_abline(slope = 1)
#   
  ##############
  #fig 1: simulated results: alpha_1 
  ##############
  mcmc_samples <- coda.samples(jags, variable.names=c("sim_alpha_1"),  n.iter=5000)
  simresultsall <- summary(mcmc_samples)    
  simresults <- data.frame(simresultsall$statistics,simresultsall$quantiles) 
  simresults$parameter<-row.names(simresults)
  simresults$param2<-substr(simresults$parameter,1,2)
  
  #hardcoded to match the number of prediction increments
  simresults$simtraitlevel <- seq(from = min(ddj_traits_subset$trait_mean), to = max(ddj_traits_subset$trait_mean), length.out = 20)
  
  ########figure 1
  basic <-  ggplot(simresults, aes(x=simtraitlevel, y = Mean)) + 
    theme_bw(base_size=16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_line() + scale_linetype_manual(values = c("solid","longdash")) + #coord_cartesian(ylim = c(-0.2, 1.6)) +
    ylab("DD")  + geom_hline(lty = 2, yintercept = 0) +
    xlab(unique(ddj_traits_subset$trait))
  
  basic_CI <- basic +
    geom_line(data = simresults, aes (x = simtraitlevel, y = X2.5.), lty = "dotted", alpha = 0.7) + #lower bound of 95%CI
    geom_line(data = simresults, aes (x = simtraitlevel, y = X97.5.), lty = "dotted", alpha = 0.7)  #upper bound of 95%CI
  
  basic_CI_points <- basic_CI + geom_point(data = ddj_traits_subset, aes(x = trait_mean, y= DD, color = trait_sp_or_g), size=4,alpha=0.5) +
                     scale_color_manual(values = c("black"))
  
  basic_CI_points_error1 <- basic_CI_points + geom_errorbar(data = ddj_traits_subset,  aes(x = trait_mean, 
                              y = DD, ymin = DD - DD_SD, ymax =  DD + DD_SD), alpha=0.2) 
  
  basic_CI_points_error1 + geom_errorbarh(data = ddj_traits_subset,  
            aes(x = trait_mean, xmin = trait_mean - trait_sd, xmax = trait_mean + trait_sd, y = DD),alpha=0.2) +
            annotate("text", label = paste("model R2 =", modelR2), x = max(ddj_traits_subset$trait_mean), y = 0.5)
  
  figname <- paste(varsubset, ".jpeg", sep = "")
  ggsave(figname, dpi = 300, height = 300, width = 300, units ="mm")
  
  print(varsubset)
  } #end variable loop
  
  
  
  
  
filter(ddj_traits, trait = "leaf_N") %>%
  
summary (lm(ddj_traits$DD[ddj_traits$trait == "leaf_N"] ~ ddj_traits$trait_mean[ddj_traits$trait == "leaf_N"]))
  
  
  