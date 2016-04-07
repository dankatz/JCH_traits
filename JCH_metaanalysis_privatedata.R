#preliminary meta-analysis including privately available data, pulling missing data from the genus
#analyzing using only the relevant traits

#set up work environment
library(ggplot2)
library(rjags)
library(dplyr)
library(metafor)
library(tidyr)

rm(list = ls())


#load in data
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/meta analysis")
mdata <- read.csv("mdata.csv")

###################################################################################
#making model for a single variable
###################################################################################


#####################
#formatting data
#####################

mdata$study <-  as.numeric(as.factor(mdata$Author))
mdata$lifestage <- NA
mdata$lifestage[mdata$LifeStage == "Seed"] <- 1
mdata$lifestage[mdata$LifeStage == "Seedling"] <- 2

mdata$dist_den <- NA
mdata$dist_den[mdata$Predictiontested == "distance"] <- 1
mdata$dist_den[mdata$Predictiontested == "density"] <- 2

mdata$OR_sd <- sqrt(mdata$OR_var)
#is there isn't species level data, take the genus level data
#leafP_area
mdata$leafP_area_mean <- mdata$leafP_area_mean_s
for(i in 1:nrow(mdata)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata$leafP_area_mean[i])){mdata$leafP_area_mean[i] <- mdata$leafP_area_mean_g[i]}
}

mdata$leafP_area_sd <- mdata$leafP_area_sd_s
for(i in 1:nrow(mdata)){
  if(is.na(mdata$leafP_area_sd[i])){mdata$leafP_area_sd[i] <- mdata$leafP_area_sd_g[i]}
}

mdata$leafP_area_sd[mdata$leafP_area_sd == 0 & !is.na(mdata$leafP_area_sd)] <- 0.001 #to prevent irrational values

#plant_height
mdata$plant_height_mean <- mdata$plant_height_mean_s
for(i in 1:nrow(mdata)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata$plant_height_mean[i])){mdata$plant_height_mean[i] <- mdata$plant_height_mean_g[i]}
}

mdata$plant_height_sd <- mdata$plant_height_sd_s
for(i in 1:nrow(mdata)){
  if(is.na(mdata$plant_height_sd[i])){mdata$plant_height_sd[i] <- mdata$plant_height_sd_g[i]}
}

mdata$plant_height_sd[mdata$plant_height_sd == 0 & !is.na(mdata$plant_height_sd)] <- 0.001 #to prevent irrational values


#datasubset for model
mdata_subset <- subset(mdata, LifeStage == "Seedling")
mdata_subset$study <-  as.numeric(as.factor(as.character(mdata_subset$Author)))

#####################
#meta-analysis
#####################

sink("model.txt")
cat("   
    model{
   
    ###### missing data
       for(i in 1:n_obs){
          variable1_mean[i] ~ dunif(variable1_mean_min, variable1_mean_max)
          variable1_sd[i] ~ dunif(variable1_sd_min, variable1_sd_max)

          variable2_mean[i] ~ dunif(variable2_mean_min, variable2_mean_max)
          variable2_sd[i] ~ dunif(variable2_sd_min, variable2_sd_max)
      }

    ###### latent vars
       for(i in 1:n_obs){
          variable1[i] ~ dnorm(variable1_mean[i], (1/variable1_sd[i]^2))
          variable1_p[i] <- max(variable1[i], 0.001) #to prevent it from going below zero

          variable2[i] ~ dnorm(variable2_mean[i], (1/variable2_sd[i]^2))
          variable2_p[i] <- max(variable2[i], 0.001) #to prevent it from going below zero
        }

    ###### model
    for(i in 1:n_obs){
      OR_mean[i] ~ dnorm(OR[i], 1/OR_var[i])

      OR[i] <- study_re[study[i]] +  #dist_den_int[dist_den[i]] +  #alpha_1 * variable1_p[i]
               alpha_1 * variable1_p[i] + alpha_2  * variable2_p[i] +
               beta_1 * percentdead[i] +
               residuals[i] 

      OR_pred[i] <- study_re[study[i]] +  alpha_1 * variable1_p[i] + alpha_2  * variable2_p[i] +
               beta_1 * percentdead[i] 
    }

    ###### priors
    for(s in 1:nstudies){study_re[s] ~ dnorm(mu, tau)}
    mu ~ dnorm(0, 0.0001)
    tau ~ dgamma(0.001, 0.001)

    for(r in 1:n_obs){residuals[r] ~ dnorm(0, tau_residuals)}
    tau_residuals  <- 1/(sigma_residuals^2)
    sigma_residuals ~ dunif(0,10)

    alpha_1 ~ dnorm(0, 0.0001)
    alpha_2 ~ dnorm(0, 0.0001)
    beta_1 ~ dnorm(0, 0.0001)

    ########## simulation for figures
       
  #######predictions for alpha_1
  sim_alpha_1[1] <- mu + alpha_1 * 0.00 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[2] <- mu + alpha_1 * 0.05 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[3] <- mu + alpha_1 * 0.10 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[4] <- mu + alpha_1 * 0.15 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[5] <- mu + alpha_1 * 0.20 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[6] <- mu + alpha_1 * 0.25 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[7] <- mu + alpha_1 * 0.30 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean
  sim_alpha_1[8] <- mu + alpha_1 * 0.35 + alpha_2  * variable2_mean_overall +  beta_1 * percentdead_overall_mean

  #######predictions for alpha_2
  sim_alpha_2[1] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 0 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[2] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 6 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[3] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 12 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[4] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 18 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[5] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 24 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[6] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 30 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[7] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 36 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[8] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 42 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[9] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 48 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[10] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 54 +  beta_1 * percentdead_overall_mean
  sim_alpha_2[11] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * 61 +  beta_1 * percentdead_overall_mean

  #######predictions for beta
  sim_beta_1[1] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.0
  sim_beta_1[2] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.1
  sim_beta_1[3] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.2
  sim_beta_1[4] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.3
  sim_beta_1[5] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.4
  sim_beta_1[6] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.5
  sim_beta_1[7] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.6
  sim_beta_1[8] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.7
  sim_beta_1[9] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.8
  sim_beta_1[10] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 0.9
  sim_beta_1[11] <- mu + alpha_1 * variable1_mean_overall + alpha_2  * variable2_mean_overall +  beta_1 * 1.0
}
    ",fill=TRUE)
sink() 

jags <- jags.model('model.txt',
                   data = list(
                     
                  #indices
                     'n_obs' = nrow(mdata_subset),
                     'study' = mdata_subset$study,
                     'nstudies' = max(mdata_subset$study),
                     
                  #Odds ratio
                    'OR_mean' = mdata_subset$OR_mean,
                    'OR_var' = mdata_subset$OR_var,
                  
                  #covariates 
                    #'null_var' = mdata_subset$rand,  #mdata_subset$rand <- runif(nrow(mdata_subset), -2,2)
                    #'lifestage' = mdata_subset$lifestage,
                    #'dist_den' = mdata_subset$dist_den,
#                     'percentdead' = mdata_subset$Nosurvivorsnear_high/(mdata_subset$No_died_near_high +
#                                                                        mdata_subset$Nosurvivorsnear_high), 
                    'percentdead' = mdata_subset$No_survivors_far_low/(mdata_subset$No_died_far_low+
                                                                       mdata_subset$No_survivors_far_low),
                    'percentdead_overall_mean' = mean(mdata_subset$No_survivors_far_low/
                                              (mdata_subset$No_died_far_low + mdata_subset$No_survivors_far_low)),

                    'variable1_mean' = mdata_subset$leafP_area_mean,
                    'variable1_sd' = mdata_subset$leafP_area_sd,
                    'variable1_mean_min' = min(mdata_subset$leafP_area_mean, na.rm = TRUE),
                    'variable1_mean_max' = max(mdata_subset$leafP_area_mean, na.rm = TRUE),
                    'variable1_sd_min' = min(mdata_subset$leafP_area_sd, na.rm = TRUE),
                    'variable1_sd_max' = max(mdata_subset$leafP_area_sd, na.rm = TRUE),
                    'variable1_mean_overall' = mean(mdata_subset$leafP_area_mean, na.rm = TRUE),

                    'variable2_mean' = mdata_subset$plant_height_mean,
                    'variable2_sd' = mdata_subset$plant_height_sd,
                    'variable2_mean_min' = min(mdata_subset$plant_height_mean, na.rm = TRUE),
                    'variable2_mean_max' = max(mdata_subset$plant_height_mean, na.rm = TRUE),
                    'variable2_sd_min' = min(mdata_subset$plant_height_sd, na.rm = TRUE),
                    'variable2_sd_max' = max(mdata_subset$plant_height_sd, na.rm = TRUE),
                    'variable2_mean_overall' = mean(mdata_subset$plant_height_mean, na.rm = TRUE)

                   ),
                   n.chains = 3,
                   n.adapt = 100)

#parameter of interest results
update(jags,n.iter=100) #update(jags,n.iter=20000) 
mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1", "alpha_2","beta_1"),  n.iter=10000)
plot(mcmc_samples)  
resultsall<-summary(mcmc_samples)    
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param2<-substr(results$parameter,1,2)

ggplot(results, aes(x = parameter, y = Mean)) + facet_wrap(~parameter, scales = "free") +
  geom_point(size=4) +   geom_errorbar(aes(ymax=X97.5., ymin = X2.5. ), width=0.15) +
  theme_bw(base_size=16) +  geom_hline(aes(x=0), lty=2) + xlab("") + ylab("coefficient estimate (mean + 95% CI)")

##############
#fig 1: simulated results: alpha_1 
##############
mcmc_samples <- coda.samples(jags, variable.names=c("sim_alpha_1"),  n.iter=10000)
simresultsall <- summary(mcmc_samples)    
simresults <- data.frame(simresultsall$statistics,simresultsall$quantiles) 
simresults$parameter<-row.names(simresults)
simresults$param2<-substr(simresults$parameter,1,2)

#extracting traitlevel
simresults$traitlevel2 <- as.numeric(substr(simresults$parameter, 13,13))
simresults$traitlevel2[simresults$traitlevel2 == 1] <- 0
simresults$traitlevel2[simresults$traitlevel2 == 2] <- 0.05
simresults$traitlevel2[simresults$traitlevel2 == 3] <- 0.10
simresults$traitlevel2[simresults$traitlevel2 == 4] <- 0.15
simresults$traitlevel2[simresults$traitlevel2 == 5] <- 0.20
simresults$traitlevel2[simresults$traitlevel2 == 6] <- 0.25
simresults$traitlevel2[simresults$traitlevel2 == 7] <- 0.30
simresults$traitlevel2[simresults$traitlevel2 == 8] <- 0.35

########figure 1
basic <-  ggplot(simresults, aes(x=traitlevel2, y = Mean)) + 
  theme_bw(base_size=16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_line() + scale_linetype_manual(values = c("solid","longdash")) + #coord_cartesian(ylim = c(-0.2, 1.6)) +
  ylab("lnOR")  + geom_hline(lty = 2) +
  xlab("leaf P / leaf area")

basic_CI <- basic +
  geom_line(data = simresults, aes (x = traitlevel2, y = X2.5.), lty = "dotted", alpha = 0.7) + #lower bound of 95%CI
  geom_line(data = simresults, aes (x = traitlevel2, y = X97.5.), lty = "dotted", alpha = 0.7)  #upper bound of 95%CI

basic_CI_points <- basic_CI + geom_point(data = mdata_subset, 
                      aes(x = leafP_area_mean, y= OR_mean, color = as.factor(study)), size=4,alpha=0.5) 

basic_CI_points_error1 <- basic_CI_points + geom_errorbar(data = mdata_subset,  aes(x = leafP_area_mean, 
       y = OR_mean, ymin = OR_mean - sqrt(OR_var), ymax =  OR_mean + sqrt(OR_var), 
       color = as.factor(study)),alpha=0.2) 

basic_CI_points_error1 + geom_errorbarh(data = mdata_subset,  aes(x = leafP_area_mean, 
        y = OR_mean, xmin = leafP_area_mean - leafP_area_sd, xmax = leafP_area_mean + leafP_area_sd, 
        color = as.factor(study)),alpha=0.2) 
ggsave("figure2.jpg", dpi = 600, height = 200, width = 300, units ="mm")



##############
#fig 2: simulated results: alpha_2 
##############
mcmc_samples <- coda.samples(jags, variable.names=c("sim_alpha_2"),  n.iter=10000)
simresultsall <- summary(mcmc_samples)    
simresults <- data.frame(simresultsall$statistics,simresultsall$quantiles) 
simresults$parameter<-row.names(simresults)
simresults$param2<-substr(simresults$parameter,1,2)

#extracting traitlevel
simresults$traitlevel2 <- substr(simresults$parameter, 13,14)
simresults$traitlevel2 <- as.numeric(gsub("[^0-9]", "", simresults$traitlevel2))
simresults$traitlevel3 <- NA
simresults$traitlevel3[simresults$traitlevel2 == 1] <- 0
simresults$traitlevel3[simresults$traitlevel2 == 2] <- 6
simresults$traitlevel3[simresults$traitlevel2 == 3] <- 12
simresults$traitlevel3[simresults$traitlevel2 == 4] <- 18
simresults$traitlevel3[simresults$traitlevel2 == 5] <- 24
simresults$traitlevel3[simresults$traitlevel2 == 6] <- 30
simresults$traitlevel3[simresults$traitlevel2 == 7] <- 36
simresults$traitlevel3[simresults$traitlevel2 == 8] <- 42
simresults$traitlevel3[simresults$traitlevel2 == 9] <- 48
simresults$traitlevel3[simresults$traitlevel2 == 10] <- 54
simresults$traitlevel3[simresults$traitlevel2 == 11] <- 61

########figure 2
basic <-  ggplot(simresults, aes(x=traitlevel3, y = Mean)) + 
  theme_bw(base_size=16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_line() + scale_linetype_manual(values = c("solid","longdash")) + coord_cartesian(xlim = c(-1, 65)) +
  ylab("lnOR")  + geom_hline(lty = 2) +
  xlab("height (m)")

basic_CI <- basic +
  geom_line(data = simresults, aes (x = traitlevel3, y = X2.5.), lty = "dotted", alpha = 0.7) + #lower bound of 95%CI
  geom_line(data = simresults, aes (x = traitlevel3, y = X97.5.), lty = "dotted", alpha = 0.7)  #upper bound of 95%CI

basic_CI_points <- basic_CI + geom_point(data = mdata_subset, 
                                         aes(x = plant_height_mean, y= OR_mean, color = as.factor(study)), size=4,alpha=0.5) 

basic_CI_points_error1 <- basic_CI_points + geom_errorbar(data = mdata_subset,  aes(x = plant_height_mean, 
                          y = OR_mean, ymin = OR_mean - sqrt(OR_var), ymax =  OR_mean + sqrt(OR_var), 
                          color = as.factor(study)),alpha=0.2) 

basic_CI_points_error1 + geom_errorbarh(data = mdata_subset,  aes(x = plant_height_mean, 
          y = OR_mean, xmin = plant_height_mean - plant_height_sd, xmax = plant_height_mean + plant_height_sd, 
          color = as.factor(study)),alpha=0.2) 
ggsave("figure2.jpg", dpi = 600, height = 200, width = 300, units ="mm")



##############
#fig 3: simulated results: beta_1 
##############
mcmc_samples <- coda.samples(jags, variable.names=c("sim_beta_1"),  n.iter=10000)
simresultsall <- summary(mcmc_samples)    
simresults <- data.frame(simresultsall$statistics,simresultsall$quantiles) 
simresults$parameter<-row.names(simresults)
simresults$param2<-substr(simresults$parameter,1,2)

#extracting traitlevel
simresults$traitlevel2 <- substr(simresults$parameter, 12,13)
simresults$traitlevel2 <- as.numeric(gsub("[^0-9]", "", simresults$traitlevel2))
simresults$traitlevel3 <- NA
simresults$traitlevel3[simresults$traitlevel2 == 1] <- 0
simresults$traitlevel3[simresults$traitlevel2 == 2] <- .1
simresults$traitlevel3[simresults$traitlevel2 == 3] <- .2
simresults$traitlevel3[simresults$traitlevel2 == 4] <- .3
simresults$traitlevel3[simresults$traitlevel2 == 5] <- .4
simresults$traitlevel3[simresults$traitlevel2 == 6] <- .5
simresults$traitlevel3[simresults$traitlevel2 == 7] <- .6
simresults$traitlevel3[simresults$traitlevel2 == 8] <- .7
simresults$traitlevel3[simresults$traitlevel2 == 9] <- .8
simresults$traitlevel3[simresults$traitlevel2 == 10] <- .9
simresults$traitlevel3[simresults$traitlevel2 == 11] <- 1

########figure 2
basic <-  ggplot(simresults, aes(x=traitlevel3, y = Mean)) + 
  theme_bw(base_size=16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_line() + scale_linetype_manual(values = c("solid","longdash")) + 
  ylab("lnOR")  + geom_hline(lty = 2) +
  xlab("background survival (proportion)")

basic_CI <- basic +
  geom_line(data = simresults, aes (x = traitlevel3, y = X2.5.), lty = "dotted", alpha = 0.7) + #lower bound of 95%CI
  geom_line(data = simresults, aes (x = traitlevel3, y = X97.5.), lty = "dotted", alpha = 0.7)  #upper bound of 95%CI

basic_CI_points <- basic_CI + geom_point(data = mdata_subset, 
     aes(x = mdata_subset$No_survivors_far_low/(mdata_subset$No_died_far_low+ mdata_subset$No_survivors_far_low),
         y= OR_mean, color = as.factor(study)), size=4,alpha=0.5) 

basic_CI_points + geom_errorbar(data = mdata_subset,  
     aes(x =  mdata_subset$No_survivors_far_low/(mdata_subset$No_died_far_low+ mdata_subset$No_survivors_far_low),
      y = OR_mean, ymin = OR_mean - sqrt(OR_var), ymax =  OR_mean + sqrt(OR_var), 
      color = as.factor(study)),alpha=0.2) 

ggsave("figure3.jpg", dpi = 600, height = 200, width = 300, units ="mm")


############
#residuals and model fit
############
#dic.jags <- (dic.samples(jags, 50000, "pD")); dic.jags #linear model: pD = 179.8, quadratic pD = 174.7, exponential = 172.7
mcmc_samples2 <- coda.samples(jags, variable.names=c("OR_pred"),  n.iter=10000)

resultsall<-summary(mcmc_samples2)    
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param2<-substr(results$parameter,1,2)
results$param7<-substr(results$parameter,1,7)
results$param12<-substr(results$parameter,1,12)

#predicted versus observed
mdata_subset$predicted <-  results$Mean[results$param2=="OR"]
summary(lm((mdata_subset$predicted ) ~ mdata_subset$OR_mean))
ggplot(mdata_subset, aes(x = predicted, y = OR_mean)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (predicted)") + ylab("OR (measured)")
ggsave("predicted_vs_observed.jpeg")





#########################################################################################
#running through each of the variables sequentially
#########################################################################################
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/meta analysis/model selection")

##########
#selecting and formatting data
##########

mdata_subset <- filter(mdata, LifeStage == "Seedling")
mdata_subset$study <-  as.numeric(as.factor(as.character(mdata_subset$Author)))

###########
#traits to loop through
###########

traitlist <- c("leaf_area", "leaf_CNratio", "leaf_lifespan", "leaf_N", "leaf_P","leafP_area", 
               "leaf_SLA", "plant_height", "plant_lifespan", "seed_n", "stem_density") #

niterations <- 50000

for(i in 1:length(traitlist)){  #start variable loop
trait <- traitlist[i] # trait <- "leafP_area"


mdata_trait <- select(mdata_subset, contains(trait))  #traits to use in this loop
names(mdata_trait) <- c("trait_mean_s", "trait_sd_s", "trait_mean_g", "trait_sd_g", "trait_n_g")
  #NOTE: this depends on the order of variables in mdata (sp vs. genus), so if that changes, this needs to be too

#if there isn't species level data for that trait, take the genus level data
mdata_trait$trait_mean <- mdata_trait$trait_mean_s
for(i in 1:nrow(mdata_trait)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata_trait$trait_mean[i])){mdata_trait$trait_mean[i] <- mdata_trait$trait_mean_g[i]}
}

mdata_trait$trait_sd <- mdata_trait$trait_sd_s
for(i in 1:nrow(mdata_trait)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata_trait$trait_sd[i])){mdata_trait$trait_sd[i] <- mdata_trait$trait_sd_g[i]}
}

mdata_trait$trait_sd[mdata_trait$trait_sd == 0 & !is.na(mdata_trait$trait_sd)] <- 0.0001 #to prevent irrational values

#####################
#meta-analysis
#####################

sink("model.txt")
cat("   
    model{
    
    ###### missing data
    for(i in 1:n_obs){
    trait_mean[i] ~ dunif(trait_mean_min, trait_mean_max)
    trait_sd[i] ~ dunif(trait_sd_min, trait_sd_max)
    }
    
    ###### latent vars
    for(i in 1:n_obs){
    trait[i] ~ dnorm(trait_mean[i], (1/trait_sd[i]^2))
    trait_p[i] <- max(trait[i], 0.01) #to prevent it from going below zero
    }
 
    ###### model
    for(i in 1:n_obs){
    OR_mean[i] ~ dnorm(OR[i], 1/OR_var[i])
    
    OR[i] <- study_re[study[i]] + 
             alpha_1 * trait_p[i] + 
             residuals[i]

    OR_pred[i] <- study_re[study[i]] + alpha_1 * trait_p[i]
    }
    
    ###### priors
    for(s in 1:nstudies){study_re[s] ~ dnorm(mu, tau)}
    mu ~ dnorm(0, 0.0001)
    tau ~ dgamma(0.001, 0.001)
    
    alpha_1 ~ dnorm(0, 0.0001)

    for(r in 1:n_obs){residuals[r] ~ dnorm(0, tau_residuals)}
    tau_residuals  <- 1/(sigma_residuals^2)
    sigma_residuals ~ dunif(0,10)
    }
    ",fill=TRUE)
sink() 

jags <- jags.model('model.txt',
                   data = list(
                     
                     #indices
                     'n_obs' = nrow(mdata_subset),
                     'study' = mdata_subset$study,
                     'nstudies' = max(mdata_subset$study),
                     
                     #Odds ratio
                     'OR_mean' = mdata_subset$OR_mean,
                     'OR_var' = mdata_subset$OR_var,
                     
                     #covariates
                     'trait_mean' = mdata_trait$trait_mean,
                     'trait_sd' = mdata_trait$trait_sd,
                     'trait_mean_min' = min(mdata_trait$trait_mean, na.rm = TRUE),
                     'trait_mean_max' = max(mdata_trait$trait_mean, na.rm = TRUE),
                     'trait_sd_min' = min(mdata_trait$trait_sd, na.rm = TRUE),
                     'trait_sd_max' = max(mdata_trait$trait_sd, na.rm = TRUE)
                   ),
                   n.chains = 3,
                   n.adapt = 100)

#parameter of interest results
update(jags,n.iter = niterations) #update(jags,n.iter=1000) 
mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1","mu"),  n.iter= niterations)
plotname <- paste(trait, "_chainhistory.pdf",sep="")
pdf(plotname);plot(mcmc_samples); dev.off()    #saving chain history into a pdf

resultsall <- summary(mcmc_samples)    
results <- data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter <- row.names(results)
results$param2 <- substr(results$parameter,1,2)
results_file_name <- paste(trait, "_results.csv",sep="")
write.csv(results, results_file_name)


#residuals and model fit
mcmc_samples2 <- coda.samples(jags, variable.names=c("OR_pred"),  n.iter= niterations)

resultsall<-summary(mcmc_samples2)    
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param2<-substr(results$parameter,1,2)
results$param7<-substr(results$parameter,1,7)
results$param12<-substr(results$parameter,1,12)

#predicted versus observed
mdata_subset$predicted <-  results$Mean[results$param2=="OR"]
model_fit <- summary(lm((mdata_subset$predicted ) ~ mdata_subset$OR_mean))
model_R2 <- paste(trait,  ", model R2 =", round(model_fit$r.squared, 2), sep = "")
ggplot(mdata_subset, aes(x = predicted, y = OR_mean)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (predicted)") + ylab("OR (measured)") + facet_wrap(~LifeStage) + ggtitle(model_R2) 
plotname <- paste(trait, "modelfit.pdf",sep="")
ggsave(plotname)

#mdata_sub$OR_mean <- mdata$OR_mean
#mdata_sub$study <- mdata$study

print(trait)
} ######end variable loop




######################################################
#looking at model results
######################################################
library(stringr)
rm(list=setdiff(ls(), c("mdata", "mdata_subset")))  #removing everthing besides survexp

traitlist <- c("leaf_area", "leaf_CNratio", "leaf_lifespan", "leaf_N", "leaf_P", "leafP_area", 
               "leaf_SLA", "plant_height", "plant_lifespan", "seed_n", "stem_density") 
filenames <- paste(traitlist, "_results.csv", sep = "")
results <- do.call("rbind", sapply(filenames, read.csv, simplify = FALSE)) 

#have to add in species and variable based on file name
results$trait <- row.names(results)
results$trait <- sub('\\..*', '', results$trait)
results$trait <- sub("_results", "", results$trait)

filter(results, parameter == "alpha_1") %>%
ggplot(aes(x = trait, y = Mean)) + 
  geom_point(size=4) +            
  geom_errorbar(aes(ymax=X97.5., ymin = X2.5. ), width=0.15) +
  theme_bw(base_size=16) + facet_wrap( ~ trait, scales = "free") +
  geom_hline(aes(x=0), lty=2) + xlab("") + ylab("coefficient estimate (mean + 95% CI)")
ggsave("results_p_leafsucker_detect_othervars.jpeg",dpi = 600, width = 12, height = 12, units = "in")



