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
#stem density
mdata$leaf_SLA_mean <- mdata$leaf_SLA_mean_s
for(i in 1:nrow(mdata)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata$leaf_SLA_mean[i])){mdata$leaf_SLA_mean[i] <- mdata$leaf_SLA_mean_g[i]}
}

mdata$leaf_SLA_sd <- mdata$leaf_SLA_sd_s
for(i in 1:nrow(mdata)){
  if(is.na(mdata$leaf_SLA_sd[i])){mdata$leaf_SLA_sd[i] <- mdata$leaf_SLA_sd_g[i]}
}

mdata$leaf_SLA_sd[mdata$leaf_SLA_sd == 0 & !is.na(mdata$leaf_SLA_sd)] <- 0.001 #to prevent irrational values

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
          variable_mean[i] ~ dunif(variable_mean_min, variable_mean_max)
          variable_sd[i] ~ dunif(variable_sd_min, variable_sd_max)
        }

    ###### latent vars
       for(i in 1:n_obs){
          variable[i] ~ dnorm(variable_mean[i], (1/variable_sd[i]^2))
          variable_p[i] <- max(variable[i], 0.001) #to prevent it from going below zero
        }

    ###### model
    for(i in 1:n_obs){
      OR_mean[i] ~ dnorm(OR[i], 1/OR_var[i])

      OR[i] <- study_re[study[i]] +  #dist_den_int[dist_den[i]] +  #alpha_1 * variable_p[i]
               alpha_1 * variable_p[i] + residuals[i]

      OR_pred[i] <- study_re[study[i]] + alpha_1 * variable_p[i]
    }

    ###### priors
    for(s in 1:nstudies){study_re[s] ~ dnorm(mu, tau)}
    mu ~ dnorm(0, 0.0001)
    tau ~ dgamma(0.001, 0.001)
#   tau  <- 1/(sigma^2)
#   sigma ~ dunif(0,10)

    for(r in 1:n_obs){residuals[r] ~ dnorm(0, tau_residuals)}
    tau_residuals  <- 1/(sigma_residuals^2)
    sigma_residuals ~ dunif(0,10)

# 
#     dist_den_int[1] ~ dnorm(0, 0.0001)
#     dist_den_int[2] ~ dnorm(0, 0.0001)
    alpha_1 ~ dnorm(0, 0.0001)

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
                    'variable_mean' = mdata_subset$leaf_SLA_mean,
                    'variable_sd' = mdata_subset$leaf_SLA_sd,
                    'variable_mean_min' = min(mdata_subset$leaf_SLA_mean, na.rm = TRUE),
                    'variable_mean_max' = max(mdata_subset$leaf_SLA_mean, na.rm = TRUE),
                    'variable_sd_min' = min(mdata_subset$leaf_SLA_sd, na.rm = TRUE),
                    'variable_sd_max' = max(mdata_subset$leaf_SLA_sd, na.rm = TRUE)
                   ),
                   n.chains = 3,
                   n.adapt = 100)

#parameter of interest results
update(jags,n.iter=100) #update(jags,n.iter=20000) 
mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1"),  n.iter=10000)
plot(mcmc_samples)  


mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1", "mu","study_re"),  n.iter=10000)
plot(mcmc_samples)  




mcmc_samples <- coda.samples(jags, variable.names=c("variable_p"),  n.iter=10000)
plot(mcmc_samples)  

#residuals and model fit
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



