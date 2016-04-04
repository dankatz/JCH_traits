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


#is there isn't species level data, take the genus level data
#stem density
mdata$leaf_photosyn_mean <- mdata$leaf_photosyn_mean_s
for(i in 1:nrow(mdata)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata$leaf_photosyn_mean[i])){mdata$leaf_photosyn_mean[i] <- mdata$leaf_photosyn_mean_g[i]}
}

mdata$leaf_photosyn_sd <- mdata$leaf_photosyn_sd_s
for(i in 1:nrow(mdata)){
  if(is.na(mdata$leaf_photosyn_sd[i])){mdata$leaf_photosyn_sd[i] <- mdata$leaf_photosyn_sd_g[i]}
}

mdata$leaf_photosyn_sd[mdata$leaf_photosyn_sd == 0 & !is.na(mdata$leaf_photosyn_sd)] <- 0.001 #to prevent irrational values

#datasubset for model
mdata_subset <- subset(mdata, LifeStage == "Seedling")
mdata_subset$study <-  as.numeric(as.factor(mdata_subset$Author))

######WHY DOES A RANDOMLY GENERATED VARIABLE COME OUT AS SIGNIFICANT A LARGE PORTION OF THE TIME???
#mdata_subset <- subset(mdata_subset, OR_mean > -5) #it isn't driven by this outlier.


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

      OR[i] <- study_re[study[i]] +  #dist_den_int[dist_den[i]] +  
               #alpha_1 * variable_p[i]
                alpha_1 * null_var[i]
    }

    ###### priors
    for(s in 1:nstudies){study_re[s] ~ dnorm(mu, tau)}
    mu ~ dnorm(0, 0.0001)
    tau ~ dgamma(0.001, 0.001)
#     tau  <- 1/(sigma^2)
#     sigma ~ dunif(0,10)

# 
#     dist_den_int[1] ~ dnorm(0, 0.0001)
#     dist_den_int[2] ~ dnorm(0, 0.0001)
    alpha_1 ~ dnorm(0, 0.0001)

    }
    ",fill=TRUE)
sink() 

rand <- runif(nrow(mdata_subset), 0,0.3)
ggplot(mdata_subset, aes(x = rand, y = OR_mean, color = as.factor(study))) + geom_point() #+ geom_smooth(method = "lm")

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
                    'null_var' = rand,
                  
                    #'lifestage' = mdata_subset$lifestage,
                    #'dist_den' = mdata_subset$dist_den,
                    'variable_mean' = mdata_subset$leaf_photosyn_mean,
                    'variable_sd' = mdata_subset$leaf_photosyn_sd,
                    'variable_mean_min' = min(mdata_subset$leaf_photosyn_mean, na.rm = TRUE),
                    'variable_mean_max' = max(mdata_subset$leaf_photosyn_mean, na.rm = TRUE),
                    'variable_sd_min' = min(mdata_subset$leaf_photosyn_sd, na.rm = TRUE),
                    'variable_sd_max' = max(mdata_subset$leaf_photosyn_sd, na.rm = TRUE)
                   ),
                   n.chains = 3,
                   n.adapt = 100)

#parameter of interest results
update(jags,n.iter=100) #update(jags,n.iter=20000) 
mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1"),  n.iter=10000)
plot(mcmc_samples)  


mcmc_samples <- coda.samples(jags, variable.names=c("dist_den_int","alpha_1", "mu","study_re"),  n.iter=10000)
plot(mcmc_samples)  




mcmc_samples <- coda.samples(jags, variable.names=c("variable_p"),  n.iter=10000)
plot(mcmc_samples)  

#residuals and model fit
#dic.jags <- (dic.samples(jags, 50000, "pD")); dic.jags #linear model: pD = 179.8, quadratic pD = 174.7, exponential = 172.7
mcmc_samples2 <- coda.samples(jags, variable.names=c("OR"),  n.iter=10000)

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

ggplot(mdata, aes(x = stem_density_mean, y = OR_mean)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  ylab("OR (measured)") + xlab("mean stem density") + facet_wrap( ~ LifeStage)

ggplot(mdata, aes(x = predicted, y = stem_density_mean, color = as.factor(study))) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (predicted)") + ylab("OR (measured)") + facet_wrap( ~ lifestage)





#########################################################################################
#running through each of the variables sequentially
#########################################################################################
setwd("Q:/Ibanez Lab/Dan Katz/JCH traits/meta analysis/model selection")

##########
#formatting data
##########

mdata$study <-  as.numeric(as.factor(mdata$Author))
mdata$lifestage <- NA
mdata$lifestage[mdata$LifeStage == "Seed"] <- 1
mdata$lifestage[mdata$LifeStage == "Seedling"] <- 2

###########
#traits to loop through
###########

traitlist <- c("leaf_area", "leaf_CNratio", "leaf_lifespan", "leaf_N", "leaf_P", "leaf_photosyn",
               "leaf_SLA", "plant_height", "plant_lifespan", "seed_n", "stem_density")

niterations <- 20000

for(i in 1:length(traitlist)){  #start variable loop
trait <- traitlist[i] # trait <- "leaf_area"

mdata_sub <- select(mdata, contains(trait))  #traits to use in this loop
names(mdata_sub) <- c("trait_mean_s", "trait_sd_s", "trait_mean_g", "trait_sd_g", "trait_n_g")
  #NOTE: this depends on the order of variables in mdata (sp vs. genus), so if that changes, this needs to be too

#if there isn't species level data for that trait, take the genus level data
mdata_sub$trait_mean <- mdata_sub$trait_mean_s
for(i in 1:nrow(mdata_sub)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata_sub$trait_mean[i])){mdata_sub$trait_mean[i] <- mdata_sub$trait_mean_g[i]}
}

mdata_sub$trait_sd <- mdata_sub$trait_sd_s
for(i in 1:nrow(mdata_sub)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata_sub$trait_sd[i])){mdata_sub$trait_sd[i] <- mdata_sub$trait_sd_g[i]}
}

mdata_sub$trait_sd[mdata_sub$trait_sd == 0 & !is.na(mdata_sub$trait_sd)] <- 0.0001 #to prevent irrational values

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
    
    OR[i] <- study_fe[study[i]] + lifestage_int[lifestage[i]]  +
    alpha_1 * trait_p[i]
    
    }
    
    ###### priors
    for(s in 1:nstudies){study_fe[s] ~ dnorm(0, 0.0001)}
    for(l in 1:2){lifestage_int[l] ~ dnorm(0, 0.0001)}
    alpha_1 ~ dnorm(0, 0.0001)
    
    }
    ",fill=TRUE)
sink() 


jags <- jags.model('model.txt',
                   data = list(
                     
                     #indices
                     'n_obs' = nrow(mdata_sub),
                     'study' = mdata$study,
                     'nstudies' = max(mdata$study),
                     
                     #Odds ratio
                     'OR_mean' = mdata$OR_mean,
                     'OR_var' = mdata$OR_var,
                     
                     #covariates
                     'lifestage' = mdata$lifestage,
                     'trait_mean' = mdata_sub$trait_mean,
                     'trait_sd' = mdata_sub$trait_sd,
                     'trait_mean_min' = min(mdata_sub$trait_mean, na.rm = TRUE),
                     'trait_mean_max' = max(mdata_sub$trait_mean, na.rm = TRUE),
                     'trait_sd_min' = min(mdata_sub$trait_sd, na.rm = TRUE),
                     'trait_sd_max' = max(mdata_sub$trait_sd, na.rm = TRUE)
                     
                   ),
                   n.chains = 3,
                   n.adapt = 100)

#parameter of interest results
update(jags,n.iter = niterations) #update(jags,n.iter=1000) 
mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1"),  n.iter= niterations)

plotname <- paste(trait, "_chainhistory.pdf",sep="")
pdf(plotname);plot(mcmc_samples); dev.off()    #saving chain history into a pdf


#residuals and model fit
mcmc_samples2 <- coda.samples(jags, variable.names=c("OR"),  n.iter= niterations)

resultsall<-summary(mcmc_samples2)    
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param2<-substr(results$parameter,1,2)
results$param7<-substr(results$parameter,1,7)
results$param12<-substr(results$parameter,1,12)

#predicted versus observed
mdata$predicted <-  results$Mean[results$param2=="OR"]
model_fit <- summary(lm((mdata$predicted ) ~ mdata$OR_mean))
model_R2 <- paste(trait,  ", model R2 =", round(model_fit$r.squared, 2), sep = "")
ggplot(mdata, aes(x = predicted, y = OR_mean)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (predicted)") + ylab("OR (measured)") + facet_wrap(~LifeStage) + ggtitle(model_R2) 
plotname <- paste(trait, "modelfit.pdf",sep="")
ggsave(plotname)

mdata_sub$OR_mean <- mdata$OR_mean
mdata_sub$study <- mdata$study

ggplot(mdata_sub, aes(x = trait_mean, y = OR_mean, color = as.factor(study))) + geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + theme_bw() + ylab("OR (measured)") + xlab(trait) + 
   theme(legend.position = "none")
plotname <- paste(trait, "trait_x_OR.pdf",sep="")
ggsave(plotname)

print(trait)
} ######end variable loop

