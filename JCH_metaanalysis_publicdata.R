#preliminary meta-analysis using the publicly available data, pulling missing data from the genus
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



#####################
#formatting data
#####################

mdata$study <-  as.numeric(as.factor(mdata$Author))
mdata$lifestage <- NA
mdata$lifestage[mdata$LifeStage == "Seed"] <- 1
mdata$lifestage[mdata$LifeStage == "Seedling"] <- 2

#is there isn't species level data, take the genus level data
#stem density
mdata$stem_density_mean <- mdata$stem_density_mean_s
for(i in 1:nrow(mdata)){  #should replace for loop with a more efficient way to do this... and make into a function
  if(is.na(mdata$stem_density_mean[i])){mdata$stem_density_mean[i] <- mdata$stem_density_mean_g[i]}
}

mdata$stem_density_sd <- mdata$stem_density_sd_s
for(i in 1:nrow(mdata)){
  if(is.na(mdata$stem_density_sd[i])){mdata$stem_density_sd[i] <- mdata$stem_density_sd_g[i]}
}

mdata$stem_density_sd[mdata$stem_density_sd == 0 & !is.na(mdata$stem_density_sd)] <- 0.001 #to prevent irrational values

#####################
#meta-analysis
#####################

sink("model.txt")
cat("   
    model{
   
    ###### missing data
       for(i in 1:n_obs){
          stem_density_mean[i] ~ dunif(stem_density_mean_min, stem_density_mean_max)
          stem_density_sd[i] ~ dunif(stem_density_sd_min, stem_density_sd_max)
        }

    ###### latent vars
       for(i in 1:n_obs){
          stem_density[i] ~ dnorm(stem_density_mean[i], (1/stem_density_sd[i]^2))
          stem_density_p[i] <- max(stem_density[i], 0.01) #to prevent it from going below zero
        }



    ###### model
    for(i in 1:n_obs){
      OR_mean[i] ~ dnorm(OR[i], 1/OR_var[i])

      OR[i] <- study_fe[study[i]] + lifestage_int[lifestage[i]]  +
                alpha_1 * stem_density_p[i]

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
                     'n_obs' = nrow(mdata),
                     'study' = mdata$study,
                     'nstudies' = max(mdata$study),
                     
                  #Odds ratio
                    'OR_mean' = mdata$OR_mean,
                    'OR_var' = mdata$OR_var,
                  
                  #covariates
                    'lifestage' = mdata$lifestage,
                    'stem_density_mean' = mdata$stem_density_mean,
                    'stem_density_sd' = mdata$stem_density_sd,
                    'stem_density_mean_min' = min(mdata$stem_density_mean, na.rm = TRUE),
                    'stem_density_mean_max' = max(mdata$stem_density_mean, na.rm = TRUE),
                    'stem_density_sd_min' = min(mdata$stem_density_sd, na.rm = TRUE),
                    'stem_density_sd_max' = max(mdata$stem_density_sd, na.rm = TRUE)
                     
                   ),
                   n.chains = 3,
                   n.adapt = 100)

#parameter of interest results
update(jags,n.iter=100) #update(jags,n.iter=1000) 
mcmc_samples <- coda.samples(jags, variable.names=c("alpha_1"),  n.iter=50000)
plot(mcmc_samples)  

#residuals and model fit
#dic.jags <- (dic.samples(jags, 50000, "pD")); dic.jags #linear model: pD = 179.8, quadratic pD = 174.7, exponential = 172.7
mcmc_samples2 <- coda.samples(jags, variable.names=c("OR"),  n.iter=50000)

resultsall<-summary(mcmc_samples2)    
results<-data.frame(resultsall$statistics,resultsall$quantiles) 
results$parameter<-row.names(results)
results$param2<-substr(results$parameter,1,2)
results$param7<-substr(results$parameter,1,7)
results$param12<-substr(results$parameter,1,12)

#predicted versus observed
mdata$predicted <-  results$Mean[results$param2=="OR"]
summary(lm((mdata$predicted ) ~ mdata$OR_mean))
ggplot(mdata, aes(x = predicted, y = OR_mean)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (predicted)") + ylab("OR (measured)")
ggsave("predicted_vs_observed.jpeg")

ggplot(mdata, aes(x = stem_density_mean, y = OR_mean)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (measured)") + ylab("mean stem density") + facet_wrap( ~ LifeStage)

ggplot(mdata, aes(x = predicted, y = stem_density_mean, color = as.factor(study))) + geom_point() + geom_smooth(method = "lm") + theme_bw() + 
  xlab("OR (predicted)") + ylab("OR (measured)") + facet_wrap( ~ lifestage)


names(mdata)

