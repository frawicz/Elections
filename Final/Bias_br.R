
### ----------------------------------------------------------
### Election polling forecasting 
### Brazilian election 2018
### 
### Fernando Rawicz
### 
### 1. Estimate Brazilian Bias
###
### ----------------------------------------------------------
rm(list = ls())

options(scipen = 9999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(shinystan)
library(stringr)
library(dplyr)
library(readxl)
library(gtools)
library(rjags)
library(runjags)
library(brms)




candidates <- c("PSL", "PT", "PDT", "PSC", "PSDB", "PSOL", "PMDB", "NOVO", "PPL", "REDE", "PSTU")

br_datafolha <- read_excel("data/pollsbrazil.xltx", sheet = "Datafolha")%>% 
  filter(Party %in% candidates)%>% 
  mutate(Partyiid = match(Party , sort(unique(Party)))) %>% 
  mutate(Bias = Poll - Result) %>% 
  mutate(Company = "Datafolha")


br_ibope <- read_excel("data/pollsbrazil.xltx", sheet = "Ibope")%>% 
  filter(Party %in% candidates)%>% 
  mutate(Partyiid = match(Party , sort(unique(Party)))) %>% 
  mutate(Bias = Poll - Result) %>% 
  mutate(Company = "Ibope")



br_bias_final <- rbind(br_datafolha, br_ibope) %>% 
  mutate(Companyiid = match(Company , sort(unique(Company))))  %>% 
  mutate(Eleciid = match(Election , sort(unique(Election)))) 





boxplot(Bias ~ Company, data=br_bias_final, main="Boxplot Bias all elections Datafolha and Ibope")


boxplot(Bias ~ Company, data=br_bias_final %>% filter(Election == "Pres"), main="Boxplot Bias President Datafolha and Ibope")



#################################################################################
set.seed(115)

mod_string = " 
model {
## sampling
for (i in 1:N){
   y[i] ~ dnorm(mu_party[PartyIndex[i]] + mu_house[HouseIndex[i]] , invsigma2)
}

## priors
for (j in 1:J){
   mu_party[j] ~ dnorm(mu, invtau2)
}

for (h in 1:H){
   mu_house[h] ~ dnorm(mu, invtau2)
}


invsigma2 ~ dgamma(3, 1)
sigma <- sqrt(pow(invsigma2, -1))

## hyperpriors

mu ~ dnorm(0, 4)
invtau2 ~ dgamma(3, 1)
tau <- sqrt(pow(invtau2, -1))

}"





params = c("mu_party", "mu_house", "mu", "invtau2", "invsigma2")



the_data <- list("y" = br_bias_final$Bias,
                 "PartyIndex" = br_bias_final$Partyiid,
                 "HouseIndex" = br_bias_final$Companyiid,
                 "N" = length(br_bias_final$Bias),
                 "J" = length(unique(br_bias_final$Partyiid)),
                 "H" = length(unique(br_bias_final$Companyiid))
                 )


mod = jags.model(textConnection(mod_string), data=the_data, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)

mod_csim = as.mcmc(do.call(rbind, mod_sim))


plot(mod_sim, ask = T)
traceplot(mod_sim, col = 1:3)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)
dic = dic.samples(mod, n.iter=1e3)


(pm_params = colMeans(mod_csim))



summary(mod_sim)


###################################################################################################################################

mod_string2 = " 
model {
## sampling
for (i in 1:N){
   y[i] ~ dnorm(mu_party[PartyIndex[i]] + mu_house[HouseIndex[i]] + mu_elec[ElecIndex[i]] , invsigma2)
}

## priors
for (j in 1:J){
   mu_party[j] ~ dnorm(mu, invtau2)
}

for (h in 1:H){
   mu_house[h] ~ dnorm(mu, invtau2)
}

for (r in 1:R){
   mu_elec[r] ~ dnorm(mu, invtau2)
}


invsigma2 ~ dgamma(3, 1)
sigma <- sqrt(pow(invsigma2, -1))

## hyperpriors

mu ~ dnorm(0, 4)
invtau2 ~ dgamma(3, 1)
tau <- sqrt(pow(invtau2, -1))

}"



set.seed(115)


params = c("mu_party", "mu_house","mu_elec", "mu", "invtau2", "invsigma2")



the_data2 <- list("y" = br_bias_final$Bias,
                 "PartyIndex" = br_bias_final$Partyiid,
                 "HouseIndex" = br_bias_final$Companyiid,
                 "ElecIndex" = br_bias_final$Eleciid,
                 "N" = length(br_bias_final$Bias),
                 "J" = length(unique(br_bias_final$Partyiid)),
                 "H" = length(unique(br_bias_final$Companyiid)),
                 "R" = length(unique(br_bias_final$Eleciid))
)


mod2 = jags.model(textConnection(mod_string2), data=the_data2, n.chains=3)
update(mod2, 1e3)

mod_sim2 = coda.samples(model=mod2,
                       variable.names=params,
                       n.iter=5e3)

mod_csim2 = as.mcmc(do.call(rbind, mod_sim2))


plot(mod_sim2, ask = T)
traceplot(mod_sim2, col = 1:3)

gelman.diag(mod_sim2)
autocorr.diag(mod_sim2)
autocorr.plot(mod_sim2)
effectiveSize(mod_sim2)
dic = dic.samples(mod2, n.iter=1e3)


(pm_params = colMeans(mod_csim2))



summary(mod_sim2)




