rm(list = ls())

options(scipen = 9999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(1234)
library(rstan)
library(shinystan)
library(stringr)
library(dplyr)
library(readxl)
library(gtools)
library(rjags)
library(runjags)
library(brms)
library(knitr)
library(DT)
library(tidyverse)
library(bsts)
library(MARSS)
library(reshape2)



nchains <- 3
niter <- 5000


election_date <- as.Date("24/01/2021", "%d/%m/%Y")

pt_polls <- read_excel("data/pollsportugal.xltx", sheet = "2021") %>% 
  mutate(days_to_election = election_date - as.Date(Date))%>% 
  mutate(iid = match(Company , sort(unique(Company)))) %>% 
  mutate("PS" = as.numeric(`PS`),
         "CHEGA" = as.numeric(`CHEGA`),
         "PCP" = as.numeric(`PCP`),
         "PSD" = as.numeric(`PSD`),
         "BE" = as.numeric(`BE`),
         "IL" = as.numeric(`IL`),
         "RIR" = as.numeric(`RIR`),
         "Ns" = as.numeric(`Ns`)
  )

party_names <- c("PS", "CHEGA", "PCP", "PSD", "BE", 
                 "IL", "RIR", "Ns" )


nParties <- length(party_names)

just_shares <- round(as.matrix(pt_polls[,party_names] / 100) * 
                       pt_polls$Sample)



just_shares2 <- just_shares %>% 
  data.frame() %>% 
  mutate(Others = `PCP` + `BE` + `IL` + `RIR`) %>% 
  dplyr::select(c("PSD", "PS", "CHEGA", "Others")) %>% 
  as.matrix()


pt_polls <- pt_polls %>%
  mutate(t= 39 - days_to_election + 1  )



houseprior <- matrix(c(rep(0, times = 20)), nrow = 5)




# forstan2 <- list(
#   # Dynamic Model
#   iid = pt_polls$iid,
#   nInst = max(pt_polls$iid),
#   nParties = 4,
#   nPolls = nrow(just_shares2),
#   nPeriods =  39,
#   date = as.numeric(pt_polls$t + 1),
#   y = just_shares2,
#   a_pred = as.numeric(c(.56, .09, .7,.28)),
#   prior_house = houseprior
# )
# 
# 
# 
# results2 <- stan(file = "combined_model2.stan", data = forstan2, warmup = 1000,
#                  iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )
# 
# 
# shinystan::launch_shinystan(results2)
# 
# 
# 
# saveRDS(results2, file = paste0("portugal_multinomial_results.RDS"))
# 
# 
# results <- readRDS(file="portugal_multinomial_results.RDS")
# 
# 



  forstanforecast <- list(
    # Dynamic Model
    iid = pt_polls$iid,
    nInst = max(pt_polls$iid),
    nParties = 4,
    nPolls = nrow(just_shares2),
    nPeriods =  39,
    nforecast = 4,
    date = as.numeric(pt_polls$t + 1),
    y = just_shares2,
    a_pred = as.numeric(c(.56, .09, .7,.28)),
    prior_house = houseprior
  )
  
  
  
  resultsforecast <- stan(file = "combined_model3forecast.stan", data = forstanforecast, warmup = 1000,
                   iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )

# stop()


  shinystan::launch_shinystan(resultsforecast)
  
  
  
  saveRDS(resultsforecast, file = paste0("portugal_multinomial_results_forecast.RDS"))
  
  
  resultsfore <- readRDS(file="portugal_multinomial_results_forecast.RDS")


  
# Forecast 10 days
  
  pt_polls10days <- pt_polls %>% 
    filter(t < 28)
  
  just_shares10days <- round(as.matrix(pt_polls10days[,party_names] / 100) * 
                               pt_polls10days$Sample)
  
  
  
  just_shares2_10days <- just_shares10days %>% 
    data.frame() %>% 
    mutate(Others = `PCP` + `BE` + `IL` + `RIR`) %>% 
    dplyr::select(c("PSD", "PS", "CHEGA", "Others")) %>% 
    as.matrix()
  
  
  
  houseprior <- matrix(c(rep(0, times = 20)), nrow = 5)
  
  
  
  
  forstanforecast10days <- list(
    # Dynamic Model
    iid = pt_polls10days$iid,
    nInst = max(pt_polls10days$iid),
    nParties = 4,
    nPolls = nrow(just_shares2_10days),
    nPeriods =  26,
    nforecast = 16,
    date = as.numeric(pt_polls10days$t + 1),
    y = just_shares2_10days,
    a_pred = as.numeric(c(.56, .09, .7,.28)),
    prior_house = houseprior
  )
  
  
  
  resultsforecast10days <- stan(file = "combined_model3forecast.stan", data = forstanforecast10days, warmup = 1000,
                          iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )
  
  
  
  saveRDS(resultsforecast10days, file = paste0("portugal_multinomial_results_forecast10days.RDS"))

  
  
  # Forecast day 1  
  
  
  
  
  
  
  pt_pollsday5 <- pt_polls %>% 
    filter(t < 5)
  
  just_sharesday5 <- round(as.matrix(pt_pollsday5[,party_names] / 100) * 
                               pt_pollsday5$Sample)
  
  
  
  just_shares2_day5 <- just_sharesday5 %>% 
    data.frame() %>% 
    mutate(Others = `PCP` + `BE` + `IL` + `RIR`) %>% 
    dplyr::select(c("PSD", "PS", "CHEGA", "Others")) %>% 
    as.matrix()
  
  
  
  houseprior <- matrix(c(rep(0, times = 20)), nrow = 5)
  
  
  
  
  forstanforecastday5 <- list(
    # Dynamic Model
    iid = pt_pollsday5$iid,
    nInst = max(pt_pollsday5$iid),
    nParties = 4,
    nPolls = nrow(just_shares2_day5),
    nPeriods =  5,
    nforecast = 36,
    date = as.numeric(pt_pollsday5$t + 1),
    y = just_shares2_day5,
    a_pred = as.numeric(c(.56, .09, .7,.28)),
    prior_house = houseprior
  )
  
  
  
  resultsforecastday5 <- stan(file = "combined_model3forecast.stan", data = forstanforecastday5, warmup = 1000,
                                iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )
  
  
  
  saveRDS(resultsforecastday5, file = paste0("portugal_multinomial_results_forecastday5.RDS"))
  
  
  
  
  shinystan(aaa)
  
  
  

