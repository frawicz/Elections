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


election_date <- as.Date("07/10/2018", "%d/%m/%Y")


br_polls <- read_excel("data/pollsbrazil.xltx")%>% 
  filter(Company == "Datafolha" | Company == "Ibope") %>% 
  # gather( key = "party", value = "share", -c("Company", "Date")
  #             ) %>% 
  mutate(days_to_election = election_date - as.Date(Date)) %>% 
  mutate(iid = match(Company , sort(unique(Company))))


party_names <- c("PSL", "PT", "PDT","PSDB", "REDE", "NOVO", "PODE", "MDB", "PATRI", "PSOL", "PSTU", 
                 "DC", "PPL", "Null", "Unknown")

nParties <- length(party_names)

just_shares <- round(as.matrix(br_polls[,party_names] / 100) * 
                       br_polls$Sample_size)


just_shares2 <- just_shares %>% 
  data.frame() %>% 
  mutate(Others = `NOVO` + `PODE` + `MDB` + `PATRI` + `PSOL` + `PSTU` + `DC` + `PPL`) %>% 
  dplyr::select(c("PSL", "PT", "PDT", "PSDB", "REDE", "Others")) %>% 
  as.matrix()


br_polls <- br_polls %>%
  mutate(t= 55 - days_to_election   )

houseprior <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0), nrow = 2)

# 
# forstan2 <- list(
#   # Dynamic Model
#   iid = br_polls$iid,
#   nInst = max(br_polls$iid),
#   nParties = 6,
#   nPolls = nrow(just_shares2),
#   nPeriods =  56, 
#   date = as.numeric(br_polls$t + 1),
#   y = just_shares2,
#   a_pred = as.numeric(c(.40, .20, .1,.1,.1,.1)),
#   prior_house = houseprior
# )
# 
# 
# 
# results2 <- stan(file = "combined_model2.stan", data = forstan2, warmup = 1000,
#                  iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )
# 
# 
# 
# saveRDS(results2, file = paste0("brazil_multinomial_results.RDS"))

# results <- readRDS(file="brazil_multinomial_results.RDS")
# 
# shinystan::launch_shinystan(results)



forstanforecast <- list(
  # Dynamic Model
  iid = br_polls$iid,
  nInst = max(br_polls$iid),
  nParties = 6,
  nPolls = nrow(just_shares2),
  nPeriods =  56, 
  nforecast = 4,
  date = as.numeric(br_polls$t + 1),
  y = just_shares2,
  a_pred = as.numeric(c(.40, .20, .1,.1,.1,.1)),
  prior_house = houseprior
)

resultsforecast <- stan(file = "combined_model3forecast.stan", data = forstanforecast, warmup = 1000,
                        iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )


# stop()

# shinystan::launch_shinystan(resultsforecast)



saveRDS(resultsforecast, file = paste0("brasil_multinomial_results_forecast.RDS"))


# resultsfore <- readRDS(file="brasil_multinomial_results_forecast.RDS")




# Forecast 10 days





br_polls10days <- br_polls %>% 
  filter(t < 43)

just_shares10days <- round(as.matrix(br_polls10days[,party_names] / 100) * 
                             br_polls10days$Sample_size)



just_shares2_10days <- just_shares10days %>% 
  data.frame() %>% 
  mutate(Others = `NOVO` + `PODE` + `MDB` + `PATRI` + `PSOL` + `PSTU` + `DC` + `PPL`) %>% 
  dplyr::select(c("PSL", "PT", "PDT", "PSDB", "REDE", "Others")) %>% 
  as.matrix()



houseprior <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0), nrow = 2)



forstanforecast10days <- list(
  # Dynamic Model
  iid = br_polls10days$iid,
  nInst = max(br_polls10days$iid),
  nParties = 6,
  nPolls = nrow(just_shares2_10days),
  nPeriods =  43, 
  nforecast = 14,
  date = as.numeric(br_polls10days$t + 1),
  y = just_shares2_10days,
  a_pred = as.numeric(c(.40, .20, .1,.1,.1,.1)),
  prior_house = houseprior
)



resultsforecast10days <- stan(file = "combined_model3forecast.stan", data = forstanforecast10days, warmup = 1000,
                              iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )



saveRDS(resultsforecast10days, file = paste0("brazil_multinomial_results_forecast10days.RDS"))



# forecast day 10



br_pollsday10 <- br_polls %>% 
  filter(t < 10)

just_sharesday10 <- round(as.matrix(br_pollsday10[,party_names] / 100) * 
                             br_pollsday10$Sample_size)



just_shares2_day10 <- just_sharesday10 %>% 
  data.frame() %>% 
  mutate(Others = `NOVO` + `PODE` + `MDB` + `PATRI` + `PSOL` + `PSTU` + `DC` + `PPL`) %>% 
  dplyr::select(c("PSL", "PT", "PDT", "PSDB", "REDE", "Others")) %>% 
  as.matrix()



houseprior <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0), nrow = 2)



forstanforecastday10 <- list(
  # Dynamic Model
  iid = br_pollsday10$iid,
  nInst = max(br_pollsday10$iid),
  nParties = 6,
  nPolls = nrow(just_shares2_day10),
  nPeriods =  10, 
  nforecast = 45,
  date = as.numeric(br_pollsday10$t + 1),
  y = just_shares2_day10,
  a_pred = as.numeric(c(.40, .20, .1,.1,.1,.1)),
  prior_house = houseprior
)



resultsforecastday10 <- stan(file = "combined_model3forecast.stan", data = forstanforecastday10, warmup = 1000,
                              iter = niter, chains = nchains, thin = 1, control = list(adapt_delta = 0.80, max_treedepth = 10) )



saveRDS(resultsforecastday10, file = paste0("brazil_multinomial_results_forecastday10.RDS"))

