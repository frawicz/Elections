rm(list = ls())

options(scipen = 9999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
library(readr)
library(stringr)
library(lubridate)
library(gridExtra)



resultsbetamultpt1 <- read_csv("data/resultspt/resultsbetamultpt1.csv")
resultsbetamultpt2 <- read_csv("data/resultspt/resultsbetamultpt2.csv")


resultspt <- rbind(resultsbetamultpt1, resultsbetamultpt2) %>% 
  mutate(Date = as.numeric(str_sub(X1, 6, -4))) %>% 
  mutate(Date = dmy("14/12/2020") + days(Date))



resultsforecastpt <- read_csv("data/resultspt/resultsforecastpt.csv") %>% filter(!str_detect(X1, "alpha"))%>% 
  mutate(Date = as.numeric(str_sub(X1, 10, -4))) %>% 
  mutate(Date = max(resultspt$Date) + days(Date)) %>% 
  filter(!str_detect(X1, "3,") & !str_detect(X1, "4,"))


resultsptfin <- resultspt %>% rbind(resultsforecastpt) #%>% 
  # mutate(Group = ifelse(str_detect(X1, "beta"), 1,2))


PSD <- resultsptfin %>% filter(str_detect(X1, ",1]")) 
PS <- resultsptfin %>% filter(str_detect(X1, ",2]"))
CHEGA <- resultsptfin %>% filter(str_detect(X1, ",3]"))
Others <- resultsptfin %>% filter(str_detect(X1, ",4]"))




p1 <- ggplot(PSD, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "PSD Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(PSD$Date[c(39)]) )+
  geom_point(aes(x= PSD$Date[c(41)], y=0.6066), colour="red", size =3)


p2 <- ggplot(PS, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "PS Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(PS$Date[c(39)]) )+
  geom_point(aes(x= PSD$Date[c(41)], y=0.1296), colour="red", size =3)


p3 <- ggplot(CHEGA, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "CHEGA Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(CHEGA$Date[c(39)]) )+
  geom_point(aes(x= CHEGA$Date[c(41)], y=0.1193), colour="red", size =3)


p4 <- ggplot(Others, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "Others Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(Others$Date[c(39)]) )+
  geom_point(aes(x= Others$Date[c(41)], y=0.1445), colour="red", size =3)


grid.arrange(p1, p2,p3,p4, nrow = 2)



