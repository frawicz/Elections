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


resultsbetamultbr1 <- read_csv("data/resultsbr/resultsbetamultbr1.csv")
resultsbetamultbr2 <- read_csv("data/resultsbr/resultsbetamultbr2.csv")
resultsbetamultbr3 <- read_csv("data/resultsbr/resultsbetamultbr3.csv")
resultsbetamultbr4 <- read_csv("data/resultsbr/resultsbetamultbr4.csv")



resultsbr <- rbind(resultsbetamultbr1, resultsbetamultbr2, resultsbetamultbr3, resultsbetamultbr4)%>% 
  mutate(Date = as.numeric(str_sub(X1, 6, -4)))  %>% 
  mutate(Date = dmy("20/08/2018") + days(Date)) 

resultsforecastbr <- read_csv("data/resultsbr/resultsforecastbr.csv") %>% filter(!str_detect(X1, "alpha"))%>% 
  mutate(Date = as.numeric(str_sub(X1, 10, -4)))%>% 
  mutate(Date = max(resultsbr$Date) + days(Date)) %>% 
  mutate(mean = as.numeric(str_replace(mean, ",", ".")),
         `2.5%` = as.numeric(str_replace(`2.5%`, ",", ".")),
         `97.5%` = as.numeric(str_replace(`97.5%`, ",", "."))
         )%>% 
  filter(!str_detect(X1, "3,") & !str_detect(X1, "4,"))


resultsbrfin <-  rbind(resultsbr, resultsforecastbr) 


PSL <- resultsbrfin %>% filter(str_detect(X1, ",1]"))
PT <- resultsbrfin %>% filter(str_detect(X1, ",2]"))
PDT <- resultsbrfin %>% filter(str_detect(X1, ",3]"))
PSDB <- resultsbrfin %>% filter(str_detect(X1, ",4]"))
REDE <- resultsbrfin %>% filter(str_detect(X1, ",5]"))
Others <- resultsbrfin %>% filter(str_detect(X1, ",6]"))



p1 <- ggplot(PSL, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "PSL Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(PSL$Date[c(56)]) )+
  geom_point(aes(x= PSL$Date[c(58)], y=0.4603), colour="red", size =3)


p2 <- ggplot(PT, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "PT Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(PT$Date[c(56)]) )+
  geom_point(aes(x= PSL$Date[c(58)], y=0.2928), colour="red", size =3)


p3 <- ggplot(PDT, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "PDT Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(PDT$Date[c(56)]) )+
  geom_point(aes(x= PSL$Date[c(58)], y=0.1247), colour="red", size =3)


p4 <- ggplot(PSDB, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "PSDB Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(PSDB$Date[c(56)]) )+
  geom_point(aes(x= PSL$Date[c(58)], y=0.0476), colour="red", size =3)


p5 <- ggplot(REDE, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "REDE Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(REDE$Date[c(56)]) )+
  geom_point(aes(x= PSL$Date[c(58)], y=0.01), colour="red", size =3)


p6 <- ggplot(Others, aes(y = mean, x = Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "grey70") +
  geom_line(aes(y = mean)) +
  labs(x = "Date",y = "%", title = "Others Valid Vote Share in %")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_vline(xintercept = as.numeric(Others$Date[c(56)]) )+
  geom_point(aes(x= PSL$Date[c(58)], y=0.0446), colour="red", size =3)


grid.arrange(p1, p2,p3,p4,p5,p6, nrow = 2)




