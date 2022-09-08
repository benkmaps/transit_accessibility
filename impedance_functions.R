library(tidyverse)
library(ggplot2)

mgaus = data.frame(t = c(25,50,75),
                   w = 0.5)
mgaus = mgaus %>% mutate(beta = -(t^2)/log(w))

nexp = data.frame(t = c(15,30,45),
                  w = 0.5)
nexp = nexp %>% mutate(beta = -log(w)/t)



times = data.frame(travel_time = 0:120)

times = times %>% 
  mutate(CUMR30 = ifelse(travel_time <= 30, 1, 0)) %>% 
  mutate(CUMR60 = ifelse(travel_time <= 60, 1, 0)) %>% 
  mutate(CUMR90 = ifelse(travel_time <= 90, 1, 0)) %>% 
  
  mutate(CUML30 = (30-travel_time) / 30 * ifelse(travel_time <= 30, 1, 0)) %>% 
  mutate(CUML60 = (60-travel_time) / 60 * ifelse(travel_time <= 60, 1, 0)) %>% 
  mutate(CUML90 = (90-travel_time) / 90 * ifelse(travel_time <= 90, 1, 0)) %>% 
  
  mutate(NEXP_1 = exp(-nexp[1,"beta"] * travel_time)) %>% 
  mutate(NEXP_2 = exp(-nexp[2,"beta"] * travel_time)) %>% 
  mutate(NEXP_3 = exp(-nexp[3,"beta"] * travel_time)) %>% 
  
  mutate(MGAUS_1 = exp(-1 * travel_time^2 / mgaus[1,"beta"])) %>% 
  mutate(MGAUS_2 = exp(-1 * travel_time^2 / mgaus[2,"beta"])) %>% 
  mutate(MGAUS_3 = exp(-1 * travel_time^2 / mgaus[3,"beta"])) %>% 
  
  pivot_longer(-c(travel_time), names_to = "fctn", values_to = "weight") %>% 
  

  
  mutate(parent_fctn = substr(fctn, 1, 4)) %>% 
  mutate(parent_fctn = factor(parent_fctn,
         ordered = T,
         levels = c("CUMR", "NEXP", "MGAU",
                    "CUML"),
         labels = c("CUMR", "NEXP", "MGAUS",
                    "CUML")))

ggplot(times, aes(x = travel_time, y = weight, col = fctn)) + geom_line(size = 2) + 
  labs(col = "Function", x = "Travel Time", y = "Weight") + facet_wrap(~parent_fctn, nrow = 2) +
  scale_color_manual(values = rep(c("#F8766D","#00BA38","#619CFF"),4))
