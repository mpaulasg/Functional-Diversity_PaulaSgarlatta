############################################################################################################
##
## Script for plotting taxonomic and functional diversity across years using Hill numbers
## 
## Code by Paula Sgarlatta, Sebastien Villeger and Camille Magneville 
##
###########################################################################################################
rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

# loading data

load(here::here("outputs/", "FD_beta_kelp_nokelp.RData") )
load(here::here("outputs/", "TD_beta_kelp_nokelp.RData") )


load(here::here("outputs/", "TD_kelp_nokelp_Hill.RData") )
load(here::here("outputs/", "FD_kelp_nokelp_Hill.RData") )


## Taxonomic dissimilarity ####

# Obtain average per year

TD_beta_k_nok <- TD_beta_kelp_nokelp %>%
  filter(Year1 == Year2, Habitat1==Habitat2) %>% 
  select(Year=Year1, q1) %>% 
  group_by(Year) %>%
  summarise( n = n(),
             q1_mean = mean(q1),
             q1_sd = sd(q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))

TD_beta_k_nok$Year <- as.numeric(TD_beta_k_nok$Year)

#plotting

plot_dissim_taxo_Hill <- ggplot(TD_beta_k_nok) +
  geom_point( aes(x=Year, y=q1_mean), stat="identity", size=2, shape=16, color="seagreen4") +
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+ q1_se), width=0.4, size=0.8, color="seagreen4") +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.25, 0.45), breaks = seq(from=0.25, to=0.45, by=0.1)  ) +
  labs(x="", y="Taxonomic dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none") 

plot_dissim_taxo_Hill


##################################### Functional dissimilarity #########################################################

# Obtain average per year 

FD_beta_k_nok <- FD_beta_kelp_nokelp %>%
  filter(Year1 == Year2, Habitat1==Habitat2) %>% 
  select(Year=Year1, q1) %>% 
  group_by(Year) %>%
  summarise( n = n(),
             q1_mean = mean(q1),
             q1_sd = sd(q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))

FD_beta_k_nok$Year <- as.numeric(FD_beta_k_nok$Year)

#plotting

plot_dissim_func_Hill <- ggplot(FD_beta_k_nok) +
  geom_point( aes(x=Year, y=q1_mean), stat="identity", size=2, shape=16, color="seagreen4") +
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+ q1_se), width=0.4, size=0.8, color="seagreen4") +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  labs(x="", y="Functional dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none") 

plot_dissim_func_Hill

## merging all plot into a single figure and saving as png ####
figure4 <- ( plot_dissim_taxo_Hill / plot_dissim_func_Hill)
ggsave(figure4, file=here::here("outputs/", "figure4_Hill.png"),
       height = 22, width = 25, unit = "cm" )
