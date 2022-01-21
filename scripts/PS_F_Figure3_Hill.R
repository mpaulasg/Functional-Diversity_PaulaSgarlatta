############################################################################################################
##
## Script for plotting taxonomic and functional diversity across habitats and years using Hill numbers
## 
## Code by Paula Sgarlatta, Sebastien Villeger and Camille Magneville 
##
###########################################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

# loading data

load(here::here("outputs/", "FD_beta_nokelp_Hill.RData") )
load(here::here("outputs/", "TD_beta_nokelp_Hill.RData") )


load(here::here("outputs/", "FD_beta_kelp_Hill.RData") )
load(here::here("outputs/", "TD_beta_kelp_Hill.RData") )


load(here::here("outputs/", "FD_beta_spatial_Hill.RData") )
load(here::here("outputs/", "TD_beta_spatial_Hill.RData") )


## taxonomic dissimilarity ####

# average dissimilarity between the kelp and no_kelp sites

#Join datasets

TD_beta_kelp$Habitat <- "Kelp"
TD_beta_nokelp$Habitat <- "No_kelp"

kelp_nokelp <- bind_rows(TD_beta_kelp,TD_beta_nokelp) 

kelp_nokelp$Year <- ifelse(kelp_nokelp$Year1 == kelp_nokelp$Year2, 
                           kelp_nokelp$Year1, NA)

taxo_diss_Hill <- kelp_nokelp  %>%
    select(Year, q0, Habitat) %>% 
  na.omit(kelp_nokelp) %>% 
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n))

#plotting

plot_dissim_taxo_Hill <- taxo_diss_Hill %>% 
  ggplot(aes(x=Year,y=q0_mean)) + 
  geom_point(aes(x=Year, y=q0_mean), stat="identity", size=3) +
  geom_errorbar( aes(x=Year, ymin=q0_mean-q0_se, ymax=q0_mean+q0_se), 
                 width=0.2, size=0.5 ) +
  labs(x="", y="Taxonomic dissimilarity (q=0)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")

plot_dissim_taxo_Hill
