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

load(here::here("outputs/", "FD_beta_kelp_nokelp.RData") )
load(here::here("outputs/", "TD_beta_kelp_nokelp.RData") )


load(here::here("outputs/", "TD_kelp_nokelp_Hill.RData") )
load(here::here("outputs/", "FD_kelp_nokelp_Hill.RData") )


load(here::here("outputs/", "FD_beta_spatial_Hill.RData") )
load(here::here("outputs/", "TD_beta_spatial_Hill.RData") )


## taxonomic dissimilarity ####

# Obtain average per year and habitat

TD_beta_k_nok <- TD_beta_kelp_nokelp %>%
  filter(Year1 == Year2, Habitat1==Habitat2) %>% 
  select(Year=Year1, Habitat=Habitat1, q1) %>% 
  group_by(Year, Habitat) %>%
  summarise( n = n(),
             q1_mean = mean(q1),
             q1_sd = sd(q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))

# Add spatial data

## Adding spatial data

TD_beta_s <- TD_beta_spatial %>% 
  filter(Habitat1==Habitat2) %>% 
  mutate(Year="2013") %>% 
  select(Year, Habitat=Habitat1, q1) %>% 
  group_by(Year, Habitat) %>%
  summarise( n = n(),
             q1_mean = mean(q1),
             q1_sd = sd(q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))


#Join datasets

TD_beta_toplot <- rbind(TD_beta_k_nok, 
                        TD_beta_s)


# color code for the 3 habitats

hab_colors <- c(Kelp="seagreen4", No_kelp="seagreen2",
                I="blue", M="lightsalmon1", O="firebrick3")
hab_shape <- c(Kelp=15, No_kelp=15,
               I=25, M=25, O=25)

#plotting

plot_dissim_taxo_Hill <- ggplot(TD_beta_toplot, 
                                mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=q1_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+ q1_se), width=0.4, size=0.8) +
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
 # scale_y_continuous( limits = c(0.2, 0.6), breaks = seq(from=0.2, to=0.6, by=0.1)  ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"), 
                     labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"), 
                    labels=c("Kelp", "No kelp", "Inshore","Midshelf", "Offshore")) + 
  #scale_shape_manual(values=hab_shape, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"),
                     #labels=c("Kelp", "No kelp", "Inshore","Midshelf", "Offshore"))+
  labs(x="", y="Taxonomic dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_dissim_taxo_Hill