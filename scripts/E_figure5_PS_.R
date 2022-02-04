################################################################################
##
## Script for plotting Kmax values
## 
## Code by Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(dplyr)
library(ggplot2)


# Load data

traits <- read.csv(here::here("data", "raw_data", "fish_traits.csv"),
                        header = T)
thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv"),
                    header = T)%>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  select(-thermal)

# Join these 2 datasets

traits_full <- inner_join(traits, thermal, by="Species", all.x=TRUE)

#Load species occurrences

nokelp <- read.csv(here::here("from_paula",
      "TemporalBRUV_species_maxN_no_kelp.csv"),header = T) %>% 
  mutate(Site1 = sub(".*_", "", Site), .before="Site") %>%
  distinct(Site1, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:108) %>% 
    filter(Abundance>0) %>% 
  select(-Abundance, -Site)

# Join with traits

nokelp_traits <- inner_join(traits_full, nokelp, by="Species", all.x=TRUE)

#Let's try to plot

#color code

thermal_col <- c(temperate="blue", tropical="red")


nokelp_kmax <- ggplot(nokelp_traits, mapping = aes(color=thermal_label)) +
  geom_point( aes( x=Site1, y=mean(Kmax)), stat="identity", size=2) +
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
plot_dissim_taxo



kelp <- read.csv(here::here("from_paula",
    "TemporalBRUV_species_maxN_kelp.csv"),header = T) %>% 
  mutate(Year = sub(".*_", "", Site), .before="Site") %>%
  distinct(Year, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:103) %>% 
  filter(Abundance>0) %>% 
  select(-Abundance, -Site)

# Join with traits

kelp_traits <- inner_join(traits_full, kelp, by="Species", all.x=TRUE)

#Let's try to plot

kelp_kmax <- ggplot(kelp_traits, mapping = aes(color=thermal_label)) +
  geom_point( aes( x=Year, y=Kmax), stat="identity", size=2) +
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")


spatial <- read.csv(here::here("from_paula",
         "SpatialUVC_species_biomass_transect.csv"),header = T) %>% 
  mutate(Habitat = sub(".*_", "", Site), .before="Site") %>%
  distinct(Habitat, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:103) %>% 
  filter(Abundance>0) %>% 
  select(-Abundance, -Site)


kelp <- read.csv(here::here("from_paula",
                            "TemporalBRUV_species_maxN_kelp.csv"),header = T) %>% 
  mutate(Year = sub(".*_", "", Site), .before="Site") %>%
  distinct(Year, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:103) %>% 
  filter(Abundance>0) %>% 
  select(-Abundance, -Site)

# Join with traits

kelp_traits <- inner_join(traits_full, kelp, by="Species", all.x=TRUE)

#Let's try to plot

kelp_kmax <- ggplot(kelp_traits, mapping = aes(color=thermal_label)) +
  geom_point( aes( x=Year, y=Kmax), stat="identity", size=2) +
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")