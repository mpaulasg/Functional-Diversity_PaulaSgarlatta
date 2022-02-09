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


nokelp_kmax <- ggplot(nokelp_traits, aes( x=Site1, y=Kmax, fill=thermal_label)) +
  geom_violin(position = "dodge", alpha=0.5, outlier.colour="transparent")+
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), binwidth = 0.01)+
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax-no kelp") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
nokelp_kmax



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

kelp_kmax <- ggplot(kelp_traits,  aes( x=Year, y=Kmax, fill=thermal_label)) +
  geom_violin(position = "dodge", alpha=0.5, outlier.colour="transparent") +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), binwidth = 0.01)+
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax - Kelp") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

kelp_kmax


## merging all plot into a single figure and saving as png ####

figure5_kmax <- ( nokelp_kmax / kelp_kmax )
ggsave(figure5_kmax, file=here::here("outputs/", "figure5_kmax.png"), 
       height = 20, width = 25, unit = "cm" )

###########################################  Spatial ###############################################################

spatial <- read.csv(here::here("data", "raw_data",
         "SpatialUVC_species_biomass_site_average.csv"),header = T) %>% 
  mutate(Habitat = sub(".*_", "", Code), .before="Code") %>%
  #distinct(Habitat, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:53) %>% 
  filter(Abundance>0) %>% 
  select(-Abundance, -Code)


# Join with traits

spatial_traits <- inner_join(traits_full, spatial, by="Species", all.x=TRUE)

#Let's try to plot

spatial_kmax <- ggplot(spatial_traits, aes(x=Habitat, y=Kmax, fill=thermal_label)) +
  geom_violin (position = "dodge") +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), binwidth = 0.01)+
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax - Spatial") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

spatial_kmax

## saving as png ####

 ggsave(spatial_kmax, file=here::here("outputs/", "figure5_kmax_spatial.png"), 
       height = 20, width = 25, unit = "cm" )

####################### Position ###########################


#Load species occurrences

nokelp_agg <- read.csv(here::here("from_paula",
                              "TemporalBRUV_species_maxN_no_kelp.csv"),header = T) %>% 
  mutate(Site1 = sub(".*_", "", Site), .before="Site") %>%
  distinct(Site1, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:108) %>% 
  filter(Abundance>0) %>% 
  select(-Site)

# Join with traits

nokelp_traits_agg <- inner_join(traits_full, nokelp_agg, by="Species", all.x=TRUE)

#Let's try to plot

#color code

thermal_col <- c(temperate="blue", tropical="red")


nokelp_agg <- ggplot(nokelp_traits_agg, aes( x=Site1, y=Abundance, colour=Agg)) +
  geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=thermal_col) + 
  #scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax-no kelp") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
nokelp_kmax



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

kelp_kmax <- ggplot(kelp_traits,  aes( x=Year, y=Kmax, fill=thermal_label)) +
  geom_violin(position = "dodge", alpha=0.5, outlier.colour="transparent") +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), binwidth = 0.01)+
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax - Kelp") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

kelp_kmax


## merging all plot into a single figure and saving as png ####

figure5_kmax <- ( nokelp_kmax / kelp_kmax )
ggsave(figure5_kmax, file=here::here("outputs/", "figure5_kmax.png"), 
       height = 20, width = 25, unit = "cm" )

###########################################  Spatial ###############################################################

spatial <- read.csv(here::here("data", "raw_data",
                               "SpatialUVC_species_biomass_site_average.csv"),header = T) %>% 
  mutate(Habitat = sub(".*_", "", Code), .before="Code") %>%
  #distinct(Habitat, .keep_all = TRUE) %>% 
  gather(Species, Abundance, 3:53) %>% 
  filter(Abundance>0) %>% 
  select(-Abundance, -Code)


# Join with traits

spatial_traits <- inner_join(traits_full, spatial, by="Species", all.x=TRUE)

#Let's try to plot

spatial_kmax <- ggplot(spatial_traits, aes(x=Habitat, y=Kmax, fill=thermal_label)) +
  geom_violin (position = "dodge") +
  geom_dotplot(binaxis='y', stackdir='center',
               position=position_dodge(1), binwidth = 0.01)+
  scale_color_manual(values=thermal_col) + 
  scale_fill_manual(values=thermal_col) + 
  #scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Kmax - Spatial") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

spatial_kmax

## saving as png ####

ggsave(spatial_kmax, file=here::here("outputs/", "figure5_kmax_spatial.png"), 
       height = 20, width = 25, unit = "cm" )
