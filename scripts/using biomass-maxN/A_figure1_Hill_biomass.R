################################################################################
##
## Script for plotting taxonomic and functional diversity across habitats and years
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(patchwork)

# loading data

load(here::here("outputs/","using biomass-maxN",  "TD_kelp_Hill.RData"))
load(here::here("outputs/","using biomass-maxN",  "FD_kelp_Hill.RData"))
load(here::here("data", "kelp_metadata.RData"))

load(here::here("outputs/", "using biomass-maxN", "TD_spatial_Hill.RData"))
load(here::here("outputs/", "using biomass-maxN", "FD_spatial_Hill.RData"))
load(here::here("data", "spatial_metadata.RData"))


## spatial trends ####

TD_spatial <- as.data.frame(TD_spatial)

Tax_spatial_Hill <- spatial_metadata %>% 
  left_join( rownames_to_column(TD_spatial, "Code"), by="Code" ) %>%
  select(Code, Habitat, FD_q0)



# mean and sd of diversity among each site for each year in each habitat type
TDspatial_Hill_toplot <- Tax_spatial_Hill %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    FD_mean = mean(FD_q0),
    FD_sd = sd(FD_q0)
  ) %>%
  mutate( FD_se = FD_sd/sqrt(n))


# color code for the 3 habitats
hab_colors <- c(Inshore= "#2C6BAA", Midshelf= "lightsalmon1", Offshore="firebrick3")

# taxonomic ----

plot_spatial_taxo_Hill <- ggplot(TDspatial_Hill_toplot) +
  geom_bar( aes(x=Habitat, y=FD_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=FD_mean-FD_se, ymax=FD_mean+FD_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,25), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none", axis.text.x = element_blank())

plot_spatial_taxo_Hill

# functional ----

FD_spatial <- as.data.frame(FD_spatial)

Fx_spatial_Hill <- spatial_metadata %>% 
  left_join( rownames_to_column(FD_spatial, "Code"), by="Code" ) %>%
  select(Code, Habitat,  FD_q1 )



# mean and sd of diversity among each site for each year in each habitat type
FDspatial_Hill_toplot <- Fx_spatial_Hill %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    q1_mean = mean(FD_q1),
    q1_sd = sd(FD_q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))
  

plot_spatial_func_Hill <- ggplot(FDspatial_Hill_toplot) +
  geom_bar( aes(x=Habitat, y=q1_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=q1_mean-q1_se, ymax=q1_mean+q1_se), width=0.1, size=0.8) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  #scale_y_continuous( limits = c(0,5), breaks = seq(from=0, to=5, by=1)  ) +
  labs(x="", y="FD Hill numbers (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_func_Hill


## temporal trends ####

# merging metadata and biodiv indices in a single table for each habitat

TD_kelp_Hill <- as.data.frame(TD_kelp_Hill)

temporal_kelp_Hill<- kelp_metadata %>% 
  left_join( rownames_to_column(TD_kelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, FD_q0) %>% 
  group_by(Year) %>% 
  summarise( 
    n = n(),
    FD_q0_mean = mean(FD_q0),
    FD_q0_sd = sd(FD_q0)
  ) %>%
  mutate( FD_q0_se = FD_q0_sd/sqrt(n)) %>% 
  mutate( Habitat = "Kelp", .before="Year" )



# taxonomic ----

plot_tempo_taxo_Hill <- ggplot(temporal_kelp_Hill, mapping = aes(color=Habitat,  fill=Habitat)) +
  geom_point (aes(x= Year, y= FD_q0_mean) , stat="identity", size=3) +
  geom_line(aes(x= Year, y= FD_q0_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=FD_q0_mean-FD_q0_se, ymax=FD_q0_mean+FD_q0_se), 
                 width=0.1, size=1) +
 scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
   labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_tempo_taxo_Hill

# functional ----

FD_kelp_Hill <- as.data.frame(FD_kelp_Hill)

Ftemporal_kelp_Hill<- kelp_metadata %>% 
  left_join( rownames_to_column(FD_kelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, FD_q1) %>% 
  group_by(Year) %>% 
  summarise( 
    n = n(),
    FD_q1_mean = mean(FD_q1),
    FD_q1_sd = sd(FD_q1)
  ) %>%
  mutate( FD_q1_se = FD_q1_sd/sqrt(n)) %>% 
  mutate( Habitat = "Kelp", .before="Year" )


plot_tempo_func_Hill <- ggplot(Ftemporal_kelp_Hill, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=FD_q1_mean), stat="identity", size=3) +
  geom_line(aes(x=Year, y=FD_q1_mean), stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=FD_q1_mean-FD_q1_se, ymax=FD_q1_mean+FD_q1_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(4,6), breaks = seq(from=0, to=6, by=1)  ) +
  scale_color_manual(values="seagreen4") +
  scale_fill_manual(values="seagreen4") + 
  labs(x="", y=" FD Hill Numbers(q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_func_Hill



## merging all plot into a single figure and saving as png ####

figure1a_Hill <- ( plot_spatial_taxo_Hill + plot_tempo_taxo_Hill ) / ( plot_spatial_func_Hill +  plot_tempo_func_Hill )
ggsave(figure1a_Hill, file=here::here("outputs/", "using biomass-maxN", "figure1a_Hill_biomass.png"),
       height = 22, width = 25, unit = "cm" )

######################## end of code ################################