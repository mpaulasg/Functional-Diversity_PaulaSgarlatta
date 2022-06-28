################################################################################
##
## Script for plotting several figures
##
## Fig. 1 - Map of study site
## 
## Fig. 2 - Species richness, biomass and functional richness
## 
## Fig. S1 - Functional Dispersion
##
## Fig. S2-S4 - Functional Identity
##
## Fig. S8 - Habitat
##
##
## Code by Paula Sgarlatta, Sebastien Villeger and Camille Magneville 
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(dplyr)
library(patchwork)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(ggplot2)
library(ozmaps)
library(cowplot)
library(ggpubr)
library(ggspatial)
library(here)

####### Fig. 1 - Map of study site

sites <- read.csv(here::here("data", "raw_data", "sites.csv"))

buffer <- 0.5

geo_bounds <- c(left = min(sites$Longitude)-buffer, 
                bottom = min(sites$Latitude)-buffer, 
                right = max(sites$Longitude)+buffer, 
                top = max(sites$Latitude)+buffer)

min_lon <- 153.05
max_lon <- 153.55
min_lat <- -30.4
max_lat <- -29.8


geo_bounds <- c(left = min_lon, bottom = min_lat, right = max_lon, top = max_lat)

Sites.grid <- expand.grid(lon_bound = c(geo_bounds[1], geo_bounds[3]), 
                          lat_bound = c(geo_bounds[2], geo_bounds[4]))

coordinates(Sites.grid) <- ~ lon_bound + lat_bound

Aus <- readOGR(dsn = "C:/Users/z5179758/Google Drive/PhD/GitHub/Functional-Diversity_PaulaSgarlatta/data/raw_data/61395_shp/australia",layer = "cstauscd_r") 

Aus_coast <- subset(Aus, FEAT_CODE != "sea")

Aus_crop <- crop(Aus_coast, extent(Sites.grid))



color_data <- c(Inshore= "#482677FF", Midshelf= "#FD6A02", 
                Offshore="#FDE725FF", Temporal= "seagreen4")



shape_data <- c(Inshore= 19, Midshelf= 19, 
                Offshore=19, Temporal= 17)

sites_plot <- ggplot()+ theme_classic() + 
  geom_polygon(data = Aus_crop, aes(x=long, y=lat, group=group), fill="grey", colour="black") +
  coord_equal(ratio = 1)+
  geom_point(data=sites, aes(x=Longitude, y=Latitude, colour = Habitat,
                             shape=Habitat, size=2)) +
  scale_color_manual(name="Habitat", values=color_data) +
  scale_shape_manual(name="Habitat", values=shape_data)+
  # annotation_scale(location = "br", plot_unit = "km")+ #width_hint = 0.5) +  #br is the location - bottom, right
  #annotation_north_arrow(location = "br", which_north = "true", 
                         #pad_x = unit(1.5, "in"), pad_y = unit(0.1, "in"),
                         #style = north_arrow_fancy_orienteering) +
  theme(panel.background = element_rect(fill = "white"),
        legend.position = c(0.83,0.6), legend.title = element_blank(),
        legend.text = element_text(size=12),
        panel.border = element_rect(colour = "black",fill = NA),
        axis.text = element_text(size = (14), colour="black"), 
        axis.title = element_blank()) + scale_size(guide = "none")+
  guides(color = guide_legend(override.aes = list(size = 3) ) ) +
 scale_y_continuous(expand = c(0,0))+
scale_x_continuous(expand = c(0,0))

sites_plot

#Plot of Australia/rectangle in the area of study

Aus <- ggplot(data = ozmap(x="country")) + 
  geom_sf(color="black", fill= "white") +
  geom_rect(xmin = 152, xmax = 154, ymin = -32, ymax = -29, 
            fill = NA, colour = "black", size = 0.5) +
  scale_fill_viridis_d(option = "plasma") +
  coord_sf(xlim = c(100.00, 160.00), ylim = c(-45.00, -10.00), expand = TRUE) +
  theme_void()

###Place both figures together

map <- ggdraw(sites_plot) +
  draw_plot(Aus, x = 0.22, y = 0.73, width = 0.28, height = 0.28)

ggsave(map, file=here::here("outputs", "Figure1.jpeg"),  
       height = 20, width = 18, unit = "cm")

################ Taxonomic and functional diversity graphs

# loading data
load(here::here("data", "spatial_metadata.RData") )
load(here::here("data", "spatial_alpha_biomass.RData") )

load(here::here("data","kelp_metadata.RData") )
load(here::here("data" ,"temporal_alpha_kelp_biomass.RData") )

## spatial trends ####

# merging metadata and biodiversity indices in a single table for each habitat type

spatial_all <- spatial_metadata %>% 
  left_join( rename (spatial, Code = Site), by="Code" ) %>%
  dplyr::select(Code, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3, Biomass)


# mean and sd of diversity among each site for each year in each habitat type
spatial_toplot <- spatial_all %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    fric_mean = mean(fric),
    fric_sd = sd(fric), 
    fdis_mean = mean(fdis),
    fdis_sd = sd(fdis), 
    fide_PC1_mean = mean(fide_PC1),
    fide_PC1_sd = sd(fide_PC1),
    fide_PC2_mean = mean(fide_PC2),
    fide_PC2_sd = sd(fide_PC2),
    fide_PC3_mean = mean(fide_PC3),
    fide_PC3_sd = sd(fide_PC3), 
    biomass_mean = mean(Biomass),
    biomass_sd = sd(Biomass)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( fric_se = fric_sd/sqrt(n)) %>% 
  mutate( fdis_se = fdis_sd/sqrt(n)) %>% 
  mutate( fide_PC1_se = fide_PC1_sd/sqrt(n)) %>%
  mutate( fide_PC2_se = fide_PC2_sd/sqrt(n)) %>%
  mutate( fide_PC3_se = fide_PC3_sd/sqrt(n)) %>% 
  mutate( biomass_se = biomass_sd/sqrt(n))
  
  spatial_toplot

# color code for the 3 habitats
hab_colors <- c(Inshore= "#482677FF", Midshelf= "#FD6A02", Offshore="#FDE725FF")

## Fig. 2 - Species richness and functional richness

# taxonomic ----

plot_spatial_taxo <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=TRic_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,30), breaks = seq(from=0, to=30, by=5)  ) +
  labs(x="", y="Species richness") +
  ggtitle(("(b)"))+
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
       legend.position = "none", axis.text.x = element_blank(),
       plot.title = element_text(size = 16, face = "bold"))
  
plot_spatial_taxo

# functional ----

plot_spatial_func <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fric_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fric_mean-fric_se, ymax=fric_mean+fric_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,0.8), breaks = seq(from=0, to=0.8, by=0.2)  ) +
  labs(x="", y="Functional richness") +
   ggtitle(("(d)"))+
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
          legend.position = "none", plot.title = element_text(size = 16, face = "bold"))
 

plot_spatial_func

# Biomass ----

plot_spatial_biomass <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=biomass_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=biomass_mean-biomass_se, ymax=biomass_mean+biomass_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  labs(x="", y="Biomass (g/125m2) ") +
  ggtitle(("(e)"))+
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none", plot.title = element_text(size = 16, face = "bold"))


plot_spatial_biomass


## temporal trends ####

# merging metadata and biodiv indices in a single table for kelp

temporal_kelp<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(temporal_alpha_kelp, "Code"), by="Code" ) %>%
  dplyr::select(Code, Site, Year, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3 )


# mean and sd of diversity among each site for each year in each habitat type
temporal_toplot <- temporal_kelp %>%
  group_by(Year, Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    fric_mean = mean(fric),
    fric_sd = sd(fric), 
    fdis_mean = mean(fdis),
    fdis_sd = sd(fdis), 
    fide_PC1_mean = mean(fide_PC1),
    fide_PC1_sd = sd(fide_PC1),
    fide_PC2_mean = mean(fide_PC2),
    fide_PC2_sd = sd(fide_PC2),
    fide_PC3_mean = mean(fide_PC3),
    fide_PC3_sd = sd(fide_PC3)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( fric_se = fric_sd/sqrt(n)) %>% 
  mutate( fdis_se = fdis_sd/sqrt(n)) %>% 
  mutate( fide_PC1_se = fide_PC1_sd/sqrt(n)) %>%
  mutate( fide_PC2_se = fide_PC2_sd/sqrt(n)) %>%
  mutate( fide_PC3_se = fide_PC3_sd/sqrt(n))
  
temporal_toplot

unique(temporal_toplot$Year)


# taxonomic ----

plot_tempo_taxo <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=TRic_mean), stat="identity", size=3) +
  geom_line(aes(x= Year, y= TRic_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,30), breaks = seq(from=0, to=30, by=5)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") + 
  labs(x="", y="Species richness") +
  ggtitle(("(a)")) + 
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none", plot.title = element_text(size = 16, face = "bold"))
  

plot_tempo_taxo
  
# functional ----

plot_tempo_func <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fric_mean), stat="identity", size=2, shape=16) +
  geom_line(aes(x= Year, y= fric_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fric_mean-fric_se, ymax=fric_mean+fric_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,0.8), breaks = seq(from=0, to=0.8, by=0.2)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") + 
  labs(x="", y="Functional richness") +
  ggtitle(("(c)")) +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none", plot.title = element_text(size = 16, face = "bold")) 
  

plot_tempo_func


## merging all plot into a single figure and saving as png ####

figure2 <- ggarrange (plot_tempo_taxo, plot_spatial_taxo, plot_tempo_func, plot_spatial_func,
  plot_spatial_biomass, ncol = 2, nrow = 3)
           

ggsave(figure2, file=here::here("outputs" , "Figure2.jpeg"),
       height = 22, width = 25, unit = "cm" )
  

######### Fig. S1 - Functional Dispersion


plot_spatial_fdis <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fdis_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fdis_mean-fdis_se, ymax=fdis_mean+fdis_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,0.8), breaks = seq(from=0, to=0.8, by=0.2)  ) +
  labs(x="", y="Functional dispersion") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_spatial_fdis

## temporal trends ####

plot_tempo_fdis <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fdis_mean), stat="identity", size=3) +
  geom_line(aes(x= Year, y= fdis_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fdis_mean-fdis_se, ymax=fdis_mean+fdis_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,0.8), breaks = seq(from=0, to=0.8, by=0.1)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")


plot_tempo_fdis

## merging all plot into a single figure and saving as png ####

figureS1_fdis <- ( plot_spatial_fdis + plot_tempo_fdis )

ggsave(figureS1_fdis, file=here::here("outputs",  "FigureS1.jpeg"),
       height = 22, width = 35, unit = "cm" )


#### Fig. S2-S3 - Functional Identity

## FIde PC1

plot_spatial_fide1 <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fide_PC1_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fide_PC1_mean-fide_PC1_se, ymax=fide_PC1_mean+fide_PC1_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  labs(x="", y="Functional Identity PC1") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_fide1

## temporal trends ####

plot_tempo_fide1 <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fide_PC1_mean), stat="identity", size=3) +
  geom_line(aes(x= Year, y= fide_PC1_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fide_PC1_mean-fide_PC1_se, ymax=fide_PC1_mean+fide_PC1_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_tempo_fide1

## merging all plot into a single figure and saving as png ####
figureS2_fide1 <- ( plot_spatial_fide1 + plot_tempo_fide1 )

ggsave(figureS2_fide1, file=here::here("outputs",  "FigureS2.jpeg"),
       height = 22, width = 35, unit = "cm" )
                     
## FIde PC2

plot_spatial_fide2 <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fide_PC2_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fide_PC2_mean-fide_PC2_se, ymax=fide_PC2_mean+fide_PC2_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  labs(x="", y="Functional Identity PC2") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_fide2

## temporal trends ####

plot_tempo_fide2 <- ggplot(temporal_toplot, 
                           mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fide_PC2_mean), stat="identity", size=3) +
  geom_line(aes(x= Year, y= fide_PC2_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fide_PC2_mean-fide_PC2_se, ymax=fide_PC2_mean+fide_PC2_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_tempo_fide2

## merging all plot into a single figure and saving as png ####
figureS3_fide2 <- ( plot_spatial_fide2 + plot_tempo_fide2 )

ggsave(figureS3_fide2, file=here::here("outputs",  "FigureS3.jpeg"),
       height = 22, width = 35, unit = "cm" )

## FIde PC3

plot_spatial_fide3 <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fide_PC3_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fide_PC3_mean-fide_PC3_se, ymax=fide_PC3_mean+fide_PC3_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  labs(x="", y="Functional Identity PC3") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_fide3

## temporal trends ####

plot_tempo_fide3 <- ggplot(temporal_toplot, 
                           mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fide_PC3_mean), stat="identity", size=3) +
  geom_line(aes(x= Year, y= fide_PC3_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fide_PC3_mean-fide_PC3_se, ymax=fide_PC3_mean+fide_PC3_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_tempo_fide3

## merging all plot into a single figure and saving as png ####

figureS4_fide3 <- ( plot_spatial_fide3 + plot_tempo_fide3 )

ggsave(figureS4_fide3, file=here::here("outputs",  "FigureS4.jpeg"),
       height = 22, width = 35, unit = "cm" )


################### Habitat


# Load data

habitat <- read.csv(here::here("data", "raw_data", "habitat_solitaries_2012.csv"))

## First graph only with coral as broad category:

habitat_toplot_full <- habitat %>% 
  group_by(Site, Habitat,Transect) %>% 
  mutate (Coral = sum(Plate_coral, Branching_coral, Columnar_coral, Colony_corals, Foliose_coral, Soft_coral)) %>% 
  dplyr::select(-Plate_coral, -Branching_coral, -Columnar_coral, -Colony_corals, -Foliose_coral, -Soft_coral) %>% 
  rename("Ecklonia radiata" = Ecklonia_radiata , "Turf & CCA" = Turf_CCA, "Other invertebrates" = Other_invertebrates,
         "Sponges & tunicates" = Sponges_tunicates, "Rock & sand" = Rock_sand) %>% 
  pivot_longer(cols= 4:10, names_to="Group") %>%
  filter(Site!="Flat_top", Site!="Look_at_me_now") %>% 
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Habitat,Transect, Group) %>%
  summarize(total=sum(value)) %>% 
  group_by(Site, Transect) %>% 
  mutate(percent_cover=(total*100/125)) %>% #125 points per transect
  dplyr::select(-total)

habitat_summary_full <-habitat_toplot_full%>%
  group_by(Habitat, Group)%>%
  summarise(Percent_cover_habitat=mean(percent_cover,na.rm=T),
            n = n(), mean = mean(percent_cover), se = sd(percent_cover/sqrt(n))) %>% 
  mutate(se) %>%
  group_by(Habitat) %>%
  arrange(desc(Group)) %>%
  mutate(
    pos = cumsum(Percent_cover_habitat),
    upper = pos + se/2,
    lower = pos - se/2
  ) %>%
  ungroup() %>%  
  filter(mean > 1)

# color code for the 3 habitats
colors_full <- c("Ecklonia radiata" = "tan3", Coral="indianred1", Macroalgae = "seagreen", 
                 "Other invertebrates" ="darkolivegreen4", "Rock & sand"  ="orange",
                 "Sponges & tunicates" = "mediumorchid1","Turf & CCA" ="red2")

## Plot figure - full

habitat_full <- ggplot(habitat_summary_full, aes(x=Habitat, y=Percent_cover_habitat, fill=Group)) +
  geom_bar(stat="identity", width = 0.7, size = 1, color = "black") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.8, width=0.1, position = "identity", color = "black")+
  scale_color_manual(values=colors_full) + 
  scale_fill_manual(values=colors_full) + 
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (18), color = "black"), axis.title = element_text(size= (18)),
        legend.key = element_rect(fill = "white"), legend.text = element_text(size=20), 
        legend.title = element_blank()) +  labs(y="Percent cover (%)", x="")

ggsave(habitat_full, file=here::here("outputs", "habitat_full.jpeg"),
       height = 16, width = 24, unit = "cm" )


### Only corals

habitat_toplot_coral <- habitat %>% 
  rename( "Plate coral" = Plate_coral, "Branching coral" = Branching_coral,
          "Columnar coral" = Columnar_coral, "Colony coral" = Colony_corals, "Foliose coral" = Foliose_coral,
          "Soft coral" = Soft_coral) %>% 
  pivot_longer(cols= 8:13, names_to="Group") %>%
  filter(Site!="Flat_top", Site!="Look_at_me_now") %>% 
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Habitat,Transect, Group) %>%
  summarize(total=sum(value)) %>% 
  group_by(Site, Transect) %>% 
  mutate(percent_cover=(total*100/125)) %>% #125 points per transect
  dplyr::select(-total)



habitat_summary_coral <-habitat_toplot_coral%>%
  group_by(Habitat, Group)%>%
  summarise(Percent_cover_habitat=mean(percent_cover,na.rm=T),
            n = n(), mean = mean(percent_cover), se = sd(percent_cover/sqrt(n))) %>% 
  mutate(se) %>%
  group_by(Habitat) %>%
  arrange(desc(Group)) %>%
  mutate(
    pos = cumsum(Percent_cover_habitat),
    upper = pos + se/2,
    lower = pos - se/2
  ) %>%
  ungroup() 

 # filter(mean > 1)

color_coral <- c("Branching coral" ="orange", "Colony coral" = "yellowgreen", "Columnar coral" ="indianred1",
                  "Foliose coral" ="cornflowerblue", "Plate coral" ="red2",
                 "Soft coral" = "darkblue")

## Plot figure

habitat_coral <- ggplot(habitat_summary_coral, aes(x=Habitat, y=Percent_cover_habitat, fill=Group)) +
  geom_bar(stat="identity", width = 0.7, size = 1, color = "black") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.8, width=0.1, position = "identity", color = "black")+
  scale_color_manual(values=color_coral) + 
  scale_fill_manual(values=color_coral) + 
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (18), color = "black"), axis.title = element_text(size= (18)),
        legend.key = element_rect(fill = "white"), legend.text = element_text(size=20), legend.title = element_blank())+
  labs(y="Percent cover (%)", x="") 

ggsave(habitat_coral, file=here::here("outputs", "habitat_coral.jpeg"),
       height = 16, width = 24, unit = "cm" )


######################## end of code #########################