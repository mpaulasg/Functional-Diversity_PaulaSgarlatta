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
library(dendextend)
library(patchwork)

# loading data
load(here::here("data", "spatial_metadata.RData") )

load(here::here("outputs/", "spatial_beta_k.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "nokelp_metadata.RData") )
load(here::here("outputs/", "temporal_beta_kelp_k.RData") )
load(here::here("outputs/", "temporal_beta_nokelp_k.RData") )


## taxonomic dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat ----
kelp_interyear_taxo_diss <- dendextend::dist_long(temporal_beta_kelp_k$taxo_diss) %>%
  left_join( select(kelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  mutate(Site = Site1 ) %>%
  select(Year, Site, taxo_diss = distance )  %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_taxo_diss <- dendextend::dist_long(temporal_beta_nokelp_k$taxo_diss) %>%
  left_join( select(nokelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  mutate(Site = Site1 ) %>%
  select(Year, Site, taxo_diss = distance )  %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----
interhab_taxo_diss <- dendextend::dist_long(spatial_beta_k$taxo_diss) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, taxo_diss = distance )

spatial_diss <- interhab_taxo_diss %>%
  group_by(Habitat1, Habitat2) %>%
  summarise(taxo_diss = mean(taxo_diss) ) %>%
  ungroup(Habitat1, Habitat2) %>%
  mutate( Habitat = paste(Habitat1, Habitat2, sep="_"), .before="Habitat1" ) %>%
  select(-Habitat1, -Habitat2) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  filter(Habitat %in%c("Midshelf_Inshore", "Offshore_Midshelf", "Offshore_Inshore") ) %>%
  droplevels("Habitat")


# merging ----
taxo_dissim_toplot <- bind_rows(kelp_interyear_taxo_diss, 
                                nokelp_interyear_taxo_diss,
                                spatial_diss) 

head(taxo_dissim_toplot)

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="seagreen4", No_Kelp="seagreen2",
                Midshelf_Inshore="mediumseagreen", Offshore_Midshelf="lightsalmon1", Offshore_Inshore="firebrick3")
hab_shape <- c(Kelp=21, No_Kelp=21,
               Midshelf_Inshore=22, Offshore_Midshelf=22, Offshore_Inshore=22)

plot_dissim_taxo <- ggplot(taxo_dissim_toplot, aes(x=Year, y=taxo_diss, group=Site) ) +
  geom_point( aes( shape= Habitat, color = Habitat, fill = Habitat), stat="identity") +
  geom_line( aes(color = Habitat), size=0.5 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_shape_manual(values=hab_shape) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.2,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="", y="Species dissimilarity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
plot_dissim_taxo



## functional dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat ----
kelp_interyear_func_diss <- dendextend::dist_long(temporal_beta_kelp_k$func_diss) %>%
  left_join( select(kelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  mutate(Site = Site1 ) %>%
  select(Year, Site, func_diss = distance )  %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_func_diss <- dendextend::dist_long(temporal_beta_nokelp_k$func_diss) %>%
  left_join( select(nokelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  mutate(Site = Site1 ) %>%
  select(Year, Site, func_diss = distance )  %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----
interhab_func_diss <- dendextend::dist_long(spatial_beta_k$func_diss) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, func_diss = distance )

spatial_diss <- interhab_func_diss %>%
  group_by(Habitat1, Habitat2) %>%
  summarise(func_diss = mean(func_diss) ) %>%
  ungroup(Habitat1, Habitat2) %>%
  mutate( Habitat = paste(Habitat1, Habitat2, sep="_"), .before="Habitat1" ) %>%
  select(-Habitat1, -Habitat2) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  filter(Habitat %in%c("Midshelf_Inshore", "Offshore_Midshelf", "Offshore_Inshore") ) %>%
  droplevels("Habitat")


# merging ----
func_dissim_toplot <- bind_rows(kelp_interyear_func_diss, 
                                nokelp_interyear_func_diss,
                                spatial_diss) 

head(func_dissim_toplot)

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="seagreen4", No_Kelp="seagreen2",
                Midshelf_Inshore="mediumseagreen", Offshore_Midshelf="lightsalmon1", Offshore_Inshore="firebrick3")
hab_shape <- c(Kelp=21, No_Kelp=21,
               Midshelf_Inshore=22, Offshore_Midshelf=22, Offshore_Inshore=22)


plot_dissim_func <- ggplot(func_dissim_toplot, aes(x=Year, y=func_diss, group=Site) ) +
  geom_point( aes( shape= Habitat, color = Habitat, fill = Habitat), stat="identity") +
  geom_line( aes(color = Habitat), size=0.5 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_shape_manual(values=hab_shape) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,0.7), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="", y="Functional dissimilarity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.background = element_blank(), legend.key=element_blank())
plot_dissim_func



#####
## merging all plot into a single figure and saving as png ####
figure2_site <- ( plot_dissim_taxo / plot_dissim_func )
ggsave(figure2_site, file=here::here("outputs/", "figure2_site_k.png"), 
       height = 12, width = 18, unit = "cm" )
