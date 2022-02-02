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

load(here::here("outputs/", "spatial_beta.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "nokelp_metadata.RData") )
load(here::here("outputs/", "temporal_beta_kelp.RData") )
load(here::here("outputs/", "temporal_beta_nokelp.RData") )


## taxonomic dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat ----
kelp_interyear_taxo_diss <- dendextend::dist_long(temporal_beta_kelp$taxo_diss) %>%
  left_join( select(kelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, distance)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             distance_mean = mean(distance),
             distance_sd = sd(distance)
  ) %>%
  mutate(distance_se = distance_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_taxo_diss <- dendextend::dist_long(temporal_beta_nokelp$taxo_diss) %>%
  left_join( select(nokelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, distance)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             distance_mean = mean(distance),
             distance_sd = sd(distance)
  ) %>%
  mutate(distance_se = distance_sd/sqrt(n))  %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----
interhab_taxo_diss <- dendextend::dist_long(spatial_beta$taxo_diss) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, distance )

spatial_diss <- interhab_taxo_diss %>%
  filter(Habitat1==Habitat2) %>%
  mutate(Habitat = Habitat1) %>% 
  group_by(Habitat) %>%
  summarise( n = n(),
             distance_mean = mean(distance),
             distance_sd = sd(distance)
  ) %>%
  mutate(distance_se = distance_sd/sqrt(n))  %>%
  mutate( Year = 2013, .before="Year" ) %>%
  droplevels("Habitat")


# merging ----
taxo_dissim_toplot <- bind_rows(kelp_interyear_taxo_diss, 
                            nokelp_interyear_taxo_diss)
                            #spatial_diss) 

head(taxo_dissim_toplot)

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="#3BB372", No_Kelp="#74E7B8",
               Inshore="#2C6BAA", Midshelf="lightsalmon1", Offshore="firebrick3")

plot_dissim_taxo <- ggplot(taxo_dissim_toplot, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point( aes( x=Year, y=distance_mean), stat="identity", size=2) +
 geom_line( aes( x=Year, y=distance_mean), stat="identity", size=1) +
  geom_errorbar( aes(x=Year, ymin=distance_mean-distance_se, ymax=distance_mean+distance_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.2,0.8), breaks = seq(from=0.2, to=0.8, by=0.2)  ) +
  labs(x="", y="Species dissimilarity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
plot_dissim_taxo



## functional dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat ----
kelp_interyear_func_diss <- dendextend::dist_long(temporal_beta_kelp$func_diss) %>%
  left_join( select(kelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, distance)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             distance_mean = mean(distance),
             distance_sd = sd(distance)
  ) %>%
  mutate(distance_se = distance_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )




# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_func_diss <- dendextend::dist_long(temporal_beta_nokelp$func_diss) %>%
  left_join( select(nokelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, distance)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             distance_mean = mean(distance),
             distance_sd = sd(distance)
  ) %>%
  mutate(distance_se = distance_sd/sqrt(n))  %>%
  mutate( Habitat = "No_Kelp", .before="Year" )


# adding dissimilarity between the 3 habitats (from 2013) ----
interhab_func_diss <- dendextend::dist_long(spatial_beta$func_diss) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, distance)

spatial_diss_func <- interhab_func_diss %>%
  filter(Habitat1==Habitat2) %>%
  mutate(Habitat = Habitat1) %>% 
  group_by(Habitat) %>%
  summarise( n = n(),
             distance_mean = mean(distance),
             distance_sd = sd(distance)
  ) %>%
  mutate(distance_se = distance_sd/sqrt(n))  %>%
  mutate( Year = 2013, .before="Year" ) %>%
  droplevels("Habitat")



# merging ----
func_dissim_toplot <- bind_rows(kelp_interyear_func_diss, 
                                nokelp_interyear_func_diss)
                                #spatial_diss_func) 

head(func_dissim_toplot)

# plotting ----


# color code for the 3 habitats

hab_colors <- c(Kelp="#3BB372", No_Kelp="#74E7B8",
                Inshore="#2C6BAA", Midshelf="lightsalmon1", Offshore="firebrick3")


plot_dissim_func <- ggplot(func_dissim_toplot, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point( aes( x=Year, y=distance_mean), stat="identity", size=2) +
  geom_line( aes( x=Year, y=distance_mean), stat="identity", size=1) +
  geom_errorbar( aes(x=Year, ymin=distance_mean-distance_se, ymax=distance_mean+distance_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0,0.6), breaks = seq(from=0, to=0.6, by=0.2)  ) +
  labs(x="", y="Functional dissimilarity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_dissim_func



#####
## merging all plot into a single figure and saving as png ####
figure2_convex_site_average <- ( plot_dissim_taxo / plot_dissim_func )
ggsave(figure2_convex_site_average, file=here::here("outputs/", "figure2_convex_site_average.png"), 
       height = 20, width = 18, unit = "cm" )
