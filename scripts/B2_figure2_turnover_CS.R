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


## taxonomic turnover ####

# values between assemblages from 2002 and later in kelp habitat, but only 
#between  the same sites ---- ----

kelp_interyear_taxo_turn <- dendextend::dist_long(temporal_beta_kelp$taxo_turn) %>%
  left_join( select(kelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, taxo_turn = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             taxo_turn_mean = mean(taxo_turn),
             taxo_turn_sd = sd(taxo_turn)
  ) %>%
  mutate( taxo_turn_se = taxo_turn_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )


# values between assemblages from 2002 and later in no kelp habitat, but only 
#between  the same sites ---- ---- ----
nokelp_interyear_taxo_turn <- dendextend::dist_long(temporal_beta_nokelp$taxo_turn) %>%
  left_join( select(nokelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, taxo_turn = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             taxo_turn_mean = mean(taxo_turn),
             taxo_turn_sd = sd(taxo_turn)
  ) %>%
  mutate( taxo_turn_se = taxo_turn_sd/sqrt(n)) %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding turnover between the 3 habitats (from 2013) ---- ----

interhab_taxo_turn <- dendextend::dist_long(spatial_beta$taxo_turn) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, taxo_turn = distance )

spatial_turn <- interhab_taxo_turn %>%
  filter(Habitat1==Habitat2) %>% 
  mutate(Habitat = Habitat1) %>% 
  group_by(Habitat) %>%
  summarise( n = n(),
             taxo_turn_mean = mean(taxo_turn),
             taxo_turn_sd = sd(taxo_turn)
  ) %>%
  mutate( taxo_turn_se = taxo_turn_sd/sqrt(n)) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  droplevels("Habitat")


# merging ----
taxo_turnim_toplot <- bind_rows(kelp_interyear_taxo_turn, 
                            nokelp_interyear_taxo_turn)
                            #spatial_turn) 

head(taxo_turnim_toplot)

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="seagreen4", No_Kelp="seagreen2",
                Inshore="#2C6BAA", Midshelf="lightsalmon1", Offshore="firebrick3")



plot_turnim_taxo <- ggplot(taxo_turnim_toplot, mapping = aes(color=Habitat)) +
  geom_point( aes(x=Year, y=taxo_turn_mean), stat="identity", size=2) +
  geom_line(aes(x=Year, y=taxo_turn_mean), stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=taxo_turn_mean-taxo_turn_se, ymax=taxo_turn_mean+taxo_turn_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors) + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0.3,0.8), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="", y="Species turnover") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
plot_turnim_taxo



## functional turnover ####

# values between assemblages from 2002 and later in kelp habitat, but only 
#between  the same sites ---- ---- ---- ----

kelp_interyear_func_turn <- dendextend::dist_long(temporal_beta_kelp$func_turn) %>%
  left_join( select(kelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1==Site2) %>% 
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, func_turn = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             func_turn_mean = mean(func_turn),
             func_turn_sd = sd(func_turn)
  ) %>%
  mutate( func_turn_se = func_turn_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat, but only 
#between  the same sites ---- ---- ---- ---- ----

nokelp_interyear_func_turn <- dendextend::dist_long(temporal_beta_nokelp$func_turn) %>%
  left_join( select(nokelp_metadata, Code, Site1=Site, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Site2=Site, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>% 
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, func_turn = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             func_turn_mean = mean(func_turn),
             func_turn_sd = sd(func_turn)
  ) %>%
  mutate( func_turn_se = func_turn_sd/sqrt(n)) %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding turnover between the 3 habitats (from 2013) ----
interhab_func_turn <- dendextend::dist_long(spatial_beta$func_turn) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, func_turn = distance )

spatial_turn <- interhab_func_turn %>%
  filter(Habitat1==Habitat2) %>% 
  mutate(Habitat = Habitat1) %>% 
  group_by(Habitat) %>%
  summarise( n = n(),
             func_turn_mean = mean(func_turn),
             func_turn_sd = sd(func_turn)
  ) %>%
   mutate( func_turn_se = func_turn_sd/sqrt(n)) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  droplevels("Habitat")


# merging ----
func_turnim_toplot <- bind_rows(kelp_interyear_func_turn, 
                                nokelp_interyear_func_turn)
                                #spatial_turn) 

head(func_turnim_toplot)

# plotting ----

plot_turnim_func <- ggplot(func_turnim_toplot, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point( aes(x=Year, y=func_turn_mean), stat="identity", size=2) +
  geom_line(aes(x=Year, y=func_turn_mean), stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=func_turn_mean-func_turn_se, ymax=func_turn_mean+func_turn_se), 
                  width=0.2, size=0.5) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors)  + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0,0.6), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="", y="Functional turnover") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.background = element_blank(), legend.key=element_blank())
plot_turnim_func



fturn <- ggplot(func_turnim_toplot, aes(x=Year, y=func_turn_mean, color=Habitat)) +
  geom_point() +
  geom_smooth(method = "lm")+
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors)  + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0,0.6), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="", y="Functional turnover") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (12),colour = "black"), axis.title = element_text(size= (14)),
        legend.title = element_text(size=12), legend.text = element_text(size=12),
        legend.background = element_blank(), legend.key=element_blank())

#####
## merging all plot into a single figure and saving as png ####
figure2_turnover <- ( plot_turnim_taxo / plot_turnim_func )
ggsave(figure2_turnover, file=here::here("outputs/", "figure2_turnover.png"), 
       height = 20, width = 18, unit = "cm" )