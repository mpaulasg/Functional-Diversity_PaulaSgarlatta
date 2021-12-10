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
load(here::here("data", "spatial_metadata.RData") )
load(here::here("outputs/", "spatial_alpha.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "nokelp_metadata.RData") )
load(here::here("outputs/", "temporal_alpha_kelp.RData") )
load(here::here("outputs/", "temporal_alpha_nokelp.RData") )


## spatial trends ####

# merging metadata and biodiv indices in a single table for each habitat type

spatial_all <- spatial_metadata %>% 
  left_join( rownames_to_column(spatial_alpha, "Code"), by="Code" ) %>%
  select(Code, Habitat, TRic=sp_richn, FRic=fric)



# mean and sd of diversity among each site for each year in each habitat type
spatial_toplot <- spatial_all %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    FRic_mean = mean(FRic),
    FRic_sd = sd(FRic)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( FRic_se = FRic_sd/sqrt(n))

spatial_toplot

# color code for the 3 habitats
hab_colors <- c(Inshore= "gold", Midshelf= "darkorange", Offshore="blue")

# taxonomic ----

plot_spatial_taxo <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=TRic_mean, color = Habitat, fill = Habitat), stat="identity") +
  geom_errorbar( aes(x=Habitat, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.5, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="Year", y="Species richness") +
  theme_bw()
plot_spatial_taxo

# functional ----

plot_spatial_func <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=FRic_mean, color = Habitat, fill = Habitat), stat="identity") +
  geom_errorbar( aes(x=Habitat, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.1, size=0.5, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="Year", y="Functional richness") +
  theme_bw()
plot_spatial_func


## temporal trends ####

# merging metadata and biodiv indices in a single table for each habitat type

temporal_kelp<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(temporal_alpha_kelp, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, FRic=fric)

temporal_nokelp<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(temporal_alpha_nokelp, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, FRic=fric)

# merging values for the 2 habitat types
temporal_all <- bind_rows(temporal_kelp, temporal_nokelp )
head(temporal_all)


# mean and sd of diversity among each site for each year in each habitat type
temporal_toplot <- temporal_all %>%
  group_by(Year, Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    FRic_mean = mean(FRic),
    FRic_sd = sd(FRic)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( FRic_se = FRic_sd/sqrt(n))
  
temporal_toplot

unique(temporal_toplot$Year)


# taxonomic ----

plot_tempo_taxo <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=TRic_mean), stat="identity", size=2, shape=22) +
  geom_errorbar( aes(x=Year, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.2) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=2)  ) +
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="Year", y="Species richness") +
  theme_bw()
plot_tempo_taxo
  
# functional ----

plot_tempo_func <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=FRic_mean), stat="identity", size=2, shape=22) +
  geom_errorbar( aes(x=Year, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.1, size=0.2) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=2)  ) +
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="Year", y="Functional richness") +
  theme_bw()
plot_tempo_func



## merging all plot into a single figure and saving as png ####
figure1 <- ( plot_spatial_taxo + plot_tempo_taxo ) / ( plot_spatial_func +  plot_tempo_func )
ggsave(figure1, file=here::here("outputs/", "figure1.png"),
       height = 12, width = 25, unit = "cm" )
  

###### [PS] Now with thermal affinity ####

# loading data

load(here::here("outputs/", "spatial_alpha_thermal.RData") )


load(here::here("outputs/", "temporal_alpha_kelp_thermal.RData") )
load(here::here("outputs/", "temporal_alpha_nokelp_thermal.RData") )


## spatial trends ####

# merging metadata and biodiv indices in a single table for each habitat type

spatial_all_thermal <- spatial_metadata %>% 
   left_join( rownames_to_column(spatial_alpha_thermal, "Code"), by="Code" ) %>%
   select(Code, Habitat, TRic=sp_richn, FRic=fric)



# mean and sd of diversity among each site for each year in each habitat type
spatial_toplot_thermal <- spatial_all_thermal %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    FRic_mean = mean(FRic),
    FRic_sd = sd(FRic)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( FRic_se = FRic_sd/sqrt(n))

spatial_toplot_thermal

# color code for the 3 habitats
hab_colors <- c(Inshore= "gold", Midshelf= "darkorange", Offshore="blue")

# taxonomic ----

plot_spatial_taxo_thermal <- ggplot(spatial_toplot_thermal) +
  geom_bar( aes(x=Habitat, y=TRic_mean, color = Habitat, fill = Habitat), stat="identity") +
  geom_errorbar( aes(x=Habitat, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.5, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="Year", y="Species richness") +
  theme_bw()
plot_spatial_taxo_thermal

# functional ----

plot_spatial_func_thermal <- ggplot(spatial_toplot_thermal) +
  geom_bar( aes(x=Habitat, y=FRic_mean, color = Habitat, fill = Habitat), stat="identity") +
  geom_errorbar( aes(x=Habitat, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.1, size=0.5, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="Year", y="Functional richness") +
  theme_bw()
plot_spatial_func_thermal


## temporal trends ####

# merging metadata and biodiv indices in a single table for each habitat type

temporal_kelp_thermal<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(temporal_alpha_kelp_thermal, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, FRic=fric)

temporal_nokelp_thermal<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(temporal_alpha_nokelp_thermal, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, FRic=fric)

# merging values for the 2 habitat types
temporal_all_thermal <- bind_rows(temporal_kelp_thermal, temporal_nokelp_thermal )
head(temporal_all_thermal)


# mean and sd of diversity among each site for each year in each habitat type
temporal_toplot_thermal <- temporal_all_thermal %>%
  group_by(Year, Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    FRic_mean = mean(FRic),
    FRic_sd = sd(FRic)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( FRic_se = FRic_sd/sqrt(n))

temporal_toplot_thermal

unique(temporal_toplot_thermal$Year)


# taxonomic ----

plot_tempo_taxo_thermal <- ggplot(temporal_toplot_thermal, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=TRic_mean), stat="identity", size=2, shape=22) +
  geom_errorbar( aes(x=Year, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.2) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=2)  ) +
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="Year", y="Species richness") +
  theme_bw()
plot_tempo_taxo_thermal

# functional ----

plot_tempo_func_thermal <- ggplot(temporal_toplot_thermal, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=FRic_mean), stat="identity", size=2, shape=22) +
  geom_errorbar( aes(x=Year, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.1, size=0.2) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=2)  ) +
  scale_y_continuous( limits = c(0,1), breaks = seq(from=0, to=1, by=0.2)  ) +
  labs(x="Year", y="Functional richness") +
  theme_bw()
plot_tempo_func_thermal



## merging all plot into a single figure and saving as png ####
figure1_thermal <- ( plot_spatial_taxo_thermal + plot_tempo_taxo_thermal ) / 
  ( plot_spatial_func_thermal +  plot_tempo_func_thermal )
ggsave(figure1_thermal, file=here::here("outputs/", "figure1_thermal.png"),
       height = 12, width = 25, unit = "cm" )
       