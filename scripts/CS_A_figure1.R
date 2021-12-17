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
hab_colors <- c(Inshore= "mediumseagreen", Midshelf= "lightsalmon1", Offshore="firebrick3")

# taxonomic ----

plot_spatial_taxo <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=TRic_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,25), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
       legend.position = "none", axis.text.x = element_blank())
plot_spatial_taxo

# functional ----

plot_spatial_func <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=FRic_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,0.35), breaks = seq(from=0, to=0.5, by=0.1)  ) +
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
          legend.position = "none")

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

# color code for kelp/no kelp

year_colors <- c(kelp= "mediumseagreen", no_kelp= "lightsalmon1")


# taxonomic ----

plot_tempo_taxo <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=TRic_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(15,37), breaks = seq(from=0, to=35, by=5)  ) +
  scale_color_manual(values=year_colors) + 
  scale_fill_manual(values=year_colors) + 
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_tempo_taxo
  
# functional ----

plot_tempo_func <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=FRic_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.3,0.8), breaks = seq(from=0, to=1, by=0.2)  ) +
  scale_color_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  scale_fill_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
           panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
           axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())
plot_tempo_func



## merging all plot into a single figure and saving as png ####
figure1 <- ( plot_spatial_taxo + plot_tempo_taxo ) / ( plot_spatial_func +  plot_tempo_func )
ggsave(figure1, file=here::here("outputs/", "figure1.png"),
       height = 22, width = 25, unit = "cm" )
  

