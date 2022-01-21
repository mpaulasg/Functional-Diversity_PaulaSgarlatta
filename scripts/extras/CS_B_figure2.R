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
  left_join( select(kelp_metadata, Code, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, taxo_diss = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             taxo_diss_mean = mean(taxo_diss),
             taxo_diss_sd = sd(taxo_diss)
  ) %>%
  mutate( taxo_diss_se = taxo_diss_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_taxo_diss <- dendextend::dist_long(temporal_beta_nokelp$taxo_diss) %>%
  left_join( select(nokelp_metadata, Code, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, taxo_diss = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             taxo_diss_mean = mean(taxo_diss),
             taxo_diss_sd = sd(taxo_diss)
  ) %>%
  mutate( taxo_diss_se = taxo_diss_sd/sqrt(n)) %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----
interhab_taxo_diss <- dendextend::dist_long(spatial_beta$taxo_diss) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, taxo_diss = distance )

spatial_diss <- interhab_taxo_diss %>%
  group_by(Habitat1, Habitat2) %>%
  summarise( n = n(),
             taxo_diss_mean = mean(taxo_diss),
             taxo_diss_sd = sd(taxo_diss)
  ) %>%
  ungroup(Habitat1, Habitat2) %>%
  mutate( taxo_diss_se = taxo_diss_sd/sqrt(n)) %>%
  mutate( Habitat = paste(Habitat1, Habitat2, sep="_"), .before="Habitat1" ) %>%
  select(-Habitat1, -Habitat2) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  filter(Habitat %in%c("Midshelf_Inshore", "Offshore_Midshelf", "Offshore_Inshore") ) %>%
  droplevels("Habitat")
typeof(spatial_diss$Year)
typeof(kelp_interyear_taxo_diss$Year)
typeof(nokelp_interyear_taxo_diss$Year)

# merging ----
taxo_dissim_toplot <- bind_rows(kelp_interyear_taxo_diss, 
                            nokelp_interyear_taxo_diss,
                            spatial_diss) 
typeof(taxo_dissim_toplot$Year)

head(taxo_dissim_toplot)

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="seagreen4", No_Kelp="seagreen2",
                Midshelf_Inshore="mediumseagreen", Offshore_Midshelf="lightsalmon1", Offshore_Inshore="firebrick3")
hab_shape <- c(Kelp=21, No_Kelp=21,
                Midshelf_Inshore=22, Offshore_Midshelf=22, Offshore_Inshore=22)


plot_dissim_taxo <- ggplot(taxo_dissim_toplot, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point(aes(x=Year, y=taxo_diss_mean, 
                shape= Habitat), stat="identity", size=3) +
  geom_errorbar( aes(x=Year, ymin=taxo_diss_mean-taxo_diss_se, ymax=taxo_diss_mean+taxo_diss_se), 
                 width=0.2, size=0.5 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_shape_manual(values=hab_shape) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.4,0.9), breaks = seq(from=0.4, to=0.9, by=0.1)  ) +
  labs(x="", y="Species dissimilarity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
plot_dissim_taxo


## functional dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat ----
kelp_interyear_func_diss <- dendextend::dist_long(temporal_beta_kelp$func_diss) %>%
  left_join( select(kelp_metadata, Code, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(kelp_metadata, Code, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, func_diss = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             func_diss_mean = mean(func_diss),
             func_diss_sd = sd(func_diss)
  ) %>%
  mutate( func_diss_se = func_diss_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_func_diss <- dendextend::dist_long(temporal_beta_nokelp$func_diss) %>%
  left_join( select(nokelp_metadata, Code, Year1=Year), by=c("rows"="Code" ) ) %>%
  left_join( select(nokelp_metadata, Code, Year2=Year), by=c("cols"="Code" ) ) %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, func_diss = distance )  %>%
  group_by(Year) %>%
  summarise( n = n(),
             func_diss_mean = mean(func_diss),
             func_diss_sd = sd(func_diss)
  ) %>%
  mutate( func_diss_se = func_diss_sd/sqrt(n)) %>%
  mutate( Habitat = "No_Kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----
interhab_func_diss <- dendextend::dist_long(spatial_beta$func_diss) %>%
  left_join( select(spatial_metadata, Code, Habitat1=Habitat), by=c("rows"="Code" ) ) %>%
  left_join( select(spatial_metadata, Code, Habitat2=Habitat), by=c("cols"="Code" ) ) %>%
  select(Habitat1, Habitat2, func_diss = distance )

spatial_diss <- interhab_func_diss %>%
  group_by(Habitat1, Habitat2) %>%
  summarise( n = n(),
             func_diss_mean = mean(func_diss),
             func_diss_sd = sd(func_diss)
  ) %>%
  ungroup(Habitat1, Habitat2) %>%
  mutate( func_diss_se = func_diss_sd/sqrt(n)) %>%
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

plot_dissim_func <- ggplot(func_dissim_toplot, mapping =aes(color=Habitat, fill=Habitat) ) +
  geom_point( aes(x=Year, y=func_diss_mean, 
                  shape= Habitat), stat="identity", size=3) +
  geom_errorbar( aes(x=Year, ymin=func_diss_mean-func_diss_se, ymax=func_diss_mean+func_diss_se), 
                  width=0.2, size=0.5 ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_Kelp", "Midshelf_Inshore", "Offshore_Midshelf", 
  "Offshore_Inshore"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_Kelp", "Midshelf_Inshore", "Offshore_Midshelf", 
 "Offshore_Inshore"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_shape_manual(values=hab_shape, name="Habitat", breaks = c("Kelp", "No_Kelp", "Midshelf_Inshore", "Offshore_Midshelf", 
  "Offshore_Inshore"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.2,0.6), breaks = seq(from=0, to=0.6, by=0.2)  ) +
  labs(x="", y="Functional dissimilarity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())
plot_dissim_func

#####
## merging all plot into a single figure and saving as png ####
figure2_diss <- ( plot_dissim_taxo / plot_dissim_func )
ggsave(figure2_diss, file=here::here("outputs/", "figure2_diss.png"), 
       height = 20, width = 18, unit = "cm" )


#######################################################################################################################

######################################## Trying with FRic #############################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(dendextend)
library(patchwork)

# loading data
load(here::here("data", "spatial_metadata.RData") )
load(here::here("outputs/", "spatial_alpha.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "nokelp_metadata.RData") )
load(here::here("outputs/", "temporal_alpha_kelp.RData") )
load(here::here("outputs/", "temporal_alpha_nokelp.RData") )


## temporal trends ####

# merging metadata and biodiv indices in a single table for each habitat

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

## Spatial trends

# merging metadata and biodiv indices in a single table for each habitat type

spatial_all <- spatial_metadata %>% 
  left_join( rownames_to_column(spatial_alpha, "Code"), by="Code" ) %>%
  select(Code, Site, Habitat, TRic=sp_richn, FRic=fric) %>% 
mutate( Year = 2013, .before="Year" )


# mean and sd of diversity among each site for each year in each habitat type
spatial_toplot <- spatial_all %>%
  group_by(Habitat, Year) %>%
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

# merging all data

fric_all <- bind_rows(temporal_toplot, spatial_toplot)
head(fric_all)


# plotting ----

# color code for the 3 habitats
hab_colors <- c(kelp="seagreen4", no_kelp="seagreen2",
                Inshore="blue", Midshelf="lightsalmon1", Offshore="firebrick3")
hab_shape <- c(kelp=21, no_kelp=21,
               Inshore=22, Midshelf=22, Offshore=22)


plot_tax <-ggplot(fric_all, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point(aes(x=Year, y=TRic_mean, shape= Habitat), stat="identity", size=2) +
  geom_errorbar( aes(x=Year, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), 
                 width=0.2, size=0.5 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_shape_manual(values=hab_shape) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(10,40), breaks = seq(from=10, to=40, by=10)  ) +
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_tax


## functional richness ####


# plotting ----

plot_fric <-ggplot(fric_all, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point(aes(x=Year, y=FRic_mean, shape= Habitat), stat="identity", size=2) +
  geom_errorbar( aes(x=Year, ymin=FRic_mean-FRic_se, ymax=FRic_mean+FRic_se), 
                 width=0.2, size=0.5 ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("kelp", "no_kelp", "Inshore", "Midshelf", "Offshore"), 
                     labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("kelp", "no_kelp", "Inshore", "Midshelf", "Offshore"), 
                    labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  scale_shape_manual(values=hab_shape, name="Habitat", breaks = c("kelp", "no_kelp", "Inshore", "Midshelf", "Offshore"), 
                     labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.3,0.8), breaks = seq(from=0.3, to=0.8, by=0.2)  ) +
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())
plot_fric



#####
## merging all plot into a single figure and saving as png ####
figure2_fric <- ( plot_tax / plot_fric )
ggsave(figure2_fric, file=here::here("outputs/", "figure2_fric.png"), 
       height = 20, width = 18, unit = "cm" )
