
################################################################################
##
## Script for merging biomass from temporal BRUV surveys 
## 
##  Sébastien Villéger
##
## NB: NOT TO BE ON FINAL GITHUB
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)


################################################################################

## temporal survey data from sites with kelp ####

# loading raw data from csv----
kelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv"))
head(kelp_metadata)  

kelp_replicates_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv"))
head(kelp_replicates_maxN)  

# merging metadata and data 
kelp <- left_join(kelp_metadata, kelp_replicates_maxN, by="Code" )
head(kelp)

# computing total biomass per site for each species ----
kelp_sites_maxN <- kelp %>% 
  pivot_longer(cols= contains("_"), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Year, species) %>%
  summarize(total=sum(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Year,sep="_"), .before="Site")

head(kelp_sites_maxN)

kelp_replicates_maxN %>% select( contains("_") ) %>% sum()
kelp_sites_maxN %>% select( contains("_") ) %>% sum()
# => same total biomass

# splitting into 2 tables ----
TemporalBRUV_kelp_metadata <- kelp_sites_maxN %>% 
  select(Code, Site, Year)

TemporalBRUV_kelp_species <- kelp_sites_maxN %>% 
  select( "Code", contains("_") )



## temporal survey data from sites with kelp ####

# loading raw data from csv----
nokelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_no_kelp.csv"))
head(nokelp_metadata)  

nokelp_replicates_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_no_kelp.csv"))
head(nokelp_replicates_maxN)  

# merging metadata and data 
nokelp <- left_join(nokelp_metadata, nokelp_replicates_maxN, by=c("Code"="Site" ))
head(nokelp)

# computing total biomass per site for each species ----
nokelp_sites_maxN <- nokelp %>% 
  pivot_longer(cols= contains("_"), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Year, species) %>%
  summarize(total=sum(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Year,sep="_"), .before="Site")

head(nokelp_sites_maxN)

nokelp_replicates_maxN %>% select( contains("_") ) %>% sum()
nokelp_sites_maxN %>% select( contains("_") ) %>% sum()
# => same total biomass

# splitting into 2 tables ----
TemporalBRUV_nokelp_metadata <- nokelp_sites_maxN %>% 
  select(Code, Site, Year)

TemporalBRUV_nokelp_species <- nokelp_sites_maxN %>% 
  select( "Code", contains("_") )





## saving as csv ----
write.csv(TemporalBRUV_nokelp_metadata, file=here::here("data", "raw_data", "TemporalBRUV_nokelp_metadata.csv"), row.names = FALSE )
write.csv(TemporalBRUV_nokelp_species, file=here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv"), row.names= FALSE)


## end of script ####
################################################################################



# metadata and data from sites (and replicates) without kelp
nokelp_metadata_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_no_kelp.csv") )
head(nokelp_metadata_all)
unique(nokelp_metadata_all$Site) # 4 sites
unique(nokelp_metadata_all$Year) # 14 years

nokelp_sp_maxN_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_no_kelp.csv") ) %>%
  column_to_rownames("Site") %>%
  as.matrix()
head(nokelp_sp_maxN_all)

# metadata and data from sites that use to have kelp and lost it:
kelp_metadata_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv") )
head(kelp_metadata_all)
unique(kelp_metadata_all$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv") ) %>%
  column_to_rownames("Site") %>%
  as.matrix()
head(kelp_sp_maxN_all) #######HERE













################################################################################
##
## Script for plotting taxonomic and functional diversity across habitats and years
##
## Fig. 1 - Species richness and functional richness
## 
## Fig. 2 - Functional Dispersion
##
## Fig. 1-3S - Functional Identity
##
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
load(here::here("outputs/", "temporal_alpha_kelp.RData") )

## spatial trends ####

# merging metadata and biodiv indices in a single table for each habitat type

spatial_all <- spatial_metadata %>% 
  left_join( rownames_to_column(spatial_alpha, "Code"), by="Code" ) %>%
  select(Code, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3)



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
    fide_PC3_sd = sd(fide_PC3)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( fric_se = fric_sd/sqrt(n)) %>% 
  mutate( fdis_se = fdis_sd/sqrt(n)) %>% 
  mutate( fide_PC1_se = fide_PC1_sd/sqrt(n)) %>%
  mutate( fide_PC2_se = fide_PC2_sd/sqrt(n)) %>%
  mutate( fide_PC3_se = fide_PC3_sd/sqrt(n))


spatial_toplot

# color code for the 3 habitats
hab_colors <- c(Inshore= "#2C6BAA", Midshelf= "lightsalmon1", Offshore="firebrick3")

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
  geom_bar( aes(x=Habitat, y=fric_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fric_mean-fric_se, ymax=fric_mean+fric_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,0.5), breaks = seq(from=0, to=0.5, by=0.1)  ) +
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_func


## temporal trends ####

# merging metadata and biodiv indices in a single table for kelp

temporal_kelp<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(temporal_alpha_kelp, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3 )




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
  scale_y_continuous( limits = c(15,30), breaks = seq(from=15, to=30, by=5)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") + 
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_tempo_taxo

# functional ----

plot_tempo_func <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fric_mean), stat="identity", size=2, shape=16) +
  geom_line(aes(x= Year, y= fric_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fric_mean-fric_se, ymax=fric_mean+fric_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") + 
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())
plot_tempo_func


## merging all plot into a single figure and saving as png ####
figure1 <- ( plot_spatial_taxo + plot_tempo_taxo ) / ( plot_spatial_func +  plot_tempo_func )

ggsave(figure1, file=here::here("outputs/", "Figure1_bis.png"),
       height = 22, width = 25, unit = "cm" )


#################################################### FDis ############################################################


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
  labs(x="", y="Functional dispersion") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_fdis

## merging all plot into a single figure and saving as png ####

figure2_fdis <- ( plot_spatial_fdis + plot_tempo_fdis )

ggsave(figure2_fdis, file=here::here("outputs/", "Figure2_fdis_bis.png"),
       height = 22, width = 35, unit = "cm" )


################################################### FIde #############################################################


plot_spatial_fide1 <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fide_PC1_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fide_PC1_mean-fide_PC1_se, ymax=fide_PC1_mean+fide_PC1_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(-0.04,0.08), breaks = seq(from=-0.04, to=0.08, by=0.02)  ) +
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
  scale_y_continuous( limits = c(-0.02,0.08), breaks = seq(from=-0.02, to=0.08, by=0.02)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="Functional Identity PC1") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_fide1

## merging all plot into a single figure and saving as png ####
figure3_fide1 <- ( plot_spatial_fide1 + plot_tempo_fide1 )
ggsave(figure3_fide1, file=here::here("outputs/", "Figure1S_fide1_bis.png"),
       height = 22, width = 35, unit = "cm" )

################################### FIde 2 #######################################################

plot_spatial_fide2 <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fide_PC2_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fide_PC2_mean-fide_PC2_se, ymax=fide_PC2_mean+fide_PC2_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  #scale_y_continuous( limits = c(0,0.01), breaks = seq(from=0, to=0.01, by=0.0025)  ) +
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
  #scale_y_continuous( limits = c(-0.03,0.09), breaks = seq(from=-0.03, to=0.09, by=0.03)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="Functional Identity PC2") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_fide2

## merging all plot into a single figure and saving as png ####
figure1_fide2 <- ( plot_spatial_fide2 + plot_tempo_fide2 )
ggsave(figure1_fide2, file=here::here("outputs/", "Figure2S_fide2_bis.png"),
       height = 22, width = 35, unit = "cm" )

############################# FIde 3

plot_spatial_fide3 <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fide_PC3_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fide_PC3_mean-fide_PC3_se, ymax=fide_PC3_mean+fide_PC3_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  #scale_y_continuous( limits = c(0,0.01), breaks = seq(from=0, to=0.01, by=0.0025)  ) +
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
  #scale_y_continuous( limits = c(-0.03,0.09), breaks = seq(from=-0.03, to=0.09, by=0.03)  ) +
  scale_color_manual(values="seagreen4") + 
  scale_fill_manual(values="seagreen4") +
  labs(x="", y="Functional Identity PC3") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_fide3

## merging all plot into a single figure and saving as png ####

figure1_fide3 <- ( plot_spatial_fide3 + plot_tempo_fide3 )
ggsave(figure1_fide3, file=here::here("outputs/", "Figure3S_fide3_bis.png"),
       height = 22, width = 35, unit = "cm" )
