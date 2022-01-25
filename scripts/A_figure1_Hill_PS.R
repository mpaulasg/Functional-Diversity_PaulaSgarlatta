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

load(here::here("outputs/", "TD_nokelp_Hill.RData") )
load(here::here("outputs/", "FD_nokelp_Hill.RData") )
load(here::here("data", "nokelp_metadata.RData"))

load(here::here("outputs/", "TD_kelp_Hill.RData"))
load(here::here("outputs/", "FD_kelp_Hill.RData"))
load(here::here("data", "kelp_metadata.RData"))

load(here::here("outputs/", "TD_spatial_Hill.RData"))
load(here::here("outputs/", "FD_spatial_Hill.RData"))
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
hab_colors <- c(Inshore= "mediumseagreen", Midshelf= "lightsalmon1", Offshore="firebrick3")

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
  select(Code, Habitat, FD_q0, FD_q1, FD_q2)



# mean and sd of diversity among each site for each year in each habitat type
FDspatial_Hill_toplot <- Fx_spatial_Hill %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    q0_mean = mean(FD_q0),
    q0_sd = sd(FD_q0)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n))
  

plot_spatial_func_Hill <- ggplot(FDspatial_Hill_toplot,
                                 mapping=aes(color = Habitat, fill = Habitat)) +
  geom_point( aes(x=Habitat, y=q0_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Habitat, ymin=q0_mean-q0_se, ymax=q0_mean+q0_se), width=0.1, size=0.8) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(3.7,4.5), breaks = seq(from=3.7, to=4.5, by=0.1)  ) +
  labs(x="", y="Hill numbers FD (q=0)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_func_Hill


## temporal trends ####

# merging metadata and biodiv indices in a single table for each habitat

TD_kelp_Hill <- as.data.frame(TD_kelp_Hill)

temporal_kelp_Hill<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(TD_kelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, Habitat, FD_q0)

TD_nokelp_Hill <- as.data.frame(TD_nokelp_Hill)

temporal_nokelp_Hill<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(TD_nokelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, Habitat, FD_q0)

# merging values for the 2 habitat types
temporal_all_Hill <- bind_rows(temporal_kelp_Hill, temporal_nokelp_Hill )


# mean and sd of diversity among each site for each year in each habitat type
temporal_toplot_Hill <- temporal_all_Hill %>%
  group_by(Year, Habitat) %>%
  summarise( 
    n = n(),
    FD_q0_mean = mean(FD_q0),
    FD_q0_sd = sd(FD_q0)
  ) %>%
  mutate( FD_q0_se = FD_q0_sd/sqrt(n))

unique(temporal_toplot_Hill$Year)

# color code for kelp/no kelp

year_colors <- c(kelp= "mediumseagreen", no_kelp= "lightsalmon1")


# taxonomic ----

plot_tempo_taxo_Hill <- ggplot(temporal_toplot_Hill, aes(x=Year, y=FD_q0_mean, fill=Habitat)) +
  geom_bar(stat="identity", position="dodge", color = "black", size=0.8) +
  geom_errorbar( aes(x=Year, ymin=FD_q0_mean-FD_q0_se, ymax=FD_q0_mean+FD_q0_se), width=0.1, size=0.8, 
                 position=position_dodge(.9), colour="black") +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
  scale_color_manual(values=year_colors) + 
  scale_fill_manual(values=year_colors) + 
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_tempo_taxo_Hill

# functional ----

FD_kelp_Hill <- as.data.frame(FD_kelp_Hill)

Ftemporal_kelp_Hill<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(FD_kelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, Habitat, FD_q0, FD_q1, FD_q2)

FD_nokelp_Hill <- as.data.frame(FD_nokelp_Hill)

Ftemporal_nokelp_Hill<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(FD_nokelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, Habitat, FD_q0, FD_q1, FD_q2)

# merging values for the 2 habitat types
Ftemporal_all_Hill <- bind_rows(Ftemporal_kelp_Hill, Ftemporal_nokelp_Hill )


# mean and sd of diversity among each site for each year in each habitat type
Ftemporal_toplot_Hill <- Ftemporal_all_Hill %>%
  group_by(Year, Habitat) %>%
  summarise( 
    n = n(),
    FD_q0_mean = mean(FD_q0),
    FD_q0_sd = sd(FD_q0)
  ) %>%
  mutate( FD_q0_se = FD_q0_sd/sqrt(n))

unique(Ftemporal_toplot_Hill$Year)

plot_tempo_func_Hill <- ggplot(Ftemporal_toplot_Hill, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=FD_q0_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin=FD_q0_mean-FD_q0_se, ymax=FD_q0_mean+FD_q0_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(4.5,6.5), breaks = seq(from=4.5, to=6.5, by=0.5)  ) +
  scale_color_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  scale_fill_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  labs(x="", y="Hill Numbers FD (q=0)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_func_Hill



## merging all plot into a single figure and saving as png ####
figure1a_Hill <- ( plot_spatial_taxo_Hill + plot_tempo_taxo_Hill ) / ( plot_spatial_func_Hill +  plot_tempo_func_Hill )
ggsave(figure1a_Hill, file=here::here("outputs/", "figure1a_Hill.png"),
       height = 22, width = 25, unit = "cm" )



######### Figure 1b ####################

plot_spatial_taxo_b <- ggplot(TDspatial_Hill_toplot) +
  geom_bar( aes(x=Habitat, y=FD_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=FD_mean-FD_se, ymax=FD_mean+FD_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,25), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_taxo_b

plot_tempo_taxo_b <- ggplot(temporal_toplot_Hill, aes(x=Year, y=FD_q0_mean, fill=Habitat)) +
  geom_bar(stat="identity", position="dodge", color = "black", size=0.8) +
  geom_errorbar( aes(x=Year, ymin=FD_q0_mean-FD_q0_se, ymax=FD_q0_mean+FD_q0_se), width=0.1, size=0.8, 
                 position=position_dodge(.9), colour="black") +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,37), breaks = seq(from=0, to=35, by=5)  ) +
  scale_color_manual(values=year_colors) + 
  scale_fill_manual(values=year_colors) + 
  labs(x="", y="") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_tempo_taxo_b


figure1b_Hill <- ( plot_spatial_taxo_b + plot_tempo_taxo_b )
ggsave(figure1b_Hill, file=here::here("outputs/", "figure1b_Hill.png"),
       height = 22, width = 25, unit = "cm" )

############################################### Figure 1c-1d ############################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

# loading data


load(here::here("outputs/", "FD_nokelp_Hill.RData") )
load(here::here("data", "nokelp_metadata.RData"))

load(here::here("outputs/", "FD_kelp_Hill.RData"))
load(here::here("data", "kelp_metadata.RData"))

load(here::here("outputs/", "FD_spatial_Hill.RData"))
load(here::here("data", "spatial_metadata.RData"))


## Functional spatial trends ####

FD_spatial <- as.data.frame(FD_spatial)

Fx_spatial_Hill <- spatial_metadata %>% 
  left_join( rownames_to_column(FD_spatial, "Code"), by="Code" ) %>%
  select(Code, Habitat, "q=0"=FD_q0, "q=1"=FD_q1, "q=2"=FD_q2)

Fx_spatial_Hill <- gather(Fx_spatial_Hill, q, qvalue, 3:5)

# mean and sd of diversity among each site in each habitat type

FDspatial_Hill_toplot_c <- Fx_spatial_Hill %>%
  group_by(Habitat, q) %>%
  summarise( 
    mean.value= mean(qvalue),
    sd.value=sd(qvalue), count=n(),
    se.mean=sd.value/sqrt(count))

# color code for the 3 habitats
hab_colors <- c(Inshore= "mediumseagreen", Midshelf= "lightsalmon1", Offshore="firebrick3")

plot_spatial_func_Hill <- ggplot(FDspatial_Hill_toplot_c) +
  geom_bar( aes(x= Habitat, y=mean.value, color = Habitat, fill = Habitat), color = "black", stat="identity", 
            size=0.8, position = position_dodge()) +
  geom_errorbar( aes( x= Habitat, ymin = mean.value - sd.value, ymax = mean.value + sd.value), 
                 width=0.1, size=0.8) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  facet_wrap(.~q, scales="free")+
  #scale_y_continuous( limits = c(3.7,4.5), breaks = seq(from=3.7, to=4.5, by=0.1)  ) +
  labs(x="", y="Hill numbers FD") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_func_Hill



# functional ----

FD_kelp_Hill <- as.data.frame(FD_kelp_Hill)

Ftemporal_kelp_Hill<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(FD_kelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, Habitat, "q=0"=FD_q0, "q=1"=FD_q1, "q=2"=FD_q2)

FD_nokelp_Hill <- as.data.frame(FD_nokelp_Hill)

Ftemporal_nokelp_Hill<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(FD_nokelp_Hill, "Code"), by="Code" ) %>%
  select(Site, Year, Habitat, "q=0"=FD_q0, "q=1"=FD_q1, "q=2"=FD_q2)

# merging values for the 2 habitat types
Ftemporal_all_Hill <- bind_rows(Ftemporal_kelp_Hill, Ftemporal_nokelp_Hill )

Ftemporal_all_Hill <- gather(Ftemporal_all_Hill, q, qvalue, 4:6)


# mean and sd of diversity among each site for each year in each habitat type
Ftemporal_toplot_Hill_c <- Ftemporal_all_Hill %>%
  group_by(Year, Habitat, q) %>%
  summarise( 
    mean.value= mean(qvalue),
    sd.value=sd(qvalue), count=n(),
    se.mean=sd.value/sqrt(count))

unique(Ftemporal_toplot_Hill_c$Year)

# color code for kelp/no kelp

year_colors <- c(kelp= "mediumseagreen", no_kelp= "lightsalmon1")

plot_tempo_func_Hill_c <- ggplot(Ftemporal_toplot_Hill_c, 
                               mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=mean.value), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin = mean.value - sd.value, ymax = mean.value + sd.value), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(4.5,6.5), breaks = seq(from=4.5, to=6.5, by=0.5)  ) +
  scale_color_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  scale_fill_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  facet_wrap(.~q, scales="free")+
  labs(x="", y="Hill Numbers FD") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_tempo_func_Hill_c



## merging all plot into a single figure and saving as png ####
figure1c_Hill <- ( plot_spatial_func_Hill + plot_tempo_func_Hill_c )
ggsave(figure1c_Hill, file=here::here("outputs/", "figure1c_Hill.png"),
       height = 22, width = 25, unit = "cm" )
