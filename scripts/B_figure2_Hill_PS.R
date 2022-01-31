############################################################################################################
##
## Script for plotting taxonomic and functional dissimilarity across habitats and years using Hill numbers
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
###########################################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)


# loading data

load(here::here("outputs/", "FD_beta_nokelp_sites.RData") )
load(here::here("outputs/", "TD_beta_nokelp_sites.RData") )


load(here::here("outputs/", "FD_beta_kelp_Hill_sites.RData") )
load(here::here("outputs/", "TD_beta_kelp_Hill_sites.RData") )


load(here::here("outputs/", "FD_beta_spatial_Hill.RData") )
load(here::here("outputs/", "TD_beta_spatial_Hill.RData") )



## taxonomic dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat, but only 
#between  the same sites ----

kelp_interyear_taxo_diss_Hill <- TD_beta_kelp_Hill_sites %>%
  filter((Year1==2002 | Year2==2002), (Site1==Site2)) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )

## Keep this one for now, but will be deleted if the one above works. 

kelp_interyear_taxo_diss_Hill <- TD_beta_kelp %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----
nokelp_interyear_taxo_diss_Hill <- TD_beta_nokelp_sites %>%
  filter((Year1==2002 | Year2==2002), (Site1 == Site2)) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "No_kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----

interhab_taxo_diss_Hill <- TD_beta_spatial %>% 
  filter(Habitat1==Habitat2) %>% 
  mutate(Habitat = Habitat1) %>% 
  group_by(Habitat) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate( q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  droplevels("Habitat")


## Keep this one for now, but will be deleted if the one above works. 
interhab_taxo_diss_Hill <- TD_beta_spatial %>% 
  group_by(Habitat1, Habitat2) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  ungroup(Habitat1, Habitat2) %>%
  mutate( q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = paste(Habitat1, Habitat2, sep="_"), .before="Habitat1" ) %>%
  select(-Habitat1, -Habitat2) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  filter(Habitat %in%c("I_M", "M_O", "I_O") ) %>%
  droplevels("Habitat")


# merging ----

kelp_interyear_taxo_diss_Hill$Year <- as.double(kelp_interyear_taxo_diss_Hill$Year)
nokelp_interyear_taxo_diss_Hill$Year <- as.double(nokelp_interyear_taxo_diss_Hill$Year)


TD_dissim_Hill <- bind_rows(kelp_interyear_taxo_diss_Hill, 
                             nokelp_interyear_taxo_diss_Hill,
                             interhab_taxo_diss_Hill) 

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="seagreen4", No_kelp="seagreen2",
                I="#2C6BAA", M="lightsalmon1", O="firebrick3")
#hab_shape <- c(Kelp=21, No_kelp=21,
               #I_M=22, M_O=22, I_O=22)


plot_dissim_taxo_Hill <- ggplot(TD_dissim_Hill, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point(aes(x=Year, y=q1_mean), stat="identity", size=3) +
  geom_line(aes(x=Year, y=q1_mean), stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+q1_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  #scale_shape_manual(values=hab_shape) + 
  scale_x_continuous( limits = c(2002, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.2,0.6), breaks = seq(from=0.2, to=0.6, by=0.1)  ) +
  labs(x="", y="Taxonomic dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (14)),
        legend.position = "none")
plot_dissim_taxo_Hill


## functional dissimilarity ####

# values between assemblages from 2002 and later in kelp habitat ----

kelp_interyear_func_diss_Hill <- FD_beta_kelp_Hill_sites %>%
  filter((Year1==2002 | Year2==2002), (Site1==Site2)) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )

## Keep this one for now, but will be deleted if the one above works. 
kelp_interyear_func_diss_Hill <- FD_beta_kelp %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "Kelp", .before="Year" )



# values between assemblages from 2002 and later in no kelp habitat ----

nokelp_interyear_func_diss_Hill <- FD_beta_nokelp_sites %>%
  filter((Year1==2002 | Year2==2002), (Site1 == Site2)) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "No_kelp", .before="Year" )


## Keep this one for now, but will be deleted if the one above works.
nokelp_interyear_func_diss_Hill <- FD_beta_nokelp %>%
  filter(Year1==2002 | Year2==2002) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, q0, q1, q2)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate(q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = "No_kelp", .before="Year" )

# adding dissimilarity between the 3 habitats (from 2013) ----

interhab_func_diss_Hill <- FD_beta_spatial %>% 
  filter(Habitat1==Habitat2) %>% 
  mutate(Habitat = Habitat1) %>% 
  group_by(Habitat) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  mutate( q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  droplevels("Habitat")


## Keep this one for now, but will be deleted if the one above works.

interhab_func_diss_Hill <- FD_beta_spatial %>% 
  group_by(Habitat1, Habitat2) %>%
  summarise( n = n(),
             q0_mean = mean(q0),
             q0_sd = sd(q0),
             q1_mean = mean(q1),
             q1_sd = sd(q1),
             q2_mean = mean(q2),
             q2_sd = sd(q2)
  ) %>%
  ungroup(Habitat1, Habitat2) %>%
  mutate( q0_se = q0_sd/sqrt(n), q1_se = q1_sd/sqrt(n), q2_se = q2_sd/sqrt(n)) %>%
  mutate( Habitat = paste(Habitat1, Habitat2, sep="_"), .before="Habitat1" ) %>%
  select(-Habitat1, -Habitat2) %>%
  mutate( Year = 2013, .before="Year" ) %>%
  filter(Habitat %in%c("I_M", "M_O", "I_O") ) %>%
  droplevels("Habitat")


# merging ----

kelp_interyear_func_diss_Hill$Year <- as.double(kelp_interyear_func_diss_Hill$Year)
nokelp_interyear_func_diss_Hill$Year <- as.double(nokelp_interyear_func_diss_Hill$Year)


FD_dissim_Hill <- bind_rows(kelp_interyear_func_diss_Hill, 
                            nokelp_interyear_func_diss_Hill,
                            interhab_func_diss_Hill) 

# plotting ----

# color code for the 3 habitats
hab_colors <- c(Kelp="#3BB372", No_kelp="#74E7B8",
                I="#2C6BAA", M="lightsalmon1", O="firebrick3")
#hab_shape <- c(Kelp=21, No_kelp=21,
               #I_M=22, M_O=22, I_O=22)

#highlight <- FD_dissim_Hill %>% 
  #filter(Year == 2002)

plot_dissim_func_Hill <- ggplot(FD_dissim_Hill, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point(aes(x=Year, y=q1_mean), stat="identity", size=2) +
  geom_line(aes(x=Year, y=q1_mean), stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+q1_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", 
    "No_kelp", "I", "M","O"),labels=c("Kelp", "No kelp", "Inshore","Midshelf",
                                      "Offshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", 
    "No_kelp", "I", "M","O"),labels=c("Kelp", "No kelp", "Inshore","Midshelf",
                                      "Offshore")) + 
  #scale_shape_manual(values=hab_shape, name="Habitat", breaks = c("Kelp", "No_kelp", "I_M", "M_O", 
                                                                 # "I_O"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_x_continuous( limits = c(2002, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,0.1), breaks = seq(from=0, to=0.1, by=0.05)  ) +
  labs(x="", y="Functional dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_dissim_func_Hill



#Probably delete this later

plot_dissim_func_Hill <- ggplot(FD_dissim_Hill, mapping = aes(color=Habitat, fill=Habitat)) +
  geom_point(aes(x=Year, y=q1_mean, 
                 shape= Habitat), stat="identity", size=3) +
  geom_point(data=highlight, aes(x=Year, y=q1_mean), stat="identity", color="red", size = 3)+
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+q1_se), 
                 width=0.2, size=0.5 ) +
    geom_errorbar( data=highlight,aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+q1_se), 
                   width=0.2, size=0.5, color="red" ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I_M", "M_O", 
  "I_O"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I_M", "M_O", 
  "I_O"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_shape_manual(values=hab_shape, name="Habitat", breaks = c("Kelp", "No_kelp", "I_M", "M_O", 
  "I_O"),labels=c("Kelp", "No kelp", "Midshelf-Inshore","Offshore-Midshelf", "Offshore-Inshore")) + 
  scale_x_continuous( limits = c(2002, 2018), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,0.1), breaks = seq(from=0, to=0.1, by=0.05)  ) +
  labs(x="", y="Functional dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_dissim_func_Hill #Need to change the legend colors

#####
## merging all plot into a single figure and saving as png ####
figure2_Hill <- ( plot_dissim_taxo_Hill / plot_dissim_func_Hill )
ggsave(figure2_Hill, file=here::here("outputs/", "figure2_Hill_site_average.png"), 
       height = 20, width = 18, unit = "cm" )


