############################################################################################################
##
## Script for plotting temporal homogenisation within each habitat
## 
## Code by Paula Sgarlatta, Sebastien Villeger and Camille Magneville 
##
###########################################################################################################
rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

# loading data

load(here::here("outputs/", "FD_beta_nokelp_Hill.RData") )
load(here::here("outputs/", "TD_beta_nokelp_Hill.RData") )


load(here::here("outputs/", "FD_beta_kelp_Hill.RData") )
load(here::here("outputs/", "TD_beta_kelp_Hill.RData") )


load(here::here("outputs/", "FD_beta_spatial_Hill.RData") )
load(here::here("outputs/", "TD_beta_spatial_Hill.RData") )


## taxonomic dissimilarity ####

# Add habitat 

TD_beta_kelp$Habitat <- "Kelp"
TD_beta_nokelp$Habitat <- "No_kelp"

# Add o column "Year" only when Year 1 and Year 2 are the same

TD_beta_kelp$Year <- ifelse(TD_beta_kelp$Year1 == TD_beta_kelp$Year2, 
                            TD_beta_kelp$Year1, NA)
TD_beta_nokelp$Year <- ifelse(TD_beta_nokelp$Year1 == TD_beta_nokelp$Year2, 
                              TD_beta_nokelp$Year1, NA)
#Omit NA's

TD_beta_kelp_done <- na.omit(TD_beta_kelp)%>% 
  select(Year, q0, q1, q2, Habitat)

TD_beta_nokelp_done <- na.omit(TD_beta_nokelp) %>% 
  select(Year, q0, q1, q2, Habitat)


## Adding spatial data

TD_beta_spatial$Habitat <- ifelse(TD_beta_spatial$Habitat1 == TD_beta_spatial$Habitat2,
                                  TD_beta_spatial$Habitat1, NA)
TD_beta_spatial$Year <- "2013"

TD_beta_spatial_done <- na.omit(TD_beta_spatial)%>% 
  select(Year, q0, q1, q2, Habitat)



#Join datasets

TD_beta_toplot <- rbind(TD_beta_kelp_done, 
                                 TD_beta_nokelp_done, TD_beta_spatial_done )



# Obtain average per year and habitat

taxo_diss_Hill <- TD_beta_toplot %>%
  select(Year, Habitat, q1) %>% 
  group_by(Year, Habitat) %>%
  summarise( n = n(),
             q1_mean = mean(q1),
             q1_sd = sd(q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))

# color code for the 3 habitats

hab_colors <- c(Kelp="seagreen4", No_kelp="seagreen2",
                I="blue", M="lightsalmon1", O="firebrick3")
hab_shape <- c(Kelp=15, No_kelp=15,
               I=25, M=25, O=25)

taxo_diss_Hill$Year <- as.numeric(taxo_diss_Hill$Year)

#plotting

plot_dissim_taxo_Hill <- ggplot(taxo_diss_Hill, 
                                mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=q1_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+ q1_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.2, 0.6), breaks = seq(from=0.2, to=0.6, by=0.1)  ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"), 
                     labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"), 
                    labels=c("Kelp", "No kelp", "Inshore","Midshelf", "Offshore")) + 
  scale_shape_manual(values=hab_shape, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"),
                     labels=c("Kelp", "No kelp", "Inshore","Midshelf", "Offshore"))+
  labs(x="", y="Taxonomic dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_dissim_taxo_Hill


## Functional dissimilarity ####

# Add habitat 

FD_beta_kelp$Habitat <- "Kelp"
FD_beta_nokelp$Habitat <- "No_kelp"

# Add o column "Year" only when Year 1 and Year 2 are the same

FD_beta_kelp$Year <- ifelse(FD_beta_kelp$Year1 == FD_beta_kelp$Year2, 
                            FD_beta_kelp$Year1, NA)
FD_beta_nokelp$Year <- ifelse(FD_beta_nokelp$Year1 == FD_beta_nokelp$Year2, 
                              FD_beta_nokelp$Year1, NA)
#Omit NA's

FD_beta_kelp_done <- na.omit(FD_beta_kelp)%>% 
  select(Year, q0, q1, q2, Habitat)

FD_beta_nokelp_done <- na.omit(FD_beta_nokelp)%>% 
  select(Year, q0, q1, q2, Habitat)

## Adding spatial data

FD_beta_spatial$Habitat <- ifelse(FD_beta_spatial$Habitat1 == FD_beta_spatial$Habitat2,
                                  FD_beta_spatial$Habitat1, NA)
FD_beta_spatial$Year <- "2013"

FD_beta_spatial_done <- na.omit(FD_beta_spatial)%>% 
  select(Year, q0, q1, q2, Habitat)


#Join datasets

FD_beta_toplot <- rbind(FD_beta_kelp_done, 
                      FD_beta_nokelp_done, FD_beta_spatial_done )

# Obtain average per year and habitat

func_diss_Hill <- FD_beta_toplot %>%
  select(Year, Habitat, q1) %>% 
  group_by(Year, Habitat) %>%
  summarise( n = n(),
             q1_mean = mean(q1),
             q1_sd = sd(q1)
  ) %>%
  mutate(q1_se = q1_sd/sqrt(n))

# color code for the 3 habitats

hab_colors <- c(Kelp="seagreen4", No_kelp="seagreen2",
                I="blue", M="lightsalmon1", O="firebrick3")
hab_shape <- c(Kelp=15, No_kelp=15,
               I=25, M=25, O=25)

func_diss_Hill$Year <- as.numeric(func_diss_Hill$Year)

#plotting

plot_dissim_func_Hill <- ggplot(func_diss_Hill, 
                                mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=q1_mean), stat="identity", size=2, shape=16) +
  geom_errorbar( aes(x=Year, ymin=q1_mean-q1_se, ymax=q1_mean+ q1_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0, 0.08), breaks = seq(from=0.02, to=0.08, by=0.02)  ) +
  scale_color_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"), 
                     labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "No_kelp", "I", "M", "O"), 
                    labels=c("Kelp", "No kelp", "Inshore", "Midshelf", "Offshore")) + 
  labs(x="", y="Functional dissimilarity (q=1)") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_dissim_func_Hill


## merging all plot into a single figure and saving as png ####
figure3_Hill <- ( plot_dissim_taxo_Hill /  plot_dissim_func_Hill )
ggsave(figure3_Hill, file=here::here("outputs/", "figure3_Hill.png"),
       height = 22, width = 25, unit = "cm" )

################################################### END OF CODE ######################################################