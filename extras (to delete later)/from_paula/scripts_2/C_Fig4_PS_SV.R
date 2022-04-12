################################################################################
##
## Script for computing functional identity
## between years for kelp and no kelp sites
## 
## Code by Paula Sgarlatta, Camille Magneville and Sébastien Villéger
##
################################################################################

rm(list=ls()) # cleaning memory


# libraries
library(reshape2)
library(tidyverse)
library(mFD)

# loading data

load(here::here("outputs", "temporal_alpha_kelp.RData") )
load(here::here("outputs", "spatial_alpha.RData") )



#######kelp - Fide on the 3 PC axes

fide3_kelp <- temporal_alpha_kelp %>% 
  select(fide_PC1,     fide_PC2, fide_PC3) 

shift3D_kelp_eucl <- dist(fide3_kelp, method = "euclidean")

shift3D_kelp_eucl_1 <- dist.to.df( list(shift3D=shift3D_kelp_eucl) ) %>% 
  as.data.frame() %>% 
  mutate(Site1=sub("_.*", "", x1), Site2=sub("_.*", "", x2),
         Year1=sub(".*_", "", x1),  Year2=sub(".*_", "", x2))


df_shift3D_kelp_toplot <- shift3D_kelp_eucl_1 %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, shift3D)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             shift3D_mean = mean(shift3D),
             shift3D_sd = sd(shift3D)
  ) %>%
  mutate(shift3D_se = shift3D_sd/sqrt(n))  %>%
  mutate( Habitat = "Kelp", .before="Year" )


shift3D_kelp_stats <- shift3D_kelp_eucl_1 %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(shift3D, Site1, Site2, Year)

#### Spatial data - Fide on the 3 PC axes

fide3_space <- spatial_alpha %>% 
  select(fide_PC1,     fide_PC2, fide_PC3) 

shift3D_space_eucl <- dist(fide3_space, method = "euclidean")


shift3D_space_eucl_1 <- dist.to.df( list(shift3D=shift3D_space_eucl) ) %>% 
  as.data.frame() %>% 
  mutate(Site1=sub(".*_", "", x1), Site2=sub(".*_", "", x2))


shift3D_space_stats <- dist.to.df( list(shift3D=shift3D_space_eucl) ) %>% 
  as.data.frame() %>% 
  mutate(Habitat1=sub(".*_", "", x1), Habitat2=sub(".*_", "", x2), Site1=sub("_.*", "", x1), 
         Site2=sub("_.*", "", x2)) %>% 
  filter(Habitat1 != Habitat2 ) %>% 
  mutate(Habitat = paste(Habitat1, Habitat2, sep="_"), Site = paste(Site1, Site2, sep="_")) %>% 
  select(shift3D, Habitat, Site)



df_shift3D_space_toplot <- shift3D_space_eucl_1 %>% 
  group_by(Site1, Site2) %>%
  summarise( n = n(),
             shift3D_mean = mean(shift3D),
             shift3D_sd = sd(shift3D)
  ) %>%
  ungroup(Site1, Site2) %>%
  mutate(shift3D_se = shift3D_sd/sqrt(n))  %>%
  mutate( Habitat = paste(Site1, Site2, sep="_"), .before="Site1" ) %>%
  select(-Site1, -Site2) %>%
  mutate( Year = 2013, .before="Year" ) %>% 
  filter(Habitat %in%c("I_M", "I_O", "M_O") ) %>%
  droplevels("Habitat")

df_shift3D_space_toplot$Year <- as.character(df_shift3D_space_toplot$Year)

# merging ----

df_shift3D <- bind_rows(df_shift3D_kelp_toplot, df_shift3D_space_toplot)

# color code for the 3 habitats

hab_colors <- c(Kelp="#3BB372",
                I_M="#2C6BAA", I_O="lightsalmon1", M_O="firebrick3")

df_shift3D$Year <- as.numeric(df_shift3D$Year)

## Plotting

plot_shift3D <- ggplot(df_shift3D, mapping = aes(color=Habitat)) +
  geom_point( aes( x=Year, y= shift3D_mean), stat="identity", size=2) +
  geom_line( aes( x=Year, y= shift3D_mean), stat="identity", size=1) +
  geom_errorbar( aes(x=Year, ymin= shift3D_mean - shift3D_se, ymax= shift3D_mean + shift3D_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors,name="Habitat", breaks = c("Kelp", "I_M", "I_O", "M_O"), 
                     labels=c("Kelp", "Inshore-Midshelf", "Inshore-Offshore", "Midshelf-Offshore")) + 
  scale_fill_manual(values=hab_colors, name="Habitat", breaks = c("Kelp", "I_M", "I_O", "M_O"), 
                    labels=c("Kelp", "Inshore-Midshelf", "Inshore-Offshore", "Midshelf-Offshore")) + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0,0.1), breaks = seq(from=0, to=0.1, by=0.02)  ) +
  labs(x="", y="Shift in Functional Identity") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_shift3D


## saving figure as png ####

ggsave(plot_shift3D, file=here::here("outputs/", "Figure4_bis.png"), 
       height = 20, width = 18, unit = "cm" )

#save for stats

save(shift3D_kelp_stats , file=here::here("outputs/", "shift3D_kelp_stats.RData") )
save(shift3D_space_stats , file=here::here("outputs/", "shift3D_space_stats.RData") )

################################################ end of script #######################################################

