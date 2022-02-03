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
load(here::here("outputs", "temporal_alpha_nokelp.RData") )
load(here::here("outputs", "temporal_alpha_kelp.RData") )


#######No kelp - Fide on the 3 PC axes

fide3_nokelp <- temporal_alpha_nokelp %>% 
  select(fide_PC1,     fide_PC2, fide_PC3) 

shift3D_nokelp_eucl <- dist(fide3_nokelp, method = "euclidean")

shift3D_nokelp_eucl_1 <- dist.to.df( list(shift3D=shift3D_nokelp_eucl) ) %>% 
  as.data.frame() %>% 
  mutate(Site1=sub("_.*", "", x1), Site2=sub("_.*", "", x2),
         Year1=sub(".*_", "", x1),  Year2=sub(".*_", "", x2))

df_shift3D_nokelp_toplot <- shift3D_nokelp_eucl_1 %>% 
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
  mutate( Habitat = "No_Kelp", .before="Year" )

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


# merging ----

df_shift3D <- bind_rows(df_shift3D_nokelp_toplot, 
                     df_shift3D_kelp_toplot)

# color code for the 3 habitats
hab_colors <- c(Kelp="#3BB372", No_Kelp="#74E7B8")

df_shift3D$Year <- as.numeric(df_shift3D$Year)

## Plotting

plot_shift3D <- ggplot(df_shift3D, mapping = aes(color=Habitat)) +
  geom_point( aes( x=Year, y= shift3D_mean), stat="identity", size=2) +
  geom_line( aes( x=Year, y= shift3D_mean), stat="identity", size=1) +
  geom_errorbar( aes(x=Year, ymin= shift3D_mean - shift3D_se, ymax= shift3D_mean + shift3D_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
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

ggsave(plot_shift3D, file=here::here("outputs/", "figure4_shift3D_SV.png"), 
       height = 20, width = 18, unit = "cm" )

################################################ end of script #######################################################

