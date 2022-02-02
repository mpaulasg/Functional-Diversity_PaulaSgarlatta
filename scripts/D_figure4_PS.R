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

# loading data
load(here::here("outputs", "temporal_alpha_nokelp.RData") )
load(here::here("outputs", "temporal_alpha_kelp.RData") )


#######No kelp

fide3_nokelp <- temporal_alpha_nokelp %>% 
  select(fide_PC3) 

fide3_nokelp_eucl <- dist(fide3_nokelp, method = "euclidean")

# Then use the mFD::dist.to.df function to ease visualizing result - don't know why is not working

fide_nokelp_eucl_1 <- dist.to.df(fide_nokelp_eucl) %>% 
  as.data.frame() %>% 
  mutate(Yearx=sub("-.*", "", x1), Yeary=sub("-.*", "", x2),
         Year1=sub(".*_", "", Yearx),  Year2=sub(".*_", "", Yeary))

### Using melt instead

df_fide3_nokelp <- melt(as.matrix(fide3_nokelp_eucl), varnames = c("SiteA", "SiteB"))

df_fide3_nokelp <- df_fide3_nokelp %>% 
  mutate(Year1=sub(".*_", "", SiteA), Year2=sub(".*_", "", SiteB), Site1=sub("_.*", "", SiteA), Site2=sub("_.*", "", SiteB)) %>% 
  select(PC3=value, Year1, Year2, Site1, Site2)

df_fide_nokelp_toplot <- df_fide3_nokelp %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, PC3)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             PC3_mean = mean(PC3),
             PC3_sd = sd(PC3)
  ) %>%
  mutate(PC3_se = PC3_sd/sqrt(n))  %>%
  mutate( Habitat = "No_Kelp", .before="Year" )


df_fide_nokelp <- cbind(df_fide1_nokelp, df_fide2_nokelp, df_fide3_nokelp) ## Not the correct one 

############# Kelp


fide3_kelp <- temporal_alpha_kelp %>% 
  select(fide_PC3) 

fide3_kelp_eucl <- dist(fide3_kelp, method = "euclidean")

df_fide3_kelp <- melt(as.matrix(fide3_kelp_eucl), varnames = c("SiteA", "SiteB"))

df_fide3_kelp <- df_fide3_kelp %>% 
  mutate(Year1=sub(".*_", "", SiteA), Year2=sub(".*_", "", SiteB), Site1=sub("_.*", "", SiteA), Site2=sub("_.*", "", SiteB)) %>% 
  select(PC3=value, Year1, Year2, Site1, Site2)

df_fide_kelp_toplot <- df_fide3_kelp %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, PC3)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             PC3_mean = mean(PC3),
             PC3_sd = sd(PC3)
  ) %>%
  mutate(PC3_se = PC3_sd/sqrt(n))  %>%
  mutate( Habitat = "Kelp", .before="Year" )

# merging ----

df_fide <- bind_rows(df_fide_nokelp_toplot, 
                     df_fide_kelp_toplot)

# color code for the 3 habitats
hab_colors <- c(Kelp="#3BB372", No_Kelp="#74E7B8")

df_fide$Year <- as.numeric(df_fide$Year)

## Plotting

plot_fide3 <- ggplot(df_fide, mapping = aes(color=Habitat)) +
  geom_point( aes( x=Year, y= PC3_mean), stat="identity", size=2) +
  geom_line( aes( x=Year, y= PC3_mean), stat="identity", size=1) +
  geom_errorbar( aes(x=Year, ymin= PC3_mean- PC3_se, ymax= PC3_mean+ PC3_se), 
                 width=0.2, size=1 ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  #scale_y_continuous( limits = c(0,0.6), breaks = seq(from=0, to=0.6, by=0.2)  ) +
  labs(x="", y="Functional Identity PC3") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())

plot_fide3


#####
## saving figure as png ####
figure2_convex_site_average <- ( plot_dissim_taxo / plot_dissim_func )
ggsave(figure2_convex_site_average, file=here::here("outputs/", "figure2_convex_site_average.png"), 
       height = 20, width = 18, unit = "cm" )

################################################ end of script #######################################################









######### No kelp

fide1_nokelp <- temporal_alpha_nokelp %>% 
  select(fide_PC1) 

fide1_nokelp_eucl <- dist(fide1_nokelp, method = "euclidean")

df_fide1_nokelp <- melt(as.matrix(fide1_nokelp_eucl), varnames = c("SiteA", "SiteB"))

df_fide1_nokelp <- df_fide1_nokelp %>% 
  mutate(Year1=sub(".*_", "", SiteA), Year2=sub(".*_", "", SiteB), Site1=sub("_.*", "", SiteA), Site2=sub("_.*", "", SiteB)) %>% 
  select(PC1=value, Year1, Year2, Site1, Site2)

fide2_nokelp <- temporal_alpha_nokelp %>% 
  select(fide_PC2) 

fide2_nokelp_eucl <- dist(fide2_nokelp, method = "euclidean")

df_fide2_nokelp <- melt(as.matrix(fide2_nokelp_eucl), varnames = c("SiteA", "SiteB"))

df_fide2_nokelp <- df_fide2_nokelp %>% 
  mutate(Year1=sub(".*_", "", SiteA), Year2=sub(".*_", "", SiteB), Site1=sub("_.*", "", SiteA), Site2=sub("_.*", "", SiteB)) %>% 
  select(PC2=value, Year1, Year2, Site1, Site2)


######### Kelp

fide1_kelp <- temporal_alpha_kelp %>% 
  select(fide_PC1) 

fide1_kelp_eucl <- dist(fide1_kelp, method = "euclidean")

df_fide1_kelp <- melt(as.matrix(fide1_kelp_eucl), varnames = c("SiteA", "SiteB"))

df_fide1_kelp <- df_fide1_kelp %>% 
  mutate(Year1=sub(".*_", "", SiteA), Year2=sub(".*_", "", SiteB), Site1=sub("_.*", "", SiteA), Site2=sub("_.*", "", SiteB)) %>% 
  select(PC1=value, Year1, Year2, Site1, Site2)

fide2_kelp <- temporal_alpha_kelp %>% 
  select(fide_PC2) 

fide2_kelp_eucl <- dist(fide2_kelp, method = "euclidean")

df_fide2_kelp <- melt(as.matrix(fide2_kelp_eucl), varnames = c("SiteA", "SiteB"))

df_fide2_kelp <- df_fide2_kelp %>% 
  mutate(Year1=sub(".*_", "", SiteA), Year2=sub(".*_", "", SiteB), Site1=sub("_.*", "", SiteA), Site2=sub("_.*", "", SiteB)) %>% 
  select(PC2=value, Year1, Year2, Site1, Site2)

