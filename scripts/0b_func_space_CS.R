################################################################################
##
## Script for computing species positions in a multidimensional space
## according to their trait values
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)

## loading ####

# load traits data ----
fish_traits <- read.csv(here::here("data", "raw_data", "fish_traits.csv"), header = T)

# load species names from surveys datasets ----
load(here::here("data", "species_allsurveys.RData") )
length(species_allsurveys) # 139 species    

# checking same species in trait and occurrences datasets ----
identical ( sort(species_allsurveys) , sort(fish_traits$Species ) ) # True

## preparing trait dataset ####

# trait values in a dataframe (species in alphabetical order) ----
sp_tr <- fish_traits %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()
head(sp_tr)

nrow(sp_tr) # 139 species

# recoding variable to match trait type ----

# looking at trait values
lapply(sp_tr, unique)

# trait type
tr_cat<-data.frame( trait_name = c("Size", "Agg", "Position", "Diet", "Kmax"),
                    trait_type = c("O","O","O", "N", "Q") )

# size as ordinal
sp_tr$Size <- factor(sp_tr$Size, 
                     levels = c("S1", "S2", "S3", "S4", "S5", "S6"),
                     ordered = TRUE)
summary(sp_tr$Size)

# aggregation as ordinal
sp_tr$Agg <- factor(sp_tr$Agg, 
                    levels = c("Solitary", "Pair", "Group"),
                    ordered = TRUE)
summary(sp_tr$Agg)

# Position as ordinal
sp_tr$Position <- factor(sp_tr$Position, 
                    levels = c("Benthic", "BenthoP", "Pelagic"),
                    ordered = TRUE)
summary(sp_tr$Position)

# diet as factor
sp_tr$Diet <- as.factor(sp_tr$Diet)
summary(sp_tr$Diet)

#Kmax as numeric

sp_tr$Kmax <- as.numeric(sp_tr$Kmax)
summary(sp_tr$Kmax)

# summary of trait data----
summary_traits <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                      sp_tr  = sp_tr)



## Computing Gower distance between species ####
sp_gower_dist <- mFD::funct.dist(sp_tr=sp_tr, tr_cat = tr_cat, 
                                   metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist) # from 0 to 1

### Compute functional spaces and their quality:

# mean absolute deviation index (as quality metric)
funct_spaces<- mFD::quality.fspaces(sp_dist = sp_gower_dist, maxdim_pcoa = 12, 
                                 deviation_weighting = "absolute", fdist_scaling = FALSE,
                                 fdendro = "average") 
funct_spaces$quality_fspaces
# => 3D space has the lowest mAD (0.055)

# illustrate quality of functional space

funct_space <- mFD::quality.fspaces.plot(
  fspaces_quality            = funct_spaces,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

funct_space

# species coordinates
sp_3D_coord<-funct_spaces$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord)

### Test correlation between traits and functional axes:

# retrieve coordinates of species:

sp_faxes_coord <- funct_spaces$"details_fspaces"$"sp_pc_coord"

# test correlation between traits and axes:
cor_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sp_tr, 
  sp_faxes_coord = sp_faxes_coord[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# get the table of correlation:

cor_tr_faxes$tr_faxes_stat

# get the plot:

cor_tr_faxes$tr_faxes_plot

## Adding thermal affinity

thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  #column_to_rownames("Species") %>% 
  select(-thermal)

#Add thermal aff to sp_faxes_coord

sp_faxes_coord <- as.data.frame(sp_faxes_coord) %>% 
  rownames_to_column("Species")

sp_faxes_coord <- inner_join(sp_faxes_coord, thermal, 
                             by="Species") %>%
  column_to_rownames("Species") %>% 
  as.matrix()

### Plot functional space:

big_plot <- mFD::funct.space.plot(sp_faxes_coord  = sp_faxes_coord[, c("PC1", "PC2", "PC3")],
                                       faxes = c("PC1", "PC2", "PC3"), name_file = NULL,
                                       faxes_nm = NULL, range_faxes = c(NA, NA),
                                       color_bg = "grey95",
                                       color_pool = "darkgreen", fill_pool = "white",
                                       shape_pool = 21, size_pool = 1,
                                       plot_ch = TRUE, color_ch = "black",
                                       fill_ch = "white", alpha_ch = 0.3,
                                       plot_vertices = TRUE, 
                                       color_vert = "blueviolet",
                                       fill_vert = "blueviolet", shape_vert = 23,
                                       size_vert = 1,
                                       plot_sp_nm = NULL, nm_size = 3, 
                                       nm_color = "black",
                                       nm_fontface = "plain",
                                       check_input = TRUE)
# Plot the graph with all pairs of axes:
big_plot$patchwork

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr, file=here::here("data/", "sp_tr.RData") )
save(tr_cat, file=here::here("data/", "tr_cat.RData") )
save(summary_traits, file=here::here("outputs/", "summary_traits.RData") )
save(sp_gower_dist, file=here::here("outputs/", "sp_gower_dist.RData") )
save(sp_3D_coord, file=here::here("outputs/", "sp_3D_coord.RData") )


##################################  end of code ######################################################################



# @@@ ADD code for plots of funct space, correl tr vs axes and save them in outputs
