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
length(species_allsurveys) # 139 species    # [PS] -> don't understand why this step. Why is species_allsurveys needed here?

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
tr_cat<-data.frame( trait_name = c("Size", "Agg", "Position", "Diet"),
                    trait_type = c("O","O","O", "N") )

# size as ordinal
sp_tr$Size <- factor(sp_tr$Size, 
                     levels = c("S2", "S3", "S4", "S5", "S6"),
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

# summary of trait data----
summary_traits <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                     sp_tr  = sp_tr)
summary_traits


## Computing Gower distance between species ####
sp_gower_dist <- mFD::funct.dist(sp_tr=sp_tr, tr_cat = tr_cat, 
                                 metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces<- mFD::quality.fspaces(sp_dist = sp_gower_dist, maxdim_pcoa = 12, 
                                    deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces$quality_fspaces
# => 3D space has the lowest mAD (0.061)

# species coordinates
sp_3D_coord<-funct_spaces$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord)



# [PS] This is the original one. Now I'll add one separate for each dataset (kelp, no kelp, spatial)

#Loading the data (if this work, this data should be saved in raw_data in CS_00_0a_fish_traits).

sp_tr_kelp <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_kelp.csv"),
                       header = T)
sp_tr_kelp <- sp_tr_kelp %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()

# from sites that never had kelp
sp_tr_nokelp <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_no_kelp.csv"),
                         header = T)

sp_tr_nokelp <-sp_tr_nokelp %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()

# traits of species from UVC surveys
sp_tr_spatial <- read.csv(here::here("from_paula", "SpatialUVC_species_traits.csv"), 
                          header = T)

sp_tr_spatial <- sp_tr_spatial %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()

####KELP 

# recoding variable to match trait type ----

# trait type
tr_cat<-data.frame( trait_name = c("Size", "Agg", "Position", "Diet"),
                    trait_type = c("O","O","O", "N") )



# size as ordinal
sp_tr_kelp$Size <- factor(sp_tr_kelp$Size, 
                          levels = c("S2", "S3", "S4", "S5", "S6"),
                          ordered = TRUE)
summary(sp_tr_kelp$Size)

# aggregation as ordinal
sp_tr_kelp$Agg <- factor(sp_tr_kelp$Agg, 
                         levels = c("Solitary", "Pair", "Group"),
                         ordered = TRUE)
summary(sp_tr_kelp$Agg)

# Position as ordinal
sp_tr_kelp$Position <- factor(sp_tr_kelp$Position, 
                              levels = c("Benthic", "BenthoP", "Pelagic"),
                              ordered = TRUE)
summary(sp_tr_kelp$Position)

# diet as factor
sp_tr_kelp$Diet <- as.factor(sp_tr_kelp$Diet)
summary(sp_tr_kelp$Diet)

# summary of trait data----
summary_traits_kelp <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                          sp_tr  = sp_tr_kelp)
summary_traits_kelp


## Computing Gower distance between species ####
sp_gower_dist_kelp <- mFD::funct.dist(sp_tr=sp_tr_kelp, tr_cat = tr_cat, 
                                      metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist_kelp) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces_kelp<- mFD::quality.fspaces(sp_dist = sp_gower_dist_kelp, maxdim_pcoa = 12, 
                                         deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces_kelp$quality_fspaces
# => 3D space has the lowest mAD (0.059)

# species coordinates
sp_3D_coord_kelp<-funct_spaces_kelp$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord_kelp)

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr_kelp, file=here::here("data/", "sp_tr_kelp.RData") )
save(summary_traits_kelp, file=here::here("outputs/", "summary_traits_kelp.RData") )
save(sp_gower_dist_kelp, file=here::here("outputs/", "sp_gower_dist_kelp.RData") )
save(sp_3D_coord_kelp, file=here::here("outputs/", "sp_3D_coord_kelp.RData") )


####NO KELP 

# size as ordinal
sp_tr_nokelp$Size <- factor(sp_tr_nokelp$Size, 
                            levels = c("S2", "S3", "S4", "S5", "S6"),
                            ordered = TRUE)
summary(sp_tr_nokelp$Size)

# aggregation as ordinal
sp_tr_nokelp$Agg <- factor(sp_tr_nokelp$Agg, 
                           levels = c("Solitary", "Pair", "Group"),
                           ordered = TRUE)
summary(sp_tr_nokelp$Agg)

# Position as ordinal
sp_tr_nokelp$Position <- factor(sp_tr_nokelp$Position, 
                                levels = c("Benthic", "BenthoP", "Pelagic"),
                                ordered = TRUE)
summary(sp_tr_nokelp$Position)

# diet as factor
sp_tr_nokelp$Diet <- as.factor(sp_tr_nokelp$Diet)
summary(sp_tr_nokelp$Diet)

# summary of trait data----
summary_traits_nokelp <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                            sp_tr  = sp_tr_nokelp)
summary_traits_nokelp


## Computing Gower distance between species ####
sp_gower_dist_nokelp <- mFD::funct.dist(sp_tr=sp_tr_nokelp, tr_cat = tr_cat, 
                                        metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist_nokelp) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces_nokelp<- mFD::quality.fspaces(sp_dist = sp_gower_dist_nokelp, maxdim_pcoa = 12, 
                                           deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces_nokelp$quality_fspaces
# => 3D space has the lowest mAD (0.058)

# species coordinates
sp_3D_coord_nokelp<-funct_spaces_nokelp$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord_nokelp)

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr_nokelp, file=here::here("data/", "sp_tr_nokelp.RData") )
save(summary_traits_nokelp, file=here::here("outputs/", "summary_traits_nokelp.RData") )
save(sp_gower_dist_nokelp, file=here::here("outputs/", "sp_gower_dist_nokelp.RData") )
save(sp_3D_coord_nokelp, file=here::here("outputs/", "sp_3D_coord_nokelp.RData") )


####SPATIAL

# size as ordinal
sp_tr_spatial$Size <- factor(sp_tr_spatial$Size, 
                             levels = c("S2", "S3", "S4", "S5", "S6"),
                             ordered = TRUE)
summary(sp_tr_spatial$Size)

# aggregation as ordinal
sp_tr_spatial$Agg <- factor(sp_tr_spatial$Agg, 
                            levels = c("Solitary", "Group"),
                            ordered = TRUE)
summary(sp_tr_spatial$Agg)

# Position as ordinal
sp_tr_spatial$Position <- factor(sp_tr_spatial$Position, 
                                 levels = c("Benthic", "BenthoP", "Pelagic"),
                                 ordered = TRUE)
summary(sp_tr_spatial$Position)

# diet as factor
sp_tr_spatial$Diet <- as.factor(sp_tr_spatial$Diet)
summary(sp_tr_spatial$Diet)

# summary of trait data----
summary_traits_spatial <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                             sp_tr  = sp_tr_spatial)
summary_traits_spatial


## Computing Gower distance between species ####
sp_gower_dist_spatial <- mFD::funct.dist(sp_tr=sp_tr_spatial, tr_cat = tr_cat, 
                                         metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist_spatial) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces_spatial<- mFD::quality.fspaces(sp_dist = sp_gower_dist_spatial, maxdim_pcoa = 12, 
                                            deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces_spatial$quality_fspaces
# => 3D space has the lowest mAD (0.052) [PS] 4D space has the lowest, should I use 4D or keep 3D to compare all together?

# species coordinates
sp_3D_coord_spatial<-funct_spaces_spatial$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord_spatial)

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr_spatial, file=here::here("data/", "sp_tr_spatial.RData") )
save(summary_traits_spatial, file=here::here("outputs/", "summary_traits_spatial.RData") )
save(sp_gower_dist_spatial, file=here::here("outputs/", "sp_gower_dist_spatial.RData") )
save(sp_3D_coord_spatial, file=here::here("outputs/", "sp_3D_coord_spatial.RData") )

# @@@ ADD code for plots of funct space, correl tr vs axes and save them in outputs



# saving ####

# trait values and trait coding dataframes ----
save(sp_tr, file=here::here("data/", "sp_tr.RData") )
save(tr_cat, file=here::here("data/", "tr_cat.RData") )
save(summary_traits, file=here::here("outputs/", "summary_traits.RData") )
save(sp_gower_dist, file=here::here("outputs/", "sp_gower_dist.RData") )
save(sp_3D_coord, file=here::here("outputs/", "sp_3D_coord.RData") )


####################### end ####


