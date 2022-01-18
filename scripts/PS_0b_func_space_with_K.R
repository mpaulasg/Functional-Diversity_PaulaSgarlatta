################################################################################
##
## Script for computing species positions in a multidimensional space
## according to their trait values
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
##
##
## [PS] Using data with Kmax as a trait. 
################################################################################ If this is going to be used, the other ones will be deleted.

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)

## loading ####

# load traits data ----
fish_traits_k <- read.csv(here::here("data", "raw_data", "fish_traits_k.csv"), header = T)

# load species names from surveys datasets ----
load(here::here("data", "species_allsurveys.RData") )
length(species_allsurveys) # 139 species    # [PS] -> don't understand why this step. Why is species_allsurveys needed here?

# checking same species in trait and occurrences datasets ----
identical ( sort(species_allsurveys) , sort(fish_traits_k$Species ) ) # True

## preparing trait dataset ####

# trait values in a dataframe (species in alphabetical order) ----
sp_tr <- fish_traits_k %>%
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

#Kmax as number

sp_tr$Kmax <- as.numeric(sp_tr$Kmax)
summary(sp_tr$Kmax)

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
# => 3D space has the lowest mAD (0.049)

# species coordinates
sp_3D_coord<-funct_spaces$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord)



# [PS] This is the original one. Now I'll add one separate for each dataset (kelp, no kelp, spatial)

#Loading the data (if this work, this data should be saved in raw_data in CS_00_0a_fish_traits).

sp_tr_kelp_k <- read.csv(here::here("data", "raw_data", "kelp_traits_K.csv"),
                       header = T)
sp_tr_kelp_k <- sp_tr_kelp_k %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()

# from sites that never had kelp
sp_tr_nokelp_k <- read.csv(here::here("data", "raw_data", "nokelp_traits_K.csv"),
                         header = T)

sp_tr_nokelp_k <-sp_tr_nokelp_k %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()

# traits of species from UVC surveys
sp_tr_spatial_k <- read.csv(here::here("data", "raw_data", "spatial_traits_K.csv"), 
                          header = T)

sp_tr_spatial_k <- sp_tr_spatial_k %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()

####KELP 

# recoding variable to match trait type ----

# trait type
tr_cat<-data.frame( trait_name = c("Size", "Agg", "Position", "Diet", "Kmax"),
                    trait_type = c("O","O","O", "N", "Q") )


# size as ordinal
sp_tr_kelp_k$Size <- factor(sp_tr_kelp_k$Size, 
                          levels = c("S2", "S3", "S4", "S5", "S6"),
                          ordered = TRUE)
summary(sp_tr_kelp_k$Size)

# aggregation as ordinal
sp_tr_kelp_k$Agg <- factor(sp_tr_kelp_k$Agg, 
                         levels = c("Solitary", "Pair", "Group"),
                         ordered = TRUE)
summary(sp_tr_kelp_k$Agg)

# Position as ordinal
sp_tr_kelp_k$Position <- factor(sp_tr_kelp_k$Position, 
                              levels = c("Benthic", "BenthoP", "Pelagic"),
                              ordered = TRUE)
summary(sp_tr_kelp_k$Position)

# diet as factor
sp_tr_kelp_k$Diet <- as.factor(sp_tr_kelp_k$Diet)
summary(sp_tr_kelp_k$Diet)

#Kmax as number

sp_tr_kelp_k$Kmax <- as.numeric(sp_tr_kelp_k$Kmax)
summary(sp_tr_kelp_k$Kmax)

# summary of trait data----
summary_traits_kelp_k <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                          sp_tr  = sp_tr_kelp_k)
summary_traits_kelp_k


## Computing Gower distance between species ####
sp_gower_dist_kelp_k <- mFD::funct.dist(sp_tr=sp_tr_kelp_k, tr_cat = tr_cat, 
                                      metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist_kelp_k) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces_kelp_k<- mFD::quality.fspaces(sp_dist = sp_gower_dist_kelp_k, maxdim_pcoa = 12, 
                                         deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces_kelp_k$quality_fspaces
# => 3D space has the lowest mAD (0.047)

# species coordinates
sp_3D_coord_kelp_k<-funct_spaces_kelp_k$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord_kelp_k)

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr_kelp_k, file=here::here("data/", "sp_tr_kelp_k.RData") )
save(summary_traits_kelp_k, file=here::here("outputs/", "summary_traits_kelp_k.RData") )
save(sp_gower_dist_kelp_k, file=here::here("outputs/", "sp_gower_dist_kelp_k.RData") )
save(sp_3D_coord_kelp_k, file=here::here("outputs/", "sp_3D_coord_kelp_k.RData") )


####NO KELP 

# size as ordinal
sp_tr_nokelp_k$Size <- factor(sp_tr_nokelp_k$Size, 
                            levels = c("S2", "S3", "S4", "S5", "S6"),
                            ordered = TRUE)
summary(sp_tr_nokelp_k$Size)

# aggregation as ordinal
sp_tr_nokelp_k$Agg <- factor(sp_tr_nokelp_k$Agg, 
                           levels = c("Solitary", "Pair", "Group"),
                           ordered = TRUE)
summary(sp_tr_nokelp_k$Agg)

# Position as ordinal
sp_tr_nokelp_k$Position <- factor(sp_tr_nokelp_k$Position, 
                                levels = c("Benthic", "BenthoP", "Pelagic"),
                                ordered = TRUE)
summary(sp_tr_nokelp_k$Position)

# diet as factor
sp_tr_nokelp_k$Diet <- as.factor(sp_tr_nokelp_k$Diet)
summary(sp_tr_nokelp_k$Diet)

#Kmax as number

sp_tr_nokelp_k$Kmax <- as.numeric(sp_tr_nokelp_k$Kmax)
summary(sp_tr_nokelp_k$Kmax)


# summary of trait data----
summary_traits_nokelp_k <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                            sp_tr  = sp_tr_nokelp_k)
summary_traits_nokelp_k


## Computing Gower distance between species ####
sp_gower_dist_nokelp_k <- mFD::funct.dist(sp_tr=sp_tr_nokelp_k, tr_cat = tr_cat, 
                                        metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist_nokelp_k) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces_nokelp_k<- mFD::quality.fspaces(sp_dist = sp_gower_dist_nokelp_k, maxdim_pcoa = 12, 
                                           deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces_nokelp_k$quality_fspaces
# => 3D space has the lowest mAD (0.047)

# species coordinates
sp_3D_coord_nokelp_k<-funct_spaces_nokelp_k$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord_nokelp_k)

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr_nokelp_k, file=here::here("data/", "sp_tr_nokelp_k.RData") )
save(summary_traits_nokelp_k, file=here::here("outputs/", "summary_traits_nokelp_k.RData") )
save(sp_gower_dist_nokelp_k, file=here::here("outputs/", "sp_gower_dist_nokelp_k.RData") )
save(sp_3D_coord_nokelp_k, file=here::here("outputs/", "sp_3D_coord_nokelp_k.RData") )


####SPATIAL

# size as ordinal
sp_tr_spatial_k$Size <- factor(sp_tr_spatial_k$Size, 
                             levels = c("S2", "S3", "S4", "S5", "S6"),
                             ordered = TRUE)
summary(sp_tr_spatial_k$Size)

# aggregation as ordinal
sp_tr_spatial_k$Agg <- factor(sp_tr_spatial_k$Agg, 
                            levels = c("Solitary", "Group"),
                            ordered = TRUE)
summary(sp_tr_spatial_k$Agg)

# Position as ordinal
sp_tr_spatial_k$Position <- factor(sp_tr_spatial_k$Position, 
                                 levels = c("Benthic", "BenthoP", "Pelagic"),
                                 ordered = TRUE)
summary(sp_tr_spatial_k$Position)

# diet as factor
sp_tr_spatial_k$Diet <- as.factor(sp_tr_spatial_k$Diet)
summary(sp_tr_spatial_k$Diet)

# summary of trait data----
summary_traits_spatial_k <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                             sp_tr  = sp_tr_spatial_k)
summary_traits_spatial_k


## Computing Gower distance between species ####
sp_gower_dist_spatial_k <- mFD::funct.dist(sp_tr=sp_tr_spatial_k, tr_cat = tr_cat, 
                                         metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist_spatial_k) # from 0 to 1


# computing PCoA-based functional spaces ----
# mean absolute deviation index (as quality metric)
funct_spaces_spatial_k <- mFD::quality.fspaces(sp_dist = sp_gower_dist_spatial_k, maxdim_pcoa = 12, 
                                            deviation_weighting = "absolute", fdist_scaling = FALSE) 
funct_spaces_spatial_k$quality_fspaces
# => 3D space has the lowest mAD (0.041)

# species coordinates
sp_3D_coord_spatial_k<-funct_spaces_spatial_k$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord_spatial_k)

# saving ####

# trait values and trait coding dataframes ----
save(sp_tr_spatial_k, file=here::here("data/", "sp_tr_spatial_k.RData") )
save(summary_traits_spatial_k, file=here::here("outputs/", "summary_traits_spatial_k.RData") )
save(sp_gower_dist_spatial_k, file=here::here("outputs/", "sp_gower_dist_spatial_k.RData") )
save(sp_3D_coord_spatial_k, file=here::here("outputs/", "sp_3D_coord_spatial_k.RData") )

# @@@ ADD code for plots of funct space, correl tr vs axes and save them in outputs



# saving ####

# trait values and trait coding dataframes ----
save(sp_tr, file=here::here("data/", "sp_tr.RData") )
save(tr_cat, file=here::here("data/", "tr_cat.RData") )
save(summary_traits, file=here::here("outputs/", "summary_traits.RData") )
save(sp_gower_dist, file=here::here("outputs/", "sp_gower_dist.RData") )
save(sp_3D_coord, file=here::here("outputs/", "sp_3D_coord.RData") )


####################### end ####


