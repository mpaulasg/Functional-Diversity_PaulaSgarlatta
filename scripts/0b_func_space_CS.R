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
length(species_allsurveys) # 119 species    

# checking same species in trait and occurrences datasets ----
identical ( sort(species_allsurveys) , sort(fish_traits$Species ) ) # True

## preparing trait dataset ####

# trait values in a dataframe (species in alphabetical order) ----
sp_tr <- fish_traits %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()
head(sp_tr)

nrow(sp_tr) # 119 species

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

# species coordinates
sp_3D_coord<-funct_spaces$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord)

#see variance explained by the PCoA axes

gower<-funct_spaces$details_deviation ## CHECK
  
# saving ####

# trait values and trait coding dataframes ----

save(sp_tr, file=here::here("data/", "sp_tr.RData") )
save(tr_cat, file=here::here("data/", "tr_cat.RData") )
save(summary_traits, file=here::here("outputs/", "summary_traits.RData") )
save(sp_gower_dist, file=here::here("outputs/", "sp_gower_dist.RData") )
save(sp_3D_coord, file=here::here("outputs/", "sp_3D_coord.RData") )
save(funct_spaces, file=here::here("outputs/", "funct_spaces.RData") )

##################################  end of code ######################################################################