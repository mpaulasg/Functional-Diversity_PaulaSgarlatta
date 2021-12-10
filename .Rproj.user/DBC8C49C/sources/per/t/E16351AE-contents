################################################################################
##
## Script for merging all traits values into a single csv file
## 
##  Sébastien Villéger
##
                                                                                ## NB: NOT TO BE ON FINAL GITHUB
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(dplyr)
library(arsenal)

# load traits data ----

# traits of species from sites surveyed through years using BRUVs 
# from sites with kelp
kelp_traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_kelp.csv"),
                        header = T)
head(kelp_traits)
names(kelp_traits) <- c("Species","Size", "Aggr", "Posi", "Diet")
nrow(kelp_traits) # 101 sp

#[PS] Adding a new datset with thermal affinity

kelp_thermal <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_only_thermal_kelp.csv"),
                         header = T)

#[PS] Check if both df have the same species

summary(comparedf(kelp_traits, kelp_thermal, by = "row.names")) #Great!

kelp_traits_with_thermal <- cbind(kelp_traits, kelp_thermal[c("Thermal_affinity_min", "Thermal_affinity_max")])

# from sites that never had kelp
nokelp_traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_no_kelp.csv"),
                          header = T)
head(nokelp_traits)
names(nokelp_traits) <- c("Species","Size", "Aggr", "Posi", "Diet")
nrow(nokelp_traits) # 106 sp

#[PS] Adding a new datset with thermal affinity

no_kelp_thermal <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_only_thermal_no_kelp.csv"),
                         header = T)

#[PS] Check if both df have the same species

summary(comparedf(nokelp_traits, no_kelp_thermal, by = "row.names")) #Great!

nokelp_traits_with_thermal <- cbind(nokelp_traits, no_kelp_thermal[c("Thermal_affinity_min", "Thermal_affinity_max")])

# traits of species from UVC surveys
spatial_traits <- read.csv(here::here("from_paula", "SpatialUVC_species_traits.csv"), 
                           header = T)
names(spatial_traits) <- c("Species", "Size", "Aggr", "Posi", "Diet")
head(spatial_traits)
nrow(spatial_traits) # 51 sp

#[PS] Adding a new datset with thermal affinity

spatial_thermal <- read.csv(here::here("from_paula", "SpatialUVC_species_traits_only_thermal.csv"),
                            header = T)

#[PS] Check if both df have the same species

summary(comparedf(spatial_traits, spatial_thermal, by = "row.names")) #Great!

spatial_traits_with_thermal <- cbind(spatial_traits, spatial_thermal[c("Thermal_affinity_min", "Thermal_affinity_max")])


# merging in a single data.frame and keeping only species present in surveys ----
fish_traits <- bind_rows( kelp_traits, nokelp_traits, spatial_traits) %>%
  distinct(Species, .keep_all = TRUE ) %>%
  as.data.frame()
head(fish_traits)

#[PS] same with thermal data

fish_traits_thermal <- bind_rows(kelp_traits_with_thermal, nokelp_traits_with_thermal, spatial_traits_with_thermal) %>%
  distinct(Species, .keep_all = TRUE ) %>%
  as.data.frame()
head(fish_traits_thermal)


nrow(fish_traits) # 142 species
nrow(fish_traits_thermal) # 142 species


sort(fish_traits$Species) #
#                                                                               

# saving as csv file
write.csv(fish_traits, file=here::here("data", "raw_data", "fish_traits.csv"), 
          row.names = FALSE )


write.csv(fish_traits_thermal, file=here::here("data", "raw_data", "fish_traits_thermal.csv"), 
          row.names = FALSE )
