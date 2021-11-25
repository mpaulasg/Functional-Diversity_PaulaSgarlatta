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


# load traits data ----

# traits of species from sites surveyed through years using BRUVs 
# from sites with kelp
kelp_traits <- read.csv(here::here("data", "TemporalBRUV_species_traits_kelp.csv"),
                        header = T)
head(kelp_traits)
names(kelp_traits) <- c("Species","Size", "Aggr", "Posi", "Diet")
nrow(kelp_traits) # 101 sp

# from sites that never had kelp
nokelp_traits <- read.csv(here::here("data", "TemporalBRUV_species_traits_no_kelp.csv"),
                          header = T)
head(nokelp_traits)
names(nokelp_traits) <- c("Species","Size", "Aggr", "Posi", "Diet")
nrow(nokelp_traits) # 106 sp

# traits of species from UVC surveys
spatial_traits <- read.csv(here::here("data", "SpatialUVC_species_traits.csv"), 
                           header = T)
names(spatial_traits) <- c("Species", "Size", "Aggr", "Posi", "Diet")
head(spatial_traits)
nrow(spatial_traits) # 51 sp

# merging in a single data.frame and keeping only species present in surveys ----
fish_traits <- bind_rows( kelp_traits, nokelp_traits, spatial_traits) %>%
  distinct(Species, .keep_all = TRUE ) %>%
  as.data.frame()
head(fish_traits)

nrow(fish_traits) # 143 species

sort(fish_traits$Species) #
#                                                                               => @@ check whether typo in "A_dussumeri"         "A_dussumieri"   

# saving as csv file
write.csv(fish_traits, file=here::here("data", "fish_traits.csv"), 
          row.names = FALSE )



