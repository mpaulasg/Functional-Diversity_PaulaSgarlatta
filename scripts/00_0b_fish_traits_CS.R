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
names(kelp_traits) <- c("Species","Size", "Agg", "Position" , "Diet")
nrow(kelp_traits) # 101 sp

# from sites that never had kelp
nokelp_traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_no_kelp.csv"),
                          header = T)
head(nokelp_traits)
names(nokelp_traits) <- c("Species","Size", "Agg", "Position" , "Diet")
nrow(nokelp_traits) # 106 sp

# traits of species from UVC surveys
spatial_traits <- read.csv(here::here("from_paula", "SpatialUVC_species_traits.csv"), 
                           header = T)
head(spatial_traits)
names(spatial_traits) <- c("Species","Size", "Agg", "Position" , "Diet")

nrow(spatial_traits) # 51 sp

#All sp with K values

k_values <- read.csv(here::here("from_paula", "fish_K_values.csv"),
                     header = T)
names(k_values) <- c("Species", "sstmean", "MaxSizeTL", "Diet", "Position", "Method", "Kmax", "Kmax_lowq", "Kmax_uppq")


#Add K values to each dataframe

kelp_traits <- inner_join(kelp_traits, k_values[ , c("Species", "Kmax")], by = c("Species"), all.x=TRUE)
  
nokelp_traits <- inner_join(nokelp_traits, k_values[ , c("Species", "Kmax")], by = c("Species"), all.x=TRUE)   

spatial_traits <- inner_join(spatial_traits, k_values[ , c("Species", "Kmax")], by = c("Species"), all.x=TRUE) 
  
# merging in a single data.frame and keeping only species present in surveys ----
fish_traits <- bind_rows( kelp_traits, nokelp_traits, spatial_traits) %>%
  distinct(Species, .keep_all = TRUE ) %>%
  as.data.frame()
head(fish_traits)

nrow(fish_traits) # 139 species

fish_traits <- as.data.frame(fish_traits[order(fish_traits$Species),])

# saving as csv file

write.csv(fish_traits, file=here::here("data", "raw_data", "fish_traits.csv"), 
          row.names = FALSE )


########################################## end of code ###############################################################
