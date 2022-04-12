
#######Not to be in final GH

rm(list=ls()) # cleaning memory


# libraries
library(tidyverse)
library(here)

spatial_thermal <- read.csv(here::here("from_paula", "SpatialUVC_species_traits_only_thermal.csv") )
kelp_thermal <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_only_thermal_kelp.csv") )


thermal_all <- bind_rows(spatial_thermal, kelp_thermal)%>%
  distinct(Species, .keep_all = TRUE ) %>%
  as.data.frame()

#Calculate thermal midpoint

thermal_all <- thermal_all %>% 
  mutate(thermal=rowMeans(thermal_all[ , c(2,3)], na.rm=TRUE)) %>% 
  select(Species, thermal)

# saving as csv file

write.csv(thermal_all, file=here::here("data", "raw_data", "thermal_all.csv"), 
          row.names = FALSE )

