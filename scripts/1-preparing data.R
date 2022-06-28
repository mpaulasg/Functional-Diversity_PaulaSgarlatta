################################################################################
##
## Script for preparing fish occurrence datasets and for computing species positions 
## in a multidimensional space according to their trait values
## 
## Code by Camille Magneville, Sebastien villeger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(dplyr)


######### Preparing data #############

## temporal survey data ####

# loading raw data from csv----


# metadata and data from sites that use to have kelp and lost it:
kelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_metadata.csv") )
unique(kelp_metadata$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()

# summary of surveys and occurrences data ----
kelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN)

# retrieve occurrence matrix:
kelp_sp_occ <- kelp_summary$asb_sp_occ

#To check how species MaxN change ovr the years:

temporal_maxn_sp <-  kelp_sp_maxN %>% 
  as.data.frame() %>% 
  rownames_to_column("Site") %>% 
  pivot_longer(!Site, names_to = "species", values_to = "Biomass")%>%  
  mutate(Year=sub(".*_", "", Site)) %>% 
  mutate(genus=sub("_.*", "", species), sp=sub(".*_", "", species)) %>% 
  mutate(Species=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  dplyr::select(-Site, -genus, -sp, -species) %>% 
  pivot_wider(names_from = Species, values_from = Biomass, values_fn= sum)


# dimensions
dim(kelp_sp_occ) # 69 assemblages * 101 species

# names of species
kelp_sp <- colnames(kelp_sp_occ) 
length(kelp_sp) # 101 sp

## This is to make a graph later

kelp_sp_data <- kelp_sp %>% 
  as.data.frame() %>% 
  rename(Species = ".") %>% 
  mutate(data_1 = "temporal")


#=> temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata <- read.csv(here::here("data", "raw_data",  "SpatialUVC_metadata_site.csv"))
head(spatial_metadata)

spatial_sp_biom <- read.csv(here::here("data", "raw_data", "SpatialUVC_species_biomass_site_average.csv")) %>%
  column_to_rownames("Code") %>% 
  as.matrix()

dim(spatial_sp_biom) # 9 assemblages * 53 species

# summary of surveys and occurrences data ----
spatial_summary <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom)

# retrieve occurrence matrix:
spatial_sp_occ <- as.data.frame(spatial_summary$asb_sp_occ)

#To check how species biomass change on each habitat:

spatial_biomass_sp <-  spatial_sp_biom %>% 
  as.data.frame() %>% 
  rownames_to_column("Site") %>% 
  pivot_longer(!Site, names_to = "species", values_to = "Biomass")%>%  
  mutate(Habitat=sub(".*_", "", Site)) %>% 
  mutate(genus=sub("_.*", "", species), sp=sub(".*_", "", species)) %>% 
  mutate(Species=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  dplyr::select(-Site, -genus, -sp, -species) %>% 
  pivot_wider(names_from = Species, values_from = Biomass, values_fn= sum)
  

# retrieve total biomass per assemblage:

spatial_biomass <- as.data.frame(spatial_summary$asb_tot_w) %>% 
  rownames_to_column("Site") %>% 
rename(Biomass = "spatial_summary$asb_tot_w")

# species names
spatial_sp <- colnames(spatial_sp_occ)
length(spatial_sp) # 53

## This is to make a graph later

spatial_sp_data <- spatial_sp %>% 
  as.data.frame() %>% 
  rename(Species = ".") %>% 
  mutate(data_2 = "spatial")


## names of species present in at least one dataset ####
sum(spatial_sp %in% kelp_sp)# 35 species shared 

species_allsurveys <- unique( c(kelp_sp ,  spatial_sp) ) 
length(species_allsurveys) # 119 species


## This is to make a graph later

species_both <- kelp_sp_data %>% 
  full_join(spatial_sp_data, by ="Species") %>% 
  replace (is.na(.), "no")


## saving dataframes #####
save(kelp_metadata, file=here::here("data", "kelp_metadata.RData") )
save(kelp_sp_occ, file=here::here("data", "kelp_sp_occ.RData") )
save(kelp_sp_maxN, file=here::here("data", "kelp_sp_maxN.RData") )


save(spatial_metadata, file=here::here("data", "spatial_metadata.RData") )
save(spatial_sp_occ, file=here::here("data","spatial_sp_occ.RData") )
save(spatial_sp_biom, file=here::here("data","spatial_sp_biom.RData") )

save(spatial_summary, file=here::here("data","spatial_summary.RData") )

save(species_allsurveys, file=here::here("data", "species_allsurveys.RData") )

write.csv(species_both, file=here::here("data", "species_both.csv"), 
           row.names = FALSE )

write.csv(spatial_biomass, file=here::here("data", "spatial_biomass.csv"), 
          row.names = FALSE )

write.csv(spatial_biomass_sp, file=here::here("data", "spatial_biomass_sp.csv"), 
          row.names = FALSE )

write.csv(temporal_maxn_sp, file=here::here("data", "temporal_maxn_sp.csv"), 
          row.names = FALSE )


######### Species position in multidimensional space #############

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

sp_faxes_coord <- funct_spaces$details_fspaces$sp_pc_coord

# species coordinates
sp_3D_coord<-funct_spaces$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord)

# saving ####

# trait values and trait coding dataframes ----

save(sp_tr, file=here::here("data/", "sp_tr.RData") )
save(tr_cat, file=here::here("data/", "tr_cat.RData") )
save(summary_traits, file=here::here("data/", "summary_traits.RData") )
save(sp_gower_dist, file=here::here("data/", "sp_gower_dist.RData") )
save(sp_3D_coord, file=here::here("data/", "sp_3D_coord.RData") )
save(funct_spaces, file=here::here("data/", "funct_spaces.RData") )
save(sp_faxes_coord, file=here::here("data/", "sp_faxes_coord.RData") )
