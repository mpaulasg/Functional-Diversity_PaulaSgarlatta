################################################################################
##
## Script for computing taxonomic and functional diversity and dissimilarity 
## between habitats using Hill numbers
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(betapart)
library(dplyr)

# loading data
load(here::here("data", "using biomass-maxN", "spatial_sp_biom.RData") )
load(here::here("outputs", "sp_3D_coord.RData") )
load(here::here("data", "using biomass-maxN", "kelp_sp_maxN.RData") )

# computing Euclidean distance between species in the 3D space
sp_dist_funct <- dist(sp_3D_coord)

## computing functional diversity based on Hill numbers  - SPATIAL ####

#Taxonomic diversity

spatial_tax_Hill   <- alpha.fd.hill (asb_sp_w = spatial_sp_biom,
                                     sp_dist = sp_dist_funct,
                                     q = 0, # Only 1 because the three values of q give same results (as is species richness)
                                     tau= "min") #This gives you the taxonomic profile

TD_spatial <- spatial_tax_Hill$asb_FD_Hill 

# Functional diversity: number of species, functional richness, dispersion and identity (along 3 axes)

spatial_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = spatial_sp_biom,
                                            sp_dist = sp_dist_funct,
                                            q = c(0, 1, 2), 
                                            tau= "mean") #This gives you functional diversity, as recommended in Chao et al 2019

FD_spatial <- spatial_alpha_FDhill$asb_FD_Hill


## computing tax/functional beta-diversity based on Hill numbers ####

#Taxonomic dissimilarity

spatial_beta_taxhill <-beta.fd.hill (asb_sp_w = spatial_sp_biom,
                                     sp_dist = sp_dist_funct,
                                     q = c(0, 1, 2),
                                     tau= "min",
                                     beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

TD_beta_spatial <- dist.to.df(spatial_beta_taxhill$beta_fd_q) %>% 
  mutate(Habitat1=sub(".*_", "", x1), Habitat2=sub(".*_", "", x2)) %>% 
  select(Habitat1, Habitat2, q0, q1, q2)

# functional dissimilarity
spatial_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = spatial_sp_biom,
                                          sp_dist = sp_dist_funct,
                                          q = c(0, 1, 2),
                                          tau= "mean",
                                          beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

FD_beta_spatial <- dist.to.df(spatial_beta_FDhill$beta_fd_q)%>% 
  mutate(Habitat1=sub(".*_", "", x1), Habitat2=sub(".*_", "", x2)) %>% 
  select(Habitat1, Habitat2, q0, q1, q2)


#################################   KELP  #########################################

#Taxonomic diversity

kelp_tax_Hill   <- alpha.fd.hill (asb_sp_w = kelp_sp_maxN,
                                  sp_dist = sp_dist_funct,
                                  q = 0, 
                                  tau= "min") #This gives you the taxonomic profile

TD_kelp_Hill <- kelp_tax_Hill$asb_FD_Hill 

# Functional diversity: number of species, functional richness, dispersion and identity (along 3 axes)

kelp_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = kelp_sp_maxN,
                                         sp_dist = sp_dist_funct,
                                         q = c(0, 1, 2), 
                                         tau= "mean") #This gives you functional diversity, as recommended in Chao et al 2019

FD_kelp_Hill <- kelp_alpha_FDhill$asb_FD_Hill


## computing tax/functional beta-diversity based on Hill numbers ####

#Taxonomic dissimilarity

kelp_beta_taxhill <-beta.fd.hill (asb_sp_w = kelp_sp_maxN,
                                  sp_dist = sp_dist_funct,
                                  q = c(0, 1, 2),
                                  tau= "min",
                                  beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

TD_beta_kelp <- dist.to.df(kelp_beta_taxhill$beta_fd_q) %>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  select(Year1, Year2, q0, q1, q2)

# functional dissimilarity
kelp_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = kelp_sp_maxN,
                                       sp_dist = sp_dist_funct,
                                       q = c(0, 1, 2),
                                       tau= "mean",
                                       beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

FD_beta_kelp <- dist.to.df(kelp_beta_FDhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  select(Year1, Year2, q0, q1, q2, x1, x2)

#For stats

FD_beta_kelp_Hill_sites <- dist.to.df(kelp_beta_FDhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  mutate(Site1=sub("_.*", "", x1),Site2=sub("_.*", "", x2) ) %>% 
  select(Year1, Year2, Site1, Site2, q0, q1, q2)

TD_beta_kelp_Hill_sites <- dist.to.df(kelp_beta_taxhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  mutate(Site1=sub("_.*", "", x1), Site2=sub("_.*", "", x2)) %>% 
  select(Year1, Year2, Site1, Site2, q0, q1, q2)



# saving ####
save(sp_dist_funct, file=here::here("outputs/", "using biomass-maxN", "sp_dist_funct_spatial_Hill.RData") )
save(TD_spatial, file=here::here("outputs/", "using biomass-maxN", "TD_spatial_Hill.RData") )
save(FD_spatial, file=here::here("outputs/", "using biomass-maxN", "FD_spatial_Hill.RData") )
save(TD_beta_spatial, file=here::here("outputs/","using biomass-maxN",  "TD_beta_spatial_Hill.RData") )
save(FD_beta_spatial, file=here::here("outputs/", "using biomass-maxN", "FD_beta_spatial_Hill.RData") )

save(TD_kelp_Hill, file=here::here("outputs/","using biomass-maxN", "TD_kelp_Hill.RData") )
save(FD_kelp_Hill, file=here::here("outputs/","using biomass-maxN", "FD_kelp_Hill.RData") )
save(TD_beta_kelp, file=here::here("outputs/", "using biomass-maxN","TD_beta_kelp_Hill.RData") )
save(FD_beta_kelp, file=here::here("outputs/", "using biomass-maxN","FD_beta_kelp_Hill.RData") )


save(FD_beta_kelp_Hill_sites, file=here::here("outputs/","using biomass-maxN", "FD_beta_kelp_Hill_sites.RData") )
save(TD_beta_kelp_Hill_sites, file=here::here("outputs/","using biomass-maxN", "TD_beta_kelp_Hill_sites.RData") )







