################################################################################
##
## Script for computing taxonomic and functional diversity and dissimilarity 
## between years using Hill numbers
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
load(here::here("data", "kelp_sp_occ.RData") )
load(here::here("data", "nokelp_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") ) #For some reason this one doesn't work
#load(here::here("outputs", "sp_3D_coord_kelp.RData") )
#load(here::here("outputs", "sp_3D_coord_nokelp.RData") )

# computing Euclidean distance between species in the 3D space

sp_dist_funct <- dist(sp_3D_coord)

#sp_dist_funct_kelp <- dist(sp_3D_coord_kelp)
#sp_dist_funct_nokelp <- dist(sp_3D_coord_nokelp)

## computing functional diversity based on Hill numbers ####

####### [PS] Same here, don't  know why it's not working for me. I need to use the specific one (in this case, sp_3D_coord_kelp/nokelp)

############################################################   KELP  #########################################################################

#Taxonomic diversity

kelp_tax_Hill   <- alpha.fd.hill (asb_sp_w = kelp_sp_occ,
                                     sp_dist = sp_dist_funct,
                                     q = 0, 
                                     tau= "min") #This gives you the taxonomic profile

TD_kelp_Hill <- kelp_tax_Hill$asb_FD_Hill 

# Functional diversity: number of species, functional richness, dispersion and identity (along 3 axes)

kelp_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = kelp_sp_occ,
                                            sp_dist = sp_dist_funct_kelp,
                                            q = c(0, 1, 2) , 
                                            tau= "mean") #This gives you functional diversity, as recommended in Chao et al 2019

FD_kelp_Hill <- kelp_alpha_FDhill$asb_FD_Hill


## computing tax/functional beta-diversity based on Hill numbers ####

#Taxonomic dissimilarity

kelp_beta_taxhill <-beta.fd.hill (asb_sp_w = kelp_sp_occ,
                                     sp_dist = sp_dist_funct_kelp,
                                     q = c(0,1, 2),
                                     tau= "min",
                                     beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

TD_beta_kelp <- dist.to.df(kelp_beta_taxhill$beta_fd_q) %>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  select(Year1, Year2, q0, q1, q2)

# functional dissimilarity
kelp_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = kelp_sp_occ,
                                          sp_dist = sp_dist_funct_kelp,
                                          q = c(0,1, 2),
                                          tau= "mean",
                                          beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

#For stats

FD_beta_kelp_Hill_sites <- dist.to.df(kelp_beta_FDhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  mutate(Site=sub("_.*", "", x1)) %>% 
  select(Year1, Year2, Site, q0, q1, q2)

FD_beta_kelp <- dist.to.df(kelp_beta_FDhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  select(Year1, Year2, q0, q1, q2)

#Saving

save(sp_dist_funct_kelp, file=here::here("outputs/", "sp_dist_funct_kelp.RData") )
save(TD_kelp_Hill, file=here::here("outputs/", "TD_kelp_Hill.RData") )
save(FD_kelp_Hill, file=here::here("outputs/", "FD_kelp_Hill.RData") )
save(TD_beta_kelp, file=here::here("outputs/", "TD_beta_kelp_Hill.RData") )
save(FD_beta_kelp, file=here::here("outputs/", "FD_beta_kelp_Hill.RData") )
save(FD_beta_kelp_Hill_sites, file=here::here("outputs/", "FD_beta_kelp_Hill_sites.RData") )


############################################################   NO KELP  #########################################################################

#Taxonomic diversity

nokelp_tax_Hill   <- alpha.fd.hill (asb_sp_w = nokelp_sp_occ,
                                  sp_dist = sp_dist_funct_nokelp,
                                  q = 2, #Recommended in Chao et al 2019
                                  tau= "min") #This gives you the taxonomic profile

TD_nokelp_Hill <- nokelp_tax_Hill$asb_FD_Hill 

# Functional diversity: number of species, functional richness, dispersion and identity (along 3 axes)

nokelp_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = nokelp_sp_occ,
                                         sp_dist = sp_dist_funct_nokelp,
                                         q = c(0, 1, 2) , 
                                         tau= "mean") #This gives you functional diversity, as recommended in Chao et al 2019

FD_nokelp_Hill <- nokelp_alpha_FDhill$asb_FD_Hill


## computing tax/functional beta-diversity based on Hill numbers ####

#Taxonomic dissimilarity

nokelp_beta_taxhill <-beta.fd.hill (asb_sp_w = nokelp_sp_occ,
                                  sp_dist = sp_dist_funct_nokelp,
                                  q = c(0,1, 2),
                                  tau= "min",
                                  beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

TD_beta_nokelp <- dist.to.df(nokelp_beta_taxhill$beta_fd_q) %>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  select(Year1, Year2, q0, q1, q2)

# functional dissimilarity
nokelp_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = nokelp_sp_occ,
                                       sp_dist = sp_dist_funct_nokelp,
                                       q = c(0,1, 2),
                                       tau= "mean",
                                       beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

#With sites - for statystical analysis

FD_beta_nokelp_sites <- dist.to.df(nokelp_beta_FDhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  mutate(Site=sub("_.*", "", x1)) %>% 
  select(Year1, Year2, Site, q0, q1, q2)


FD_beta_nokelp <- dist.to.df(nokelp_beta_FDhill$beta_fd_q)%>% 
  mutate(Year1=sub(".*_", "", x1), Year2=sub(".*_", "", x2)) %>% 
  select(Year1, Year2, q0, q1, q2)

#From Seb

tau_mean<-mean(sp_dist_funct)
tau_mean # 0.39
test_occ <- spatial_sp_occ[1:2,]
sp_test<-names(which(apply(test_occ,2,max)>0))
dist_test<-as.matrix(sp_dist_funct)[sp_test,sp_test]
test_occ[,sp_test]
unique_test<-names(which(apply(test_occ,2,sum)==1))
dist_unique_test<-dist_test[unique_test,unique_test]
summary(as.dist(dist_unique_test))   # all distances between species unique to different assemblages < mean(dist) 


# summary
lapply(spatial_beta_FDhill, summary)
# => beta q0 null because all most distinct pairs of species are in all assemblages
# => beta q=1 still very low, similar composition 

# saving ####

save(sp_dist_funct_nokelp, file=here::here("outputs/", "sp_dist_funct_nokelp.RData") )
save(TD_nokelp_Hill, file=here::here("outputs/", "TD_nokelp_Hill.RData") )
save(FD_nokelp_Hill, file=here::here("outputs/", "FD_nokelp_Hill.RData") )
save(TD_beta_nokelp, file=here::here("outputs/", "TD_beta_nokelp_Hill.RData") )
save(FD_beta_nokelp, file=here::here("outputs/", "FD_beta_nokelp_Hill.RData") )
save(FD_beta_nokelp_sites, file=here::here("outputs/", "FD_beta_nokelp_sites.RData") )

