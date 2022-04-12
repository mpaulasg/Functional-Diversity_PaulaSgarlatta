################################################################################
##
## Script for kelp spatial FD analysis (Paula Sgarlatta)
##
## Start: 25 Feb 2021 - End: 07 Ap 2021
## Code by Camille Magneville
##
##1st part: Indices per habitat (for graphs)
##
## 1/ First load and transform dataset
## 2/ Then, basic FD analysis
## 3/ Then, compute alpha indices (FDiv, FRic, FEve, FIde)##
##
## 2nd part: Taxonomic analysis (richness and taxonomic
## beta diversity)
##
################################################################################

rm(list=ls()) # cleaning memory


#### 1 - Load (and transform) datasets ####
####  Data per site ####


# load sp*tr data:
sp_tr <- read.csv("data/SpatialUVC_species_traits.csv")

# transform columns type:
sp_tr$Size <- as.ordered(sp_tr$Size)
ordered(sp_tr$Size)
sp_tr$Aggregation <- as.ordered(sp_tr$Aggregation)
ordered(sp_tr$Aggregation)
# I have changed the order of the groups with Solitary < Pair < Group ...
# ... instead of Group < Pair < Solitary:
sp_tr[, "Aggregation"] <- factor(sp_tr[, "Aggregation"], levels = c("Solitary", "Pair", "Group"),
                         ordered = TRUE)
ordered(sp_tr$Aggregation)
sp_tr$Position <- as.ordered(sp_tr$Position)
ordered(sp_tr$Position)
sp_tr$Diet <- as.factor(sp_tr$Diet)

# add species as row names:
rownames(sp_tr) <- sp_tr$species_id
sp_tr <- tibble::column_to_rownames(sp_tr, var = "Species")

# load asb data:
asb_sp_w <- read.csv("data/SpatialUVC_species_biomass_site_average.csv")

# group between habitats:
asb_sp_w$hab <- c("Inshore", "Inshore", "Inshore", "Midshelf", "Midshelf", 
                  "Midshelf", "Offshore", "Offshore", "Offshore")

# create a new asb dataframe for habitat (and not sites):
asb_sp_w_hab <- as.data.frame(matrix(ncol = ncol(asb_sp_w), nrow = 3))
colnames(asb_sp_w_hab) <- colnames(asb_sp_w)
asb_sp_w_hab$site_id <- c("Inshore", "Midshelf", "Offshore")

rownames(asb_sp_w_hab) <- asb_sp_w_hab$site_id


# fill this new df:
for (c in colnames(asb_sp_w[, - c(1, ncol(asb_sp_w))])) {
  hab_I <- mean(asb_sp_w[which(asb_sp_w$hab == "Inshore"), c])
  asb_sp_w_hab["Inshore", c] <- hab_I
  hab_M <- mean(asb_sp_w[which(asb_sp_w$hab == "Midshelf"), c])
  asb_sp_w_hab["Midshelf", c] <- hab_M
  hab_O <- mean(asb_sp_w[which(asb_sp_w$hab == "Offshore"), c])
  asb_sp_w_hab["Offshore", c] <- hab_O
}
# remove Site and site_id :
asb_sp_w_hab <- asb_sp_w_hab[, - c(1, ncol(asb_sp_w_hab))]

# remove hab column:
asb_sp_w_hab <- asb_sp_w_hab[, - c(ncol(asb_sp_w_hab))]

# remove the period column (added because contained in the colnames of asb_sp_w cf l54):

asb_sp_w_hab<- asb_sp_w_hab[, -ncol(asb_sp_w_hab)]

# load trait types data:
sp_tr_cat <- read.csv("data/SpatialUVC_traits_type.csv")



################################################################################



#### 2 - Basic FD analysis ####


### Summary of the data: 

library(remotes)
#remotes::install_github("CmlMagneville/mFD")
library(mFD)

## summary for sp*tr data:

sp_tr_summ <- mFD::sp.tr.summary(tr_cat = sp_tr_cat, 
                                 sp_tr  = sp_tr)
# repartition of traits:
sp_tr_summ$tr_summary_list


## summary for asb data:

asb_sp_w_hab <- as.matrix(asb_sp_w_hab)

asb_sp_summ <- mFD::asb.sp.summary(asb_sp_w = asb_sp_w_hab) 

# retrieve occurrence matrix:
asb_sp_occ <- asb_sp_summ$asb_sp_occ

# total biomass of each species:
asb_sp_summ$sp_tot_w

# total biomass for each asb (ie habitat):
asb_sp_summ$asb_tot_w

# species richness per habitat:
asb_sp_summ$asb_sp_richn

# list of species present in each assemblage (ie habitat):
asb_sp_summ$asb_sp_nm


### Compute functional distance:

# compute functional distance based on traits:
sp_dist <- mFD::funct.dist(
  sp_tr  = sp_tr,
  tr_cat = sp_tr_cat,
  metric = "gower",
  scale_euclid  = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA  = TRUE)

# is there a lot of distances = 0? (group into FEs or not?)
sp_dist2 <- as.matrix(sp_dist)
inds <- which(sp_dist2 == min(sp_dist2), arr.ind=TRUE)
inds
# it is ok: only between P_spilurus and D_punctulatus...
# ... We don't need to gather into FEs
sp_dist

### Compute functional spaces and their quality:

# compute the quality of functional spaces:
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist, 
  maxdim_pcoa         = 10, 
  deviation_weighting = "absolute", 
  fdist_scaling       = FALSE, 
  fdendro             = "average")
# ... same warning message saying that some distance equal 0 but we have ...
# ... checked and it is ok!

# retrieve table of quality metric:
fspaces_quality$quality_fspaces
# We can see that the functional space with the lowest quality metric is ...
# ... the one with four dimensions (pcoa_4d), we will thus continue with ...
# ... this functional space!

# illustrate quality of functional space (not useful but can be used as ...
# ... Supplementary Figure for instance):
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")



### Test correlation between traits and functional axes:

# retrieve coordinates of species:
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

# test correlation between traits and axes:
cor_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sp_tr, 
  sp_faxes_coord = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# get the table of correlation:
cor_tr_faxes$tr_faxes_stat

# get the plot:
cor_tr_faxes$tr_faxes_plot

# we can see for instance that...
# ... PC1 is driven by Aggregation with Solitary species on the right and ...
# ... Grouped species on the left, by Diet, thermal affin min and max.



### Plot functional space:

big_plot <- mFD::funct.space.plot(sp_faxes_coord  = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
                                  faxes = NULL, name_file = NULL,
                                  faxes_nm = NULL, range_faxes = c(NA, NA),
                                  color_bg = "grey95",
                                  color_pool = "darkgreen", fill_pool = "white",
                                  shape_pool = 21, size_pool = 1,
                                  plot_ch = TRUE, color_ch = "black",
                                  fill_ch = "white", alpha_ch = 0.3,
                                  plot_vertices = TRUE, 
                                  color_vert = "blueviolet",
                                  fill_vert = "blueviolet", shape_vert = 23,
                                  size_vert = 1,
                                  plot_sp_nm = NULL, nm_size = 3, 
                                  nm_color = "black",
                                  nm_fontface = "plain",
                                  check_input = TRUE)
# Plot the graph with all pairs of axes:
big_plot$patchwork



################################################################################

#### 3 - Compute functional alpha diversity indices ####

### Compute indices:
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = asb_sp_w_hab,
  ind_vect         = c("fdis", "fspe", "fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)


fd_ind_values <- alpha_fd_indices$functional_diversity_indices

#write.csv(fd_ind_values, file = "data/fd_ind_habitat.csv")


################################################################################



#### 4 - Compute beta diversity indices ####


### Compute beta diversity indices:

beta_fd_indices <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = asb_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

beta_fd_indices$pairasb_fbd_indices

# Jordan Gacutan code ammend (18/08/21)

indices = c("Total dissimilarity", "Turnover", "Nestedness")
base = c("Inshore-Midshelf", "Inshore-Offshore", "Midshelf-Offshore")

library (purrr)
library (dplyr)
library (reshape2)
  
temp = map2_dfc(beta_fd_indices$pairasb_fbd_indices, indices, function(x,y){
  mat.x = as.matrix(x)
  diag(mat.x) = NA 
    
  diss.df = melt(mat.x) %>% as_tibble() %>% tidyr::drop_na()
  
  diss.df_fin = distinct(diss.df, value, .keep_all = T)
  names(diss.df_fin) = c("X1", "X2", paste0(y))
  
  output = diss.df_fin %>% select(y)
  
  return(output) 
    
})

fin_table = bind_cols(Pairs = base, temp)

#write.csv(fin_table, file = "data/Fbeta_diversity_spatial_site.csv")


#########################################################################
## 
##
##Functional diversity indices per site/transect (not between habitats)
##
##
#########################################################################

#### Load (and transform) datasets ####


# Use sp*tr data from before:

sp_tr 

# load asb data:
asb_sp_w_s <- read.csv("data/SpatialUVC_species_biomass_site_average.csv")

sp_transect <- read.csv("data/SpatialUVC_species_biomass_transect.csv")

# create a new asb dataframe for transect/site:
asb_sp_w_hab_s <- as.data.frame(asb_sp_w_s)

sp_transect_habitat <- as.data.frame(sp_transect)



# add transects as row names:
rownames(asb_sp_w_hab_s) <- asb_sp_w_hab_s$species_id
asb_sp_w_hab_s <- tibble::column_to_rownames(asb_sp_w_hab_s, var = "Sites")

rownames(sp_transect_habitat) <- sp_transect$species_id
sp_transect_habitat <- tibble::column_to_rownames(sp_transect_habitat, var = "Sites")


# Use trait types data from before:
sp_tr_cat


################################################################################



#### 2 - Basic FD analysis ####


### Summary of the data: 

library(remotes)
#remotes::install_github("CmlMagneville/mFD")
library(mFD)

## summary for sp*tr data:

sp_tr_summ_s <- mFD::sp.tr.summary(tr_cat = sp_tr_cat, 
                                 sp_tr  = sp_tr)

sp_tr_summ_transect <- mFD::sp.tr.summary(tr_cat = sp_tr_cat, 
                                   sp_tr  = sp_tr)

# repartition of traits:
sp_tr_summ_s$tr_summary_list

sp_tr_summ_transect$tr_summary_list


## summary for asb data:

asb_sp_w_hab_s <- as.matrix(asb_sp_w_hab_s)

asb_sp_summ_s <- mFD::asb.sp.summary(asb_sp_w = asb_sp_w_hab_s)

sp_transect_habitat <- as.matrix(sp_transect_habitat)

sp_summ_transect <- mFD::asb.sp.summary(asb_sp_w = sp_transect_habitat)

# retrieve occurrence matrix:
asb_sp_occ_s <- asb_sp_summ_s$asb_sp_occ

sp_occ_transect <- sp_summ_transect$asb_sp_occ

# total biomass of each species:
asb_sp_summ_s$sp_tot_w

sp_summ_transect$sp_tot_w

# total biomass for each asb (ie Site):
asb_sp_summ_s$asb_tot_w

sp_summ_transect$asb_tot_w

#write.csv(asb_sp_summ_s$asb_tot_w, "data/biomass_site.csv")

# species richness per site:
asb_sp_summ_s$asb_sp_richn

sp_summ_transect$asb_sp_richn

# list of species present in each assemblage (ie site):
asb_sp_summ_s$asb_sp_nm

sp_summ_transect$asb_sp_nm


### Compute functional distance:

 #compute functional distance based on traits:
sp_dist <- mFD::funct.dist(
   sp_tr  = sp_tr,
   tr_cat = sp_tr_cat,
   metric = "gower",
   scale_euclid  = "scale_center",
   ordinal_var = "classic",
   weight_type = "equal",
   stop_if_NA  = TRUE)


# is there a lot of distances = 0? (group into FEs or not?)
#No need to check, same than last time

### Compute functional spaces and their quality:

 #compute the quality of functional spaces:
 fspaces_quality <- mFD::quality.fspaces(
   sp_dist             = sp_dist, 
   maxdim_pcoa         = 10, 
   deviation_weighting = "absolute", 
   fdist_scaling       = FALSE, 
   fdendro             = "average")
# # ... same warning message saying that some distance equal 0 but we have ...
# # ... checked and it is ok!
 
# # retrieve table of quality metric:
 fspaces_quality$quality_fspaces

# ### Test correlation between traits and functional axes:
 
# # retrieve coordinates of species:
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"
 
# # test correlation between traits and axes:
 cor_tr_faxes <- mFD::traits.faxes.cor(
   sp_tr          = sp_tr, 
   sp_faxes_coord = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")], 
   plot           = FALSE)

## No need to run this as is the same then the example below (per habitat)

################################################################################


#### 3 - Compute functional alpha diversity indices ####

### Compute indices:
 
#Eliminating rows that have less than 4 sp
 
 sp_transect_habitat <- sp_transect_habitat[-c(2,10), ]
 
 sp_transect_habitat <- as.matrix(sp_transect_habitat)
 
 alpha_fd_indices_site <- mFD::alpha.fd.multidim(
   sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
   asb_sp_w         = asb_sp_w_hab_s,
   ind_vect         = c("fdis", "fric", "fspe"),
   scaling          = TRUE,
   check_input      = TRUE,
   details_returned = TRUE)
 
 fd_values_site <- alpha_fd_indices_site$functional_diversity_indices
 
 #add habitat
 
 fd_values_site$habitat <- (c("I","I", "I","M","M", "M",
                              "O", "O","O"))
 

 alpha_fd_indices_transect <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = sp_transect_habitat,
  ind_vect         = c("fdis", "fric", "fspe"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_values_transect <- alpha_fd_indices_transect$functional_diversity_indices

#add habitat

fd_values_transect$habitat <- (c("I","I", "I", "I", "I", "I", 
 "I", "I", "I", "I", "M","M", "M", "M", "M", "M", 
 "M","M", "M", "M", "M", "M", "O", "O","O", "O", "O", "O",
 "O", "O","O", "O", "O", "O"))

#write.csv(fd_values_transect, file="data/Diversity_spatial_transect.csv")

################################################################################



#### 4 - Compute functional beta diversity indices ####


### Compute beta diversity indices:

beta_fd_indices_s <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = asb_sp_occ_s,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

beta_fd_indices_s$pairasb_fbd_indices$jac_diss

#Eliminating rows that have less than 4 sp

sp_occ_transect <- sp_occ_transect[-c(2,10), ]

sp_occ_transect <- as.matrix(sp_occ_transect)

beta_fd_indices_transect <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_occ       = sp_occ_transect,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

beta_fd_indices_transect$pairasb_fbd_indices$jac_diss




################################################################
###############################################################
#2nd part - Taxonomic analysis

#rm(list=ls()) # cleaning memory

#Load packages

library (tidyr)
library(dplyr)
library(stringr)
library (forcats)
library(rsq)
library(margins)
library(betapart)
library(vegan)
library(tidyselect)
library(reshape)


### load all data

species_site <- asb_sp_occ_s
species_transect <- sp_occ_transect

# Sample metadata and explanatory variables ----
meta_data_site <- read.csv(here::here("data", "SpatialUVC_metadata_site.csv")) %>% 
  mutate(Replicate = Code) %>% 
  mutate(Site = str_sub(Replicate, end = -3)) %>% 
  `rownames<-`(.$Code)

meta_data_transect <- meta_data <- read.csv(here::here("data", "SpatialUVC_metadata_transect.csv")) %>% 
  mutate(Replicate = Code) %>% 
  mutate(Site = str_sub(Replicate, end = -3)) %>% 
  `rownames<-`(.$Code)

#### Species dissimilarity (beta diversity) between habitats ####

# calculate taxonomic beta diversity - jaccard dissimilarity and 
#turnover 

totalBETA_tax_spatial_site <- beta.pair(species_site, index.family = "jaccard")
totaldis_tax_spatial_site <- as.matrix(totalBETA_tax_spatial_site$beta.jac)
diag(totaldis_tax_spatial_site) <- NA
turnover_tax_spatial_site <- as.matrix(totalBETA_tax_spatial_site$beta.jtu)
diag(turnover_tax_spatial_site) <- NA


totalBETA_tax_spatial_transect <- beta.pair(species_transect, index.family = "jaccard")
turnover_tax_spatial_transect <- as.matrix(totalBETA_tax_spatial_transect$beta.jtu)
diag(turnover_tax_spatial_transect) <- NA
totaldis_tax_spatial_transect <- as.matrix(totalBETA_tax_spatial_transect$beta.jac)
diag(totaldis_tax_spatial_transect) <- NA

