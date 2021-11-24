################################################################################
##
## Script for importing all data from csv files
## 
## Code by Camille Magneville, Sébastien villéger and Paula Sgarlatta
##
## 1/ load and preparing trait datasets - kelp and no kelp sites
## 2/ 
################################################################################

#### 1 - Load (and transform) datasets ####

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(mFD)

## 1 trait data ####

# load traits data ----

# species from sites with kelp
traits_kelp <- read.csv("data/TemporalBRUV_species_traits_kelp.csv",
                        header = TRUE, row.names = 1)
head(traits_kelp)

#Traits for sp presents in sites that never had kelp
traits_no_kelp <- read.csv("data/TemporalBRUV_species_traits_no_kelp.csv",
                           header = TRUE, row.names = 1)
head(traits_no_kelp)

# merging in a single data.frame and trait named renamed
sp_tr <- as.data.frame( rbind( traits_kelp, traits_no_kelp) )
names(sp_tr) <- c("Size", "Aggr", "Posi", "Diet")
head(sp_tr)

# recoding variable to match trait type ---

# looking at trait values
lapply(sp_tr, unique)

# trait type
tr_cat<-data.frame( trait_name = c("Size", "Aggr", "Posi", "Diet"),
                    trait_type = c("O","O","O","N") )

# size as ordinal
sp_tr$Size <- factor(sp_tr$Size, 
                     levels = c("S2", "S3", "S4", "S5", "S6", "S7"),
                     ordered = TRUE)
summary(sp_tr$Size)

# aggregation as ordinal
sp_tr$Aggr <- factor(sp_tr$Aggr, 
                    levels = c("Solitary", "Pair", "Group"),
                    ordered = TRUE)
summary(sp_tr$Aggr)

# Position as ordinal
sp_tr$Posi <- factor(sp_tr$Posi, 
                    levels = c("Benthic", "BenthoP", "Pelagic"),
                    ordered = TRUE)
summary(sp_tr$Posi)

# diet as factor
sp_tr$Diet <- as.factor(sp_tr$Diet)
summary(sp_tr$Diet)

summary(sp_tr)

####################### # => trait data ready ####

## 2 species occurrences data ####

# load temporal maxN data - only sites that use to have kelp and lost it:
kelp <- read.csv("data/TemporalBRUV_species_maxN_kelp.csv")
rownames(kelp) <- kelp$species_id
kelp <- tibble::column_to_rownames(kelp, var = "Site")

# load trait types data - same for both datasets (kelp and no kelp):
traits_type <- read.csv("data/TemporalBRUV_traits_type.csv")

#Loading data from sites that never had kelp

no_kelp <- read.csv("data/TemporalBRUV_species_maxN_no_kelp.csv")

rownames(no_kelp) <- no_kelp$species_id
no_kelp <- tibble::column_to_rownames(no_kelp, var = "Site")




library(remotes)
#remotes::install_github("CmlMagneville/mFD")
library(mFD)

## summary for sp*tr data:

sp_tr_summ_kelp <- mFD::sp.tr.summary(tr_cat = traits_type, 
                                 sp_tr  = traits_kelp)
sp_tr_summ_no_kelp <- mFD::sp.tr.summary(tr_cat = traits_type, 
                                 sp_tr  = traits_no_kelp)

# repartition of traits:

sp_tr_summ_kelp$tr_summary_list

sp_tr_summ_no_kelp$tr_summary_list


## summary for asb data:

kelp <- as.matrix(kelp)

no_kelp <- as.matrix(no_kelp)

asb_sp_summ_kelp <- mFD::asb.sp.summary(asb_sp_w = kelp)

asb_sp_summ_no_kelp <- mFD::asb.sp.summary(asb_sp_w = no_kelp)

# retrieve occurrence matrix:
asb_sp_occ_kelp <- asb_sp_summ_kelp$asb_sp_occ

asb_sp_occ_no_kelp <- asb_sp_summ_no_kelp$asb_sp_occ

# total maxN of each species:
asb_sp_summ_kelp$sp_tot_w

asb_sp_summ_no_kelp$sp_tot_w

# total maxN for each asb (ie replicate):
asb_sp_summ_kelp$asb_tot_w

asb_sp_summ_no_kelp$asb_tot_w

# species richness per replicate:
asb_sp_summ_kelp$asb_sp_richn

asb_sp_summ_no_kelp$asb_sp_richn

# list of species present in each assemblage (ie replicate):
asb_sp_summ_kelp$asb_sp_nm

asb_sp_summ_no_kelp$asb_sp_nm


### Compute functional distance:

# compute functional distance based on traits:
sp_dist_kelp <- mFD::funct.dist(
  sp_tr  = traits_kelp,
  tr_cat = traits_type,
  metric = "gower",
  scale_euclid  = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA  = TRUE)

# is there a lot of distances = 0? (group into FEs or not?)
sp_dist_kelp2 <- as.matrix(sp_dist_kelp)
inds_kelp <- which(sp_dist_kelp2 == min(sp_dist_kelp2), arr.ind=TRUE)
inds_kelp

##PS: Not sure about this, can we double-check again? Do we 
## need to gather these into FEntities?


sp_dist_no_kelp <- mFD::funct.dist(
  sp_tr  = traits_no_kelp,
  tr_cat = traits_type,
  metric = "gower",
  scale_euclid  = "scale_center",
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA  = TRUE)


sp_dist_no_kelp2 <- as.matrix(sp_dist_no_kelp)
inds_no_kelp <- which(sp_dist_no_kelp2 == min(sp_dist_no_kelp2), arr.ind=TRUE)
inds_no_kelp
##PS: Not sure about this, can we double-check again? Do we 
## need to gather these into FEntities?


### Compute functional spaces and their quality:

# compute the quality of functional spaces:
fspaces_quality_kelp <- mFD::quality.fspaces(
  sp_dist             = sp_dist_kelp, 
  maxdim_pcoa         = 10, 
  deviation_weighting = "absolute", 
  fdist_scaling       = FALSE, 
  fdendro             = "average")
# ... same warning message saying that some distance equal 0 but we have ...
# ... checked and it is ok! [PS] please double-check!

fspaces_quality_no_kelp <- mFD::quality.fspaces(
  sp_dist             = sp_dist_no_kelp, 
  maxdim_pcoa         = 10, 
  deviation_weighting = "absolute", 
  fdist_scaling       = FALSE, 
  fdendro             = "average")

# ... same warning message saying that some distance equal 0 but we have ...
# ... checked and it is ok! [PS] please double-check!

# retrieve table of quality metric:
fspaces_quality_kelp$quality_fspaces

#[PS]Should we use 3 dimensions?

fspaces_quality_no_kelp$quality_fspaces

#[PS]Should we use 3 dimensions?

# illustrate quality of functional space (not useful but can be used as ...
# ... Supplementary Figure for instance):
funct_space_kelp <- mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_kelp,
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

funct_space_kelp

##[PS] Is this correct? It looks diferent from the first one I obtained

funct_space_no_kelp <- mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality_no_kelp,
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

funct_space_no_kelp

### Test correlation between traits and functional axes:

# retrieve coordinates of species:

sp_faxes_coord_kelp <- fspaces_quality_kelp$"details_fspaces"$"sp_pc_coord"

sp_faxes_coord_no_kelp <- fspaces_quality_no_kelp$"details_fspaces"$"sp_pc_coord"

# test correlation between traits and axes:
cor_tr_faxes_kelp <- mFD::traits.faxes.cor(
  sp_tr          = traits_kelp, 
  sp_faxes_coord = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

cor_tr_faxes_no_kelp <- mFD::traits.faxes.cor(
  sp_tr          = traits_no_kelp, 
  sp_faxes_coord = sp_faxes_coord_no_kelp[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# get the table of correlation:

cor_tr_faxes_kelp$tr_faxes_stat

cor_tr_faxes_no_kelp$tr_faxes_stat

# get the plot:

cor_tr_faxes_kelp$tr_faxes_plot

cor_tr_faxes_no_kelp$tr_faxes_plot


### Plot functional space:

big_plot_kelp <- mFD::funct.space.plot(sp_faxes_coord  = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")],
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
big_plot_kelp$patchwork

big_plot_no_kelp <- mFD::funct.space.plot(sp_faxes_coord  = sp_faxes_coord_no_kelp[, c("PC1", "PC2", "PC3")],
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
big_plot_no_kelp$patchwork

################################################################################


#### 3 - Compute and plot functional alpha diversity indices ####

### Compute indices:

alpha_fd_indices_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")],
  asb_sp_w         = kelp,
  ind_vect         = c("fdis", "fspe", "fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# get indices table:
fd_temporal_kelp <- alpha_fd_indices_kelp$functional_diversity_indices

#add year as column 

library(dplyr)
library(stringr)

# Sample metadata and explanatory variables ----

meta_data_kelp <- read.csv(here::here("data", "TemporalBRUV_species_metadata_kelp.csv")) %>% 
  mutate(Replicate = Code) %>% 
  mutate(Site = str_sub(Replicate, end = -7)) %>% 
  `rownames<-`(.$Code)

meta_data_kelp$Richness <- asb_sp_summ_kelp$asb_sp_richn

fd_temporal_kelp$year <- meta_data_kelp$Year


alpha_fd_indices_no_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_no_kelp[, c("PC1", "PC2", "PC3")],
  asb_sp_w         = no_kelp,
  ind_vect         = c("fdis", "fspe", "fric"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# get indices table:
fd_temporal_no_kelp <- alpha_fd_indices_no_kelp$functional_diversity_indices

#add year as column 

# Sample metadata and explanatory variables ----

meta_data_no_kelp <- read.csv(here::here("data", "TemporalBRUV_species_metadata_no_kelp.csv")) %>% 
  mutate(Replicate = Code) %>% 
  mutate(Site = str_sub(Replicate, end = -7)) %>% 
  `rownames<-`(.$Code)

meta_data_no_kelp$Richness <- asb_sp_summ_no_kelp$asb_sp_richn

fd_temporal_no_kelp$year <- meta_data_no_kelp$Year

### Compute beta diversity indices:

beta_fd_indices_kelp <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_kelp[, c("PC1", "PC2", "PC3")],
  asb_sp_occ       = asb_sp_occ_kelp,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

beta_fd_indices_kelp$pairasb_fbd_indices

# Obtaining matrix of each component of beta diversity

totaldis_func_kelp <- as.matrix(beta_fd_indices_kelp$pairasb_fbd_indices$jac_diss)
diag(totaldis_func_kelp) <- NA
turnover_func_kelp <- as.matrix(beta_fd_indices_kelp$pairasb_fbd_indices$jac_turn)
diag(turnover_func_kelp) <- NA


beta_fd_indices_no_kelp <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord_no_kelp[, c("PC1", "PC2", "PC3")],
  asb_sp_occ       = asb_sp_occ_no_kelp,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

beta_fd_indices_no_kelp$pairasb_fbd_indices

# Obtaining matrix of each component of beta diversity

totaldis_func_no_kelp <- as.matrix(beta_fd_indices_no_kelp$pairasb_fbd_indices$jac_diss)
diag(totaldis_func_no_kelp) <- NA
turnover_func_no_kelp <- as.matrix(beta_fd_indices_no_kelp$pairasb_fbd_indices$jac_turn)
diag(turnover_func_no_kelp) <- NA

#############################################################
###Taxonomic diversity

#rm(list=ls()) # cleaning memory

#Load packages

library(betapart)
library(vegan)


### load all data

species_kelp <- asb_sp_occ_kelp

#    make sure species and site data are in the same order
rownames(species_kelp) == rownames(meta_data_kelp)
meta_data_kelp <- meta_data_kelp[rownames(species_kelp), ]

#### Species dissimilarity (beta diversity) between years ####

# calculate taxonomic beta diversity - jaccard dissimilarity and 
#turnover components

totalBETA_tax_kelp <- beta.pair(species_kelp, index.family = "jaccard")
totaldis_tax_kelp <- as.matrix(totalBETA_tax_kelp$beta.jac)
diag(totaldis_tax_kelp) <- NA
turnover_tax_kelp <- as.matrix(totalBETA_tax_kelp$beta.jtu)
diag(turnover_tax_kelp) <- NA

species_no_kelp <- asb_sp_occ_no_kelp

#    make sure species and site data are in the same order
rownames(species_no_kelp) == rownames(meta_data_no_kelp)
meta_data_no_kelp <- meta_data_no_kelp[rownames(species_no_kelp), ]

#### Species dissimilarity (beta diversity) between years ####

# calculate taxonomic beta diversity - jaccard dissimilarity and 
#turnover components

totalBETA_tax_no_kelp <- beta.pair(species_no_kelp, index.family = "jaccard")
totaldis_tax_no_kelp <- as.matrix(totalBETA_tax_no_kelp$beta.jac)
diag(totaldis_tax_no_kelp) <- NA
turnover_tax_no_kelp <- as.matrix(totalBETA_tax_no_kelp$beta.jtu)
diag(turnover_tax_no_kelp) <- NA
