################################################################################
##
## Script for plotting changes in functional richness across habitats and over years
## in a multidimensional space
## 
## 
## Code by Camille Magneville, Sebastien Villeger and Paula Sgarlatta
##
################################################################################


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(patchwork)
library(mFD)

# loading data
load(here::here("outputs", "sp_3D_coord.RData") )

load(here::here("data", "spatial_metadata.RData") )
load(here::here("data", "using biomass-maxN", "spatial_sp_biom.RData") )
load(here::here("outputs/", "using biomass-maxN", "spatial_fd_biomass.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "using biomass-maxN", "kelp_sp_maxN.RData") )
load(here::here("outputs/", "using biomass-maxN", "temporal_fd_kelp_biomass.RData") )


## settings ####

# vertices of all fe in 4D ----
pool_vert_nm <- spatial_fd$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(sp_3D_coord)
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]

# loading thermal affinity data:
thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  #column_to_rownames("Species") %>% 
  select(-thermal)

#Colors for plots

thermal_aff_colors <- c(tropical = "lightsalmon1", temperate = "#2C6BAA")

## temporal kelp ####

# computing occurrences of species in each year (we will use 2002-2009-2018)

kelp_years_sp_occ <- rbind( 
  y2002 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2002"),"Code"],],2,max ),
  y2009 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2009"),"Code"],],2,max ),
  y2018 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2018"),"Code"],],2,max ))  


# # Retrieve species coordinates matrix for year 2002:
kelp_2002_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2002") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2002_occ_2 <- kelp_2002_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2002_occ_2, thermal, by= "Species")

kelp_2002_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()


## Compute FRic values #### 

# compute FRic for all habitats  ---
 Fric <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
                          asb_sp_w = kelp_2002_occ_thermal,
                          ind_vect = c("fric"),
                          scaling = TRUE,
                          details_returned = TRUE)

## Because there are not enough tropical species, I'll only do it in 2 axes

sp_2d_coord <- sp_3D_coord[, c("PC1", "PC2")]

sp_filter_2002 <- mFD::sp.filter(asb_nm          = c("y2002"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final <- sp_filter_2002$`species coordinates`[, c("PC1", "PC2")]

sp_thermal <- sp_2d_coord_final %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

# compute FRic for all habitats  ---
Fric <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                          asb_sp_w = kelp_2002_occ_thermal,
                          ind_vect = c("fric"),
                          scaling = TRUE,
                          details_returned = TRUE)


## plotting  ####

# list to store ggplot
ggplot_2002 <- list()

# pairs of axes
pairs_axes <- list(c(1,2))

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC2"), 
                              color_bg = "grey95")
  
  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                        sp_coord2D = sp_3D_coord,
                        vertices_nD = pool_vert_nm,
                        plot_pool = FALSE,
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
  
 
  # species present in trop:
  sp_trop <- sp_thermal$Species[which(sp_thermal$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp <- sp_thermal$Species[which(sp_thermal$thermal_label == "temperate")]
  
  # vertices in trop:
  vert_trop <- Fric$details$asb_vert_nm$tropical
  
  # vertices in temp:
  vert_temp <- Fric$details$asb_vert_nm$temperate
  
    # plot convex hull of assemblage but not species
  ggplot_z2 <-fric.plot(ggplot_bg = ggplot_z, 
                        asb_sp_coord2D = list(asb1 = sp_2d_coord_final[sp_trop, xy], 
                                              asb2 = sp_2d_coord_final[sp_temp, xy]),
                        asb_vertices_nD = list(asb1 = vert_trop, 
                                               asb2 = vert_temp),
                        plot_sp = TRUE,
                        color_sp = thermal_aff_colors,
                        fill_sp = c(asb1 = "white", asb2 = "white"),
                        size_sp = c(asb1 = 1, asb2 = 1),
                        shape_sp = c(asb1 = 16, asb2 = 16),
                        color_vert = thermal_aff_colors,
                        fill_vert = thermal_aff_colors,
                        size_vert = c(asb1 = 4, asb2 = 4),
                        shape_vert = c(asb1 = 16, asb2 = 16),
                        alpha_ch = c(asb1 = 0, asb2 = 0),
                        color_ch = c(asb1 = NA, asb2 =" red"),
                        fill_ch = c(asb1 = NA, asb2 = NA))

  
# legend and title
if (z==1) {
    ggplot_z2  <- ggplot_z2  + labs(title = "2002" )+
    theme(plot.title = element_text(size = 20, color = "red"))
}
  
  # ggplot stored in list
  ggplot_2002[[z]] <- ggplot_z2
  
}# end of z

################ 2009

# # Retrieve species coordinates matrix for year 2009:

kelp_2009_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2009") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2009_occ_2 <- kelp_2009_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2009_occ_2, thermal, by= "Species")

kelp_2009_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2009 <- mFD::sp.filter(asb_nm          = c("y2009"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2009 <- sp_filter_2009$`species coordinates`[, c("PC1", "PC2")]

sp_thermal_2009 <- sp_2d_coord_final_2009 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric_2009 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                          asb_sp_w = kelp_2009_occ_thermal,
                          ind_vect = c("fric"),
                          scaling = TRUE,
                          details_returned = TRUE)

## plotting  ####

# list to store ggplot
ggplot_2009 <- list()

# pairs of axes
pairs_axes <- list( c(1,2) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
    # species present in trop:
  sp_trop_2009 <- sp_thermal_2009$Species[which(sp_thermal_2009$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp_2009 <- sp_thermal_2009$Species[which(sp_thermal_2009$thermal_label == "temperate")]
  
  # vertices in trop:
  vert_trop_2009 <- Fric_2009$details$asb_vert_nm$tropical
  
  # vertices in temp:
  vert_temp_2009 <- Fric_2009$details$asb_vert_nm$temperate
  
  # plot convex hull of assemblage but not species
  ggplot_z_2009 <-fric.plot(ggplot_bg = ggplot_z, 
                        asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2009[sp_trop_2009, xy], 
                                              asb2 = sp_2d_coord_final_2009[sp_temp_2009, xy]),
                        asb_vertices_nD = list(asb1 = vert_trop_2009, 
                                               asb2 = vert_temp_2009),
                        plot_sp = TRUE,
                        color_sp = thermal_aff_colors,
                        fill_sp = c(asb1 = "white", asb2 = "white"),
                        size_sp = c(asb1 = 1, asb2 = 1),
                        shape_sp = c(asb1 = 16, asb2 = 16),
                        color_vert = thermal_aff_colors,
                        fill_vert = thermal_aff_colors,
                        size_vert = c(asb1 = 4, asb2 = 4),
                        shape_vert = c(asb1 = 16, asb2 = 16),
                        alpha_ch = c(asb1 = 0, asb2 = 0),
                        color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
                        fill_ch = c(asb1 = NA, asb2 = NA))
  
    # legend and title
  if (z==1) {
    ggplot_z_2009 <- ggplot_z_2009  + labs(title = "2009" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # ggplot stored in list
  ggplot_2009[[z]] <- ggplot_z_2009
  
}# end of z

################ 2018

# # Retrieve species coordinates matrix for year 2018:

kelp_2018_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2018") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2018_occ_2 <- kelp_2018_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2018_occ_2, thermal, by= "Species")

kelp_2018_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2018 <- mFD::sp.filter(asb_nm          = c("y2018"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2018 <- sp_filter_2018$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2018 <- sp_2d_coord_final_2018 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric_2018 <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
                               asb_sp_w = kelp_2018_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)

## plotting  ####

# list to store ggplot
ggplot_2018 <- list()

# pairs of axes
pairs_axes <- list( c(1,2) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
    # species present in trop:
  sp_trop_2018 <- sp_thermal_2018$Species[which(sp_thermal_2018$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp_2018 <- sp_thermal_2018$Species[which(sp_thermal_2018$thermal_label == "temperate")]
  
  # vertices in trop:
  vert_trop_2018 <- Fric_2018$details$asb_vert_nm$tropical
  
  # vertices in temp:
  vert_temp_2018 <- Fric_2018$details$asb_vert_nm$temperate
  
  # plot convex hull of assemblage but not species
  ggplot_z_2018 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2018[sp_trop_2018, xy], 
                                                  asb2 = sp_2d_coord_final_2018[sp_temp_2018, xy]),
                            asb_vertices_nD = list(asb1 = vert_trop_2018, 
                                                   asb2 = vert_temp_2018),
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 4, asb2 = 4),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = NA, asb2 ="#00C19A"),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2018 <- ggplot_z_2018  + labs(title = "2018" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # ggplot stored in list
  ggplot_2018[[z]] <- ggplot_z_2018
  
}# end of z
  
  
######################## SPATIAL ####

# computing occurrences of species in each habitat
hab_sp_occ <- rbind( 
  Inshore = apply(spatial_sp_biom [spatial_metadata[which(spatial_metadata$Habitat=="Inshore"),"Code"],],2,max ),
  Midshelf = apply(spatial_sp_biom [spatial_metadata[which(spatial_metadata$Habitat=="Midshelf"),"Code"],],2,max ),
  Offshore = apply(spatial_sp_biom [spatial_metadata[which(spatial_metadata$Habitat=="Offshore"),"Code"],],2,max )
)  

########## INSHORE

# # Retrieve species coordinates matrix for inshore:

inshore_occ <- hab_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "Inshore") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

inshore_occ_2 <- inshore_occ %>% 
as.data.frame() %>% 
  gather(Species, Abundance, 1:53)

add_thermal <- left_join(inshore_occ_2, thermal, by= "Species")

inshore_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_inshore <- mFD::sp.filter(asb_nm          = c("Inshore"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = hab_sp_occ)

sp_2d_coord_final_inshore <- sp_filter_inshore$`species coordinates`[, c("PC1", "PC2")]

sp_thermal_inshore <- sp_2d_coord_final_inshore %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric_inshore <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = inshore_occ,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)

#### Because there are only 2 tropical species, I cannot compute tropical species


## plotting  ####

# list to store ggplot
ggplot_inshore <- list()

# pairs of axes
pairs_axes <- list( c(1,2) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
   # species present in trop:
  #sp_trop_inshore <- sp_thermal_inshore$Species[which(sp_thermal_inshore$thermal_label == "tropical")]
  
  # species present in temp:
  #sp_temp_inshore <- sp_thermal_inshore$Species[which(sp_thermal_inshore$thermal_label == "temperate")]
  
    # vertices in trop:
  #vert_trop_inshore <- Fric_inshore$details$asb_vert_nm$tropical
  
  # vertices in temp:
  #vert_temp_inshore <- Fric_inshore$details$asb_vert_nm$temperate
  
  # vertices
  vert_inshore <- Fric_inshore$details$asb_vert_nm$Inshore
  
  # plot convex hull of assemblage but not species
  ggplot_z_inshore <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = list("Inshore" = sp_2d_coord_final_inshore),
                            asb_vertices_nD = list("Inshore"= vert_inshore),
                            plot_sp = TRUE,
                            color_sp = c("Inshore" = "#2C6BAA"),
                            fill_sp = c("Inshore"= "white"),
                            size_sp = c("Inshore" = 1),
                            shape_sp = c("Inshore" = 16),
                            color_vert = c("Inshore" = "#2C6BAA"),
                            fill_vert = c("Inshore" = "#2C6BAA"),
                            size_vert = c("Inshore"= 4),
                            shape_vert = c("Inshore" = 16),
                            alpha_ch = c("Inshore" = 0),
                            color_ch = c("Inshore" = "#2C6BAA"),
                            fill_ch = c("Inshore" = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_inshore <- ggplot_z_inshore  + labs(title = "Inshore" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # ggplot stored in list
  ggplot_inshore[[z]] <- ggplot_z_inshore
  
}# end of z


################## MIDSHELF

# # Retrieve species coordinates matrix for midshelf:

midshelf_occ <- hab_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "Midshelf") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

midshelf_occ_2 <- midshelf_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:53)

add_thermal <- left_join(midshelf_occ_2, thermal, by= "Species")

midshelf_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_midshelf <- mFD::sp.filter(asb_nm          = c("Midshelf"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = hab_sp_occ)

sp_2d_coord_final_midshelf <- sp_filter_midshelf$`species coordinates`[, c("PC1", "PC2")]

sp_thermal_midshelf <- sp_2d_coord_final_midshelf %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric_midshelf <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
                               asb_sp_w = midshelf_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)

## plotting  ####

# list to store ggplot
ggplot_midshelf <- list()

# pairs of axes
pairs_axes <- list( c(1,2) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # species present in trop:
  sp_trop_midshelf <- sp_thermal_midshelf$Species[which(sp_thermal_midshelf$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp_midshelf <- sp_thermal_midshelf$Species[which(sp_thermal_midshelf$thermal_label == "temperate")]
  
  # vertices in trop:
  vert_trop_midshelf <- Fric_midshelf$details$asb_vert_nm$tropical
  
  # vertices in temp:
  vert_temp_midshelf <- Fric_midshelf$details$asb_vert_nm$temperate
  
  #Vertices
  
  vert_midshelf <- Fric_midshelf$details$asb_vert_nm$Midshelf
  
  # plot convex hull of assemblage but not species
  ggplot_z_midshelf <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = list(asb1 = sp_2d_coord_final_midshelf[sp_trop_midshelf, xy], 
                                                  asb2 = sp_2d_coord_final_midshelf[sp_temp_midshelf, xy]),
                            asb_vertices_nD = list(asb1 = vert_trop_midshelf, 
                                                   asb2 = vert_temp_midshelf),
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 4, asb2 = 4),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = NA, asb2 ="lightsalmon1"),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_midshelf <- ggplot_z_midshelf  + labs(title = "Midshelf" )+
      theme(plot.title = element_text(size = 20, color = "lightsalmon1"))
  }
  
  # ggplot stored in list
  ggplot_midshelf[[z]] <- ggplot_z_midshelf
  
}# end of z


################## OFFSHORE

# # Retrieve species coordinates matrix for midshelf:

offshore_occ <- hab_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "Offshore") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

offshore_occ_2 <- offshore_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:53)

add_thermal <- left_join(offshore_occ_2, thermal, by= "Species")

offshore_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_offshore <- mFD::sp.filter(asb_nm          = c("Offshore"), 
                                     sp_faxes_coord = sp_3D_coord, 
                                     asb_sp_w       = hab_sp_occ)

sp_2d_coord_final_offshore <- sp_filter_offshore$`species coordinates`[, c("PC1", "PC2")]

sp_thermal_offshore <- sp_2d_coord_final_offshore %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric_offshore <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
                                   asb_sp_w = offshore_occ_thermal,
                                   ind_vect = c("fric"),
                                   scaling = TRUE,
                                   details_returned = TRUE)

## plotting  ####

# list to store ggplot
ggplot_offshore <- list()

# pairs of axes
pairs_axes <- list( c(1,2) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # species present in trop:
  sp_trop_offshore <- sp_thermal_offshore$Species[which(sp_thermal_offshore$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp_offshore <- sp_thermal_offshore$Species[which(sp_thermal_offshore$thermal_label == "temperate")]
  
  # vertices in trop:
  vert_trop_offshore <- Fric_offshore$details$asb_vert_nm$tropical
  
  # vertices in temp:
  vert_temp_offshore <- Fric_offshore$details$asb_vert_nm$temperate
  
  #Vertices
  
  vert_offshore <- Fric_offshore$details$asb_vert_nm$Offshore
  
  # plot convex hull of assemblage but not species
  ggplot_z_offshore <-fric.plot(ggplot_bg = ggplot_z, 
                                asb_sp_coord2D = list(asb1 = sp_2d_coord_final_offshore[sp_trop_offshore, xy], 
                                                      asb2 = sp_2d_coord_final_offshore[sp_temp_offshore, xy]),
                                asb_vertices_nD = list(asb1 = vert_trop_offshore, 
                                                       asb2 = vert_temp_offshore),
                                plot_sp = TRUE,
                                color_sp = thermal_aff_colors,
                                fill_sp = c(asb1 = "white", asb2 = "white"),
                                size_sp = c(asb1 = 1, asb2 = 1),
                                shape_sp = c(asb1 = 16, asb2 = 16),
                                color_vert = thermal_aff_colors,
                                fill_vert = thermal_aff_colors,
                                size_vert = c(asb1 = 4, asb2 = 4),
                                shape_vert = c(asb1 = 16, asb2 = 16),
                                alpha_ch = c(asb1 = 0, asb2 = 0),
                                color_ch = c(asb1 = NA, asb2 ="firebrick3"),
                                fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_offshore <- ggplot_z_offshore  + labs(title = "Offshore" )+
      theme(plot.title = element_text(size = 20, color = "firebrick3"))
  }
  
  # ggplot stored in list
  ggplot_offshore[[z]] <- ggplot_z_offshore
  
}# end of z

## merging all plots into a single figure and saving as png ####

figure9 <- ( ggplot_2002[[1]] + ggplot_2009[[1]] +  ggplot_2018[[1]] )/
   ( ggplot_inshore[[1]]  + ggplot_midshelf[[1]] +  ggplot_offshore[[1]] )


ggsave(figure9, file=here::here("outputs/", "using biomass-maxN", "Figure9_biomass.png"),
       height = 16, width = 24, unit = "cm" )


#################################################################################################################

############ Without tropical/temperate
## settings ####

## temporal kelp ####

# computing occurrences of species in each year (we will use 2002-2010-2018)

kelp_years_sp_occ <- rbind( 
  y2002 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2002"),"Code"],],2,max ),
  y2010 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2010"),"Code"],],2,max ),
  y2018 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2018"),"Code"],],2,max ))  


sp_filter_2002 <- mFD::sp.filter(asb_nm          = c("y2002"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_faxes_coord_3D <- sp_filter_2002$`species coordinates`[, c("PC1", "PC2")]

# Set faxes limits:
# set range of axes if c(NA, NA):
range_sp_coord  <- range(sp_3D_coord)
range_faxes_lim <- range_sp_coord + c(-1, 1)*(range_sp_coord[2] - 
                                                range_sp_coord[1]) * 0.05

# coordinates of all species ----
pool_coord<-spatial_fd$details$sp_faxes_coord[,c("PC1", "PC2")] 

# vertices of all fe in 4D ----
pool_vert_nm<-spatial_fd$details$pool_vert_nm

# Retrieve the background plot:
ggplot_bg <- mFD::background.plot(
  range_faxes = range_faxes_lim, 
  faxes_nm    = c("PC 1", "PC 2"), 
  color_bg    = "grey90") 


# convex hull of all species pool

ggplot_z<-pool.plot(ggplot_bg=ggplot_bg,
                    sp_coord2D=pool_coord,
                    vertices_nD=pool_vert_nm,
                    plot_pool=FALSE,
                    color_ch=NA, fill_ch="white", alpha_ch=1)

# Retrieve vertices names:
vert_nm <- vertices(sp_faxes_coord_3D, 
                    order_2D = TRUE, check_input = TRUE)

# compute FRic for kelp 2002  ---

ggplot_fric <- mFD::fric.plot(
  ggplot_bg       = ggplot_z,
  asb_sp_coord2D  = list(y2002 = sp_faxes_coord_3D),
  asb_vertices_nD = list(y2002 = vert_nm),
  plot_sp         = TRUE,
  color_ch        = c("y2002" = "red"), 
  fill_ch         = c("y2002" = "white"), 
  alpha_ch        = c("y2002" = 0.3),
  size_sp = c("y2002" = 1),
  shape_sp = c("y2002" = 16),
  color_sp = c("y2002" = "red"),
  fill_sp = c("y2002" = "red"),
  size_vert = c("y2002" = 1),
  color_vert = c("y2002" = "red"),
  fill_vert = c("y2002" = "red"),
  shape_vert = c("y2002" = 16))
ggplot_fric

ggsave(ggplot_fric, file=here::here("outputs/", "kelp_2002.png"),
       height = 16, width = 24, unit = "cm" )
