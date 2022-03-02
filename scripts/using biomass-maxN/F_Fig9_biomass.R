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
library(ggpubr)

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

thermal_aff_colors <- c(temperate = "#2C6BAA", tropical = "lightsalmon1")

## temporal kelp ####

# computing occurrences of species in each year 

kelp_years_sp_occ <- rbind( 
  y2002 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2002"),"Code"],],2,max ),
  y2003 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2003"),"Code"],],2,max ),
  y2004 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2004"),"Code"],],2,max ),
  y2005 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2005"),"Code"],],2,max ),
  y2006 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2006"),"Code"],],2,max ),
  y2007 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2007"),"Code"],],2,max ),
  y2008 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2008"),"Code"],],2,max ),
  y2009 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2009"),"Code"],],2,max ),
  y2010 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2010"),"Code"],],2,max ),
  y2011 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2011"),"Code"],],2,max ),
  y2013 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2013"),"Code"],],2,max ),
  y2015 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2015"),"Code"],],2,max ),
  y2018 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2018"),"Code"],],2,max ))  

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

sp_2d_coord <- sp_3D_coord[, c("PC1", "PC2")]

##################### 2002 #################################

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

sp_filter_2002 <- mFD::sp.filter(asb_nm          = c("y2002"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2002 <- sp_filter_2002$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2002 <- sp_2d_coord_final_2002 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2002 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2002_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2002$functional_diversity_indices
fd_details <- Fric_2002$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2002$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2002$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Epinephelus_fasciatus","Gymnothorax_thyrsoidea","Thalassoma_lutescens", 
                 "Austrolabrus_maculatus", "Brachaelurus_waddi" ,"Centropogon_australis","Orectolobus_halei",
                "Scorpis_lineolata", "Seriola_lalandi" , "Trachinops_taeniatus" )

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2002 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2002 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2002 <- ggplot_z_2002  + labs(title = "2002" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2002 <- ggplot_z_2002 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2002 <- ggplot_z_2002  + labs(title = "2002" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2002[[z]] <- ggplot_z_2002
  
}# end of z









# # # Retrieve species coordinates matrix for year 2002:
# kelp_2002_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2002") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2002_occ_2 <- kelp_2002_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2002_occ_2, thermal, by= "Species")
# 
# kelp_2002_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
#  Fric <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
#                           asb_sp_w = kelp_2002_occ_thermal,
#                           ind_vect = c("fric"),
#                           scaling = TRUE,
#                           details_returned = TRUE)
# 
# ## Because there are not enough tropical species, I'll only do it in 2 axes
# 
# sp_filter_2002 <- mFD::sp.filter(asb_nm          = c("y2002"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final <- sp_filter_2002$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal <- sp_2d_coord_final %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# # compute FRic for all habitats  ---
# Fric <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                           asb_sp_w = kelp_2002_occ_thermal,
#                           ind_vect = c("fric"),
#                           scaling = TRUE,
#                           details_returned = TRUE)
# 
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2002 <- list()
# 
# # pairs of axes
# pairs_axes <- list(c(1,2))
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop <- sp_thermal$Species[which(sp_thermal$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp <- sp_thermal$Species[which(sp_thermal$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop <- Fric$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp <- Fric$details$asb_vert_nm$temperate
#   
#     # plot convex hull of assemblage but not species
#   ggplot_z2 <-fric.plot(ggplot_bg = ggplot_z, 
#                         asb_sp_coord2D = list(asb1 = sp_2d_coord_final[sp_trop, xy], 
#                                               asb2 = sp_2d_coord_final[sp_temp, xy]),
#                         asb_vertices_nD = list(asb1 = vert_trop, 
#                                                asb2 = vert_temp),
#                         plot_sp = TRUE,
#                         color_sp = thermal_aff_colors,
#                         fill_sp = c(asb1 = "white", asb2 = "white"),
#                         size_sp = c(asb1 = 1, asb2 = 1),
#                         shape_sp = c(asb1 = 16, asb2 = 16),
#                         color_vert = thermal_aff_colors,
#                         fill_vert = thermal_aff_colors,
#                         size_vert = c(asb1 = 4, asb2 = 4),
#                         shape_vert = c(asb1 = 16, asb2 = 16),
#                         alpha_ch = c(asb1 = 0, asb2 = 0),
#                         color_ch = c(asb1 = NA, asb2 =" red"),
#                         fill_ch = c(asb1 = NA, asb2 = NA))
# 
#   
# # legend and title
# if (z==1) {
#     ggplot_z2  <- ggplot_z2  + labs(title = "2002" )+
#     theme(plot.title = element_text(size = 20, color = "red"))
# }
#   
#   # ggplot stored in list
#   ggplot_2002[[z]] <- ggplot_z2
#   
# }# end of z

################ 2003

# # Retrieve species coordinates matrix for year 2003:

kelp_2003_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2003") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2003_occ_2 <- kelp_2003_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2003_occ_2, thermal, by= "Species")

kelp_2003_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2003 <- mFD::sp.filter(asb_nm          = c("y2003"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2003 <- sp_filter_2003$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2003 <- sp_2d_coord_final_2003 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2003 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2003_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2003$functional_diversity_indices
fd_details <- Fric_2003$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2003$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2003$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi" ,
                 "Chromis_hypsilepis" , "Orectolobus_maculatus",
                 "Pseudocaranx_dentex","Scorpis_lineolata", 
                 "Nelusetta_ayraud", "Siganus_fuscescens",
                 "Thalassoma_lunare", "Thalassoma_lutescens" )

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2003 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2003 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2003 <- ggplot_z_2003  + labs(title = "2003" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2003 <- ggplot_z_2003 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2003 <- ggplot_z_2003  + labs(title = "2003" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2003[[z]] <- ggplot_z_2003
  
}# end of z



# # Retrieve species coordinates matrix for year 2003:

# kelp_2003_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2003") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2003_occ_2 <- kelp_2003_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2003_occ_2, thermal, by= "Species")
# 
# kelp_2003_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2003 <- mFD::sp.filter(asb_nm          = c("y2003"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2003 <- sp_filter_2003$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2003 <- sp_2d_coord_final_2003 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2003 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2003_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2003 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2003 <- sp_thermal_2003$Species[which(sp_thermal_2003$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2003 <- sp_thermal_2003$Species[which(sp_thermal_2003$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2003 <- Fric_2003$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2003 <- Fric_2003$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2003 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2003[sp_trop_2003, xy], 
#                                                   asb2 = sp_2d_coord_final_2003[sp_temp_2003, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2003, 
#                                                    asb2 = vert_temp_2003),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2003 <- ggplot_z_2003  + labs(title = "2003" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2003[[z]] <- ggplot_z_2003
#   
# }# end of z

################ 2004

# # Retrieve species coordinates matrix for year 2004:

kelp_2004_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2004") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2004_occ_2 <- kelp_2004_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2004_occ_2, thermal, by= "Species")

kelp_2004_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2004 <- mFD::sp.filter(asb_nm          = c("y2004"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2004 <- sp_filter_2004$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2004 <- sp_2d_coord_final_2004 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2004 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2004_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2004$functional_diversity_indices
fd_details <- Fric_2004$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2004$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2004$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi" ,
                 "Orectolobus_maculatus","Sarda_australis",
                 "Scorpis_lineolata", "Trachinops_taeniatus" ,
                 "Chaetodon_flavirostris", "Nelusetta_ayraud", 
                 "Paracaesio_xanthura" ,"Sufflamen_chrysopterum",
                 "Thalassoma_lunare", "Thalassoma_lutescens" )

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2004 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2004 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2004 <- ggplot_z_2004  + labs(title = "2004" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2004 <- ggplot_z_2004 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2004 <- ggplot_z_2004  + labs(title = "2004" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2004[[z]] <- ggplot_z_2004
  
}# end of z











# kelp_2004_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2004") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2004_occ_2 <- kelp_2004_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2004_occ_2, thermal, by= "Species")
# 
# kelp_2004_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2004 <- mFD::sp.filter(asb_nm          = c("y2004"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2004 <- sp_filter_2004$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2004 <- sp_2d_coord_final_2004 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2004 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2004_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2004 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2004 <- sp_thermal_2004$Species[which(sp_thermal_2004$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2004 <- sp_thermal_2004$Species[which(sp_thermal_2004$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2004 <- Fric_2004$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2004 <- Fric_2004$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2004 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2004[sp_trop_2004, xy], 
#                                                   asb2 = sp_2d_coord_final_2004[sp_temp_2004, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2004, 
#                                                    asb2 = vert_temp_2004),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2004 <- ggplot_z_2004  + labs(title = "2004" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2004[[z]] <- ggplot_z_2004
#   
# }# end of z

################ 2005

# # Retrieve species coordinates matrix for year 2005:

kelp_2005_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2005") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2005_occ_2 <- kelp_2005_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2005_occ_2, thermal, by= "Species")

kelp_2005_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2005 <- mFD::sp.filter(asb_nm          = c("y2005"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2005 <- sp_filter_2005$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2005 <- sp_2d_coord_final_2005 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2005 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2005_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2005$functional_diversity_indices
fd_details <- Fric_2005$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2005$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2005$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi" ,
                 "Orectolobus_maculatus","Pseudocaranx_dentex",
                 "Scorpis_lineolata", "Trachinops_taeniatus" ,
                 "Lutjanus_russellii", "Nelusetta_ayraud", 
                 "Paracaesio_xanthura" , "Siganus_fuscescens",
                 "Thalassoma_lunare" )

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2005 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2005 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2005 <- ggplot_z_2005  + labs(title = "2005" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2005 <- ggplot_z_2005 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2005 <- ggplot_z_2005  + labs(title = "2005" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2005[[z]] <- ggplot_z_2005
  
}# end of z











# kelp_2005_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2005") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2005_occ_2 <- kelp_2005_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2005_occ_2, thermal, by= "Species")
# 
# kelp_2005_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2005 <- mFD::sp.filter(asb_nm          = c("y2005"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2005 <- sp_filter_2005$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2005 <- sp_2d_coord_final_2005 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2005 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2005_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2005 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2005 <- sp_thermal_2005$Species[which(sp_thermal_2005$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2005 <- sp_thermal_2005$Species[which(sp_thermal_2005$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2005 <- Fric_2005$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2005 <- Fric_2005$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2005 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2005[sp_trop_2005, xy], 
#                                                   asb2 = sp_2d_coord_final_2005[sp_temp_2005, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2005, 
#                                                    asb2 = vert_temp_2005),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2005 <- ggplot_z_2005  + labs(title = "2005" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2005[[z]] <- ggplot_z_2005
#   
# }# end of z

################ 2006

# # Retrieve species coordinates matrix for year 2006:

kelp_2006_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2006") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2006_occ_2 <- kelp_2006_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2006_occ_2, thermal, by= "Species")

kelp_2006_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2006 <- mFD::sp.filter(asb_nm          = c("y2006"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2006 <- sp_filter_2006$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2006 <- sp_2d_coord_final_2006 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2006 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2006_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2006$functional_diversity_indices
fd_details <- Fric_2006$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2006$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2006$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Brachaelurus_waddi" ,"Eupetrichthys_angustipes",
                 "Orectolobus_halei", "Parapercis_stricticeps",
                 "Pseudocaranx_dentex","Abudefduf_bengalensis",
                 "Scorpis_lineolata", "Trachinops_taeniatus" ,
                 "Abudefduf_bengalensis", "Nelusetta_ayraud", 
                 "Seriola_rivoliana", "Siganus_fuscescens",
                 "Stethojulis_interrupta", "Thalassoma_lunare" )

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2006 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2006 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2006 <- ggplot_z_2006  + labs(title = "2006" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2006 <- ggplot_z_2006 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2006 <- ggplot_z_2006  + labs(title = "2006" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2006[[z]] <- ggplot_z_2006
  
}# end of z














# kelp_2006_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2006") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2006_occ_2 <- kelp_2006_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2006_occ_2, thermal, by= "Species")
# 
# kelp_2006_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2006 <- mFD::sp.filter(asb_nm          = c("y2006"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2006 <- sp_filter_2006$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2006 <- sp_2d_coord_final_2006 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2006 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2006_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2006 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2006 <- sp_thermal_2006$Species[which(sp_thermal_2006$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2006 <- sp_thermal_2006$Species[which(sp_thermal_2006$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2006 <- Fric_2006$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2006 <- Fric_2006$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2006 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2006[sp_trop_2006, xy], 
#                                                   asb2 = sp_2d_coord_final_2006[sp_temp_2006, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2006, 
#                                                    asb2 = vert_temp_2006),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2006 <- ggplot_z_2006  + labs(title = "2006" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2006[[z]] <- ggplot_z_2006
#   
# }# end of z

################ 2007

# # Retrieve species coordinates matrix for year 2007:

kelp_2007_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2007") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2007_occ_2 <- kelp_2007_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2007_occ_2, thermal, by= "Species")

kelp_2007_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2007 <- mFD::sp.filter(asb_nm          = c("y2007"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2007 <- sp_filter_2007$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2007 <- sp_2d_coord_final_2007 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2007 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2007_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2007$functional_diversity_indices
fd_details <- Fric_2007$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2007$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2007$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi" ,
                 "Carcharias_taurus", "Kyphosus_sydneyanus",
                 "Orectolobus_maculatus", "Scorpis_lineolata",
                 "Trachinops_taeniatus","Aulostomus_chinensis",
                 "Chaetodon_flavirostris", "Echeneis_naucrates",
                 "Lutjanus_russellii", "Pterocaesio_digramma",
                 "Thalassoma_lunare")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2007 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2007 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2007 <- ggplot_z_2007  + labs(title = "2007" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2007 <- ggplot_z_2007 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2007 <- ggplot_z_2007  + labs(title = "2007" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2007[[z]] <- ggplot_z_2007
  
}# end of z










# kelp_2007_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2007") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2007_occ_2 <- kelp_2007_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2007_occ_2, thermal, by= "Species")
# 
# kelp_2007_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2007 <- mFD::sp.filter(asb_nm          = c("y2007"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2007 <- sp_filter_2007$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2007 <- sp_2d_coord_final_2007 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2007 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2007_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2007 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2007 <- sp_thermal_2007$Species[which(sp_thermal_2007$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2007 <- sp_thermal_2007$Species[which(sp_thermal_2007$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2007 <- Fric_2007$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2007 <- Fric_2007$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2007 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2007[sp_trop_2007, xy], 
#                                                   asb2 = sp_2d_coord_final_2007[sp_temp_2007, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2007, 
#                                                    asb2 = vert_temp_2007),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2007 <- ggplot_z_2007  + labs(title = "2007" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2007[[z]] <- ggplot_z_2007
#   
# }# end of z


################ 2008

# # Retrieve species coordinates matrix for year 2008:

kelp_2008_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2008") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2008_occ_2 <- kelp_2008_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2008_occ_2, thermal, by= "Species")

kelp_2008_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2008 <- mFD::sp.filter(asb_nm          = c("y2008"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2008 <- sp_filter_2008$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2008 <- sp_2d_coord_final_2008 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2008 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2008_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2008$functional_diversity_indices
fd_details <- Fric_2008$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2008$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2008$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi" ,
                 "Orectolobus_maculatus", "Scorpis_lineolata",
                 "Seriola_dumerili","Trachinops_taeniatus",
                 "Lutjanus_russellii", "Nelusetta_ayraud" ,
                 "Synodus_jaculum")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2008 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2008 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "red", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2008 <- ggplot_z_2008  + labs(title = "2008" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2008 <- ggplot_z_2008 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2008 <- ggplot_z_2008  + labs(title = "2008" )+
      theme(plot.title = element_text(size = 20, color = "red"))
  }
  # ggplot stored in list
  ggplot_2008[[z]] <- ggplot_z_2008
  
}# end of z




# kelp_2008_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2008") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2008_occ_2 <- kelp_2008_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2008_occ_2, thermal, by= "Species")
# 
# kelp_2008_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2008 <- mFD::sp.filter(asb_nm          = c("y2008"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2008 <- sp_filter_2008$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2008 <- sp_2d_coord_final_2008 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2008 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2008_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2008 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2008 <- sp_thermal_2008$Species[which(sp_thermal_2008$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2008 <- sp_thermal_2008$Species[which(sp_thermal_2008$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2008 <- Fric_2008$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2008 <- Fric_2008$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2008 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2008[sp_trop_2008, xy], 
#                                                   asb2 = sp_2d_coord_final_2008[sp_temp_2008, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2008, 
#                                                    asb2 = vert_temp_2008),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2008 <- ggplot_z_2008  + labs(title = "2008" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2008[[z]] <- ggplot_z_2008
#   
# }# end of z
# 


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

sp_2d_coord_final_2009 <- sp_filter_2009$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2009 <- sp_2d_coord_final_2009 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2009 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2009_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2009$functional_diversity_indices
fd_details <- Fric_2009$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2009$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2009$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Epinephelus_fasciatus", "Thalassoma_lunare",
                 "Thalassoma_lutescens", "Brachaelurus_waddi",
                 "Chaetodon_guentheri", "Kyphosus_sydneyanus",
                 "Orectolobus_maculatus", "Parma_unifasciata",
                 "Pseudocaranx_dentex", "Scorpis_lineolata",
                 "Trachinops_taeniatus")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2009 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2009 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "#2C6BAA", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2009 <- ggplot_z_2009  + labs(title = "2009" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2009 <- ggplot_z_2009 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2009 <- ggplot_z_2009  + labs(title = "2009" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  # ggplot stored in list
  ggplot_2009[[z]] <- ggplot_z_2009
  
}# end of z






# kelp_2009_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2009") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2009_occ_2 <- kelp_2009_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2009_occ_2, thermal, by= "Species")
# 
# kelp_2009_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2009 <- mFD::sp.filter(asb_nm          = c("y2009"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2009 <- sp_filter_2009$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2009 <- sp_2d_coord_final_2009 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2009 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                           asb_sp_w = kelp_2009_occ_thermal,
#                           ind_vect = c("fric"),
#                           scaling = TRUE,
#                           details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2009 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#     # species present in trop:
#   sp_trop_2009 <- sp_thermal_2009$Species[which(sp_thermal_2009$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2009 <- sp_thermal_2009$Species[which(sp_thermal_2009$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2009 <- Fric_2009$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2009 <- Fric_2009$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2009 <-fric.plot(ggplot_bg = ggplot_z, 
#                         asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2009[sp_trop_2009, xy], 
#                                               asb2 = sp_2d_coord_final_2009[sp_temp_2009, xy]),
#                         asb_vertices_nD = list(asb1 = vert_trop_2009, 
#                                                asb2 = vert_temp_2009),
#                         plot_sp = TRUE,
#                         color_sp = thermal_aff_colors,
#                         fill_sp = c(asb1 = "white", asb2 = "white"),
#                         size_sp = c(asb1 = 1, asb2 = 1),
#                         shape_sp = c(asb1 = 16, asb2 = 16),
#                         color_vert = thermal_aff_colors,
#                         fill_vert = thermal_aff_colors,
#                         size_vert = c(asb1 = 4, asb2 = 4),
#                         shape_vert = c(asb1 = 16, asb2 = 16),
#                         alpha_ch = c(asb1 = 0, asb2 = 0),
#                         color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                         fill_ch = c(asb1 = NA, asb2 = NA))
#   
#     # legend and title
#   if (z==1) {
#     ggplot_z_2009 <- ggplot_z_2009  + labs(title = "2009" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2009[[z]] <- ggplot_z_2009
#   
# }# end of z
# 


################ 2010

# # Retrieve species coordinates matrix for year 2010:

kelp_2010_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2010") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2010_occ_2 <- kelp_2010_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2010_occ_2, thermal, by= "Species")

kelp_2010_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2010 <- mFD::sp.filter(asb_nm          = c("y2010"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2010 <- sp_filter_2010$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2010 <- sp_2d_coord_final_2010 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2010 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2010_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2010$functional_diversity_indices
fd_details <- Fric_2010$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2010$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2010$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi",
                 "Orectolobus_maculatus","Parapercis_ramsayi",
                 "Scorpis_lineolata","Seriola_lalandi",
                 "Trachinops_taeniatus", "Aetobatus_narinari" ,
                 "Labroides_dimidiatus", "Paracaesio_xanthura",
                 "Platax_teira")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2010 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2010 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "#2C6BAA", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2010 <- ggplot_z_2010  + labs(title = "2010" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2010 <- ggplot_z_2010 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2010 <- ggplot_z_2010  + labs(title = "2010" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  # ggplot stored in list
  ggplot_2010[[z]] <- ggplot_z_2010
  
}# end of z









# kelp_2010_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2010") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2010_occ_2 <- kelp_2010_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2010_occ_2, thermal, by= "Species")
# 
# kelp_2010_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2010 <- mFD::sp.filter(asb_nm          = c("y2010"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2010 <- sp_filter_2010$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2010 <- sp_2d_coord_final_2010 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2010 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2010_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2010 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2010 <- sp_thermal_2010$Species[which(sp_thermal_2010$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2010 <- sp_thermal_2010$Species[which(sp_thermal_2010$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2010 <- Fric_2010$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2010 <- Fric_2010$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2010 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2010[sp_trop_2010, xy], 
#                                                   asb2 = sp_2d_coord_final_2010[sp_temp_2010, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2010, 
#                                                    asb2 = vert_temp_2010),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2010 <- ggplot_z_2010  + labs(title = "2010" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2010[[z]] <- ggplot_z_2010
#   
# }# end of z

################ 2011

# # Retrieve species coordinates matrix for year 2011:

kelp_2011_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2011") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2011_occ_2 <- kelp_2011_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2011_occ_2, thermal, by= "Species")

kelp_2011_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2011 <- mFD::sp.filter(asb_nm          = c("y2011"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2011 <- sp_filter_2011$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2011 <- sp_2d_coord_final_2011 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2011 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2011_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2011$functional_diversity_indices
fd_details <- Fric_2011$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2011$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2011$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi",
                 "Kyphosus_sydneyanus", "Orectolobus_halei",
                 "Pseudocaranx_dentex","Scorpis_lineolata",
                 "Trachinops_taeniatus", "Aulostomus_chinensis" ,
                 "Thalassoma_lunare", "Thalassoma_lutescens")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2011 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2011 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "#2C6BAA", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2011 <- ggplot_z_2011  + labs(title = "2011" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2011 <- ggplot_z_2011 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2011 <- ggplot_z_2011  + labs(title = "2011" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  # ggplot stored in list
  ggplot_2011[[z]] <- ggplot_z_2011
  
}# end of z





# kelp_2011_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2011") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2011_occ_2 <- kelp_2011_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2011_occ_2, thermal, by= "Species")
# 
# kelp_2011_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2011 <- mFD::sp.filter(asb_nm          = c("y2011"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2011 <- sp_filter_2011$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2011 <- sp_2d_coord_final_2011 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2011 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2011_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2011 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2011 <- sp_thermal_2011$Species[which(sp_thermal_2011$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2011 <- sp_thermal_2011$Species[which(sp_thermal_2011$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2011 <- Fric_2011$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2011 <- Fric_2011$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2011 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2011[sp_trop_2011, xy], 
#                                                   asb2 = sp_2d_coord_final_2011[sp_temp_2011, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2011, 
#                                                    asb2 = vert_temp_2011),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2011 <- ggplot_z_2011  + labs(title = "2011" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2011[[z]] <- ggplot_z_2011
#   
# }# end of z


############# 2013

# # Retrieve species coordinates matrix for year 2013:

kelp_2013_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2013") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2013_occ_2 <- kelp_2013_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2013_occ_2, thermal, by= "Species")

kelp_2013_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2013 <- mFD::sp.filter(asb_nm          = c("y2013"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2013 <- sp_filter_2013$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2013 <- sp_2d_coord_final_2013 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2013 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2013_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2013$functional_diversity_indices
fd_details <- Fric_2013$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2013$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2013$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi",
                 "Orectolobus_halei","Parapercis_ramsayi",
                 "Scorpis_lineolata","Seriola_hippos",
                 "Trachinops_taeniatus", 
                "Labroides_dimidiatus",
                 "Lutjanus_russellii" , "Siganus_fuscescens")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2013 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2013 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "#2C6BAA", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2013 <- ggplot_z_2013  + labs(title = "2013" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2013 <- ggplot_z_2013 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2013 <- ggplot_z_2013  + labs(title = "2013" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  # ggplot stored in list
  ggplot_2013[[z]] <- ggplot_z_2013
  
}# end of z








# kelp_2013_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2013") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2013_occ_2 <- kelp_2013_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2013_occ_2, thermal, by= "Species")
# 
# kelp_2013_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2013 <- mFD::sp.filter(asb_nm          = c("y2013"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2013 <- sp_filter_2013$`species coordinates`[, c("PC1", "PC2", "PC3")]
# 
# sp_thermal_2013 <- sp_2d_coord_final_2013 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2013 <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
#                                asb_sp_w = kelp_2013_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2013 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2013 <- sp_thermal_2013$Species[which(sp_thermal_2013$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2013 <- sp_thermal_2013$Species[which(sp_thermal_2013$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2013 <- Fric_2013$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2013 <- Fric_2013$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2013 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2013[sp_trop_2013, xy], 
#                                                   asb2 = sp_2d_coord_final_2013[sp_temp_2013, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2013, 
#                                                    asb2 = vert_temp_2013),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#00C19A"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2013 <- ggplot_z_2013  + labs(title = "2013" )+
#       theme(plot.title = element_text(size = 20, color = "#00C19A"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2013[[z]] <- ggplot_z_2013
#   
# }# end of z


################ 2015

# # Retrieve species coordinates matrix for year 2015:

kelp_2015_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2015") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2015_occ_2 <- kelp_2015_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2015_occ_2, thermal, by= "Species")

kelp_2015_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2015 <- mFD::sp.filter(asb_nm          = c("y2015"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2015 <- sp_filter_2015$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2015 <- sp_2d_coord_final_2015 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2015 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2015_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2015$functional_diversity_indices
fd_details <- Fric_2015$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2015$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2015$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi",
                 "Orectolobus_halei","Parapercis_stricticeps",
                 "Scorpis_lineolata","Seriola_dumerili",
                 "Trachinops_taeniatus", 
                 "Labroides_dimidiatus","Lethrinus_nebulosus",
                 "Lutjanus_russellii" , "Plectorhinchus_flavomaculatus",
                 "Pterocaesio_digramma")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2015 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2015 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "#2C6BAA", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2015 <- ggplot_z_2015  + labs(title = "2015" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2015 <- ggplot_z_2015 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2015 <- ggplot_z_2015  + labs(title = "2015" )+
      theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
  }
  # ggplot stored in list
  ggplot_2015[[z]] <- ggplot_z_2015
  
}# end of z








# kelp_2015_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2015") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2015_occ_2 <- kelp_2015_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2015_occ_2, thermal, by= "Species")
# 
# kelp_2015_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2015 <- mFD::sp.filter(asb_nm          = c("y2015"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2015 <- sp_filter_2015$`species coordinates`[, c("PC1", "PC2")]
# 
# sp_thermal_2015 <- sp_2d_coord_final_2015 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# 
# ## Compute FRic values #### 
# 
# # compute FRic for all habitats  ---
# Fric_2015 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2015_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_2015 <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_2015 <- sp_thermal_2015$Species[which(sp_thermal_2015$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_2015 <- sp_thermal_2015$Species[which(sp_thermal_2015$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_2015 <- Fric_2015$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_2015 <- Fric_2015$details$asb_vert_nm$temperate
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2015 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2015[sp_trop_2015, xy], 
#                                                   asb2 = sp_2d_coord_final_2015[sp_temp_2015, xy]),
#                             asb_vertices_nD = list(asb1 = vert_trop_2015, 
#                                                    asb2 = vert_temp_2015),
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2015 <- ggplot_z_2015  + labs(title = "2015" )+
#       theme(plot.title = element_text(size = 20, color = "#2C6BAA"))
#   }
#   
#   # ggplot stored in list
#   ggplot_2015[[z]] <- ggplot_z_2015
#   
# }# end of z

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

Fric_2018 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2018_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2018$functional_diversity_indices
fd_details <- Fric_2018$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2018$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2018$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Brachaelurus_waddi",
                 "Orectolobus_maculatus","Parapercis_stricticeps",
                 "Pseudocaranx_dentex", "Scorpis_lineolata",
                 "Trachinops_taeniatus", 
                 "Lutjanus_russellii" , "Paracaesio_xanthura",
                 "Stethojulis_interrupta",
                 "Thalassoma_lunare")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_2018 <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_2018 <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "#00C19A", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2018 <- ggplot_z_2018  + labs(title = "2018" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2018 <- ggplot_z_2018 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_2018 <- ggplot_z_2018  + labs(title = "2018" )+
      theme(plot.title = element_text(size = 20, color = "#00C19A"))
  }
  # ggplot stored in list
  ggplot_2018[[z]] <- ggplot_z_2018
  
}# end of z







# kelp_2018_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
#   as.data.frame() %>%
#   rownames_to_column("Sites") %>%
#   filter(Sites == "y2018") %>%
#   column_to_rownames("Sites") %>%
#   as.matrix()
# 
# kelp_2018_occ_2 <- kelp_2018_occ %>% 
#   as.data.frame() %>% 
#   gather(Species, Abundance, 1:101)
# 
# add_thermal <- left_join(kelp_2018_occ_2, thermal, by= "Species")
# 
# kelp_2018_occ_thermal <- add_thermal %>% 
#   spread(Species, Abundance, fill=0) %>% 
#   column_to_rownames("thermal_label") %>% 
#   as.matrix()
# 
# sp_filter_2018 <- mFD::sp.filter(asb_nm          = c("y2018"), 
#                                  sp_faxes_coord = sp_3D_coord, 
#                                  asb_sp_w       = kelp_years_sp_occ)
# 
# sp_2d_coord_final_2018 <- sp_filter_2018$`species coordinates`[, c("PC1", "PC2", "PC3")]
# 
# sp_thermal_2018 <- sp_2d_coord_final_2018 %>% 
#   as.data.frame() %>% 
#   rownames_to_column("Species") %>% 
#   left_join(add_thermal, by="Species")
# 
# ## Compute FRic values #### 
# 
# Fric_2018 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
#                                asb_sp_w = kelp_2018_occ_thermal,
#                                ind_vect = c("fric"),
#                                scaling = TRUE,
#                                details_returned = TRUE)
# 
# 
# ## retrieve names of main input:
# asb_fd_ind <- Fric_2018$functional_diversity_indices
# fd_details <- Fric_2018$details
# 
# ### Prepare data for plotting:
# 
# ## get coordinates of species:
# sp_faxes_coord <- fd_details$sp_faxes_coord
# 
# ## get number of dimensions in input:
# nb_dim <- ncol(sp_faxes_coord)
# 
# ## Define arguments
# 
# faxes               = NULL
# faxes_nm            = NULL
# range_faxes         = c(NA, NA)
# plot_asb_nm <- c("temperate", "tropical")
# plot_sp_nm <- c("Austrolabrus_maculatus","Brachaelurus_waddi", "Orectolobus_maculatus",
#                 "Parapercis_stricticeps","Pseudocaranx_dentex","Scorpis_lineolata", "Trachinops_taeniatus",
#                 "Choerodon_venustus", "Lutjanus_adetii", "Lutjanus_fulviflamma","Lutjanus_russellii",
#                 "Paracaesio_xanthura", "Parupeneus_multifasciatus", "Siganus_fuscescens", "Stethojulis_interrupta",
#                 "Thalassoma_lunare", "Thalassoma_lutescens")
# 
# ## define arguments values and prepare data for plotting:
# 
# # give faxes identity if faxes set to NULL:
# if (is.null(faxes)) {
#   faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
# }
# 
# # give faxes names if faxes set to NULL:
# if (is.null(faxes_nm)) {
#   faxes_nm <- faxes
# }
# names(faxes_nm) <- faxes
# 
# # get number of axes:
# nb_faxes <- length(faxes)
# 
# # get combinations of axes on plot:
# axes_plot <- utils::combn(faxes, 2)
# plot_nb   <- ncol(axes_plot)
# 
# # set range of axes if c(NA, NA):
# if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
#   range_sp_coord  <- range(sp_faxes_coord)
#   range_faxes <- range_sp_coord +
#     c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
# }
# 
# # create a dataframe with species coordinates and option (vertices + label)
#  sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")
# 
# # if some species names to be plotted, adding a character variable to sp_faxes_coord:
# 
# if (! is.null(plot_sp_nm)) {
#   sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
# }
#  
#  # get vertices of the convex hull of the species pool:
#  vert_pool <- fd_details$pool_vert_nm
#  
#  # retrieve names and weights of species present in each assemblage:
#  
#  # get names of assemblages:
#  pool <- "pool"
#  asb1 <- plot_asb_nm[1]
#  nm_asb <- asb1
#  asb2 <- plot_asb_nm[2]
#  nm_asb <- paste(nm_asb, asb2, sep = "_")
# 
#  sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
#  sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))
# 
#  ## plotting  ####
# 
# # list to store ggplot
# ggplot_2018 <- list()
# 
# # pairs of axes
# 
# for (z in (1:plot_nb)) {
#   
#   # names of axes
#   xy_z <- axes_plot[1:2, z]
#   
#   # get species coordinates along the 2 axes:
#   sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
#   colnames(sp_coord_xy) <- c("x", "y")
#   
#   # list with dataframes for plot:
#   asb_sp_coord2D_k <- list()
#   asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
#   vertices_nD_k <- list()
#   vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
#   asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
#   vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_2018 <-fric.plot(ggplot_bg = ggplot_z, 
#                             asb_sp_coord2D = asb_sp_coord2D_k,
#                             asb_vertices_nD = vertices_nD_k,
#                             plot_sp = TRUE,
#                             color_sp = thermal_aff_colors,
#                             fill_sp = c(asb1 = "white", asb2 = "white"),
#                             size_sp = c(asb1 = 1, asb2 = 1),
#                             shape_sp = c(asb1 = 16, asb2 = 16),
#                             color_vert = thermal_aff_colors,
#                             fill_vert = thermal_aff_colors,
#                             size_vert = c(asb1 = 4, asb2 = 4),
#                             shape_vert = c(asb1 = 16, asb2 = 16),
#                             alpha_ch = c(asb1 = 0, asb2 = 0),
#                             color_ch = c(asb1 = "#00C19A", asb2 =NA),
#                             fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_2018 <- ggplot_z_2018  + labs(title = "2018" )+
#       theme(plot.title = element_text(size = 20, color = "#00C19A"))
#   }
#   
#   # add species names if needed:
#   if (! is.null(plot_sp_nm)) {
#     x <- NULL
#     y <- NULL
#     ggplot_z_2018 <- ggplot_z_2018 +
#       ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
#                                ggplot2::aes_string(x = xy_z[1],
#                                                    y = xy_z[2],
#                                                    label = "label"),
#                                size = 3, colour= "black",
#                                fontface = "plain",
#                                max.overlaps = Inf,
#                                box.padding = grid::unit(2, 'lines'),
#                                force = 5,
#                                arrow = grid::arrow(length = grid::unit(0.02,
#                                                                        'npc')),
#                                segment.color = "black")
#   }
#   
#                             
#                             
# # legend and title
#   if (z==1) {
#     ggplot_z_2018 <- ggplot_z_2018  + labs(title = "2018" )+
#       theme(plot.title = element_text(size = 20, color = "#00C19A"))
#   }
#   # ggplot stored in list
#   ggplot_2018[[z]] <- ggplot_z_2018
#   
# }# end of z
#   


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

## retrieve names of main input:
asb_fd_ind <- Fric_inshore$functional_diversity_indices
fd_details <- Fric_inshore$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices 
vert <- Fric_inshore$details$asb_vert_nm$Inshore
vert

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("Inshore")
plot_sp_nm <- c( "Austrolabrus_maculatus", "Achoerodus_viridis" ,
                 "Girella_elevata","Pempheris_affinis",
                 "Parma_oligolepis", "Scorpis_lineolata")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_inshore <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_inshore <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = c(asb1 = "#2C6BAA"),
                            fill_sp = c(asb1 = "white"),
                            size_sp = c(asb1 = 1),
                            shape_sp = c(asb1 = 16),
                            color_vert = c(asb1 = "#2C6BAA"),
                            fill_vert = c(asb1 = "#2C6BAA"),
                            size_vert = c(asb1 = 4),
                            shape_vert = c(asb1 = 16),
                            alpha_ch = c(asb1 = 0),
                            color_ch = c(asb1 = "#2C6BAA"),
                            fill_ch = c(asb1 = NA))
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_inshore <- ggplot_z_inshore +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
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
Fric_midshelf <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                                   asb_sp_w = midshelf_occ_thermal,
                                   ind_vect = c("fric"),
                                   scaling = TRUE,
                                   details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_midshelf$functional_diversity_indices
fd_details <- Fric_midshelf$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


## Vertices to see which species to plot:

# vertices in trop:
vert_trop_midshelf <- Fric_midshelf$details$asb_vert_nm$tropical
vert_trop_midshelf

# vertices in temp:
vert_temp_midshelf <- Fric_midshelf$details$asb_vert_nm$temperate
vert_temp_midshelf

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c("Aulostomus_chinensis",  "Halichoeres_margaritaceus", 
                "Siganus_fuscescens", "Stegastes_apicalis",
                "Thalassoma_lunare","Thalassoma_lutescens",
                "Achoerodus_viridis", "Kyphosus_sydneyanus",
                "Parma_unifasciata", "Pempheris_affinis",
                 "Scorpis_lineolata")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)


# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))


## plotting  ####

# list to store ggplot
ggplot_midshelf <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_midshelf <-fric.plot(ggplot_bg = ggplot_z, 
                                asb_sp_coord2D = asb_sp_coord2D_k,
                                asb_vertices_nD = vertices_nD_k,
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
                                color_ch = c(asb1 = "lightsalmon1", asb2 =NA),
                                fill_ch = c(asb1 = NA, asb2 = NA))
  
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_midshelf <- ggplot_z_midshelf +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
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
Fric_offshore <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                                   asb_sp_w = offshore_occ_thermal,
                                   ind_vect = c("fric"),
                                   scaling = TRUE,
                                   details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_offshore$functional_diversity_indices
fd_details <- Fric_offshore$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


## Vertices to see which species to plot:

# vertices in trop:
vert_trop_offshore <- Fric_offshore$details$asb_vert_nm$tropical
vert_trop_offshore

# vertices in temp:
vert_temp_offshore <- Fric_offshore$details$asb_vert_nm$temperate
vert_temp_offshore

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c("Aulostomus_chinensis", "Caesio_caerulaurea", "Halichoeres_margaritaceus",  "Leptojulis_cyanopleura",
                 "Stegastes_gascoynei", "Thalassoma_lunare","Halichoeres_margaritaceus",
                "Leptojulis_cyanopleura", "Labroides_dimidiatus", "Mulloidichthys_vanicolensis", 
                "Naso_unicornis", "Plectorhinchus_flavomaculatus", "Paracanthurus_hepatus", 
                "Achoerodus_viridis", "Kyphosus_sydneyanus", "Pempheris_affinis",
                "Pseudocaranx_dentex","Scorpis_lineolata")

## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to sp_faxes_coord:

if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))


## plotting  ####

# list to store ggplot
ggplot_offshore <- list()

# pairs of axes

for (z in (1:plot_nb)) {
  
  # names of axes
  xy_z <- axes_plot[1:2, z]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_z])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  
  # plot convex hull of assemblage but not species
  ggplot_z_offshore <-fric.plot(ggplot_bg = ggplot_z, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
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
                            color_ch = c(asb1 = "firebrick3", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_offshore <- ggplot_z_offshore +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 3, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  
  
  # legend and title
  if (z==1) {
    ggplot_z_offshore <- ggplot_z_offshore  + labs(title = "Offshore" )+
      theme(plot.title = element_text(size = 20, color = "firebrick3"))
  }
  # ggplot stored in list
  ggplot_offshore[[z]] <- ggplot_z_offshore
  
}# end of z


# 
# ## plotting  ####
# 
# # list to store ggplot
# ggplot_offshore <- list()
# 
# # pairs of axes
# pairs_axes <- list( c(1,2) )
# 
# for (z in 1:length(pairs_axes)) {
#   
#   # names of axes   
#   xy <- pairs_axes[[z]]
#   
#   # species present in trop:
#   sp_trop_offshore <- sp_thermal_offshore$Species[which(sp_thermal_offshore$thermal_label == "tropical")]
#   
#   # species present in temp:
#   sp_temp_offshore <- sp_thermal_offshore$Species[which(sp_thermal_offshore$thermal_label == "temperate")]
#   
#   # vertices in trop:
#   vert_trop_offshore <- Fric_offshore$details$asb_vert_nm$tropical
#   
#   # vertices in temp:
#   vert_temp_offshore <- Fric_offshore$details$asb_vert_nm$temperate
#   
#   #Vertices
#   
#   vert_offshore <- Fric_offshore$details$asb_vert_nm$Offshore
#   
#   # plot convex hull of assemblage but not species
#   ggplot_z_offshore <-fric.plot(ggplot_bg = ggplot_z, 
#                                 asb_sp_coord2D = list(asb1 = sp_2d_coord_final_offshore[sp_trop_offshore, xy], 
#                                                       asb2 = sp_2d_coord_final_offshore[sp_temp_offshore, xy]),
#                                 asb_vertices_nD = list(asb1 = vert_trop_offshore, 
#                                                        asb2 = vert_temp_offshore),
#                                 plot_sp = TRUE,
#                                 color_sp = thermal_aff_colors,
#                                 fill_sp = c(asb1 = "white", asb2 = "white"),
#                                 size_sp = c(asb1 = 1, asb2 = 1),
#                                 shape_sp = c(asb1 = 16, asb2 = 16),
#                                 color_vert = thermal_aff_colors,
#                                 fill_vert = thermal_aff_colors,
#                                 size_vert = c(asb1 = 4, asb2 = 4),
#                                 shape_vert = c(asb1 = 16, asb2 = 16),
#                                 alpha_ch = c(asb1 = 0, asb2 = 0),
#                                 color_ch = c(asb1 = NA, asb2 ="firebrick3"),
#                                 fill_ch = c(asb1 = NA, asb2 = NA))
#   
#   # legend and title
#   if (z==1) {
#     ggplot_z_offshore <- ggplot_z_offshore  + labs(title = "Offshore" )+
#       theme(plot.title = element_text(size = 20, color = "firebrick3"))
#   }
#   
#   # ggplot stored in list
#   ggplot_offshore[[z]] <- ggplot_z_offshore
#   
# }# end of z

## merging all plots into a single figure and saving as png ####

figure9A <- (ggplot_2002[[1]] + ggplot_2009[[1]] + ggplot_2018[[1]])
               
figure9B <-(ggplot_inshore[[1]] + ggplot_midshelf[[1]] + ggplot_offshore[[1]])

###SAVE EACH FIG INDIVIDUALLY


ggsave(figure9, file=here::here("outputs/", "using biomass-maxN", "Figure9_biomass.png"),
       height = 16, width = 24, unit = "cm" )


figure4S <- ( ggplot_2002[[1]] + ggplot_2003[[1]] +  ggplot_2004[[1]] +  ggplot_2005[[1]] + ggplot_2006[[1]] +  ggplot_2007[[1]] ) /
  ( ggplot_2008[[1]] + ggplot_2009[[1]] +  ggplot_2011[[1]] + ggplot_2013[[1]] + ggplot_2015[[1]] +  ggplot_2018[[1]] )

ggsave(figure4S , file=here::here("outputs/", "using biomass-maxN", "Figure4s_biomass.png"),
       height = 30, width = 30, unit = "cm" )

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
