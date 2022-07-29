################################################################################
##
## Script for plotting changes in the multidimensional functional space across 
## habitats and over years in a multidimensional space
##
##  * changes in space and time (only 2002, 2009 and 2018) - Figure 4
##
##  * changes in time (all years) - Figure S7
## 
## 
## Code by Paula Sgarlatta, Camille Magneville and Sebastien Villeger 
##
################################################################################


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(dplyr)
library(here)
library(patchwork)
library(mFD)
library(ggpubr)
library(fishualize)

# loading data

load(here::here("data", "sp_3D_coord.RData") )

load(here::here("data", "spatial_metadata.RData") )
load(here::here("data", "spatial_sp_biom.RData") )
load(here::here("data", "spatial_fd_biomass.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "kelp_sp_maxN.RData") )
load(here::here("data", "temporal_fd_kelp_biomass.RData") )

# loading thermal affinity data

thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  select(-thermal)

sp_3D_coord_bg <- sp_3D_coord #This one is just for the background polygon

# First, shorten species names 

sp_3D_coord <- sp_3D_coord %>% 
  as.data.frame() %>% 
  rownames_to_column("Species_2") %>% 
  mutate(genus=sub("_.*", "", Species_2), sp=sub(".*_", "", Species_2)) %>% 
  mutate(Species=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  select(-genus, -sp, -Species_2) %>% 
  column_to_rownames("Species")

sp_3D_coord <- as.matrix(sp_3D_coord)

spatial_sp_biom <- spatial_sp_biom %>% 
  as.data.frame() %>% 
  rownames_to_column("Site") %>% 
  pivot_longer(!Site, names_to="species", values_to="count") %>% 
  mutate(genus=sub("_.*", "", species), sp=sub(".*_", "", species)) %>% 
  mutate(Species=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  dplyr::select(-genus, -sp, -species) %>% 
  pivot_wider(names_from = Species, values_from=count) %>% 
  column_to_rownames("Site")

kelp_sp_maxN_1 <- kelp_sp_maxN %>% 
  as.data.frame() %>% 
  rownames_to_column("Site") %>% 
  pivot_longer(!Site, names_to="species", values_to="count") %>% 
  mutate(genus=sub("_.*", "", species), sp=sub(".*_", "", species)) %>% 
  mutate(Species=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  dplyr::select(-genus, -sp, -species) %>% 
  pivot_wider(names_from = Species, values_from=count) %>% 
  column_to_rownames("Site")

thermal <- thermal %>% 
  mutate(genus=sub("_.*", "", Species), sp=sub(".*_", "", Species)) %>% 
  mutate(Species=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  select(-genus, -sp) 

## settings ####

# vertices of all fe in 4D ----

pool_vert_nm <- spatial_fd$details$pool_vert_nm 

# range of axes
range_faxes_coord <- range(sp_3D_coord_bg)
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]

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
  y2014 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2014"),"Code"],],2,max ),
  y2015 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2015"),"Code"],],2,max ),
  y2018 = apply(kelp_sp_maxN [kelp_metadata[which(kelp_metadata$Year=="2018"),"Code"],],2,max ))  

# background with axes range set + title - different options for graphs

ggplot_z_full <- background.plot(range_faxes = range_axes,
                            faxes_nm = c("PC1", "PC2"), 
                            color_bg = "grey95")

ggplot_z_PC1 <-background.plot(range_faxes = range_axes,
                               faxes_nm = c("PC1", ""), 
                               color_bg = "grey95")


ggplot_z_PC2 <-background.plot(range_faxes = range_axes,
                               faxes_nm = c("", "PC2"), 
                               color_bg = "grey95")

ggplot_z_empty <-background.plot(range_faxes = range_axes,
                               faxes_nm = c("", ""), 
                               color_bg = "grey95")


# convex hull of global species pool

ggplot_z_full <- pool.plot(ggplot_bg = ggplot_z_full,
                      sp_coord2D = sp_3D_coord_bg,
                      vertices_nD = pool_vert_nm,
                      plot_pool = FALSE,
                      color_ch = "NA", fill_ch = "white", alpha_ch = 1)

ggplot_z_PC1 <- pool.plot(ggplot_bg = ggplot_z_PC1,
                           sp_coord2D = sp_3D_coord_bg,
                           vertices_nD = pool_vert_nm,
                           plot_pool = FALSE,
                           color_ch = "NA", fill_ch = "white", alpha_ch = 1)


ggplot_z_PC2 <- pool.plot(ggplot_bg = ggplot_z_PC2,
                          sp_coord2D = sp_3D_coord_bg,
                          vertices_nD = pool_vert_nm,
                          plot_pool = FALSE,
                          color_ch = "NA", fill_ch = "white", alpha_ch = 1)

ggplot_z_empty <- pool.plot(ggplot_bg = ggplot_z_empty,
                          sp_coord2D = sp_3D_coord_bg,
                          vertices_nD = pool_vert_nm,
                          plot_pool = FALSE,
                          color_ch = "NA", fill_ch = "white", alpha_ch = 1)




sp_2d_coord <- sp_3D_coord[, c("PC1", "PC2")]

##################### 2002 #################################

# # Retrieve species coordinates matrix for year 2002:

kelp_2002_occ <- kelp_years_sp_occ %>% 
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2002") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2002_occ_2 <- kelp_2002_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- dplyr::full_join( kelp_2002_occ_2, thermal,by= "Species") #not working!!!

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
plot_sp_nm <- c( "T. taeniatus", "E. fasciatus","G. thyrsoidea","T. lutescens", 
                 "B. waddi" ,"O. halei", "A. maculatus", "S. lalandi" )

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
nm_asb <- paste(nm_asb, asb2, sep = ". ")

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
  # plot convex hull of assemblage 
  ggplot_z_2002 <-fric.plot(ggplot_bg = ggplot_z_PC2, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
    # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    ggplot_z_2002 <- ggplot_z_2002 +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = 4, colour= "black",
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
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  # ggplot stored in list
  ggplot_2002[[z]] <- ggplot_z_2002
  
}# end of z

################### Now without species names 

plot_sp_nm <- NULL

# plot convex hull of assemblage 
ggplot_2002_no_sp <-fric.plot(ggplot_bg = ggplot_z_PC2, 
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_vertices_nD = vertices_nD_k,
                          plot_sp = TRUE,
                          color_sp = thermal_aff_colors,
                          fill_sp = c(asb1 = "white", asb2 = "white"),
                          size_sp = c(asb1 = 1, asb2 = 1),
                          shape_sp = c(asb1 = 16, asb2 = 16),
                          color_vert = thermal_aff_colors,
                          fill_vert = thermal_aff_colors,
                          size_vert = c(asb1 = 1, asb2 = 1),
                          shape_vert = c(asb1 = 16, asb2 = 16),
                          alpha_ch = c(asb1 = 0, asb2 = 0),
                          color_ch = c(asb1 = "seagreen4", asb2 =NA),
                          fill_ch = c(asb1 = NA, asb2 = NA))

ggplot_2002_no_sp <- ggplot_2002_no_sp  + labs(title = "2002" )+
  theme(plot.title = element_text(size = 20, color = "seagreen4"))

################ 2003

# # Retrieve species coordinates matrix for year 2003:

kelp_2003_occ <- kelp_years_sp_occ %>% 
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

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = ". ")

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
  ggplot_z_2003 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2003 <- ggplot_z_2003  + labs(title = "2003" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  
     # ggplot stored in list
  ggplot_2003[[z]] <- ggplot_z_2003
  
}# end of z

################ 2004

# # Retrieve species coordinates matrix for year 2004:

kelp_2004_occ <- kelp_years_sp_occ %>% 
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
plot_sp_nm <- NULL

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

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = ". ")

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
  ggplot_z_2004 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2004 <- ggplot_z_2004  + labs(title = "2004" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  
  # ggplot stored in list
  ggplot_2004[[z]] <- ggplot_z_2004
  
}# end of z

################ 2005

# # Retrieve species coordinates matrix for year 2005:

kelp_2005_occ <- kelp_years_sp_occ %>% 
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
plot_sp_nm <- NULL

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

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = ". ")

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
  ggplot_z_2005 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2005 <- ggplot_z_2005  + labs(title = "2005" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  
  # ggplot stored in list
  ggplot_2005[[z]] <- ggplot_z_2005
  
}# end of z

################ 2006

# # Retrieve species coordinates matrix for year 2006:

kelp_2006_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2006 <-fric.plot(ggplot_bg = ggplot_z_PC2, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2006 <- ggplot_z_2006  + labs(title = "2006" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  
  # ggplot stored in list
  ggplot_2006[[z]] <- ggplot_z_2006
  
}# end of z

################ 2007

# # Retrieve species coordinates matrix for year 2007:

kelp_2007_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2007 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2007 <- ggplot_z_2007  + labs(title = "2007" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  
  # ggplot stored in list
  ggplot_2007[[z]] <- ggplot_z_2007
  
}# end of z

################ 2008

# # Retrieve species coordinates matrix for year 2008:

kelp_2008_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2008 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2008 <- ggplot_z_2008  + labs(title = "2008" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  
  # ggplot stored in list
  ggplot_2008[[z]] <- ggplot_z_2008
  
}# end of z

################ 2009

# # Retrieve species coordinates matrix for year 2009:

kelp_2009_occ <- kelp_years_sp_occ %>% 
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
plot_sp_nm <- c( "E. fasciatus", "T. lunare",
                 "T. lutescens", "B.waddi",
                 "K. sydneyanus",
                 "P. unifasciata",
                 "P. dentex", 
                 "T. taeniatus")

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
  ggplot_z_2009 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2009 <- ggplot_z_2009  + labs(title = "2009" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
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
                               size = 4, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
  # ggplot stored in list
  ggplot_2009[[z]] <- ggplot_z_2009
  
}# end of z

######### Now without species names 

plot_sp_nm <- NULL

# plot convex hull of assemblage but not species
ggplot_2009_no_sp <-fric.plot(ggplot_bg = ggplot_z_empty, 
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_vertices_nD = vertices_nD_k,
                          plot_sp = TRUE,
                          color_sp = thermal_aff_colors,
                          fill_sp = c(asb1 = "white", asb2 = "white"),
                          size_sp = c(asb1 = 1, asb2 = 1),
                          shape_sp = c(asb1 = 16, asb2 = 16),
                          color_vert = thermal_aff_colors,
                          fill_vert = thermal_aff_colors,
                          size_vert = c(asb1 = 1, asb2 = 1),
                          shape_vert = c(asb1 = 16, asb2 = 16),
                          alpha_ch = c(asb1 = 0, asb2 = 0),
                          color_ch = c(asb1 = "seagreen4", asb2 =NA),
                          fill_ch = c(asb1 = NA, asb2 = NA))

# legend and title

ggplot_2009_no_sp <- ggplot_2009_no_sp  + labs(title = "2009" )+
    theme(plot.title = element_text(size = 20, color = "seagreen4"))

################ 2010

# # Retrieve species coordinates matrix for year 2010:

kelp_2010_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2010 <-fric.plot(ggplot_bg = ggplot_z_full, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2010 <- ggplot_z_2010  + labs(title = "2010" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  # ggplot stored in list
  ggplot_2010[[z]] <- ggplot_z_2010
  
}# end of z

################ 2011

# # Retrieve species coordinates matrix for year 2011:

kelp_2011_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2011 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2011 <- ggplot_z_2011  + labs(title = "2011" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  # ggplot stored in list
  ggplot_2011[[z]] <- ggplot_z_2011
  
}# end of z

############# 2013

# # Retrieve species coordinates matrix for year 2013:

kelp_2013_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2013 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2013 <- ggplot_z_2013  + labs(title = "2013" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  # ggplot stored in list
  ggplot_2013[[z]] <- ggplot_z_2013
  
}# end of z


############# 2014

# # Retrieve species coordinates matrix for year 2014:

kelp_2014_occ <- kelp_years_sp_occ %>% 
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2014") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

kelp_2014_occ_2 <- kelp_2014_occ %>% 
  as.data.frame() %>% 
  gather(Species, Abundance, 1:101)

add_thermal <- left_join(kelp_2014_occ_2, thermal, by= "Species")

kelp_2014_occ_thermal <- add_thermal %>% 
  spread(Species, Abundance, fill=0) %>% 
  column_to_rownames("thermal_label") %>% 
  as.matrix()

sp_filter_2014 <- mFD::sp.filter(asb_nm          = c("y2014"), 
                                 sp_faxes_coord = sp_3D_coord, 
                                 asb_sp_w       = kelp_years_sp_occ)

sp_2d_coord_final_2014 <- sp_filter_2014$`species coordinates`[, c("PC1", "PC2", "PC3")]

sp_thermal_2014 <- sp_2d_coord_final_2014 %>% 
  as.data.frame() %>% 
  rownames_to_column("Species") %>% 
  left_join(add_thermal, by="Species")

## Compute FRic values #### 

Fric_2014 <- alpha.fd.multidim(sp_faxes_coord = sp_2d_coord,
                               asb_sp_w = kelp_2014_occ_thermal,
                               ind_vect = c("fric"),
                               scaling = TRUE,
                               details_returned = TRUE)


## retrieve names of main input:
asb_fd_ind <- Fric_2014$functional_diversity_indices
fd_details <- Fric_2014$details

### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric_2014$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric_2014$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")

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
ggplot_2014 <- list()

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
  ggplot_z_2014 <-fric.plot(ggplot_bg = ggplot_z_PC1, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2014 <- ggplot_z_2014  + labs(title = "2014" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  # ggplot stored in list
  ggplot_2014[[z]] <- ggplot_z_2014
  
}# end of z

################ 2015

# # Retrieve species coordinates matrix for year 2015:

kelp_2015_occ <- kelp_years_sp_occ %>% 
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
  ggplot_z_2015 <-fric.plot(ggplot_bg = ggplot_z_full, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2015 <- ggplot_z_2015  + labs(title = "2015" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
  }
  # ggplot stored in list
  ggplot_2015[[z]] <- ggplot_z_2015
  
}# end of z

################ 2018

# # Retrieve species coordinates matrix for year 2018:

kelp_2018_occ <- kelp_years_sp_occ %>% 
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
plot_sp_nm <- c( "A. maculatus", "B. waddi",
                 "O. maculatus","P. stricticeps",
                 "L. russellii" , "P. xanthura",
                 "S. interrupta",
                 "T. lunare")

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
  ggplot_z_2018 <-fric.plot(ggplot_bg = ggplot_z_empty, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "seagreen4", asb2 =NA),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # legend and title
  if (z==1) {
    ggplot_z_2018 <- ggplot_z_2018  + labs(title = "2018" )+
      theme(plot.title = element_text(size = 20, color = "seagreen4"))
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
                               size = 4, colour= "black",
                               fontface = "plain",
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = "black")
  }
  
# ggplot stored in list
  ggplot_2018[[z]] <- ggplot_z_2018
  
}# end of z


############## Now without species names

plot_sp_nm <- NULL

# plot convex hull of assemblage but not species
ggplot_2018_no_sp <-fric.plot(ggplot_bg = ggplot_z_PC1, 
                          asb_sp_coord2D = asb_sp_coord2D_k,
                          asb_vertices_nD = vertices_nD_k,
                          plot_sp = TRUE,
                          color_sp = thermal_aff_colors,
                          fill_sp = c(asb1 = "white", asb2 = "white"),
                          size_sp = c(asb1 = 1, asb2 = 1),
                          shape_sp = c(asb1 = 16, asb2 = 16),
                          color_vert = thermal_aff_colors,
                          fill_vert = thermal_aff_colors,
                          size_vert = c(asb1 = 1, asb2 = 1),
                          shape_vert = c(asb1 = 16, asb2 = 16),
                          alpha_ch = c(asb1 = 0, asb2 = 0),
                          color_ch = c(asb1 = "seagreen4", asb2 =NA),
                          fill_ch = c(asb1 = NA, asb2 = NA))

# legend and title

ggplot_2018_no_sp <- ggplot_2018_no_sp  + labs(title = "2018" )+
    theme(plot.title = element_text(size = 20, color = "seagreen4"))

######################## SPATIAL ####

# computing occurrences of species in each habitat
hab_sp_occ <- rbind( 
  Inshore = apply(spatial_sp_biom [spatial_metadata[which(spatial_metadata$Habitat=="Inshore"),"Code"],],2,max ),
  Midshelf = apply(spatial_sp_biom [spatial_metadata[which(spatial_metadata$Habitat=="Midshelf"),"Code"],],2,max ),
  Offshore = apply(spatial_sp_biom [spatial_metadata[which(spatial_metadata$Habitat=="Offshore"),"Code"],],2,max )
)  

########## INSHORE

# # Retrieve species coordinates matrix for inshore:

inshore_occ <- hab_sp_occ %>% 
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
plot_sp_nm <- c( "A. maculatus", "A. viridis" ,
                 "G. elevata","P. affinis",
                 "P. oligolepis")

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
  ggplot_z_inshore <-fric.plot(ggplot_bg = ggplot_z_full, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = c(asb1 = "#2C6BAA"),
                            fill_sp = c(asb1 = "white"),
                            size_sp = c(asb1 = 1),
                            shape_sp = c(asb1 = 16),
                            color_vert = c(asb1 = "#2C6BAA"),
                            fill_vert = c(asb1 = "#2C6BAA"),
                            size_vert = c(asb1 = 1),
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
                               size = 4, colour= "black",
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

midshelf_occ <- hab_sp_occ %>% 
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
plot_sp_nm <- c("A. chinensis",  "H. margaritaceus", 
                "S. fuscescens", 
                "T. lunare","T. lutescens",
                "A. viridis", "K. sydneyanus",
                "P. unifasciata")

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
  ggplot_z_midshelf <-fric.plot(ggplot_bg = ggplot_z_PC1, 
                                asb_sp_coord2D = asb_sp_coord2D_k,
                                asb_vertices_nD = vertices_nD_k,
                                plot_sp = TRUE,
                                color_sp = thermal_aff_colors,
                                fill_sp = c(asb1 = "white", asb2 = "white"),
                                size_sp = c(asb1 = 1, asb2 = 1),
                                shape_sp = c(asb1 = 16, asb2 = 16),
                                color_vert = thermal_aff_colors,
                                fill_vert = thermal_aff_colors,
                                size_vert = c(asb1 = 1, asb2 = 1),
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
                               size = 4, colour= "black",
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

offshore_occ <- hab_sp_occ %>% 
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
plot_sp_nm <- c("C. caerulaurea",
                "H. margaritaceus",
                "L. dimidiatus", "M. vanicolensis", 
                "N. unicornis", "P. flavomaculatus", "P. hepatus", 
                "A. viridis")

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
  ggplot_z_offshore <-fric.plot(ggplot_bg = ggplot_z_PC1, 
                            asb_sp_coord2D = asb_sp_coord2D_k,
                            asb_vertices_nD = vertices_nD_k,
                            plot_sp = TRUE,
                            color_sp = thermal_aff_colors,
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = thermal_aff_colors,
                            fill_vert = thermal_aff_colors,
                            size_vert = c(asb1 = 1, asb2 = 1),
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
                               size = 4, colour= "black",
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
      theme(plot.title = element_text(size = 20, color = "firebrick3")) +
      add_fishape(family = "Acanthuridae",
                  option = "Naso_unicornis",
                  xmin =  0.25 ,xmax = 0.35, ymin = -0.20, ymax = -0.35,
                  fill = "black",
                  alpha = 1)
  }
  # ggplot stored in list
  ggplot_offshore[[z]] <- ggplot_z_offshore
}# end of z

## merging all plots into the different figures and saving as jpeg ####

figure4 <- (ggplot_2002[[1]] + ggplot_2009[[1]] + ggplot_2018[[1]])/(ggplot_inshore[[1]] + ggplot_midshelf[[1]] + ggplot_offshore[[1]])

ggsave(figure4, file=here::here("outputs", "Figure4.jpeg"),
       height = 25, width = 45, unit = "cm" )

figure7s_a <-  (ggplot_2002_no_sp | ggplot_2003[[1]] | ggplot_2004[[1]] |  ggplot_2005[[1]])/
              (ggplot_2006[[1]] |  ggplot_2007[[1]] | ggplot_2008[[1]] | ggplot_2009_no_sp )/
  (ggplot_2010[[1]] | ggplot_2011[[1]] | ggplot_2013[[1]] | ggplot_2014[[1]])/
 (ggplot_2015[[1]] |  ggplot_2018_no_sp)

ggsave(figure7s_a , file=here::here("outputs", "Figure7s.jpeg"),
       height = 30, width = 30, unit = "cm" )

################################### end of code ##########################################