################################################################################
##
## Script for plotting filling of functional space by fish assemblages
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
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
load(here::here("data", "spatial_sp_occ.RData") )
load(here::here("outputs/", "spatial_fd.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "kelp_sp_occ.RData") )
load(here::here("data", "nokelp_metadata.RData") )

load(here::here("data", "nokelp_sp_occ.RData") )


## settings ####

## Prepare thermal affinity data

thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  column_to_rownames("Species") %>% 
  select(-thermal)

## temporal kelp ####

# computing occurrences of species in each habitat
kelp_years_sp_occ <- rbind( 
  y2002 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2002"),"Code"],],2,max ),
  y2008 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2008"),"Code"],],2,max ),
  y2013 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2013"),"Code"],],2,max ),
  y2018 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2018"),"Code"],],2,max )
)  

# compute FRic for kelp  ---
kelp_years_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_3D_coord, 
                                         asb_sp_w = kelp_years_sp_occ,
                                         ind_vect = c("fric"), 
                                         scaling = TRUE, 
                                         details_returned = TRUE
)


kelp_years_multidimFD$details


## plotting  ####

## get number of dimensions in input:
nb_dim <- ncol(sp_3D_coord)

## define arguments:
faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("y2002", "y2008", "y2013", "y2018")
plot_sp_nm <- NULL



## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_3D_coord)[1:min(c(4, nb_dim))]
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
  range_sp_coord  <- range(sp_3D_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)...
# ... if required:
sp_faxes_coord_plot <- data.frame(sp_3D_coord, label = "")

# if some species names to be plotted, adding a character variable to ...
# ... sp_faxes_coord:
if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
pool_vert_nm<-spatial_fd$details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")
asb3 <- plot_asb_nm[3]
nm_asb <- paste(nm_asb, asb3, sep = "_")
asb4 <- plot_asb_nm[4]
nm_asb <- paste(nm_asb, asb4, sep = "_")


sp_asb1 <- names(which(kelp_years_multidimFD$details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(kelp_years_multidimFD$details$asb_sp_occ[asb2, ] == 1))
sp_asb3 <- names(which(kelp_years_multidimFD$details$asb_sp_occ[asb3, ] == 1))
sp_asb4 <- names(which(kelp_years_multidimFD$details$asb_sp_occ[asb4, ] == 1))

### Plot FRic:

# list to store ggplot
panels_fric_kelp <- list()

# loop on combinations:
for (k in (1:plot_nb)) {
  
  # names of axes
  xy_k <- axes_plot[1:2, k]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- kelp_years_multidimFD$details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- kelp_years_multidimFD$details$asb_vert_nm[[asb2]]
  asb_sp_coord2D_k[["asb3"]] <- sp_coord_xy[sp_asb3, ]
  vertices_nD_k[["asb3"]] <- kelp_years_multidimFD$details$asb_vert_nm[[asb3]]
  asb_sp_coord2D_k[["asb4"]] <- sp_coord_xy[sp_asb4, ]
  vertices_nD_k[["asb4"]] <- kelp_years_multidimFD$details$asb_vert_nm[[asb4]]
  
  # background = axes defined by range of values and names as specified:
  plot_k <- mFD::background.plot(range_faxes, faxes_nm = xy_k, color_bg = "grey95")
  
  # add species pool:
  plot_k <- mFD::pool.plot(ggplot_bg = plot_k,
                           sp_coord2D = sp_coord_xy,
                           vertices_nD = pool_vert_nm,
                           plot_pool = TRUE,
                           color_ch = NA,
                           fill_ch = "white",
                           alpha_ch = 1,
                           shape_pool = 3,
                           size_pool = 0.7,
                           color_pool = "grey50",
                           fill_pool = NA,
                           shape_vert = 3,
                           size_vert = 0.7,
                           color_vert = "grey50",
                           fill_vert = NA)
  
  # plot 2D convex hulls and points for the 4 assemblages:
  plot_k <- mFD::fric.plot(ggplot_bg = plot_k,
                           asb_sp_coord2D = asb_sp_coord2D_k,
                           asb_vertices_nD = vertices_nD_k,
                           plot_sp = TRUE,
                           color_ch = c(asb1 = "green3", asb2 = "blue3",
                                        asb3 = "gold1", asb4 = "red"),
                           fill_ch = c(asb1 = "green3", asb2 = "blue3",
                                       asb3 = "gold1", asb4 = "red"),
                           alpha_ch = c(asb1 = 0.5, asb2 = 0.3, asb3 = 0.2, asb4 = 0.1),
                           shape_sp = c(asb1 = 15, asb2 = 17, asb3 = 16, asb4 = 12),
                           size_sp = c(asb1 = 1, asb2 = 1, asb3 = 1, asb4 = 1),
                           color_sp = c(asb1 = "green3", asb2 = "blue3",
                                        asb3 = "gold1", asb4 = "red"),
                           fill_sp = c(asb1 = "green3", asb2 = "blue3",
                                       asb3 = "gold1", asb4 = "red"),
                           shape_vert = c(asb1 = 15, asb2 = 17, asb3 = 16, asb4 = 12),
                           size_vert = c(asb1 = 1, asb2 = 1, asb3 = 1, asb4 = 1),
                           color_vert = c(asb1 = "green3", asb2 = "blue3",
                                          asb3 = "gold1", asb4 = "red"),
                           fill_vert = c(asb1 = "green3", asb2 = "blue3",
                                         asb3 = "gold1", asb4 = "red"))
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    plot_k <- plot_k +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_k[1],
                                                   y = xy_k[2],
                                                   label = "label"),
                               size = size_sp_nm, colour= color_sp_nm,
                               fontface = fontface_sp_nm,
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = color_sp_nm)
  }
  
  # save plot in a list:
  panels_fric_kelp[[k]] <- plot_k
  
}
  
if (k==1) {
  plot_k <- plot_k + 
    ggtitle("Temporal - No Kelp")
} # end of k


################ temporal no kelp ##############################

# computing occurrences of species in each habitat
nokelp_years_sp_occ <- rbind( 
  y2002 = apply(nokelp_sp_occ [nokelp_metadata[which(nokelp_metadata$Year=="2002"),"Code"],],2,max ),
  y2008 = apply(nokelp_sp_occ [nokelp_metadata[which(nokelp_metadata$Year=="2008"),"Code"],],2,max ),
  y2013 = apply(nokelp_sp_occ [nokelp_metadata[which(nokelp_metadata$Year=="2013"),"Code"],],2,max ),
  y2018 = apply(nokelp_sp_occ [nokelp_metadata[which(nokelp_metadata$Year=="2018"),"Code"],],2,max )
)  

# compute FRic for all habitats  ---
nokelp_years_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_3D_coord, 
                                           asb_sp_w = nokelp_years_sp_occ,
                                           ind_vect = c("fric"), 
                                           scaling = TRUE, 
                                           details_returned = TRUE
)


## plotting  ####

## get number of dimensions in input:
nb_dim <- ncol(sp_3D_coord)

## define arguments:
faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("y2002", "y2008", "y2013", "y2018")
plot_sp_nm <- NULL


## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_3D_coord)[1:min(c(4, nb_dim))]
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
  range_sp_coord  <- range(sp_3D_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)...
# ... if required:
sp_faxes_coord_plot <- data.frame(sp_3D_coord, label = "")

# if some species names to be plotted, adding a character variable to ...
# ... sp_faxes_coord:
if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
pool_vert_nm<-spatial_fd$details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")
asb3 <- plot_asb_nm[3]
nm_asb <- paste(nm_asb, asb3, sep = "_")
asb4 <- plot_asb_nm[4]
nm_asb <- paste(nm_asb, asb4, sep = "_")


sp_asb1 <- names(which(nokelp_years_multidimFD$details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(nokelp_years_multidimFD$details$asb_sp_occ[asb2, ] == 1))
sp_asb3 <- names(which(nokelp_years_multidimFD$details$asb_sp_occ[asb3, ] == 1))
sp_asb4 <- names(which(nokelp_years_multidimFD$details$asb_sp_occ[asb4, ] == 1))

### Plot FRic:

# list to store ggplot
panels_fric_nokelp <- list()

# loop on combinations:
for (k in (1:plot_nb)) {
  
  # names of axes
  xy_k <- axes_plot[1:2, k]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- nokelp_years_multidimFD$details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- nokelp_years_multidimFD$details$asb_vert_nm[[asb2]]
  asb_sp_coord2D_k[["asb3"]] <- sp_coord_xy[sp_asb3, ]
  vertices_nD_k[["asb3"]] <- nokelp_years_multidimFD$details$asb_vert_nm[[asb3]]
  asb_sp_coord2D_k[["asb4"]] <- sp_coord_xy[sp_asb4, ]
  vertices_nD_k[["asb4"]] <- nokelp_years_multidimFD$details$asb_vert_nm[[asb4]]
  
  # background = axes defined by range of values and names as specified:
  plot_k <- mFD::background.plot(range_faxes, faxes_nm = xy_k, color_bg = "grey95")
  
  # add species pool:
  plot_k <- mFD::pool.plot(ggplot_bg = plot_k,
                           sp_coord2D = sp_coord_xy,
                           vertices_nD = pool_vert_nm,
                           plot_pool = TRUE,
                           color_ch = NA,
                           fill_ch = "white",
                           alpha_ch = 1,
                           shape_pool = 3,
                           size_pool = 0.7,
                           color_pool = "grey50",
                           fill_pool = NA,
                           shape_vert = 3,
                           size_vert = 0.7,
                           color_vert = "grey50",
                           fill_vert = NA)
  
  # plot 2D convex hulls and points for the 4 assemblages:
  plot_k <- mFD::fric.plot(ggplot_bg = plot_k,
                           asb_sp_coord2D = asb_sp_coord2D_k,
                           asb_vertices_nD = vertices_nD_k,
                           plot_sp = TRUE,
                           color_ch = c(asb1 = "green3", asb2 = "blue3",
                                        asb3 = "gold1", asb4 = "red"),
                           fill_ch = c(asb1 = "green3", asb2 = "blue3",
                                       asb3 = "gold1", asb4 = "red"),
                           alpha_ch = c(asb1 = 0.5, asb2 = 0.3, asb3 = 0.2, asb4 = 0.1),
                           shape_sp = c(asb1 = 15, asb2 = 17, asb3 = 16, asb4 = 12),
                           size_sp = c(asb1 = 1, asb2 = 1, asb3 = 1, asb4 = 1),
                           color_sp = c(asb1 = "green3", asb2 = "blue3",
                                        asb3 = "gold1", asb4 = "red"),
                           fill_sp = c(asb1 = "green3", asb2 = "blue3",
                                       asb3 = "gold1", asb4 = "red"),
                           shape_vert = c(asb1 = 15, asb2 = 17, asb3 = 16, asb4 = 12),
                           size_vert = c(asb1 = 1, asb2 = 1, asb3 = 1, asb4 = 1),
                           color_vert = c(asb1 = "green3", asb2 = "blue3",
                                          asb3 = "gold1", asb4 = "red"),
                           fill_vert = c(asb1 = "green3", asb2 = "blue3",
                                         asb3 = "gold1", asb4 = "red"))
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    plot_k <- plot_k +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_k[1],
                                                   y = xy_k[2],
                                                   label = "label"),
                               size = size_sp_nm, colour= color_sp_nm,
                               fontface = fontface_sp_nm,
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = color_sp_nm)
  }
  
  # save plot in a list:
  panels_fric_nokelp[[k]] <- plot_k
  
}

if (k==1) {
  plot_k <- plot_k + 
    ggtitle("Temporal - No Kelp")
} # end of k



################### spatial ##########################

# computing occurrences of species in each habitat
hab_sp_occ <- rbind( 
  Inshore = apply(spatial_sp_occ [spatial_metadata[which(spatial_metadata$Habitat=="Inshore"),"Code"],],2,max ),
  Midshelf = apply(spatial_sp_occ [spatial_metadata[which(spatial_metadata$Habitat=="Midshelf"),"Code"],],2,max ),
  Offshore = apply(spatial_sp_occ [spatial_metadata[which(spatial_metadata$Habitat=="Offshore"),"Code"],],2,max )
)  

# compute FRic for all habitats  ---
hab_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_3D_coord, 
                                  asb_sp_w = hab_sp_occ,
                                  ind_vect = c("fric"), 
                                  scaling = TRUE, 
                                  details_returned = TRUE
)

# color code for the 3 habitats
hab_colors <- c(Inshore="gold1", Midshelf="red2", Offshore="blue3")

## plotting  ####

## get number of dimensions in input:
nb_dim <- ncol(sp_3D_coord)

## define arguments:
faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("Inshore", "Midshelf", "Offshore")
plot_sp_nm <- NULL


## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_3D_coord)[1:min(c(4, nb_dim))]
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
  range_sp_coord  <- range(sp_3D_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)...
# ... if required:
sp_faxes_coord_plot <- data.frame(sp_3D_coord, label = "")

# if some species names to be plotted, adding a character variable to ...
# ... sp_faxes_coord:
if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
pool_vert_nm<-spatial_fd$details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")
asb3 <- plot_asb_nm[3]
nm_asb <- paste(nm_asb, asb3, sep = "_")


sp_asb1 <- names(which(hab_multidimFD$details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(hab_multidimFD$details$asb_sp_occ[asb2, ] == 1))
sp_asb3 <- names(which(hab_multidimFD$details$asb_sp_occ[asb3, ] == 1))

### Plot FRic:

# list to store ggplot
panels_fric_spatial <- list()

# loop on combinations:
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
  vertices_nD_k[["asb1"]] <- hab_multidimFD$details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- hab_multidimFD$details$asb_vert_nm[[asb2]]
  asb_sp_coord2D_k[["asb3"]] <- sp_coord_xy[sp_asb3, ]
  vertices_nD_k[["asb3"]] <- hab_multidimFD$details$asb_vert_nm[[asb3]]
  
  
  # background = axes defined by range of values and names as specified:
  plot_k <- mFD::background.plot(range_faxes, faxes_nm = xy_z, color_bg = "grey95")
  
  # add species pool:
  plot_k <- mFD::pool.plot(ggplot_bg = plot_k,
                           sp_coord2D = sp_coord_xy,
                           vertices_nD = pool_vert_nm,
                           plot_pool = TRUE,
                           color_ch = NA,
                           fill_ch = "white",
                           alpha_ch = 1,
                           shape_pool = 3,
                           size_pool = 0.7,
                           color_pool = "grey50",
                           fill_pool = NA,
                           shape_vert = 3,
                           size_vert = 0.7,
                           color_vert = "grey50",
                           fill_vert = NA)
  
  # plot 2D convex hulls and points for the 3 assemblages:
  plot_k <- mFD::fric.plot(ggplot_bg = plot_k,
                           asb_sp_coord2D = asb_sp_coord2D_k,
                           asb_vertices_nD = vertices_nD_k,
                           plot_sp = TRUE,
                           color_ch = c(asb1 = "green3", asb2 = "blue3",
                                        asb3 = "gold1"),
                           fill_ch = c(asb1 = "green3", asb2 = "blue3",
                                       asb3 = "gold1"),
                           alpha_ch = c(asb1 = 0.5, asb2 = 0.3, asb3 = 0.2),
                           shape_sp = c(asb1 = 15, asb2 = 17, asb3 = 16),
                           size_sp = c(asb1 = 1, asb2 = 1, asb3 = 1),
                           color_sp = c(asb1 = "green3", asb2 = "blue3",
                                        asb3 = "gold1"),
                           fill_sp = c(asb1 = "green3", asb2 = "blue3",
                                       asb3 = "gold1"),
                           shape_vert = c(asb1 = 15, asb2 = 17, asb3 = 16),
                           size_vert = c(asb1 = 1, asb2 = 1, asb3 = 1),
                           color_vert = c(asb1 = "green3", asb2 = "blue3",
                                          asb3 = "gold1"),
                           fill_vert = c(asb1 = "green3", asb2 = "blue3",
                                         asb3 = "gold1"))
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    plot_k <- plot_k +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_z[1],
                                                   y = xy_z[2],
                                                   label = "label"),
                               size = size_sp_nm, colour= color_sp_nm,
                               fontface = fontface_sp_nm,
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = color_sp_nm)
  }
  
  # save plot in a list:
  panels_fric_spatial[[k]] <- plot_k
  
}
#end of z

                        
  # legend
  if (k==1) {
    nlabels_z <- length(hab_colors)
    plot_k <- plot_k + 
      geom_text(aes(x = rep(range_axes[1]*0.9, nlabels_z ),
                    y = range_axes[2]*(1.05-0.15*(1:nlabels_z)),
                    label = names(hab_colors),
                    color = hab_colors ),
                hjust = "left",
                show.legend = FALSE) +
      ggtitle("Spatial")
  }
  
  # end of z



## merging all plots into a single figure and saving as png ####
figure3 <- ( panels_fric_nokelp[[1]] +  panels_fric_kelp[[1]] + panels_fric_spatial[[1]] ) / ( 
  
  panels_fric_nokelp[[2]] +  panels_fric_kelp[[2]] + panels_fric_spatial[[2]] )


ggsave(figure3, file=here::here("outputs/", "figure3.png"),
       height = 16, width = 24, unit = "cm" )
