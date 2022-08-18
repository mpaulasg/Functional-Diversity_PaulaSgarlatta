################################################################################
##
## Script for plotting several figures/tables
##
##  Fig. S1 - relation between traits and PcoAs
##
##  Table S5 - Correlation of traits and PCoAs
##
##  Fig. 3 a-b - Functional space with thermal affinity information    
##
## Code by Camille Magneville, Paula Sgarlatta and Sebastien Villeger 
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(patchwork)
library(ggplot2)
library(mFD)


# loading data

load(here::here("data", "funct_spaces.RData") )
load(here::here("data", "sp_faxes_coord.RData") )
load(here::here("data", "sp_tr.RData") )

load(here::here("data", "spatial_fd_biomass.RData"))


## Illustrate quality of functional space

qual_space <- mFD::quality.fspaces.plot(
  fspaces_quality            = funct_spaces,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


plot(qual_space) #Won't include this as a figure

## Table with mAD values

mAD_values <- funct_spaces$quality_fspaces

### Test correlation between traits and functional axes:

#Change to full traits names

sp_tr <- sp_tr %>% 
  rename("Maximum reported length" = "Size","Common aggregation" = "Agg",
         "Vertical position" = "Position")

cor_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sp_tr, 
  sp_faxes_coord = sp_faxes_coord[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# get the table of correlation:

corr_table <- as.data.frame(cor_tr_faxes$tr_faxes_stat)  # Table S5

# get the plot:

plot_corr <- cor_tr_faxes$tr_faxes_plot


## Save figures

ggsave(plot_corr, file=here::here("outputs", "FigureS1.jpeg"),
       height = 16, width = 30, unit = "cm" )

write.csv(corr_table, file=here::here("outputs", "Correlation_traits_table_S5.csv"), 
          row.names = FALSE)


#### Functional space with thermal affinity information - Fig. 3

#Loading thermal affinity data

thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  dplyr::select(-thermal)


# Add thermal aff to sp_faxes_coord
sp_3D_coord <- spatial_fd$details$sp_faxes_coord

sp_faxes_coord <- as.data.frame(sp_3D_coord) %>% 
  rownames_to_column("Species")

sp_faxes_coord <- inner_join(sp_faxes_coord, thermal, 
                             by="Species") 

# Change class of thermal_affinity column: character to factor

sp_faxes_coord$thermal_label <- as.factor(sp_faxes_coord$thermal_label)
class(sp_faxes_coord$thermal_label) # ok


# Create a new assemblage*species df with assemblages being either tropical or temperate

asb_sp <- sp_faxes_coord[, c(1, 5)]

asb_sp_new <- asb_sp %>% 
  add_column(present = as.numeric(1)) %>% 
  pivot_wider(names_from = thermal_label, values_from = present)

asb_sp_new[is.na(asb_sp_new)] <- as.numeric(0)

asb_sp_new2 <- t(asb_sp_new)

colnames(asb_sp_new2) <- asb_sp_new2[1, ]

asb_sp_new2 <- asb_sp_new2[-1, ]


# create a dataframe that will contain the same values as asb_sp_new2 because (I don't knwo why) 
# I cannot convert character into numeric

asb_sp_new3 <- as.data.frame(matrix(nrow = nrow(asb_sp_new2), ncol = ncol(asb_sp_new2)))#(s)
colnames(asb_sp_new3) <- colnames(asb_sp_new2)
rownames(asb_sp_new3) <- rownames(asb_sp_new2)

for (i in (1:nrow(asb_sp_new3))) {
  for (j in (1:ncol(asb_sp_new3))) {
    asb_sp_new3[i, j] <- as.numeric(asb_sp_new2[i, j])
  }
}

asb_sp_new3 <- as.matrix(asb_sp_new3)


## settings ####

# vertices of all fe in 4D ----
pool_vert_nm <- spatial_fd$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(sp_3D_coord[,1:3])
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord, 
                                         asb_sp_w = asb_sp_new3,
                                         ind_vect = c("fric"), 
                                         scaling = TRUE, 
                                         details_returned = TRUE)

# color code for thermal affinity

thermal_aff_colors <- c(tropical = "firebrick1", temperate = "#2C6BAA")


## plotting  ####

# list to store ggplot
ggplot_pc <- list()

# pairs of axes
pairs_axes <- list(c(1,2), c(1,3), c(2, 3))

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                            faxes_nm = paste0("PC", xy), 
                            color_bg = "grey95")

  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                      sp_coord2D = sp_3D_coord[,xy],
                      vertices_nD = pool_vert_nm,
                      plot_pool = FALSE,
                      color_ch = "black", fill_ch = "white", alpha_ch = 1)
  
  

    # species present in trop:
    sp_trop <- sp_faxes_coord$Species[which(sp_faxes_coord$thermal_label == "tropical")]
    
    # species present in temp:
    sp_temp <- sp_faxes_coord$Species[which(sp_faxes_coord$thermal_label == "temperate")]
    
    # vertices in trop:
    vert_trop <- Fric$details$asb_vert_nm$tropical
    
    # vertices in temp:
    vert_temp <- Fric$details$asb_vert_nm$temperate
    
    
    # plot convex hull of assemblage but not species
    ggplot_z2 <-fric.plot(ggplot_bg = ggplot_z, 
                         asb_sp_coord2D = list(asb1 = sp_3D_coord[sp_trop, xy], 
                                               asb2 = sp_3D_coord[sp_temp, xy]),
                         asb_vertices_nD = list(asb1 = vert_trop, 
                                                asb2 = vert_temp),
                         plot_sp = TRUE,
                         color_sp = thermal_aff_colors,
                         fill_sp = c(asb1 = "white", asb2 = "white"),
                         size_sp = c(asb1 = 3, asb2 = 3),
                         shape_sp = c(asb1 = 16, asb2 = 16),
                         color_vert = thermal_aff_colors,
                         fill_vert = thermal_aff_colors,
                         size_vert = c(asb1 = 3, asb2 = 3),
                         shape_vert = c(asb1 = 16, asb2 = 16),
                         alpha_ch = c(asb1 = 0, asb2 = 0),
                         color_ch = c(asb1 = NA, asb2 = NA),
                         fill_ch = c(asb1 = NA, asb2 = NA))
  
 
  # ggplot stored in list
  ggplot_pc[[z]] <- ggplot_z2
  
  
}# end of z

## Compute Caption:

# plot white basic window:
plot_caption <- ggplot2::ggplot(data.frame(x = range_axes, 
                                           y = range_axes),
                                ggplot2::aes(x = x, y = y)) +
  ggplot2::scale_x_continuous(limits = range_axes, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = range_axes, expand = c(0, 0)) +
  ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
  ggplot2::geom_rect(xmin = range_axes[1], xmax = range_axes[2],
                     ymin = range_axes[1], ymax = range_axes[2],
                     fill = "white", colour ="white")

## merging all plots into a single figure and saving as png ####

figure3 <- ggplot_pc [[1]] + ggplot_pc [[2]]

figureS2 <- (ggplot_pc [[1]] + ggplot_pc [[2]]) + 
                        ggplot_pc [[3]]

ggsave(figure3, file=here::here("outputs",  "Figure3.jpeg"),
       height = 16, width = 25, unit = "cm" )

ggsave(figureS2, file=here::here("outputs",  "FigureS2.jpeg"),
       height = 12, width = 32, unit = "cm" )


############################## end of code ############################################


## Delete this???


#try adding names to understand which species are the ones with extreme values

#Loading thermal affinity data

thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  dplyr::select(-thermal)


# Add thermal aff to sp_faxes_coord
sp_3D_coord <- spatial_fd$details$sp_faxes_coord

sp_faxes_coord <- as.data.frame(sp_3D_coord) %>% 
  rownames_to_column("Species")

sp_faxes_coord <- inner_join(sp_faxes_coord, thermal, 
                             by="Species") 

# Change class of thermal_affinity column: character to factor

sp_faxes_coord$thermal_label <- as.factor(sp_faxes_coord$thermal_label)
class(sp_faxes_coord$thermal_label) # ok


# Create a new assemblage*species df with assemblages being either tropical or temperate

asb_sp <- sp_faxes_coord[, c(1, 5)]

asb_sp_new <- asb_sp %>% 
  add_column(present = as.numeric(1)) %>% 
  pivot_wider(names_from = thermal_label, values_from = present)

asb_sp_new[is.na(asb_sp_new)] <- as.numeric(0)

asb_sp_new2 <- t(asb_sp_new)

colnames(asb_sp_new2) <- asb_sp_new2[1, ]

asb_sp_new2 <- asb_sp_new2[-1, ]


# create a dataframe that will contain the same values as asb_sp_new2 because (I don't knwo why) 
# I cannot convert character into numeric

asb_sp_new3 <- as.data.frame(matrix(nrow = nrow(asb_sp_new2), ncol = ncol(asb_sp_new2)))#(s)
colnames(asb_sp_new3) <- colnames(asb_sp_new2)
rownames(asb_sp_new3) <- rownames(asb_sp_new2)

for (i in (1:nrow(asb_sp_new3))) {
  for (j in (1:ncol(asb_sp_new3))) {
    asb_sp_new3[i, j] <- as.numeric(asb_sp_new2[i, j])
  }
}

asb_sp_new3 <- as.matrix(asb_sp_new3)


## settings ####

# vertices of all fe in 4D ----
pool_vert_nm <- spatial_fd$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(sp_3D_coord[,1:3])
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord, 
                          asb_sp_w = asb_sp_new3,
                          ind_vect = c("fric"), 
                          scaling = TRUE, 
                          details_returned = TRUE)

## retrieve names of main input:
asb_fd_ind <- Fric$functional_diversity_indices
fd_details <- Fric$details

# ## get coordinates of species:
# sp_faxes_coord <- fd_details$sp_faxes_coord

# color code for thermal affinity

thermal_aff_colors <- c(tropical = "firebrick1", temperate = "#2C6BAA")

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)


#Check vertices to choose which species to plot

# vertices in trop:
vert_trop <- Fric$details$asb_vert_nm$tropical
vert_trop

# vertices in temp:
vert_temp <- Fric$details$asb_vert_nm$temperate
vert_temp

## Define arguments

faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("temperate", "tropical")
plot_sp_nm <- c( "Centropogon_australis" ,   "Chaetodon_guentheri" ,     "Chrysophrys_auratus"  ,    "Dasyatis_brevicaudata"  , 
                 "Enoplosus_armatus" ,       "Epinephelus_daemelii" ,    "Eupetrichthys_angustipes", "Heterodontus_galeatus",   
                  "Kyphosus_sydneyanus",      "Orectolobus_maculatus",    "Parma_oligolepis" ,        "Sarda_australis" ,        
                "Scorpaena_cardinalis" ,    "Scorpis_lineolata"  ,      "Trachinops_taeniatus" ,
                "Abudefduf_bengalensis"  ,   "Aetobatus_narinari" ,       "Aluterus_monoceros" ,       "Caesio_caerulaurea"  ,     
                 "Chaetodon_kleinii"  ,       "Epinephelus_fasciatus"  ,   "Halichoeres_margaritaceus", "Leptojulis_cyanopleura" ,  
                "Lethrinus_nebulosus" ,      "Lutjanus_russellii"   ,     "Nelusetta_ayraud"  ,        "Pterocaesio_digramma",     
                 "Seriola_rivoliana"  ,       "Stegastes_apicalis" ,       "Stethojulis_interrupta" ,   "Synodus_jaculum" ,         
               "Thalassoma_jansenii"  )

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
nm_asb <- paste(nm_asb, asb2, sep = " ")

sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))

## plotting  ####

# list to store ggplot
ggplot_pc <- list()

# pairs of axes
pairs_axes <- list(c(1,2), c(1,3), c(2, 3))

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy_z <- pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                              faxes_nm = paste0("PC", xy_z), 
                              color_bg = "grey95")
  
  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                        sp_coord2D = sp_3D_coord[,xy_z],
                        vertices_nD = pool_vert_nm,
                        plot_pool = FALSE,
                        color_ch = "black", fill_ch = "white", alpha_ch = 1)
  
  
  
  # species present in trop:
  sp_trop <- sp_faxes_coord$Species[which(sp_faxes_coord$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp <- sp_faxes_coord$Species[which(sp_faxes_coord$thermal_label == "temperate")]
  
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
  ggplot_z2 <-fric.plot(ggplot_bg = ggplot_z, 
                        asb_sp_coord2D = list(asb1 = sp_3D_coord[sp_trop, xy], 
                                              asb2 = sp_3D_coord[sp_temp, xy]),
                        asb_vertices_nD = list(asb1 = vert_trop, 
                                               asb2 = vert_temp),
                        plot_sp = TRUE,
                        color_sp = thermal_aff_colors,
                        fill_sp = c(asb1 = "white", asb2 = "white"),
                        size_sp = c(asb1 = 3, asb2 = 3),
                        shape_sp = c(asb1 = 16, asb2 = 16),
                        color_vert = thermal_aff_colors,
                        fill_vert = thermal_aff_colors,
                        size_vert = c(asb1 = 3, asb2 = 3),
                        shape_vert = c(asb1 = 16, asb2 = 16),
                        alpha_ch = c(asb1 = 0, asb2 = 0),
                        color_ch = c(asb1 = NA, asb2 = NA),
                        fill_ch = c(asb1 = NA, asb2 = NA))
  
  
  # ggplot stored in list
  ggplot_pc[[z]] <- ggplot_z2
  
  
}# end of z

## Compute Caption:

# plot white basic window:
plot_caption <- ggplot2::ggplot(data.frame(x = range_axes, 
                                           y = range_axes),
                                ggplot2::aes(x = x, y = y)) +
  ggplot2::scale_x_continuous(limits = range_axes, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = range_axes, expand = c(0, 0)) +
  ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
  ggplot2::geom_rect(xmin = range_axes[1], xmax = range_axes[2],
                     ymin = range_axes[1], ymax = range_axes[2],
                     fill = "white", colour ="white")


## merging all plots into a single figure and saving as png ####
figure <- panels.to.patchwork(ggplot_pc, plot_caption = plot_caption)

figure3 <- ggplot_pc [[1]]

figure3S <- (ggplot_pc [[1]] + ggplot_pc [[2]]) + 
  ggplot_pc [[3]]


ggsave(figure3, file=here::here("outputs",  "Figure3.jpeg"),
       height = 16, width = 16, unit = "cm" )

ggsave(figure3S, file=here::here("outputs",  "figure3S.jpeg"),
       height = 16, width = 32, unit = "cm" )

