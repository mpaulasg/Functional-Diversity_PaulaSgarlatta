################################################################################
##
## Script for plotting:
##
##  *quality of functional space and relation between traits and PcoAs
##
##  *functional space comparing spatial/temporal data (part of Figure 1)
##
##  *functional space with thermal affinity information (Figure 5)
##   
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


plot(qual_space)

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


ggsave(qual_space, file=here::here("outputs", "Figure_qual_space_extra.jpeg"),
       height = 20, width = 50, unit = "cm" )

ggsave(plot_corr, file=here::here("outputs", "FigureS6.jpeg"),
       height = 16, width = 30, unit = "cm" )

write.csv(corr_table, file=here::here("outputs", "Correlation_traits_table_S5.csv"), 
          row.names = FALSE)

write.csv(mAD_values, file=here::here("outputs", "mAD_values.csv"), 
          row.names = FALSE)



##### Functional space comparing spatial/temporal data

#Loading data from both spatial/temporal

species_both <- read.csv(here::here("data", "species_both.csv")) %>% 
  mutate(type_data= if_else(data_1 == "temporal" & data_2 == "spatial", "both",
                            if_else(data_1 == "temporal" & data_2 == "no", "temporal","spatial"))) %>% 
  dplyr::select(-data_1, -data_2)


# Add data type to sp_faxes_coord
sp_3D_coord <- spatial_fd$details$sp_faxes_coord

sp_faxes_coord <- as.data.frame(sp_3D_coord) %>% 
  rownames_to_column("Species")

sp_faxes_coord <- inner_join(sp_faxes_coord, species_both, 
                             by="Species") 

# Change class of data type column: character to factor

sp_faxes_coord$type_data <- as.factor(sp_faxes_coord$type_data)
class(sp_faxes_coord$type_data) # ok

# Create a new assemblage*species df with assemblages being either spatial, temporal or both

asb_sp <- sp_faxes_coord[, c(1, 5)]

asb_sp_new <- asb_sp %>% 
  add_column(present = as.numeric(1)) %>% 
  pivot_wider(names_from = type_data, values_from = present)

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

# color code for data type

data_colors <- c(spatial = "#FDE725FF", temporal = "seagreen4", both = "#80471C")


## plotting  ####

  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC2"), 
                              color_bg = "grey95")
  
  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                        sp_coord2D = sp_3D_coord,
                        vertices_nD = pool_vert_nm,
                        plot_pool = FALSE,
                        color_ch = "black", fill_ch = "white", alpha_ch = 1)
  
  # get names of assemblages:
  pool <- "pool"
  plot_asb_nm <- c("spatial", "temporal", "both")
  asb1 <- plot_asb_nm[1]
  nm_asb <- asb1
  asb2 <- plot_asb_nm[2]
  nm_asb <- paste(nm_asb, asb2, sep = "_")
  asb3 <- plot_asb_nm[3]
  nm_asb <- paste(nm_asb, asb3, sep = "_")
  
  
  sp_asb1 <- names(which(Fric$details$asb_sp_occ[asb1, ] == 1))
  sp_asb2 <- names(which(Fric$details$asb_sp_occ[asb2, ] == 1))
  sp_asb3 <- names(which(Fric$details$asb_sp_occ[asb3, ] == 1))
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_3D_coord[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- Fric$details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_3D_coord[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- Fric$details$asb_vert_nm[[asb2]]
  asb_sp_coord2D_k[["asb3"]] <- sp_3D_coord[sp_asb3, ]
  vertices_nD_k[["asb3"]] <- Fric$details$asb_vert_nm[[asb3]]
  
  # plot convex hull of assemblage but not species
  ggplot_z2 <-fric.plot(ggplot_bg = ggplot_z, 
                        asb_sp_coord2D = asb_sp_coord2D_k,
                        asb_vertices_nD = vertices_nD_k,
                        plot_sp = TRUE,
                        color_sp = data_colors,
                        fill_sp = c(asb1 = "white", asb2 = "white", asb3="white"),
                        size_sp = c(asb1 = 3, asb2 = 3, asb3=3),
                        shape_sp = c(asb1 = 16, asb2 = 16, asb=16),
                        color_vert = data_colors,
                        fill_vert = data_colors,
                        size_vert = c(asb1 = 3, asb2 = 3, asb3=3),
                        shape_vert = c(asb1 = 16, asb2 = 16, asb=16),
                        alpha_ch = c(asb1 = 0, asb2 = 0, asb3=0),
                        color_ch = c(asb1 = NA, asb2 = NA, asb3=NA),
                        fill_ch = c(asb1 = NA, asb2 = NA, asb3=NA))
  
  
  
## Save figure
  
ggsave(ggplot_z2, file=here::here("outputs",  "Figure1b.jpeg"),
       height = 16, width = 16, unit = "cm" )


###################### Functional space with thermal affinity information

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
figure <- panels.to.patchwork(ggplot_pc, plot_caption = plot_caption)

figure3 <- ggplot_pc [[1]]


ggsave(figure3, file=here::here("outputs",  "Figure3.jpeg"),
       height = 16, width = 16, unit = "cm" )

############################## end of code ############################################
