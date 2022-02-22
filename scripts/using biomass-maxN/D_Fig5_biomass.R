################################################################################
##
## Script for plotting functional space with thermal affinity information
## 
## Code by Camille Magneville, Paula Sgarlatta and Sebastien Villeger 
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(patchwork)
library(mFD)


# loading data
load(here::here("outputs", "using biomass-maxN", "spatial_fd_biomass.RData"))

# loading thermal affinity data:
thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  #column_to_rownames("Species") %>% 
  select(-thermal)


# Add thermal aff to sp_faxes_coord:
sp_3D_coord <- spatial_fd$details$sp_faxes_coord

sp_faxes_coord <- as.data.frame(sp_3D_coord) %>% 
  rownames_to_column("Species")

sp_faxes_coord <- inner_join(sp_faxes_coord, thermal, 
                             by="Species") 

# Change class of thermal_affinity column: character to factor:
sp_faxes_coord$thermal_label <- as.factor(sp_faxes_coord$thermal_label)
class(sp_faxes_coord$thermal_label) # ok


# Create a new assemblage*species df with assemblages being either tropical ...
# ... or temperate:

asb_sp <- sp_faxes_coord[, c(1, 5)]

asb_sp_new <- asb_sp %>% 
  add_column(present = as.numeric(1)) %>% 
  pivot_wider(names_from = thermal_label, values_from = present)

asb_sp_new[is.na(asb_sp_new)] <- as.numeric(0)

asb_sp_new2 <- t(asb_sp_new)

colnames(asb_sp_new2) <- asb_sp_new2[1, ]

asb_sp_new2 <- asb_sp_new2[-1, ]


# create a dataframe that will contain the same values as asb_sp_new2 ...
# ... because (I don't knwo why) I can not convert character into numeric:
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

thermal_aff_colors <- c(tropical = "lightsalmon1", temperate = "#2C6BAA")


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
                         size_sp = c(asb1 = 1, asb2 = 1),
                         shape_sp = c(asb1 = 16, asb2 = 16),
                         color_vert = thermal_aff_colors,
                         fill_vert = thermal_aff_colors,
                         size_vert = c(asb1 = 4, asb2 = 4),
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


# complete with points info:
# plot_caption <- plot_caption +
  # # ggplot2::geom_point(x = range_axes[1] + spread_faxes*0.125,
  #                     y = range_axes[2] - spread_faxes*0.20,
  #                     fill = thermal_aff_colors[[1]], 
  #                     color = thermal_aff_colors[[1]],
  #                     shape = 16,
  #                     size = 4) + 
  
  # ggplot2::geom_text(x = range_axes[1] + spread_faxes*0.45,
  #                    y = range_axes[2] - spread_faxes*0.20,
  #                    label = "thermal affinity = tropical",
  #                    colour = thermal_aff_colors[[1]], size = 4) + 
  
 
  
  # #ggplot2::geom_text(x = range_axes[1] + spread_faxes*0.50,
  #                    y = range_axes[2] - spread_faxes*0.50,
  #                    label = "thermal affinity = temperate", 
  #                    colour = thermal_aff_colors[[2]], size = 4) 


## merging all plots into a single figure and saving as png ####
figure <- panels.to.patchwork(ggplot_pc, plot_caption = plot_caption)


ggsave(figure, file=here::here("outputs/", "using biomass-maxN",  "Figure5_biomass.png"),
       height = 16, width = 24, unit = "cm" )
