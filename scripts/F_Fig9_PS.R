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
load(here::here("data", "spatial_sp_occ.RData") )
load(here::here("outputs/", "spatial_fd.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "kelp_sp_occ.RData") )
load(here::here("outputs/", "temporal_fd_kelp.RData") )


## settings ####

# vertices of all fe in 4D ----
pool_vert_nm <- spatial_fd$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(sp_3D_coord)
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]


## temporal kelp ####

# computing occurrences of species in each year (we will use 2002-2010-2018)

kelp_years_sp_occ <- rbind( 
  y2002 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2002"),"Code"],],2,max ),
  y2010 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2010"),"Code"],],2,max ),
  y2018 = apply(kelp_sp_occ [kelp_metadata[which(kelp_metadata$Year=="2018"),"Code"],],2,max ))  


# # Retrieve species coordinates matrix for year 2002:
kelp_2002_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2002") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

# loading thermal affinity data:
thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  #column_to_rownames("Species") %>% 
  select(-thermal)

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

## Because there are not enouh tropical species, I'll only do it in 2 axes

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
ggplot_pc <- list()

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
                        color_sp = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                        fill_sp = c(asb1 = "white", asb2 = "white"),
                        size_sp = c(asb1 = 1, asb2 = 1),
                        shape_sp = c(asb1 = 16, asb2 = 16),
                        color_vert = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                        fill_vert = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                        size_vert = c(asb1 = 4, asb2 = 4),
                        shape_vert = c(asb1 = 16, asb2 = 16),
                        alpha_ch = c(asb1 = 0, asb2 = 0),
                        color_ch = c(asb1 = NA, asb2 =" red"),
                        fill_ch = c(asb1 = NA, asb2 = NA))
  
  # ggplot stored in list
  ggplot_pc[[z]] <- ggplot_z2
  
  
}# end of z


################ 2010

# # Retrieve species coordinates matrix for year 2010:

kelp_2010_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2010") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

# loading thermal affinity data:
thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  #column_to_rownames("Species") %>% 
  select(-thermal)

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

# compute FRic for all habitats  ---
Fric_2010 <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord,
                          asb_sp_w = kelp_2010_occ_thermal,
                          ind_vect = c("fric"),
                          scaling = TRUE,
                          details_returned = TRUE)

## plotting  ####

# list to store ggplot
ggplot_2010 <- list()

# pairs of axes
pairs_axes <- list( c(1,2), c(1,3) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                              faxes_nm = paste0("PC", xy), 
                              color_bg = "grey95")
  
  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                        sp_coord2D = sp_3D_coord,
                        vertices_nD = pool_vert_nm,
                        plot_pool = FALSE,
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
  
  
  
  # species present in trop:
  sp_trop_2010 <- sp_thermal_2010$Species[which(sp_thermal_2010$thermal_label == "tropical")]
  
  # species present in temp:
  sp_temp_2010 <- sp_thermal_2010$Species[which(sp_thermal_2010$thermal_label == "temperate")]
  
  # vertices in trop:
  vert_trop_2010 <- Fric_2010$details$asb_vert_nm$tropical
  
  # vertices in temp:
  vert_temp_2010 <- Fric_2010$details$asb_vert_nm$temperate
  
  # plot convex hull of assemblage but not species
  ggplot_z_2010 <-fric.plot(ggplot_bg = ggplot_z, 
                        asb_sp_coord2D = list(asb1 = sp_2d_coord_final_2010[sp_trop_2010, xy], 
                                              asb2 = sp_2d_coord_final_2010[sp_temp_2010, xy]),
                        asb_vertices_nD = list(asb1 = vert_trop_2010, 
                                               asb2 = vert_temp_2010),
                        plot_sp = TRUE,
                        color_sp = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                        fill_sp = c(asb1 = "white", asb2 = "white"),
                        size_sp = c(asb1 = 1, asb2 = 1),
                        shape_sp = c(asb1 = 16, asb2 = 16),
                        color_vert = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                        fill_vert = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                        size_vert = c(asb1 = 4, asb2 = 4),
                        shape_vert = c(asb1 = 16, asb2 = 16),
                        alpha_ch = c(asb1 = 0, asb2 = 0),
                        color_ch = c(asb1 = NA, asb2 ="#2C6BAA"),
                        fill_ch = c(asb1 = NA, asb2 = NA))
  
  # ggplot stored in list
  ggplot_2010[[z]] <- ggplot_z_2010
  
  
}# end of z



################ 2018

# # Retrieve species coordinates matrix for year 2018:

kelp_2018_occ <- kelp_years_sp_occ %>% # I'm sure there is an easier way to do this...
  as.data.frame() %>%
  rownames_to_column("Sites") %>%
  filter(Sites == "y2018") %>%
  column_to_rownames("Sites") %>%
  as.matrix()

# loading thermal affinity data:
thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(thermal_label= if_else(thermal>"23", "tropical", "temperate")) %>%   
  #column_to_rownames("Species") %>% 
  select(-thermal)

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
pairs_axes <- list( c(1,2), c(1,3) )

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                              faxes_nm = paste0("PC", xy), 
                              color_bg = "grey95")
  
  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                        sp_coord2D = sp_3D_coord,
                        vertices_nD = pool_vert_nm,
                        plot_pool = FALSE,
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
  
  
  
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
                            color_sp = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                            fill_sp = c(asb1 = "white", asb2 = "white"),
                            size_sp = c(asb1 = 1, asb2 = 1),
                            shape_sp = c(asb1 = 16, asb2 = 16),
                            color_vert = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                            fill_vert = c(asb1 = "lightsalmon1", asb2 = "#2C6BAA"),
                            size_vert = c(asb1 = 4, asb2 = 4),
                            shape_vert = c(asb1 = 16, asb2 = 16),
                            alpha_ch = c(asb1 = 0, asb2 = 0),
                            color_ch = c(asb1 = "yellow", asb2 ="green"),
                            fill_ch = c(asb1 = NA, asb2 = NA))
  
  # ggplot stored in list
  ggplot_2018[[z]] <- ggplot_z_2018
  
  ## merging all plots into a single figure and saving as png ####
  figure9 <- ( ggplot_z2 + ggplot_2010[[1]] +  ggplot_2018[[1]] )
    
#######################HERE!! Change color of convex hull 
  
  
}# end of z


## spatial ####

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

# FIde
# fide_hab<- spatial_metadata %>%
#   select(Habitat) %>%
#   bind_cols( select(spatial_fd$functional_diversity_indices, fide_PC1, fide_PC2 , fide_PC3) )


# color code for the 3 habitats
hab_colors <- c(Inshore="gold1", Midshelf="red2", Offshore="blue3")


## plotting  ####

# list to store ggplot
ggplot_spatial<-list()

# pairs of axes
pairs_axes<-list( c(1,2), c(1,3) )

for ( a in 1:length(pairs_axes) ) {
  
  # names of axes   
  xy<-pairs_axes[[a]]
  
  # background with axes range set + title
  ggplot_a<-background.plot(range_faxes=range_axes,
                            faxes_nm=paste0("PC", xy), 
                            color_bg="grey95")
  
  # convex hull of species pool
  ggplot_a<-pool.plot(ggplot_bg=ggplot_a,
                      sp_coord2D=pool_coord[,xy],
                      vertices_nD=pool_vert_nm,
                      plot_pool=FALSE,
                      color_ch=NA, fill_ch="white", alpha_ch=1)
  
  # loop on habitats
  for (b in row.names(hab_sp_occ) ) {
    
    # color given habitat
    col_b<-as.character(hab_colors[b])
    
    # species present in v
    sp_b<-names(which(hab_multidimFD$details$asb_sp_occ[b,]==1))
    
    # plot convex hull of assemblage but not species
    ggplot_a<-fric.plot( ggplot_bg=ggplot_a, 
                         asb_sp_coord2D=list(vv=pool_coord[sp_b,xy]),
                         asb_vertices_nD=list(vv=hab_multidimFD$details$asb_vert_nm[[b]]),
                         plot_sp = FALSE,
                         # size_sp = c(vv=2),
                         # size_vert = c(vv=2),
                         # shape_sp = c(vv=4),
                         # shape_vert = c(vv=7),
                         #color_sp = c(vv=col_b),
                         color_ch = c(vv=col_b),
                         #color_vert = c(vv=col_b),
                         #fill_vert = c(vv=col_b),
                         #fill_sp = c(vv=col_b),
                         fill_ch= c(vv=col_b),
                         alpha_ch = c(vv=0.1))
    
  }# end of b
  
  # average position of species (in the 3 sites) of each habitat
  ggplot_a <- ggplot_a +
    # geom_point(data=fide_hab, aes_string( x=paste0("fide_PC",xy[1]), y=paste0("fide_PC",xy[2]),
    #                            color = "Habitat", shape = "Habitat"), size=0.9, show.legend = FALSE) +
    scale_colour_manual( values = hab_colors )
  
  
  # legend
  if (a==1) {
    nlabels_b <- length(hab_colors)
    ggplot_a <- ggplot_a + 
      geom_text(aes(x = rep(range_axes[1]*0.85, nlabels_b ),
                    y = range_axes[2]*(1.05-0.15*(1:nlabels_b)),
                    label = names(hab_colors) ) ,
                color = hab_colors ,
                hjust = "left",
                show.legend = FALSE) +
      ggtitle("Spatial")
  } #end of legend
  
  # ggplot stored in list
  ggplot_spatial[[a]] <- ggplot_a
  
  
}# end of a



## merging all plots into a single figure and saving as png ####
figure3 <- ( ggplot_temporal_kelp[[1]] + ggplot_spatial[[1]] ) / ( 
  
  ggplot_temporal_kelp[[2]] + ggplot_spatial[[2]])


ggsave(figure3, file=here::here("outputs/", "Figure3_bis.png"),
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
