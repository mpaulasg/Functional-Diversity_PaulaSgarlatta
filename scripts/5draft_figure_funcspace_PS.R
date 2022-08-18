################################################################################
##
## Script for plotting changes in the multidimensional functional space across 
## habitats and over years in a multidimensional space
##
##  Figure 3 - Changes in space and time (only 2002, 2009 and 2018)
##
##  Figure S5 - changes in time (all years) and in PC3.
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

# loading data
load(here::here("data", "sp_3D_coord.RData") )
load(here::here("data", "spatial_metadata.RData") )
load(here::here("data", "spatial_sp_biom.RData") )
load(here::here("data", "spatial_fd_biomass.RData") )

load(here::here("data", "kelp_metadata.RData") )
load(here::here("data", "kelp_sp_maxN.RData") )
load(here::here("data", "temporal_fd_kelp_biomass.RData") )


# loading thermal affinity data
sp_thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(affinity= if_else(thermal>"23", "tropical", "temperate")) %>%   
  dplyr::select(-thermal)


## preparing  ###

# names of Functional axes to plot
axes_nm <- c("PC1", "PC2")

# coordinates of species along those axes
sp_coord_2D <- sp_3D_coord[,axes_nm]

# color code for thermal affinity
col_temp <- "#2C6BAA"
col_trop <- "firebrick1"

# table with species coordinates, full species names (as in other objects) 
#  shortened species names (for displaying) and thermal affinity

tab_sp <- sp_coord_2D %>% 
  as.data.frame() %>% 
  rownames_to_column("Species_name") %>% 
  mutate(genus=sub("_.*", "", Species_name), sp=sub(".*_", "", Species_name)) %>% 
  mutate(Species_code=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  dplyr::select(-genus, -sp) %>%
  left_join(sp_thermal, by=c("Species_name"="Species")) %>%
  mutate(col_affinity = case_when(affinity=="temperate" ~ col_temp,
                                  affinity=="tropical" ~ col_trop) )

head(tab_sp)

dim(sp_coord_2D)
dim(kelp_sp_maxN)
dim(spatial_sp_biom)

# names of temperate species
temperate_sp <- tab_sp %>%
  filter(affinity=="temperate") %>%
  pull(Species_name)

# range of axes given values along selected PC
range_faxes_coord <- range(sp_coord_2D)
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]


###### TEMPORAL ####

##### All years for SM

# names of years to illustrate

# years_all <-c("y2002", "y2003", "y2004", "y2005", "y2006", "y2007", "y2008", "y2009", "y2010",
#                "y2011", "y2013", "y2014", "y2015","y2018")

years_all <-c("2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010",
              "2011", "2013", "2014", "2015","2018")

# computing average biomass of species in all the years
temp_sp_biom <- kelp_sp_maxN %>%
  as.data.frame %>%
  rownames_to_column("Code") %>%
  left_join( dplyr::select(kelp_metadata, Code, Year )  ) %>%
  dplyr::select(-Code) %>%
  group_by(Year) %>%
  summarise( across(.cols = everything(), .fns = mean, .names=NULL ) ) 

#Adding a letter because is confusing with a number

#temp_sp_biom$Year<- paste0("y", temp_sp_biom$Year) 

temp_sp_biom <-  temp_sp_biom %>% 
  column_to_rownames("Year") %>%
  as.matrix()

range( apply(temp_sp_biom, 1, max) ) # all species present in at least one year


# compute FRic and FIde for all habitats in 2D (not for values, only for plotting in 2D)
temp_mFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D, 
                           asb_sp_w = temp_sp_biom,
                           ind_vect = c("fric", "fide"), 
                           scaling = TRUE, details_returned = TRUE)

# vertices of the global pool in 2D ----
pool_vert_nm_2D<-temp_mFD$details$pool_vert_nm
temp_mFD$details$asb_vert_nm$y2018


# compute only with temperate species (not for values, only for plotting in 2D)
temp_temperate_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D[temperate_sp,], 
                                            asb_sp_w = temp_sp_biom[, colnames(temp_sp_biom) %in% temperate_sp],
                                            ind_vect = c("fric"), 
                                            scaling = TRUE, details_returned = TRUE)

# list to store panels
temporal_ggplot<-list()

#loop on years

for (b in years_all) {
  
  # background
  ggplot_b <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC2"), 
                              color_bg = "grey95")
  
  # convex hull of global species pool
  
  ggplot_b <- pool.plot(ggplot_bg = ggplot_b,
                        sp_coord2D = sp_coord_2D,
                        vertices_nD = pool_vert_nm_2D,
                        plot_pool = FALSE,
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
  # adding title
  ggplot_b <- ggplot_b + labs(title = b)
  
  # plot convex hull of all species and of only temperate species
  ggplot_b <-fric.plot( ggplot_bg=ggplot_b,
                       asb_sp_coord2D=list(all=sp_coord_2D),
                       asb_vertices_nD= list(all=temp_mFD$details$asb_vert_nm[[b]]),
                       plot_sp = FALSE,
                       color_ch=c(all="black"),
                       fill_ch=c(all=NA),
                       alpha_ch=c(all=0.1))
  
  
  
  # plot species weights by thermal affinity using plot.fide showing mean value
  
  # species present in v with temperate or tropical affinity
  sp_b<-row.names(temp_mFD$details$asb_sp_faxes_coord[[b]])
  temp_sp_b<-tab_sp %>%
    filter(Species_name %in% sp_b & affinity=="temperate")
  temp_sp_b$Species_name
  trop_sp_b<-tab_sp %>%
    filter(Species_name %in% sp_b & affinity=="tropical")
  trop_sp_b$Species_name
  
  ggplot_b<- fide.plot(ggplot_bg=ggplot_b,
                       asb_sp_coord2D=list(temp=temp_mFD$details$asb_sp_faxes_coord[[b]][temp_sp_b$Species_name,],
                                           trop=temp_mFD$details$asb_sp_faxes_coord[[b]][trop_sp_b$Species_name,]),
                       asb_sp_relatw=list(temp=temp_mFD$details$asb_sp_relatw[b,temp_sp_b$Species_name],
                                          trop=temp_mFD$details$asb_sp_relatw[b,trop_sp_b$Species_name]),
                       asb_fide_coord2D=list(temp=temp_mFD$functional_diversity_indices[b,paste0("fide_", axes_nm)],
                                             trop=temp_mFD$functional_diversity_indices[b,paste0("fide_", axes_nm)]),
                       plot_sp = TRUE,
                       shape_sp = c(temp=21, trop=21),
                       color_sp =c(temp=col_temp, trop=col_trop),
                       fill_sp = c(temp=col_temp, trop=col_trop),
                       color_fide= c(temp="white", trop="white"),
                       fill_fide= c(temp=NA, trop=NA),
                       shape_fide = c(temp=1, trop=1),
                       size_fide = c(temp=0.0001, trop=0.0001),
                       linetype_segment = c(temp=1, trop=1),
                       color_segment = c(temp=NA, trop=NA),
                       width_segment = c(temp=0, trop=0)
  )
  # storing ggplot in list
  temporal_ggplot[[b]]<-ggplot_b
  
  
}# end of b


#### ONLY 3 YEARS (FOR MAIN PAPER)

# names of years to illustrate

years_names<-c("2002", "2009", "2018")

# list with names of species to display for each year

sp_code_temporal <- list( y2002 = c( "T. taeniatus", "S. lineolata", "A. strigatus", "B. waddi",
                                     "A. maculatus", "P. unifasciata", "C. picta"),
                        y2009= c( "T. taeniatus","A. strigatus", "B. waddi",
                                   "S. lineolata", "P. unifasciata",
                                  "C. picta"),
                       y2018 = c("T. taeniatus", "A. strigatus","B. waddi", "A. maculatus",
                                 "S. lineolata",
                                   "P. unifasciata", "L. russellii",
                                "S. interrupta", "P. xanthura", "S. fuscescens","C. picta"))


# list to store panels
temporal_ggplot<-list()

#loop on years

for (b in years_names) {
  
  # background
  ggplot_b <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC2"), 
                              color_bg = "grey95")
  
  # convex hull of global species pool and all species as small grey crosses
  ggplot_b <- pool.plot(ggplot_bg = ggplot_b,
                        sp_coord2D = sp_coord_2D,
                        vertices_nD = pool_vert_nm_2D,
                        plot_pool = TRUE,
                        shape_pool = 3, size_pool = 0.8, color_pool = "grey50", 
                        shape_vert = 3, size_vert = 1,   color_vert = "grey50",
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
  
  
  
  
  # plot convex hull of all species and of only temperate species
  ggplot_b<-fric.plot( ggplot_bg=ggplot_b,
                       asb_sp_coord2D=list(all=sp_coord_2D),
                       asb_vertices_nD= list(all=temp_mFD$details$asb_vert_nm[[b]]),
                       plot_sp = FALSE,
                       color_ch=c(all="black"),
                       fill_ch=c(all=NA),
                       alpha_ch=c(all=0.1))
  
  ggplot_b<- fide.plot(ggplot_bg=ggplot_b,
                       asb_sp_coord2D=list(temp=temp_mFD$details$asb_sp_faxes_coord[[b]][temp_sp_b$Species_name,],
                                           trop=temp_mFD$details$asb_sp_faxes_coord[[b]][trop_sp_b$Species_name,]),
                       asb_sp_relatw=list(temp=temp_mFD$details$asb_sp_relatw[b,temp_sp_b$Species_name],
                                          trop=temp_mFD$details$asb_sp_relatw[b,trop_sp_b$Species_name]),
                       asb_fide_coord2D=list(temp=temp_mFD$functional_diversity_indices[b,paste0("fide_", axes_nm)],
                                             trop=temp_mFD$functional_diversity_indices[b,paste0("fide_", axes_nm)]),
                       plot_sp = TRUE,
                       shape_sp = c(temp=21, trop=21),
                       color_sp =c(temp=col_temp, trop=col_trop),
                       fill_sp = c(temp=col_temp, trop=col_trop),
                       color_fide= c(temp="white", trop="white"),
                       fill_fide= c(temp=NA, trop=NA),
                       shape_fide = c(temp=1, trop=1),
                       size_fide = c(temp=0.0001, trop=0.0001),
                       linetype_segment = c(temp=1, trop=1),
                       color_segment = c(temp=NA, trop=NA),
                       width_segment = c(temp=0, trop=0)
  )
  
  #codes of species to display
  sp_codes_temporal <- tab_sp %>%
    filter(Species_code %in% sp_code_temporal[[b]])

  ggplot_b <- ggplot_b +
    ggrepel::geom_text_repel(data =  sp_codes_temporal,
                             mapping = ggplot2::aes_string(x = axes_nm[1],
                                                           y = axes_nm[2],
                                                           label = "Species_code"
                             ),
                             colour= sp_codes_temporal$col_affinity,
                             size = 4,
                             fontface = "plain",
                             segment.color = "black",
                             max.overlaps = Inf,
                             direction = "x",
                             box.padding = grid::unit(1, 'lines'),
                             force = 5,
                             arrow = grid::arrow(length = grid::unit(0.02, 'npc')),
                             show.legend = FALSE
  )

   # storing ggplot in list
  temporal_ggplot[[b]]<-ggplot_b
  
  
    ## THIS IS NOT WORKING - FIX IT 
  
}# end of b

## Now plotting PC1 vs PC3

# names of Functional axes to plot

axes_nm <- c("PC1", "PC3")

# coordinates of species along those axes
sp_coord_2D <- sp_3D_coord[,axes_nm]

# table with species coordinates, full species names (as in other objects) 
#  shortened species names (for displaying) and thermal affinity

tab_sp <- sp_coord_2D %>% 
  as.data.frame() %>% 
  rownames_to_column("Species_name") %>% 
  mutate(genus=sub("_.*", "", Species_name), sp=sub(".*_", "", Species_name)) %>% 
  mutate(Species_code=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  dplyr::select(-genus, -sp) %>%
  left_join(sp_thermal, by=c("Species_name"="Species")) %>%
  mutate(col_affinity = case_when(affinity=="temperate" ~ col_temp,
                                  affinity=="tropical" ~ col_trop) )

head(tab_sp)

# names of temperate species
temperate_sp <- tab_sp %>%
  filter(affinity=="temperate") %>%
  pull(Species_name)


# compute FRic and FIde for all habitats in 2D (not for values, only for plotting in 2D)
temp_mFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D, 
                            asb_sp_w = temp_sp_biom,
                            ind_vect = c("fric", "fide"), 
                            scaling = TRUE, details_returned = TRUE)

# vertices of the global pool in 2D ----
pool_vert_nm_2D<-temp_mFD$details$pool_vert_nm
temp_mFD$details$asb_vert_nm$y2018


# compute only with temperate species (not for values, only for plotting in 2D)
temp_temperate_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D[temperate_sp,], 
                                             asb_sp_w = temp_sp_biom[, colnames(temp_sp_biom) %in% temperate_sp],
                                             ind_vect = c("fric"), 
                                             scaling = TRUE, details_returned = TRUE)

# list to store panels
temporal_ggplot_2<-list()

#loop on years

for (c in years_names) {
  
  # background
  ggplot_c <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC3"), 
                              color_bg = "grey95")
  
  # convex hull of global species pool and all species as small grey crosses
  ggplot_c <- pool.plot(ggplot_bg = ggplot_c,
                        sp_coord2D = sp_coord_2D,
                        vertices_nD = pool_vert_nm_2D,
                        plot_pool = TRUE,
                        shape_pool = 3, size_pool = 0.8, color_pool = "grey50", 
                        shape_vert = 3, size_vert = 1,   color_vert = "grey50",
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
  
  
  # plot convex hull of all species and of only temperate species
  ggplot_c<-fric.plot( ggplot_bg=ggplot_c,
                       asb_sp_coord2D=list(all=sp_coord_2D),
                       asb_vertices_nD= list(all=temp_mFD$details$asb_vert_nm[[c]]),
                       plot_sp = FALSE,
                       color_ch=c(all="black"),
                       fill_ch=c(all=NA),
                       alpha_ch=c(all=0.1))
  
  
  
  # plot species weights by thermal affinity using plot.fide showing mean value
  
  # species present in v with temperate or tropical affinity
  sp_c<-row.names(temp_mFD$details$asb_sp_faxes_coord[[c]])
  temp_sp_c<-tab_sp %>%
    filter(Species_name %in% sp_c & affinity=="temperate")
  temp_sp_c$Species_name
  trop_sp_c<-tab_sp %>%
    filter(Species_name %in% sp_c & affinity=="tropical")
  trop_sp_c$Species_name
  
  ggplot_c<- fide.plot(ggplot_bg=ggplot_c,
                       asb_sp_coord2D=list(temp=temp_mFD$details$asb_sp_faxes_coord[[c]][temp_sp_c$Species_name,],
                                           trop=temp_mFD$details$asb_sp_faxes_coord[[c]][trop_sp_c$Species_name,]),
                       asb_sp_relatw=list(temp=temp_mFD$details$asb_sp_relatw[c,temp_sp_c$Species_name],
                                          trop=temp_mFD$details$asb_sp_relatw[c,trop_sp_c$Species_name]),
                       asb_fide_coord2D=list(temp=temp_mFD$functional_diversity_indices[c,paste0("fide_", axes_nm_2)],
                                             trop=temp_mFD$functional_diversity_indices[c,paste0("fide_", axes_nm_2)]),
                       plot_sp = TRUE,
                       shape_sp = c(temp=21, trop=21),
                       color_sp =c(temp=col_temp, trop=col_trop),
                       fill_sp = c(temp=col_temp, trop=col_trop),
                       color_fide= c(temp="white", trop="white"),
                       fill_fide= c(temp=NA, trop=NA),
                       shape_fide = c(temp=1, trop=1),
                       size_fide = c(temp=0.0001, trop=0.0001),
                       linetype_segment = c(temp=1, trop=1),
                       color_segment = c(temp=NA, trop=NA),
                       width_segment = c(temp=0, trop=0)
  )
  
  # storing ggplot in list
  temporal_ggplot_2[[c]]<-ggplot_c
  
  
}# end of c  ## NOT WORKING - FIX IT


##### SPATIAL #####

axes_nm <- c("PC1", "PC2")
  

# names of habitats to illustrate
hab_names<-c("Inshore", "Midshelf", "Offshore")

# coordinates of species along those axes
sp_coord_2D <- sp_3D_coord[,axes_nm]

# computing average biomass of species in the 3 habitats
hab_sp_biom <- spatial_sp_biom %>%
  as.data.frame %>%
  rownames_to_column("Code") %>%
  left_join( dplyr::select(spatial_metadata, Code, Habitat )  ) %>%
  dplyr::select(-Code) %>%
  group_by(Habitat) %>%
  summarise( across(.cols = everything(), .fns = mean, .names=NULL ) ) %>%
  column_to_rownames("Habitat") %>%
 as.matrix()

range( apply(hab_sp_biom, 1, max) ) # all species present in at least one of the 3 habitats


# compute FRic and FIde for all habitats in 2D (not for values, only for plotting in 2D)
hab_mFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D, 
                                     asb_sp_w = hab_sp_biom,
                                     ind_vect = c("fric", "fide"), 
                                     scaling = TRUE, details_returned = TRUE)

# vertices of the global pool in 2D ----
pool_vert_nm_2D<-hab_mFD$details$pool_vert_nm


# compute only with temperate species (not for values, only for plotting in 2D)
hab_temperate_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D[temperate_sp,], 
                                  asb_sp_w = hab_sp_biom[, colnames(hab_sp_biom) %in% temperate_sp],
                                  ind_vect = c("fric"), 
                                  scaling = TRUE, details_returned = TRUE)


hab_temperate_multidimFD$details$asb_vert_nm$Offshore

hab_temperate_multidimFD$details$asb_vert_nm$Inshore

# list with names of species to display for each habitat    

sp_code_spatial<- list( Inshore = c(  "A. viridis", "M. fuscus", "S. lineolata"),
 Midshelf= c("A. viridis", "S. lineolata", "P. microlepidotus"),
 Offshore = c(  "M. vanicolensis", "C. caerulaurea", "P. microlepidotus",
              "N. unicornis", "M. fuscus", "C. kleinii", "P. hepatus",  "A. viridis", "S. lineolata")
 )

# list to store panels
hab_ggplot<-list()

#loop on habitats
for (v in hab_names) {

  # background
  ggplot_v <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC2"), 
                              color_bg = "grey95")
  
  # adding title
  ggplot_v <- ggplot_v + labs(title = v )
  
  # convex hull of global species pool and all species as small grey crosses
  ggplot_v <- pool.plot(ggplot_bg = ggplot_v,
                        sp_coord2D = sp_coord_2D,
                        vertices_nD = pool_vert_nm_2D,
                        plot_pool = TRUE,
                        shape_pool = 3, size_pool = 0.8, color_pool = "grey50", 
                        shape_vert = 3, size_vert = 1,   color_vert = "grey50",
                        color_ch = "NA", fill_ch = "white", alpha_ch = 1)
    
  
  
  
  # plot convex hull of all species and of only temperate species
  ggplot_v<-fric.plot( ggplot_bg=ggplot_v,
                       asb_sp_coord2D=list(all=sp_coord_2D),
                       asb_vertices_nD= list(all=hab_mFD$details$asb_vert_nm[[v]]),
                       plot_sp = FALSE,
                       color_ch=c(all="black"),
                       fill_ch=c(all=NA),
                       alpha_ch=c(all=0.1))
  
  
  
  # plot species weights by thermal affinity using plot.fide showing mean value
  
  # species present in v with temperate or tropical affinity
  sp_v<-row.names(hab_mFD$details$asb_sp_faxes_coord[[v]])
  temp_sp_v<-tab_sp %>%
    filter(Species_name %in% sp_v & affinity=="temperate")
  trop_sp_v<-tab_sp %>%
    filter(Species_name %in% sp_v & affinity=="tropical")
  
  ggplot_v<- fide.plot(ggplot_bg=ggplot_v,
                       asb_sp_coord2D=list(temp=hab_mFD$details$asb_sp_faxes_coord[[v]][temp_sp_v$Species_name,],
                                           trop=hab_mFD$details$asb_sp_faxes_coord[[v]][trop_sp_v$Species_name,]),
                       asb_sp_relatw=list(temp=hab_mFD$details$asb_sp_relatw[v,temp_sp_v$Species_name],
                                          trop=hab_mFD$details$asb_sp_relatw[v,trop_sp_v$Species_name]),
                       asb_fide_coord2D=list(temp=hab_mFD$functional_diversity_indices[v,paste0("fide_", axes_nm)],
                                             trop=hab_mFD$functional_diversity_indices[v,paste0("fide_", axes_nm)]),
                       plot_sp = TRUE,
                       shape_sp = c(temp=21, trop=21),
                       color_sp =c(temp=col_temp, trop=col_trop),
                       fill_sp = c(temp=col_temp, trop=col_trop),
                       color_fide= c(temp="white", trop="white"),
                       fill_fide= c(temp=NA, trop=NA),
                       shape_fide = c(temp=1, trop=1),
                       size_fide = c(temp=0.0001, trop=0.0001),
                       linetype_segment = c(temp=1, trop=1),
                       color_segment = c(temp=NA, trop=NA),
                       width_segment = c(temp=0, trop=0)
                       )
  
# codes of species to display
  # sp_codes <- tab_sp %>%
  #   filter(Species_code %in% sp_code_spatial[[v]])
  # 
  # ggplot_v <- ggplot_v +
  #   ggrepel::geom_text_repel(data =  sp_codes,
  #         # @: replace Species_code with Species_name if full names provided in sp_code_spatial
  #                           mapping = ggplot2::aes_string(x = axes_nm[1],
  #                                                y = axes_nm[2],
  #                                                label = "Species_code"
  #                                                ),
  #         colour= sp_codes$col_affinity,
  #                            size = 4,
  #                            fontface = "plain",
  #                            segment.color = "black",
  #                            max.overlaps = Inf,
  #                            box.padding = grid::unit(2, 'lines'),
  #                            force = 5,
  #                            direction = "x",
  #                            arrow = grid::arrow(length = grid::unit(0.02, 'npc')),
  #                            show.legend = FALSE
  #                            )

  # storing ggplot in list
  hab_ggplot[[v]]<-ggplot_v
  
  
 }# end of v


# arranging panels

figureS5a <- (temporal_ggplot[[1]] | temporal_ggplot[[2]] | temporal_ggplot[[3]]|temporal_ggplot[[4]])/
              (temporal_ggplot[[5]] | temporal_ggplot[[6]] | temporal_ggplot[[7]]| temporal_ggplot[[8]])/
              (temporal_ggplot[[9]] | temporal_ggplot[[10]]| temporal_ggplot[[11]] | temporal_ggplot[[12]])/ 
                          (temporal_ggplot[[13]] | temporal_ggplot[[14]])

figure3_years <- temporal_ggplot[[1]] + temporal_ggplot[[2]] + temporal_ggplot[[3]]

figureS5b <- temporal_ggplot_2[[1]] + temporal_ggplot_2[[2]] + temporal_ggplot_2[[3]]

figure3_habitat <- (hab_ggplot[[1]] + hab_ggplot[[2]] + hab_ggplot[[3]])

figureS6 <- (hab_ggplot[[1]] + hab_ggplot[[2]] + hab_ggplot[[3]])

# saving as jepg ----

ggsave(figureS5a, file=here::here("outputs", "FigureS5a.jpeg"),
       height = 30, width = 30, unit = "cm" )

ggsave(figureS5b, file=here::here("outputs", "FigureS5b.jpeg"),
       height = 20, width = 30, unit = "cm" )

ggsave(figure3_years, file=here::here("outputs", "Figure3_years.jpeg"),
       height = 20, width = 30, unit = "cm" )

ggsave(figure3_habitat, file=here::here("outputs", "Figure3_habitat.jpeg"),
       height = 20, width = 30, unit = "cm" )

ggsave(figure3_habitat, file=here::here("outputs", "Figure3_habitat_empty.jpeg"),
       height = 20, width = 30, unit = "cm" )

ggsave(figureS6, file=here::here("outputs", "FigureS6.jpeg"),
       height = 20, width = 30, unit = "cm" )





