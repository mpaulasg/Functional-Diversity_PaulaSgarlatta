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

# names of habitats to illustrate
years_names<-c("y2002", "y2009", "y2018")

# computing average biomass of species in the 3 habitats
temp_sp_biom <- kelp_sp_maxN %>%
  as.data.frame %>%
  rownames_to_column("Code") %>%
  left_join( dplyr::select(kelp_metadata, Code, Year )  ) %>%
  dplyr::select(-Code) %>%
  group_by(Year) %>%
  summarise( across(.cols = everything(), .fns = mean, .names=NULL ) ) 

#Adding a letter becase is confusing with a number

temp_sp_biom$Year<- paste0("y", temp_sp_biom$Year) 

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


# compute only with temperate species (not for values, only for plotting in 2D)
temp_temperate_multidimFD<-alpha.fd.multidim(sp_faxes_coord = sp_coord_2D[temperate_sp,], 
                                            asb_sp_w = temp_sp_biom[, colnames(temp_sp_biom) %in% temperate_sp],
                                            ind_vect = c("fric"), 
                                            scaling = TRUE, details_returned = TRUE)




# list with names of species to display for each habitat    # @ Paula: would be better to provide full names and to pick codes from the table "tab_sp"

sp_code_temporal <- list( y2002 = c( "T. taeniatus", "E. fasciatus","G. thyrsoidea","T. lutescens", 
                 "B. waddi" ,"O. halei", "A. maculatus", "S. lalandi" ),
                        y2009= c( "E. fasciatus", "T. lunare",
                                  "T. lutescens", "B.waddi",
                                  "K. sydneyanus",
                                  "P. unifasciata",
                                  "P. dentex", 
                                  "T. taeniatus"),
                       y2018 = c( "A. maculatus", "B. waddi",
                                      "O. maculatus","P. stricticeps",
                                      "L. russellii" , "P. xanthura",
                                      "S. interrupta",
                                      "T. lunare"))



# list to store panels
temporal_ggplot<-list()

#loop on habitats
for (b in years_names) {
  
  # background
  ggplot_b <- background.plot(range_faxes = range_axes,
                              faxes_nm = c("PC1", "PC2"), 
                              color_bg = "grey95")
  
  # adding title
  ggplot_b <- ggplot_b + labs(title = b )
  
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
  
  
  
  # plot species weights by thermal affinity using plot.fide showing mean value
  
  # species present in v with temperate or tropical affinity
  sp_b<-row.names(temp_mFD$details$asb_sp_faxes_coord[[b]])
  temp_sp_b<-tab_sp %>%
    filter(Species_name %in% sp_b & affinity=="temperate")
  trop_sp_b<-tab_sp %>%
    filter(Species_name %in% sp_b & affinity=="tropical")
  
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
  
  # codes of species to display
  sp_codes_temporal <- tab_sp %>%
    filter(Species_code %in% sp_code_temporal[[b]])
  
  ggplot_b <- ggplot_b +
    ggrepel::geom_text_repel(data =  sp_codes_temporal, 
                             # @: replace Species_code with Species_name if full names provided in sp_code_spatial
                             mapping = ggplot2::aes_string(x = axes_nm[1],
                                                           y = axes_nm[2],
                                                           label = "Species_code"
                             ),
                             colour= sp_codes_temporal$col_affinity,
                             size = 4,
                             fontface = "plain",
                             segment.color = "black",
                             max.overlaps = Inf,
                             box.padding = grid::unit(2, 'lines'),
                             force = 5,
                             arrow = grid::arrow(length = grid::unit(0.02, 'npc')),
                             show.legend = FALSE
    )
  
   # storing ggplot in list
  temporal_ggplot[[b]]<-ggplot_b
  
  
}# end of b

##### SPATIAL #####
  

# names of habitats to illustrate
hab_names<-c("Inshore", "Midshelf", "Offshore")

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




# list with names of species to display for each habitat    # @ Paula: would be better to provide full names and to pick codes from the table "tab_sp"

sp_code_spatial<- list( Inshore = c( "A. maculatus", "A. viridis" , "G. elevata", "P. affinis", "P. oligolepis"),
 Midshelf= c("A. chinensis",  "H. margaritaceus",  "S. fuscescens",  "T. lunare", 
             "T. lutescens", "A. viridis", "K. sydneyanus", "P. unifasciata"),
 Offshore = c("C. caerulaurea","H. margaritaceus",  "L. dimidiatus", "M. vanicolensis",  
              "N. unicornis", "P. flavomaculatus", "P. hepatus",  "A. viridis")
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
  sp_codes <- tab_sp %>%
    filter(Species_code %in% sp_code_spatial[[v]])
  
  ggplot_v <- ggplot_v +
    ggrepel::geom_text_repel(data =  sp_codes, 
          # @: replace Species_code with Species_name if full names provided in sp_code_spatial
                            mapping = ggplot2::aes_string(x = axes_nm[1],
                                                 y = axes_nm[2],
                                                 label = "Species_code"
                                                 ),
          colour= sp_codes$col_affinity,
                             size = 4,
                             fontface = "plain",
                             segment.color = "black",
                             max.overlaps = Inf,
                             box.padding = grid::unit(2, 'lines'),
                             force = 5,
                             arrow = grid::arrow(length = grid::unit(0.02, 'npc')),
                             show.legend = FALSE
                             )

  # legend and title
  if (v==3) {
    ggplot_v <- ggplot_v  +
      theme(plot.title = element_text(size = 20, color = "#FDE725FF")) +
      add_fishape(family = "Acanthuridae",
                  option = "Naso_unicornis",
                  xmin =  0.25 ,xmax = 0.35, ymin = -0.20, ymax = -0.35,
                  fill = "black",
                  alpha = 1)
      } # Check this one, don't know why it's not working

  # storing ggplot in list
  hab_ggplot[[v]]<-ggplot_v
  
  
 }# end of v


# arranging panels
figure4 <-  (temporal_ggplot[[1]] + temporal_ggplot[[2]] + temporal_ggplot[[3]])/
  (hab_ggplot[[1]] + hab_ggplot[[2]] + hab_ggplot[[3]])



# saving as jepg ----
ggsave(figure4, file=here::here("outputs", "Figure4_v2.jpeg"),
       height = 30, width = 30, unit = "cm" )


