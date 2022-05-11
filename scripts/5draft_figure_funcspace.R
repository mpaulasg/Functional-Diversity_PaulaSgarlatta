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


# loading thermal affinity data
sp_thermal <- read.csv(here::here("data", "raw_data", "thermal_all.csv")) %>% 
  mutate(affinity= if_else(thermal>"23", "tropical", "temperate")) %>%   
  select(-thermal)

  
## preparing ###

 # names of Functional axes to plot
 axes_nm <- c("PC1", "PC2")
 
 # coordinates of species along those axes
 sp_coord_2D <- sp_3D_coord[,axes_nm]
 
 # color code for thermal affinity
col_temp <- "blue3"
col_trop <- "orange3"
 
# table with species coordinates, full species names (as in other objects) 
#  shortened species names (for displaying) and thermal affinity

tab_sp <- sp_coord_2D %>% 
  as.data.frame() %>% 
  rownames_to_column("Species_name") %>% 
  mutate(genus=sub("_.*", "", Species_name), sp=sub(".*_", "", Species_name)) %>% 
  mutate(Species_code=paste(substr(genus, 1, 1), sp, sep = ". ")) %>% 
  select(-genus, -sp) %>%
  left_join(sp_thermal, by=c("Species_name"="Species")) %>%
  mutate(col_affinity = case_when(affinity=="temperate" ~ col_temp,
                                  affinity=="tropical" ~ col_trop) )

head(tab_sp)

dim(sp_coord_2D)
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


## habitats ####

# names of habitats to illustrate
hab_names<-c("Inshore", "Midshelf", "Offshore")

# computing average biomass of species in the 3 habitats
hab_sp_biom <- spatial_sp_biom %>%
  as.data.frame %>%
  rownames_to_column("Code") %>%
  left_join( select(spatial_metadata, Code, Habitat )  ) %>%
  select(-Code) %>%
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
                       asb_sp_coord2D=list(all=sp_coord_2D, temp=sp_coord_2D[temperate_sp,]),
                       
                       asb_vertices_nD=list(all=hab_mFD$details$asb_vert_nm[[v]],
                                            temp=hab_temperate_multidimFD$details$asb_vert_nm[[v]]),
                       plot_sp = FALSE,
                       color_ch=c(all="black", temp=col_temp),
                       fill_ch=c(all=NA, temp=NA),
                       alpha_ch=c(all=0.1, temp=0.1)
                       )
  
  
  
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


  # storing ggplot in list
  hab_ggplot[[v]]<-ggplot_v
  
  
 }# end of v


# arranging panels
figure4 <- # (ggplot_2002[[1]] + ggplot_2009[[1]] + ggplot_2018[[1]])/
  (hab_ggplot[[1]] + hab_ggplot[[2]] + hab_ggplot[[3]])



# saving as jepg ----
ggsave(figure4, file=here::here("outputs", "Figure4_seb.jpeg"),
       height = 25, width = 45, unit = "cm" )


