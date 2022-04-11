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
load(here::here("outputs/", "temporal_fd_kelp.RData") )

load(here::here("data", "nokelp_metadata.RData") )
load(here::here("data", "nokelp_sp_occ.RData") )
load(here::here("outputs/", "temporal_fd_nokelp.RData") )


## settings ####

# coordinates of all species ----
pool_coord<-spatial_fd$details$sp_faxes_coord 

# vertices of all fe in 4D ----
pool_vert_nm<-spatial_fd$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(pool_coord[,1:3])
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1


names(spatial_fd$details$asb_vert_nm)

## temporal kelp ####


# color code for the 4 years
years_colors <- c(y2002="green3", y2008="blue3", y2013="grey50", y2018="red3")

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

# FIde for kelp sites
kelp_fide_years<- kelp_metadata %>%
  mutate(Year=paste0("y",Year)) %>%
  select(Year) %>%
  bind_cols( select(temporal_fd_kelp$functional_diversity_indices, fide_PC1, fide_PC2 , fide_PC3) ) %>%
  filter( Year %in% names(years_colors) )


## plotting  ####

# list to store ggplot
ggplot_temporal_kelp<-list()

# pairs of axes
pairs_axes<-list( c(1,2), c(1,3) )

for ( z in 1:length(pairs_axes) ) {
  
  # names of axes   
  xy<-pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z<-background.plot(range_faxes=range_axes,
                            faxes_nm=paste0("PC", xy), 
                            color_bg="grey95")
  
  # convex hull of species pool
  ggplot_z<-pool.plot(ggplot_bg=ggplot_z,
                      sp_coord2D=pool_coord[,xy],
                      vertices_nD=pool_vert_nm,
                      plot_pool=FALSE,
                      color_ch=NA, fill_ch="white", alpha_ch=1)
  
  # loop on years
  for (h in row.names(kelp_years_sp_occ) ) {
    
    # color given habitat
    col_h<-as.character(years_colors[h])
    
    # species present in v
    sp_h<-names(which(kelp_years_multidimFD$details$asb_sp_occ[h,]==1))
    
    # plot convex hull of assemblage but not species
    ggplot_z<-fric.plot( ggplot_bg=ggplot_z, 
                         asb_sp_coord2D=list(vv=pool_coord[sp_h,xy]),
                         asb_vertices_nD=list(vv=kelp_years_multidimFD$details$asb_vert_nm[[h]]),
                         plot_sp = FALSE,
                         color_ch=c(vv=col_h),
                         fill_ch=c(vv=col_h),
                         alpha_ch=c(vv=0.1)
    )
  }# end of h
  
  # average position of species (in the 4 sites) for each year
  ggplot_z <- ggplot_z +
    geom_point(data=kelp_fide_years, aes_string( x=paste0("fide_PC",xy[1]), y=paste0("fide_PC",xy[2]),
                                          color = "Year", shape = "Year"), size=0.9, show.legend = FALSE) +
    scale_colour_manual( values = years_colors )
  
  
  # legend and title
  if (z==1) {
    nlabels_y <- length(years_colors)
    ggplot_z <- ggplot_z + 
      geom_text(aes(x = rep(range_axes[1]*0.95, nlabels_y ),
                    y = range_axes[2]*(1.05-0.15*(1:nlabels_y)),
                    label = substr(names(years_colors),2,5) ),
                    color = years_colors,
                hjust = "left",
                show.legend = FALSE) +
      ggtitle("Temporal - Kelp")
  }
  
  # ggplot stored in list
  ggplot_temporal_kelp[[z]]<-ggplot_z
  
}# end of z



## temporal no kelp ####

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

# FIde for nokelp sites
nokelp_fide_years<- nokelp_metadata %>%
  mutate(Year=paste0("y",Year)) %>%
  select(Year) %>%
  bind_cols( select(temporal_fd_nokelp$functional_diversity_indices, fide_PC1, fide_PC2 , fide_PC3) ) %>%
  filter( Year %in% names(years_colors) )



## plotting  ####

# list to store ggplot
ggplot_temporal_nokelp<-list()

# pairs of axes
pairs_axes<-list( c(1,2), c(1,3) )

for ( z in 1:length(pairs_axes) ) {
  
  # names of axes   
  xy<-pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z<-background.plot(range_faxes=range_axes,
                            faxes_nm=paste0("PC", xy), 
                            color_bg="grey95")
  
  # convex hull of species pool
  ggplot_z<-pool.plot(ggplot_bg=ggplot_z,
                      sp_coord2D=pool_coord[,xy],
                      vertices_nD=pool_vert_nm,
                      plot_pool=FALSE,
                      color_ch=NA, fill_ch="white", alpha_ch=1)
  
  # loop on years
  for (h in row.names(nokelp_years_sp_occ) ) {
    
    # color given habitat
    col_h<-as.character(years_colors[h])
    
    # species present in v
    sp_h<-names(which(nokelp_years_multidimFD$details$asb_sp_occ[h,]==1))
    
    # plot convex hull of assemblage but not species
    ggplot_z<-fric.plot( ggplot_bg=ggplot_z, 
                         asb_sp_coord2D=list(vv=pool_coord[sp_h,xy]),
                         asb_vertices_nD=list(vv=nokelp_years_multidimFD$details$asb_vert_nm[[h]]),
                         plot_sp = FALSE,
                         color_ch=c(vv=col_h),
                         fill_ch=c(vv=col_h),
                         alpha_ch=c(vv=0.1)
    )
    
  }# end of v

  # average position of species (in the 4 sites) for each year
  ggplot_z <- ggplot_z +
    geom_point(data=nokelp_fide_years, aes_string( x=paste0("fide_PC",xy[1]), y=paste0("fide_PC",xy[2]),
                                                 color = "Year", shape = "Year"), size=0.9, show.legend = FALSE) +
    scale_colour_manual( values = years_colors )

  # title  
  if (z==1) {
    ggplot_z <- ggplot_z + 
      ggtitle("Temporal - No Kelp")
  }
  
  # ggplot stored in list
  ggplot_temporal_nokelp[[z]]<-ggplot_z
  
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
fide_hab<- spatial_metadata %>%
  select(Habitat) %>%
  bind_cols( select(spatial_fd$functional_diversity_indices, fide_PC1, fide_PC2 , fide_PC3) )


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
                         color_ch = c(vv=col_b),
                         fill_ch= c(vv=col_b),
                         alpha_ch = c(vv=0.1))
    
  }# end of b
  
  # average position of species (in the 3 sites) of each habitat
  ggplot_a <- ggplot_a +
    geom_point(data=fide_hab, aes_string( x=paste0("fide_PC",xy[1]), y=paste0("fide_PC",xy[2]),
                               color = "Habitat", shape = "Habitat"), size=0.9, show.legend = FALSE) +
    scale_colour_manual( values = hab_colors )
  
  
  # legend
  if (a==1) {
    nlabels_b <- length(hab_colors)
    ggplot_a <- ggplot_a + 
      geom_text(aes(x = rep(range_axes[1]*0.95, nlabels_b ),
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
figure3 <- ( ggplot_temporal_nokelp[[1]] +  ggplot_temporal_kelp[[1]] + ggplot_spatial[[1]] ) / ( 
  
  ggplot_temporal_nokelp[[2]] +  ggplot_temporal_kelp[[2]] + ggplot_spatial[[2]])


ggsave(figure3, file=here::here("outputs/", "Figure3.png"),
       height = 16, width = 24, unit = "cm" )
