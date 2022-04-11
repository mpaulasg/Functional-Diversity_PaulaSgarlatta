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



load(here::here("data", "nokelp_metadata.RData") )
load(here::here("data", "nokelp_sp_occ.RData") )
load(here::here("outputs/", "temporal_fd_nokelp.RData") )



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



########### For statistics

## temporal survey data ####

# loading raw data from csv----

# metadata and data from sites (and replicates) without kelp
nokelp_metadata_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_no_kelp.csv") )
head(nokelp_metadata_all)
unique(nokelp_metadata_all$Site) # 4 sites
unique(nokelp_metadata_all$Year) # 14 years

nokelp_sp_maxN_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_no_kelp.csv") ) %>%
  column_to_rownames("Site") %>%
  as.matrix()
head(nokelp_sp_maxN_all)

# metadata and data from sites that use to have kelp and lost it:
kelp_metadata_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv") )
head(kelp_metadata_all)
unique(kelp_metadata_all$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv") ) %>%
  column_to_rownames("Site") %>%
  as.matrix()
head(kelp_sp_maxN_all)


# summary of surveys and occurrences data ----
kelp_summary_all <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN_all)
nokelp_summary_all <- mFD::asb.sp.summary(asb_sp_w = nokelp_sp_maxN_all)

# retrieve occurrence matrix:
kelp_sp_occ_all <- kelp_summary_all$asb_sp_occ
nokelp_sp_occ_all <- nokelp_summary_all$asb_sp_occ


# dimensions
dim(kelp_sp_occ_all) # 204 assemblages * 101 species
dim(nokelp_sp_occ_all) # 167 assemblages * 106 species

############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata_all <- read.csv(here::here("from_paula",  "SpatialUVC_metadata_transect.csv"))
head(spatial_metadata_all)

spatial_sp_biom_all <- read.csv(here::here("from_paula", "SpatialUVC_species_biomass_transect.csv")) %>%
  column_to_rownames("Site") %>% 
  as.matrix()

dim(spatial_sp_biom_all) # 36 assemblages * 51 species

# summary of surveys and occurrences data ----
spatial_summary_all <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom_all)

# retrieve occurrence matrix:
spatial_sp_occ_all <- spatial_summary_all$asb_sp_occ

## saving dataframes #####
save(kelp_metadata_all, file=here::here("data", "kelp_metadata_all.RData") )
save(kelp_sp_occ_all, file=here::here("data", "kelp_sp_occ_all.RData") )
save(kelp_summary_all, file=here::here("data", "kelp_summary_all.RData") )

save(nokelp_metadata_all, file=here::here("data", "nokelp_metadata_all.RData") )
save(nokelp_sp_occ_all, file=here::here("data", "nokelp_sp_occ_all.RData") )
save(nokelp_summary_all, file=here::here("data", "nokelp_summary_all.RData") )

save(spatial_metadata_all, file=here::here("data", "spatial_metadata_all.RData") )
save(spatial_sp_occ_all, file=here::here("data", "spatial_sp_occ_all.RData") )
save(spatial_summary_all, file=here::here("data", "spatial_summary_all.RData") )

############################################### end of script ########################################################


# metadata and data from sites without kelp
nokelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_metadata.csv") )
head(nokelp_metadata)
unique(nokelp_metadata$Site) # 4 sites
unique(nokelp_metadata$Year) # 14 years

nokelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(nokelp_sp_maxN)


nokelp <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv") )

nokelp_summary <- mFD::asb.sp.summary(asb_sp_w = nokelp_sp_maxN)

nokelp_sp_occ <- nokelp_summary$asb_sp_occ

dim(nokelp_sp_occ) # 56 assemblages * 106 species

nokelp_sp <- colnames(nokelp_sp_occ) 
length(nokelp_sp) # 106 sp

sum(kelp_sp %in% nokelp_sp)# 83 species shared

temporal_sp <- unique(c(kelp_sp, nokelp_sp))

length(temporal_sp) # 124 unique species ##[PS] This is not giving unique species, not sure what is exactly giving?

temporal_sp_kelp <- as.data.frame(setdiff(kelp_sp, nokelp_sp)) #[PS] Sp only in kelp

colnames(temporal_sp_kelp) <- "Species"

temporal_sp_nokelp <- as.data.frame(setdiff(nokelp_sp, kelp_sp)) #[PS] Sp only in no kelp

colnames(temporal_sp_nokelp) <- "Species"

temporal_sp_unique <- bind_rows(temporal_sp_kelp, temporal_sp_nokelp)

temporal_sp_unique <- temporal_sp_unique[order(temporal_sp_unique$Species),] 

length(temporal_sp_unique) #41 unique species - THIS IS CORRECT

save(nokelp_metadata, file=here::here("data", "nokelp_metadata.RData") )
save(nokelp_sp_occ, file=here::here("data", "nokelp_sp_occ.RData") )
save(nokelp_summary, file=here::here("data", "nokelp_summary.RData") )

# from sites that never had kelp
nokelp_traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_no_kelp.csv"),
                          header = T)
head(nokelp_traits)
names(nokelp_traits) <- c("Species","Size", "Agg", "Position" , "Diet")
nrow(nokelp_traits) # 106 sp

write.csv(nokelp_traits, file=here::here("data", "raw_data", "nokelp_traits.csv"), 
          row.names = FALSE )


## temporal survey data from sites without kelp ####

# loading raw data from csv----
nokelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_no_kelp.csv"))
head(nokelp_metadata)  

nokelp_replicates_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_no_kelp.csv"))
head(nokelp_replicates_maxN)  

# merging metadata and data 
nokelp <- left_join(nokelp_metadata, nokelp_replicates_maxN, by=c("Code"="Site") )
head(nokelp)

# computing total biomass per site for each species ----
nokelp_sites_maxN <- nokelp %>% 
  pivot_longer(cols= contains("_"), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Year, species) %>%
  summarize(total=sum(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Year,sep="_"), .before="Site")

head(nokelp_sites_maxN)

nokelp_replicates_maxN %>% select( contains("_") ) %>% sum()
nokelp_sites_maxN %>% select( contains("_") ) %>% sum()
# => same total biomass

# splitting into 2 tables ----
TemporalBRUV_nokelp_metadata <- nokelp_sites_maxN %>% 
  select(Code, Site, Year)

TemporalBRUV_nokelp_species <- nokelp_sites_maxN %>% 
  select( "Code", contains("_") )

## saving as csv ----
write.csv(TemporalBRUV_nokelp_metadata, file=here::here("data", "raw_data", "TemporalBRUV_nokelp_metadata.csv"), row.names = FALSE )
write.csv(TemporalBRUV_nokelp_species, file=here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv"), row.names = FALSE)


nokelp_thermal <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_only_thermal_no_kelp.csv") )

load(here::here("data", "nokelp_sp_occ.RData") )


## computing taxonomic and functional diversity for no kelp sites ####

# number of species, functional richness, dispersion and identity (along 3 axes)

temporal_fd_nokelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = nokelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_nokelp <- temporal_fd_nokelp$functional_diversity_indices

save(temporal_fd_nokelp, file=here::here("outputs/", "temporal_fd_nokelp.RData") )
save(temporal_alpha_nokelp, file=here::here("outputs/", "temporal_alpha_nokelp.RData") )

## For statistics

# loading data
load(here::here("data", "kelp_sp_occ_all.RData") )
load(here::here("data", "nokelp_sp_occ_all.RData") )
load(here::here("outputs", "sp_3D_coord.RData") ) 
load(here::here("data", "spatial_sp_occ_all.RData") )
load(here::here("outputs", "sp_3D_coord.RData") )

## computing taxonomic and functional diversity for no kelp sites ####

# number of species, functional richness, dispersion and identity (along 3 axes)

temporal_fd_nokelp_all <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = nokelp_sp_occ_all,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_nokelp_all <- temporal_fd_nokelp_all$functional_diversity_indices

## computing taxonomic and functional diversity for kelp sites ####

temporal_fd_kelp_all <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = kelp_sp_occ_all,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp_all <- temporal_fd_kelp_all$functional_diversity_indices

## computing taxonomic and functional diversity for spatial data ####

spatial_fd_all <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = spatial_sp_occ_all,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

spatial_alpha_all <- spatial_fd_all$functional_diversity_indices

# saving ####

# trait values and trait coding dataframes ----
save(temporal_fd_nokelp_all, file=here::here("outputs/", "temporal_fd_nokelp_all.RData") )
save(temporal_alpha_nokelp_all, file=here::here("outputs/", "temporal_alpha_nokelp_all.RData") )

save(temporal_fd_kelp_all, file=here::here("outputs/", "temporal_fd_kelp_all.RData") )
save(temporal_alpha_kelp_all, file=here::here("outputs/", "temporal_alpha_kelp_all.RData") )

save(spatial_fd_all, file=here::here("outputs/", "spatial_fd_all.RData") )
save(spatial_alpha_all, file=here::here("outputs/", "spatial_alpha_all.RData") )

load(here::here("data", "nokelp_metadata.RData") )
load(here::here("outputs/", "temporal_alpha_nokelp.RData") )


temporal_nokelp<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(temporal_alpha_nokelp, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3)

# merging values for the 2 habitat types
temporal_all <- bind_rows(temporal_kelp, temporal_nokelp )
head(temporal_all)

## end of no kelp ####


load(here::here("outputs", "temporal_alpha_nokelp.RData") )

nokelp_stats <- temporal_alpha_nokelp %>%
  rownames_to_column(var = "Site1") %>% 
  mutate(Site=sub("_.*", "", Site1),Year=sub(".*_", "", Site1), .before="Site1")%>% 
  select(-Site1)


#######No kelp - Fide on the 3 PC axes

fide3_nokelp <- temporal_alpha_nokelp %>% 
  select(fide_PC1,     fide_PC2, fide_PC3) 

shift3D_nokelp_eucl <- dist(fide3_nokelp, method = "euclidean")

shift3D_nokelp_eucl_1 <- dist.to.df( list(shift3D=shift3D_nokelp_eucl) ) %>% 
  as.data.frame() %>% 
  mutate(Site1=sub("_.*", "", x1), Site2=sub("_.*", "", x2),
         Year1=sub(".*_", "", x1),  Year2=sub(".*_", "", x2))


df_shift3D_nokelp_toplot <- shift3D_nokelp_eucl_1 %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>%
  select(Year, shift3D)  %>%
  group_by(Year) %>%
  summarise( n = n(),
             shift3D_mean = mean(shift3D),
             shift3D_sd = sd(shift3D)
  ) %>%
  mutate(shift3D_se = shift3D_sd/sqrt(n))  %>%
  mutate( Habitat = "No_Kelp", .before="Year" )


shift3D_nokelp_stats <- shift3D_nokelp_eucl_1 %>% 
  filter(Year1==2002 | Year2==2002) %>%
  filter(Site1 == Site2) %>%
  mutate(Year = pmax ( Year1, Year2 ) ) %>% 
  select(shift3D, Site1, Site2, Year)



