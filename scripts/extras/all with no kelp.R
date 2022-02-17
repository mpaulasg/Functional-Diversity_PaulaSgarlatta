################################################################################
##
## Script for merging all traits values into a single csv file
## 
##  Sébastien Villéger
##
## NB: NOT TO BE ON FINAL GITHUB
################################################################################


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(dplyr)
library(arsenal)

# load traits data ----

# traits of species from sites surveyed through years using BRUVs 
# from sites with kelp
kelp_traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_kelp.csv"),
                        header = T)
head(kelp_traits)
names(kelp_traits) <- c("Species","Size", "Agg", "Position" , "Diet")
nrow(kelp_traits) # 101 sp

# from sites that never had kelp
nokelp_traits <- read.csv(here::here("from_paula", "TemporalBRUV_species_traits_no_kelp.csv"),
                          header = T)
head(nokelp_traits)
names(nokelp_traits) <- c("Species","Size", "Agg", "Position" , "Diet")
nrow(nokelp_traits) # 106 sp

# traits of species from UVC surveys
spatial_traits <- read.csv(here::here("from_paula", "SpatialUVC_species_traits.csv"), 
                           header = T)
head(spatial_traits)
names(spatial_traits) <- c("Species","Size", "Agg", "Position" , "Diet")

nrow(spatial_traits) # 53 sp

#All sp with K values

k_values <- read.csv(here::here("from_paula", "fish_K_values.csv"),
                     header = T)
names(k_values) <- c("Species", "sstmean", "MaxSizeTL", "Diet", "Position", "Method", "Kmax", "Kmax_lowq", "Kmax_uppq")


#Add K values to each dataframe

kelp_traits <- inner_join(kelp_traits, k_values[ , c("Species", "Kmax")], by = c("Species"), all.x=TRUE)

nokelp_traits <- inner_join(nokelp_traits, k_values[ , c("Species", "Kmax")], by = c("Species"), all.x=TRUE)   

spatial_traits <- inner_join(spatial_traits, k_values[ , c("Species", "Kmax")], by = c("Species"), all.x=TRUE) 

# merging in a single data.frame and keeping only species present in surveys ----
fish_traits <- bind_rows( kelp_traits, nokelp_traits, spatial_traits) %>%
  distinct(Species, .keep_all = TRUE ) %>%
  as.data.frame()
head(fish_traits)

nrow(fish_traits) # 139 species

fish_traits <- as.data.frame(fish_traits[order(fish_traits$Species),])

# saving as csv file

write.csv(fish_traits, file=here::here("scripts", "extras", "fish_traits_extra.csv"), 
          row.names = FALSE )

write.csv(kelp_traits, file=here::here("scripts", "extras", "kelp_traits_extra.csv"), 
          row.names = FALSE )

write.csv(nokelp_traits, file=here::here("scripts", "extras", "nokelp_traits_extra.csv"), 
          row.names = FALSE )

write.csv(spatial_traits, file=here::here("scripts", "extras", "spatial_traits_extra.csv"), 
          row.names = FALSE )


########################################## end of code ###############################################################


###############################################################################
##
## Script for merging biomass from temporal BRUV surveys 
## 
##  Sébastien Villéger
##
## NB: NOT TO BE ON FINAL GITHUB
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

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
write.csv(TemporalBRUV_nokelp_metadata, file=here::here("scripts", "extras", "TemporalBRUV_nokelp_metadata_extra.csv"), row.names = FALSE )
write.csv(TemporalBRUV_nokelp_species, file=here::here("scripts", "extras", "TemporalBRUV_nokelp_species_extra.csv"), row.names = FALSE)


## end of no kelp ####
################################################################################

## temporal survey data from sites with kelp ####
# loading raw data from csv----
kelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv"))
head(kelp_metadata)  

kelp_replicates_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv"))
head(kelp_replicates_maxN)  

# merging metadata and data 
kelp <- left_join(kelp_metadata, kelp_replicates_maxN, by="Code" )
head(kelp)

# computing total biomass per site for each species ----
kelp_sites_maxN <- kelp %>% 
  pivot_longer(cols= contains("_"), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Year, species) %>%
  summarize(total=sum(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Year,sep="_"), .before="Site")

head(kelp_sites_maxN)

kelp_replicates_maxN %>% select( contains("_") ) %>% sum()
kelp_sites_maxN %>% select( contains("_") ) %>% sum()
# => same total biomass

# splitting into 2 tables ----
TemporalBRUV_kelp_metadata <- kelp_sites_maxN %>% 
  select(Code, Site, Year)

TemporalBRUV_kelp_species <- kelp_sites_maxN %>% 
  select( "Code", contains("_") )

## saving as csv ----
write.csv(TemporalBRUV_kelp_metadata, file=here::here("scripts", "extras", "TemporalBRUV_kelp_metadata_extra.csv"), row.names = FALSE )
write.csv(TemporalBRUV_kelp_species, file=here::here("scripts", "extras", "TemporalBRUV_kelp_species_extra.csv"), row.names= FALSE)


################################################################################
##
## Script for preparing fish occurrence datasets
## 
## Code by Camille Magneville, Sébastien villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(dplyr)

## temporal survey data ####

# loading raw data from csv----

# metadata and data from sites without kelp
nokelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_metadata.csv") )
head(nokelp_metadata)
unique(nokelp_metadata$Site) # 4 sites
unique(nokelp_metadata$Year) # 14 years

nokelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(nokelp_sp_maxN)

# metadata and data from sites that use to have kelp and lost it:
kelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_metadata.csv") )
head(kelp_metadata)
unique(kelp_metadata$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(kelp_sp_maxN)

### Joining both datasets

kelp <- read.csv(here::here("scripts", "extras", "TemporalBRUV_kelp_species_extra.csv") )
nokelp <- read.csv(here::here("scripts", "extras", "TemporalBRUV_nokelp_species_extra.csv") )


# summary of surveys and occurrences data ----
kelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN)
nokelp_summary <- mFD::asb.sp.summary(asb_sp_w = nokelp_sp_maxN)

# retrieve occurrence matrix:
kelp_sp_occ <- kelp_summary$asb_sp_occ
nokelp_sp_occ <- nokelp_summary$asb_sp_occ


# dimensions
dim(kelp_sp_occ) # 69 assemblages * 101 species
dim(nokelp_sp_occ) # 56 assemblages * 106 species

# names of species
kelp_sp <- colnames(kelp_sp_occ) 
length(kelp_sp) # 101 sp

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


############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata <- read.csv(here::here("data", "raw_data",  "SpatialUVC_metadata_site.csv"))
head(spatial_metadata)

spatial_sp_biom <- read.csv(here::here("data", "raw_data", "SpatialUVC_species_biomass_site_average.csv")) %>%
  column_to_rownames("Code") %>% 
  as.matrix()

dim(spatial_sp_biom) # 9 assemblages * 53 species

# summary of surveys and occurrences data ----

spatial_summary <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom)

# retrieve occurrence matrix:
spatial_sp_occ <- spatial_summary$asb_sp_occ

# species names
spatial_sp <- colnames(spatial_sp_occ)
length(spatial_sp) # 53

## names of species present in at least one dataset ####
sum(spatial_sp %in% temporal_sp)# 36 species shared 

species_allsurveys <- unique( c(temporal_sp ,  spatial_sp) ) 
length(species_allsurveys) # 139 species

## saving dataframes #####

save(kelp_metadata, file=here::here("scripts", "extras", "kelp_metadata_extra.RData") )
save(kelp_sp_occ, file=here::here("scripts", "extras", "kelp_sp_occ_extra.RData") )
save(kelp_summary, file=here::here("scripts", "extras", "kelp_summary_extra.RData") )

save(nokelp_metadata, file=here::here("scripts", "extras", "nokelp_metadata_extra.RData") )
save(nokelp_sp_occ, file=here::here("scripts", "extras", "nokelp_sp_occ_extra.RData") )
save(nokelp_summary, file=here::here("scripts", "extras", "nokelp_summary_extra.RData") )

save(spatial_metadata, file=here::here("scripts", "extras", "spatial_metadata_extra.RData") )
save(spatial_sp_occ, file=here::here("scripts", "extras", "spatial_sp_occ_extra.RData") )

save(spatial_summary, file=here::here("scripts", "extras", "spatial_summary_extra.RData") )

save(species_allsurveys, file=here::here("scripts", "extras", "species_allsurveys_extra.RData") )

################### end of code ##############################################################


################################################################################
##
## Script for computing species positions in a multidimensional space
## according to their trait values
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)

## loading ####

# load traits data ----
fish_traits <- read.csv(here::here("scripts", "extras", "fish_traits_extra.csv"), header = T)

# load species names from surveys datasets ----
load(here::here("scripts", "extras", "species_allsurveys_extra.RData") )
length(species_allsurveys) # 139 species    

# checking same species in trait and occurrences datasets ----
identical ( sort(species_allsurveys) , sort(fish_traits$Species ) ) # True

## preparing trait dataset ####

# trait values in a dataframe (species in alphabetical order) ----
sp_tr <- fish_traits %>%
  arrange("Species") %>%
  column_to_rownames("Species") %>%
  as.data.frame()
head(sp_tr)

nrow(sp_tr) # 139 species

# recoding variable to match trait type ----
# looking at trait values
lapply(sp_tr, unique)

# trait type
tr_cat<-data.frame( trait_name = c("Size", "Agg", "Position", "Diet", "Kmax"),
                    trait_type = c("O","O","O", "N", "Q") )

# size as ordinal
sp_tr$Size <- factor(sp_tr$Size, 
                     levels = c("S1", "S2", "S3", "S4", "S5", "S6"),
                     ordered = TRUE)
summary(sp_tr$Size)

# aggregation as ordinal
sp_tr$Agg <- factor(sp_tr$Agg, 
                    levels = c("Solitary", "Pair", "Group"),
                    ordered = TRUE)
summary(sp_tr$Agg)

# Position as ordinal
sp_tr$Position <- factor(sp_tr$Position, 
                         levels = c("Benthic", "BenthoP", "Pelagic"),
                         ordered = TRUE)
summary(sp_tr$Position)

# diet as factor
sp_tr$Diet <- as.factor(sp_tr$Diet)
summary(sp_tr$Diet)

#Kmax as numeric

sp_tr$Kmax <- as.numeric(sp_tr$Kmax)
summary(sp_tr$Kmax)

# summary of trait data----
summary_traits <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                     sp_tr  = sp_tr)


## Computing Gower distance between species ####
sp_gower_dist <- mFD::funct.dist(sp_tr=sp_tr, tr_cat = tr_cat, 
                                 metric="gower")
# => no need to compute FE since all indices are not sensitive to redundant species
range(sp_gower_dist) # from 0 to 1

### Compute functional spaces and their quality:

# mean absolute deviation index (as quality metric)
funct_spaces<- mFD::quality.fspaces(sp_dist = sp_gower_dist, maxdim_pcoa = 12, 
                                    deviation_weighting = "absolute", fdist_scaling = FALSE,
                                    fdendro = "average") 

funct_spaces$quality_fspaces

# => 3D space has the lowest mAD (0.055)

# species coordinates
sp_3D_coord<-funct_spaces$details_fspaces$sp_pc_coord[,1:3]
summary(sp_3D_coord)
# saving ####

# trait values and trait coding dataframes ----

save(sp_tr, file=here::here("scripts", "extras", "sp_tr_extra.RData") )
save(tr_cat, file=here::here("scripts", "extras", "tr_cat_extra.RData") )
save(summary_traits, file=here::here("scripts", "extras", "summary_traits_extra.RData") )
save(sp_gower_dist, file=here::here("scripts", "extras", "sp_gower_dist_extra.RData") )
save(sp_3D_coord, file=here::here("scripts", "extras", "sp_3D_coord_extra.RData") )
save(funct_spaces, file=here::here("scripts", "extras", "funct_spaces_extra.RData") )
save(sp_faxes_coord, file=here::here("scripts", "extras", "sp_faxes_coord_extra.RData") )

##################################  end of code ######################################################################

################################################################################
##
## Script for computing taxonomic and functional diversity between years for kelp and no kelp sites,
##
## and also for spatial data
## 


# loading data
load(here::here("scripts", "extras", "kelp_sp_occ_extra.RData") )
load(here::here("scripts", "extras", "nokelp_sp_occ_extra.RData") )
load(here::here("scripts", "extras", "spatial_sp_occ_extra.RData") )
load(here::here("scripts", "extras", "sp_3D_coord_extra.RData") )

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

## computing taxonomic and functional diversity for kelp sites ####
temporal_fd_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = kelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp <- temporal_fd_kelp$functional_diversity_indices

## computing taxonomic and functional diversity for spatial data ####

spatial_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = spatial_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)
spatial_alpha <- spatial_fd$functional_diversity_indices

# saving ####

# trait values and trait coding dataframes ----
save(temporal_fd_nokelp, file=here::here("scripts", "extras", "temporal_fd_nokelp_extra.RData") )
save(temporal_alpha_nokelp, file=here::here("scripts", "extras", "temporal_alpha_nokelp_extra.RData") )

save(temporal_fd_kelp, file=here::here("scripts", "extras", "temporal_fd_kelp_extra.RData") )
save(temporal_alpha_kelp, file=here::here("scripts", "extras", "temporal_alpha_kelp_extra.RData") )

save(spatial_fd, file=here::here("scripts", "extras", "spatial_fd_extra.RData") )
save(spatial_alpha, file=here::here("scripts", "extras", "spatial_alpha_extra.RData") )


#################################### end of script ####################################################################

################################################################################
##
## Script for plotting taxonomic and functional diversity across habitats and years
##
## Fig. 1 - Species richness and functional richness
## 
## Fig. 2 - Functional Dispersion
##
## Fig. 1-3S - Functional Identity
##
##
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(patchwork)

# loading data
load(here::here("scripts", "extras", "spatial_metadata_extra.RData") )
load(here::here("scripts", "extras", "spatial_alpha_extra.RData") )

load(here::here("scripts", "extras", "kelp_metadata_extra.RData") )
load(here::here("scripts", "extras", "nokelp_metadata_extra.RData") )
load(here::here("scripts", "extras", "temporal_alpha_kelp_extra.RData") )
load(here::here("scripts", "extras", "temporal_alpha_nokelp_extra.RData") )


## spatial trends ####

# merging metadata and biodiv indices in a single table for each habitat type

spatial_all <- spatial_metadata %>% 
  left_join( rownames_to_column(spatial_alpha, "Code"), by="Code" ) %>%
  select(Code, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3)



# mean and sd of diversity among each site for each year in each habitat type
spatial_toplot <- spatial_all %>%
  group_by(Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    fric_mean = mean(fric),
    fric_sd = sd(fric), 
    fdis_mean = mean(fdis),
    fdis_sd = sd(fdis), 
    fide_PC1_mean = mean(fide_PC1),
    fide_PC1_sd = sd(fide_PC1),
    fide_PC2_mean = mean(fide_PC2),
    fide_PC2_sd = sd(fide_PC2),
    fide_PC3_mean = mean(fide_PC3),
    fide_PC3_sd = sd(fide_PC3)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( fric_se = fric_sd/sqrt(n)) %>% 
  mutate( fdis_se = fdis_sd/sqrt(n)) %>% 
  mutate( fide_PC1_se = fide_PC1_sd/sqrt(n)) %>%
  mutate( fide_PC2_se = fide_PC2_sd/sqrt(n)) %>%
  mutate( fide_PC3_se = fide_PC3_sd/sqrt(n))


spatial_toplot

# color code for the 3 habitats
hab_colors <- c(Inshore= "#2C6BAA", Midshelf= "lightsalmon1", Offshore="firebrick3")

# taxonomic ----

plot_spatial_taxo <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=TRic_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  scale_y_continuous( limits = c(0,25), breaks = seq(from=0, to=35, by=5)  ) +
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none", axis.text.x = element_blank())
plot_spatial_taxo

# functional ----

plot_spatial_func <- ggplot(spatial_toplot) +
  geom_bar( aes(x=Habitat, y=fric_mean, color = Habitat, fill = Habitat), stat="identity", color = "black", size=0.8) +
  geom_errorbar( aes(x=Habitat, ymin=fric_mean-fric_se, ymax=fric_mean+fric_se), width=0.1, size=0.8, colour="black" ) +
  scale_color_manual(values=hab_colors) + 
  scale_fill_manual(values=hab_colors) + 
  #scale_y_continuous( limits = c(0,0.4), breaks = seq(from=0, to=0.4, by=0.1)  ) +
  scale_y_continuous( limits = c(0,0.5), breaks = seq(from=0, to=0.5, by=0.1)  ) +
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")

plot_spatial_func

## temporal trends ####

# merging metadata and biodiv indices in a single table for each habitat

temporal_kelp<- kelp_metadata %>% 
  mutate(Habitat="kelp") %>%
  left_join( rownames_to_column(temporal_alpha_kelp, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3 )

temporal_nokelp<- nokelp_metadata %>% 
  mutate(Habitat="no_kelp") %>%
  left_join( rownames_to_column(temporal_alpha_nokelp, "Code"), by="Code" ) %>%
  select(Code, Site, Year, Habitat, TRic=sp_richn, fric, fdis, fide_PC1, fide_PC2, fide_PC3)

# merging values for the 2 habitat types
temporal_all <- bind_rows(temporal_kelp, temporal_nokelp )
head(temporal_all)

write.csv(temporal_all, file=here::here("scripts", "extras", "temporal_all.csv"),
                                       row.names = FALSE)

# mean and sd of diversity among each site for each year in each habitat type

temporal_toplot <- temporal_all %>%
    group_by(Year, Habitat) %>%
  summarise( 
    n = n(),
    TRic_mean = mean(TRic),
    TRic_sd = sd(TRic),
    fric_mean = mean(fric),
    fric_sd = sd(fric), 
    fdis_mean = mean(fdis),
    fdis_sd = sd(fdis), 
    fide_PC1_mean = mean(fide_PC1),
    fide_PC1_sd = sd(fide_PC1),
    fide_PC2_mean = mean(fide_PC2),
    fide_PC2_sd = sd(fide_PC2),
    fide_PC3_mean = mean(fide_PC3),
    fide_PC3_sd = sd(fide_PC3)
  ) %>%
  mutate( TRic_se = TRic_sd/sqrt(n))  %>%
  mutate( fric_se = fric_sd/sqrt(n)) %>% 
  mutate( fdis_se = fdis_sd/sqrt(n)) %>% 
  mutate( fide_PC1_se = fide_PC1_sd/sqrt(n)) %>%
  mutate( fide_PC2_se = fide_PC2_sd/sqrt(n)) %>%
  mutate( fide_PC3_se = fide_PC3_sd/sqrt(n))
temporal_toplot

unique(temporal_toplot$Year)

# color code for kelp/no kelp

year_colors <- c(kelp= "seagreen4", no_kelp= "seagreen2")


# taxonomic ----

plot_tempo_taxo <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=TRic_mean), stat="identity", size=3) +
  geom_line(aes(x= Year, y= TRic_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=TRic_mean-TRic_se, ymax=TRic_mean+TRic_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(15,40), breaks = seq(from=15, to=40, by=5)  ) +
  scale_color_manual(values=year_colors) + 
  scale_fill_manual(values=year_colors) + 
  #scale_y_continuous( limits = c(15,30), breaks = seq(from=15, to=30, by=5)  ) +
  #scale_color_manual(values="seagreen4") + 
  #scale_fill_manual(values="seagreen4") + 
  labs(x="", y="Species richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_tempo_taxo

# functional ----

plot_tempo_func <- ggplot(temporal_toplot, 
                          mapping=aes(color = Habitat, fill = Habitat) ) +
  geom_point( aes(x=Year, y=fric_mean), stat="identity", size=2, shape=16) +
  geom_line(aes(x= Year, y= fric_mean) , stat="identity", size=1)+
  geom_errorbar( aes(x=Year, ymin=fric_mean-fric_se, ymax=fric_mean+fric_se), width=0.4, size=0.8) +
  scale_x_continuous( limits = c(2001, 2019), breaks = seq(from=2002, to=2018, by=4)  ) +
  scale_y_continuous( limits = c(0.4,0.8), breaks = seq(from=0.4, to=0.8, by=0.2)  ) +
  scale_color_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  scale_fill_manual(values=year_colors, name="Habitat", breaks = c("kelp", "no_kelp"), labels=c("Kelp", "No kelp")) + 
  #scale_color_manual(values="seagreen4") + 
  #scale_fill_manual(values="seagreen4") + 
  labs(x="", y="Functional richness") +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.title = element_text(size=14), legend.text = element_text(size=14),
        legend.background = element_blank(), legend.key=element_blank())
plot_tempo_func


## merging all plot into a single figure and saving as png ####
figure1 <- ( plot_spatial_taxo + plot_tempo_taxo ) / ( plot_spatial_func +  plot_tempo_func )

ggsave(figure1, file=here::here("outputs/", "Figure1.png"),
              height = 22, width = 25, unit = "cm" )
       
       