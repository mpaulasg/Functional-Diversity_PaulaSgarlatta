################################################################################
##
## Script for spatial FD analysis
##
## ## Code taken from Boulanger et al. 2021 paper: 
#  "Environmental DNA metabarcoding reveals and unpacks a 
#  biodiversity conservation paradox in Mediterranean marine reserves."  
#  https://github.com/eboulanger/MEDeDNA--reserves
#
#  Modified by Paula Sgarlatta
##
## Script to run:
#' * richness analyses 
#' * dissimilarity analyses
#' 
################################################################################

rm(list=ls()) # cleaning memory

##Load packages

library(vegan)
library(glmmTMB)
library(car)

load(here::here("outputs", "FD_beta_kelp_Hill_sites.RData") )

## Preparing data for stats

beta_stats <- FD_beta_kelp_Hill_sites %>%
  filter(Year1==2002) %>% 
  mutate(Year=paste(Year1, Year2, sep="-")) %>% 
  select(Year, Site, q0, q1, q2)


##GLMM

FBeta_Hill_spatial <- glmmTMB(q2 ~ Year + (1|Site), data = beta_stats, family = beta_family())

summary(FBeta_Hill_spatial)

mod.res <- simulateResiduals(FBeta_Hill_spatial)
plot(mod.res)


mod<-lm(q2 ~ Year,
        data = beta_stats)

summary(mod)

Anova(mod)

library(DHARMa)

mod.res <- simulateResiduals(mod)
plot(mod.res)

## Pairwise comparisons

library(emmeans)

emmeans(mod, pairwise ~ Year)


#ANOVA for S

s_aov_site <- aov(sp_richn~habitat ,data=fd_values_site)

summary(s_aov_site)

#Check  normality of residuals visually via 
#histogram and QQ-plot

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(s_aov_site$residuals)

# QQ-plot

qqPlot(s_aov_site$residuals,
       id = FALSE # id = FALSE to remove point identification
)

#Homocedasticity

plot(s_aov_site)

####[PS] Probably better to do this analysis with transect data

s_aov_transect <- aov(sp_richn~habitat ,data=fd_values_transect)

summary(s_aov_transect)

#Not significant 

#Check  normality of residuals visually via 
#histogram and QQ-plot

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(s_aov_transect$residuals)

# QQ-plot

qqPlot(s_aov_transect$residuals,
       id = FALSE # id = FALSE to remove point identification
)

#Homocedasticity

plot(s_aov_transect)

#[PS] Residuals  don't look great - transformation?

# distance-based redundancy analysis
# with jaccard dissimilarity

fric_dbrda_spatial <- capscale(alpha_fd_indices_site$functional_diversity_indices$fric ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(fric_dbrda_spatial, by = "terms", permutations = 9999)

#Significcant - p=0.0152*

meta_data_for_functional <- meta_data_transect[-c(2,10), ]

fric_dbrda_transect <- capscale(alpha_fd_indices_transect$functional_diversity_indices$fric ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(fric_dbrda_transect, by = "terms", permutations = 9999)

#significative - p=0.0074**


fdis_dbrda_spatial <- capscale(alpha_fd_indices_site$functional_diversity_indices$fdis ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(fdis_dbrda_spatial, by = "terms", permutations = 9999)

#Not significant-p=0.207

fdis_dbrda_transect <- capscale(alpha_fd_indices_transect$functional_diversity_indices$fdis ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(fdis_dbrda_transect, by = "terms", permutations = 9999)

#Not significant-p=0.08

fspe_dbrda_spatial <- capscale(alpha_fd_indices_site$functional_diversity_indices$fspe ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(fspe_dbrda_spatial, by = "terms", permutations = 9999)

#Not significant-p=0.909

fspe_dbrda_transect <- capscale(alpha_fd_indices_transect$functional_diversity_indices$fspe ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(fspe_dbrda_transect, by = "terms", permutations = 9999)

#Not significant-p=0.494

##Beta diversity

tax_diss_dbrda_spatial  <- capscale(totaldis_tax_spatial_site ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(tax_diss_dbrda_spatial, by = "terms", permutations = 9999)

#Significant-p=0.0025

tax_diss_dbrda_transect  <- capscale(totaldis_tax_spatial_transect ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(tax_diss_dbrda_transect, by = "terms", permutations = 9999)

#Significant-p=1e-04****

tax_turn_dbrda_spatial  <- capscale(turnover_tax_spatial_site ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(tax_turn_dbrda_spatial, by = "terms", permutations = 9999)

#Significant-p=0.0032

tax_turn_dbrda_transect  <- capscale(turnover_tax_spatial_transect ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(tax_turn_dbrda_transect, by = "terms", permutations = 9999)

#Significant-p=1e-04****

###Functional beta diversity

func_diss_dbrda_spatial<- capscale(beta_fd_indices_s$pairasb_fbd_indices$jac_diss ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(func_diss_dbrda_spatial, by = "terms", permutations = 9999)

#Significant-p=0.0034****

func_diss_dbrda_transect<- capscale(beta_fd_indices_transect$pairasb_fbd_indices$jac_diss ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(func_diss_dbrda_transect, by = "terms", permutations = 9999)

#Significant-p=1e-04****

func_turn_dbrda_spatial<- capscale(beta_fd_indices_s$pairasb_fbd_indices$jac_turn ~ Habitat, data = meta_data_site, distance = "jaccard")
anova(func_turn_dbrda_spatial, by = "terms", permutations = 9999)

#Significant - p=0.0234

func_turn_dbrda_transect<- capscale(beta_fd_indices_transect$pairasb_fbd_indices$jac_turn ~ Habitat, data = meta_data_for_functional, distance = "jaccard")
anova(func_turn_dbrda_transect, by = "terms", permutations = 9999)

#Not significant - p=0.1014

