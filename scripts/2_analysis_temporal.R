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
#
################################################################################

#rm(list=ls()) # cleaning memory

##Load packages


library(ggplot2)
library(car)
library(rstatix)

#ANOVA for S

s_aov_year_kelp <- aov(sp_richn~year ,data=fd_temporal_kelp)

summary(s_aov_year_kelp)

#Check  normality of residuals visually via 
#histogram and QQ-plot

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(s_aov_year_kelp$residuals)

# QQ-plot

qqPlot(s_aov_year_kelp$residuals,
       id = FALSE # id = FALSE to remove point identification
)

#Homocedasticity

plot(s_aov_year_kelp)

####[PS] Assumptions met!

TukeyHSD(s_aov_year_kelp)  ## check


s_aov_year_no_kelp <- aov(sp_richn~year ,data=fd_temporal_no_kelp)

summary(s_aov_year_no_kelp)

#Check  normality of residuals visually via 
#histogram and QQ-plot

par(mfrow = c(1, 2)) # combine plots

# histogram
hist(s_aov_year_no_kelp$residuals)

# QQ-plot

qqPlot(s_aov_year_no_kelp$residuals,
       id = FALSE # id = FALSE to remove point identification
)

#Homocedasticity

plot(s_aov_year_no_kelp)

####[PS] Assumptions met!

TukeyHSD(s_aov_year_no_kelp)  ## check

#Significative - p < 2.2e-16

#Check assumptions

qqnorm(residuals(s_aov_year))

# distance-based redundancy analysis
# with jaccard dissimilarity

fric_dbrda_kelp <- capscale( fd_temporal_kelp$fric~ Year, data = meta_data_kelp, distance = "jaccard")
anova(fric_dbrda_kelp, by = "terms", permutations = 9999)

#Significant - p<1e-04

fric_dbrda_no_kelp <- capscale( fd_temporal_no_kelp$fric~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(fric_dbrda_no_kelp, by = "terms", permutations = 9999)

#Significant - p<1e-04

##GRAPH

fdis_dbrda_kelp <- capscale( fd_temporal_kelp$fdis~ Year, data = meta_data_kelp, distance = "jaccard")
anova(fdis_dbrda_kelp, by = "terms", permutations = 9999)

#Not significant - 0.1906

fdis_dbrda_no_kelp <- capscale( fd_temporal_no_kelp$fdis~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(fdis_dbrda_no_kelp, by = "terms", permutations = 9999)

# Significant - p=0.0027**


fspe_dbrda_kelp <- capscale( fd_temporal_kelp$fspe~ Year, data = meta_data_kelp, distance = "jaccard")
anova(fspe_dbrda_kelp, by = "terms", permutations = 9999)

#Significant - p=0.0057**

fspe_dbrda_no_kelp <- capscale( fd_temporal_no_kelp$fspe~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(fspe_dbrda_no_kelp, by = "terms", permutations = 9999)

#Significant - p=0.01*

###Taxonomic beta diversity

tax_diss_dbrda_kelp  <- capscale(totaldis_tax_kelp ~ Year, data = meta_data_kelp, distance = "jaccard")
anova(tax_diss_dbrda_kelp, by = "terms", permutations = 9999)

tax_diss_dbrda_no_kelp  <- capscale(totaldis_tax_no_kelp ~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(tax_diss_dbrda_no_kelp, by = "terms", permutations = 9999)

#Significant - p=1e-04

tax_turn_dbrda_kelp  <- capscale(turnover_tax_kelp ~ Year, data = meta_data_kelp, distance = "jaccard")
anova(tax_turn_dbrda_kelp, by = "terms", permutations = 9999)

#Significant - p<1e-04

tax_turn_dbrda_no_kelp  <- capscale(turnover_tax_no_kelp ~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(tax_turn_dbrda_no_kelp, by = "terms", permutations = 9999)

#Significant - p=0.0021**

#####Functional beta diversity

func_diss_dbrda_kelp<- capscale(totaldis_func_kelp ~ Year, data = meta_data_kelp, distance = "jaccard")
anova(func_diss_dbrda_kelp, by = "terms", permutations = 9999)

#Significant - p<1e-04

func_diss_dbrda_no_kelp<- capscale(totaldis_func_no_kelp ~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(func_diss_dbrda_no_kelp, by = "terms", permutations = 9999)

#Significant - p<1e-04

func_turn_dbrda_kelp<- capscale(turnover_func_kelp ~ Year, data = meta_data_kelp, distance = "jaccard")
anova(func_turn_dbrda_kelp, by = "terms", permutations = 9999)

#Not significant - p=0.08

func_turn_dbrda_no_kelp<- capscale(turnover_func_no_kelp ~ Year, data = meta_data_no_kelp, distance = "jaccard")
anova(func_turn_dbrda_no_kelp, by = "terms", permutations = 9999)

## Significant - p=6e-04***

