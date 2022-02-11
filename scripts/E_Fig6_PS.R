################################################################################
##
## Script for plotting quality of functional space and relation between traits and 
##
## PcoAs
## 
## Code by Paula Sgarlatta, Camille Magneville and Sebastien Villeger 
##
################################################################################


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(ggplot2)


## loading ####

load(here::here("outputs", "funct_spaces.RData") )
load(here::here("outputs", "sp_faxes_coord.RData") )
load(here::here("data", "sp_tr.RData") )



# illustrate quality of functional space

qual_space <- mFD::quality.fspaces.plot(
  fspaces_quality            = funct_spaces,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")


plot(qual_space)

### Test correlation between traits and functional axes:

cor_tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = sp_tr, 
  sp_faxes_coord = sp_faxes_coord[, c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

# get the table of correlation:

corr_table <- as.data.frame(cor_tr_faxes$tr_faxes_stat)  

# get the plot:

plot_corr <- cor_tr_faxes$tr_faxes_plot


ggsave(qual_space, file=here::here("outputs/", "Figure4S.png"),
       height = 16, width = 30, unit = "cm" )

ggsave(plot_corr, file=here::here("outputs/", "Figure5S.png"),
       height = 16, width = 30, unit = "cm" )

write.csv(corr_table, file=here::here("outputs/", "Correlation_traits_table.csv"), 
          row.names = FALSE)


########################################## end of code ###############################################################
