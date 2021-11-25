################################################################################
##
## Script for importing all data from csv files
## 
## Code by Camille Magneville, Sébastien villéger and Paula Sgarlatta
##
## 1/ load and preparing trait datasets - kelp and no kelp sites
## 2/ 
################################################################################

#### 1 - Load (and transform) datasets ####

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)

## 1 trait data ####

# load traits data ----

# species from sites with kelp
traits_kelp <- read.csv(here::here("data/TemporalBRUV_species_traits_kelp.csv"),
                        header = TRUE, row.names = 1)
head(traits_kelp)
nrow(traits_kelp)

#Traits for sp presents in sites that never had kelp
traits_no_kelp <- read.csv(here::here("data/TemporalBRUV_species_traits_no_kelp.csv"),
                           header = TRUE, row.names = 1)
head(traits_no_kelp)

# merging in a single data.frame and trait named shortened (4 characters) ----
sp_tr <- as.data.frame( rbind( traits_kelp, traits_no_kelp) )
names(sp_tr) <- c("Size", "Aggr", "Posi", "Diet")
head(sp_tr)

nrow(sp_tr)
unique(row.names(sp_tr) )

# recoding variable to match trait type ---

# looking at trait values
lapply(sp_tr, unique)

# trait type
tr_cat<-data.frame( trait_name = c("Size", "Aggr", "Posi", "Diet"),
                    trait_type = c("O","O","O","N") )

# size as ordinal
sp_tr$Size <- factor(sp_tr$Size, 
                     levels = c("S2", "S3", "S4", "S5", "S6", "S7"),
                     ordered = TRUE)
summary(sp_tr$Size)

# aggregation as ordinal
sp_tr$Aggr <- factor(sp_tr$Aggr, 
                    levels = c("Solitary", "Pair", "Group"),
                    ordered = TRUE)
summary(sp_tr$Aggr)

# Position as ordinal
sp_tr$Posi <- factor(sp_tr$Posi, 
                    levels = c("Benthic", "BenthoP", "Pelagic"),
                    ordered = TRUE)
summary(sp_tr$Posi)

# diet as factor
sp_tr$Diet <- as.factor(sp_tr$Diet)
summary(sp_tr$Diet)

# summary of trait data----
summary_traits <- mFD::sp.tr.summary(tr_cat = tr_cat, 
                                      sp_tr  = sp_tr)
summary_traits

# saving trait values and trait coding dataframes ----
save(sp_tr, file=here::here("data/", "sp_tr.RData") )
save(tr_cat, file=here::here("data/", "tr_cat.RData") )
save(summary_traits, file=here::here("outputs/", "summary_traits.RData") )


####################### # => trait data ready ####




