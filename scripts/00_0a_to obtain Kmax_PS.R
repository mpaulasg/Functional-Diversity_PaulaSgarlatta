rm(list=ls()) # cleaning memory
# libraries
library(tidyverse)
library(here)
library(dplyr)
library(xgboost)
library(rfishprod)


#Load data

for_K <- read.csv(here::here("from_paula", "fish_for_K.csv"),
                 header = T)

#To predict Kmax, you need diferent parameters, which are:

#mean sea surface temperature(sstmean)
# Diet: herbivores/detritivores (HerDet), herbivores/macroalgivores(HerMac), omnivores(Omnivr), 
#planktivores(Plktiv), sessile invertivores(InvSes), mobile invertivores(InvMob),  
#fish and cephalopod predators (FisCep).
#Position: Pelagic reef dwelling (PelgDw), Pelagic reef associated (PelgAs), Bentho-pelagic reef dwelling (BtPlDw), 
#Benthic reef associated (BnthAs)or dwelling (BnthDw). -> See Morais & Bellwood 2018 for more details.
#Method: 

# Then you need to be sure that the parameters in your data match perfectly with db (data from the package). The
#parameters needs to be written exactly as in db

for_K <-tidytrait (for_K, db) 

#Awesome! Now, we will use the formula from Morais & Bellwood 2018

fmod <- formula (~ sstmean + MaxSizeTL + Diet + Position + Method)

# Predicting Kmax, the standardised VBGF parameter (Recommendation: use 100s to 1000s iterations) # It takes several 
#minutes: in my computer, ~ 7 min

datagr <- predKmax (for_K, 
                    dataset = db,
                    fmod = fmod,
                    niter = 1000,
                    return = 'pred')

# apply the prediction

datagr <- datagr$pred

#Now you can see Kmax and the quantiles attached to your data.

#Save data

write.csv(datagr, file=here::here("from_paula", "fish_K_values.csv"), 
          row.names = FALSE )
