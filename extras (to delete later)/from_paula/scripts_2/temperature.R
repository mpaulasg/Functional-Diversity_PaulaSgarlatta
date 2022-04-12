################################################################################
##
## Script to obtain temperature
## 
## Code by Paula Sgarlatta
##
################################################################################


rm(list=ls()) # cleaning memory
# libraries
library(tidyverse)
library(here)
library(dplyr)
library(lubridate)


#Load data

temp <- read.csv(here::here("data", "raw_data",  "temperature.csv"))

#Date is in USA format as is the only way lubridate works (at least for me)

temp_august <- temp %>% 
  mutate(Month = month(Date)) %>% 
  mutate(Year = year(Date)) %>% 
  mutate(Day = day(Date)) %>% 
  filter(Year %in% c("2002", "2003", "2004", "2005", "2006","2007",
    "2008", "2009", "2010", "2011","2013","2014", "2015", "2018")) %>% 
  filter(Month %in% "8") %>% 
  mutate(average_midshelf= rowMeans(cbind(NWMA, Groper), na.rm=T))
  
unique(temp_august$Year) #14 years

temp_toplot <- temp_august %>%
  group_by(Year) %>%
  summarise( 
    n = n(),
    temp_mean = mean(average_midshelf),
    temp_sd = sd(average_midshelf)) %>% 
  mutate( temp_se = temp_sd/sqrt(n))

plot_temp <- ggplot(temp_toplot, aes(x= Year, y= temp_mean)) +
  geom_line() +
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks.y = element_blank(), 
        panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(size = (14),colour = "black"), axis.title = element_text(size= (16)),
        legend.position = "none")
plot_temp

