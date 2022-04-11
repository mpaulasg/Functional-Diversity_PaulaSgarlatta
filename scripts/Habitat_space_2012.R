################################################################################
##
## Script for analysign habitat data from the cross-shelf gradient
##
## 
## Code by Paula Sgarlatta
##
################################################################################


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)

# Load data

habitat <- read.csv(here::here("data", "using biomass-maxN", "habitat_solitaries_2012.csv"))

habitat_toplot <- habitat %>% 
  pivot_longer(cols= 4:10, names_to="Group") %>%
  filter(Site!="Flat_top", Site!="Look_at_me_now") %>% 
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Habitat,Transect, Group) %>%
  summarize(total=sum(value)) %>% 
  group_by(Site, Transect) %>% 
  mutate(percent_cover=(total*100/125)) %>% #125 points per transect
  select(-total)
  


habitat_summary <-habitat_toplot%>%
  group_by(Habitat, Group)%>%
  summarise(Percent_cover_habitat=mean(percent_cover,na.rm=T),
            n = n(), mean = mean(percent_cover), se = sd(percent_cover/sqrt(n))) %>% 
  mutate(se) %>%
  group_by(Habitat) %>%
  arrange(desc(Group)) %>%
  mutate(
    pos = cumsum(Percent_cover_habitat),
    upper = pos + se/2,
    lower = pos - se/2
  ) %>%
  ungroup() %>% 
  filter(mean > 1)



####(FIGURE 1) ----

habitat_v3 <- ggplot(habitat_summary, aes(x=Habitat, y=Percent_cover_habitat, fill=Group)) +
    geom_bar(stat="identity", width = 0.7, size = 1, color = "black") + 
    geom_errorbar(aes(ymin = lower, ymax = upper), size = 0.8, width=0.1, position = "identity", color = "black")+
  theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), axis.ticks = element_blank(), 
                  panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
                  axis.text = element_text(size = (18), color = "black"), axis.title = element_text(size= (18)),
                  legend.key = element_rect(fill = "white"), legend.text = element_text(size=20))+
  labs(y="Percent cover (%)", x="") + scale_fill_discrete(name = "", 
labels = c("Coral", "Ecklonia radiata", "Macroalgae", "Other invertebrates",
           "Rock & sand", "Sponges & tunicates", "Turf & CCA"))

## Error bars not working

ggsave(habitat_v3, file=here::here("outputs/", "habitat_stack.png"),
       height = 16, width = 24, unit = "cm" )


############# Using 2019 data

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(tidytext)

# Load data

habitat <- read.csv(here::here("data", "using biomass-maxN", "PQs_cross_shelf_to_plot_2019.csv"))

habitat_toplot <- habitat %>% 
  pivot_longer(cols= 6:13, names_to="Group") %>%
  filter(Site!="Whitmore") %>% 
  mutate(value=as.numeric(value)) %>%
  group_by(PQ, Transect, Group, Site, Habitat) %>%
  summarize(total=sum(value)) %>% 
  group_by(PQ, Transect, Group, Site, Habitat) %>% 
  mutate(percent_cover=(total*100/25)) %>% #25 points per PQ
  select(-total)


####(FIGURE X-1) ----
habitat_summary <-habitat_toplot%>%
  group_by(Site, Habitat, Group)%>%
  summarise(Percent_cover_site=mean(percent_cover,na.rm=T),
            stdev=sd(percent_cover))

habitat_plot <- ggplot()+
  geom_jitter(data=habitat_toplot, aes(x= reorder_within(Group, -percent_cover, Site, fun=mean), #note here that I included -Percent_cover to reorder descending
                                       y = percent_cover, colour=Group), width=0.15, alpha=0.4)+
  facet_wrap(.~Site, scales="free", nrow=3)+
  geom_errorbar(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site, ymin=Percent_cover_site-stdev,ymax=Percent_cover_site+stdev),fill = "black",width=.15)+
  geom_point(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site), size=4.5)+
  geom_point(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site, colour=Group), size=4)+
  theme_bw()+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

ggsave(habitat_plot, file=here::here("outputs/", "habitat_sites_2019.png"),
       height = 16, width = 24, unit = "cm" )

####(FIGURE X-2) ----
habitat_summary_habitat <-habitat_toplot%>%
  group_by(Habitat, Group)%>%
  summarise(Percent_cover_habitat=mean(percent_cover,na.rm=T),
            stdev=sd(percent_cover))

habitat_plot_hab <- ggplot()+
  geom_jitter(data=habitat_toplot, aes(x= reorder_within(Group, -percent_cover, Habitat, fun=mean), #note here that I included -Percent_cover to reorder descending
                                       y = percent_cover, colour=Group), width=0.15, alpha=0.4)+
  facet_wrap(.~Habitat, scales="free")+
  geom_errorbar(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Habitat, fun=mean), y = Percent_cover_habitat, ymin=Percent_cover_habitat-stdev,ymax=Percent_cover_habitat+stdev),fill = "black",width=.15)+
  geom_point(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Habitat, fun=mean), y = Percent_cover_habitat), size=4.5)+
  geom_point(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Habitat, fun=mean), y = Percent_cover_habitat, colour=Group), size=4)+
  theme_bw()+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

ggsave(habitat_plot_hab, file=here::here("outputs/", "habitat_cross-shelf_2019.png"),
       height = 16, width = 24, unit = "cm" )

####(FIGURE X-3) ----

habitat_v3<- ggplot(habitat_summary_habitat, aes(x=Habitat, y=Percent_cover_habitat, fill=Group)) +
  geom_bar(stat="identity", 
           position="stack",
           colour="black") + 
  theme_bw()
# geom_errorbar(aes(ymin=Percent_cover_habitat-stdev,ymax=Percent_cover_habitat+stdev), 
#               size=0.5, width=.25, stat = "identity",
#               position = "identity")

## Error bars not working

ggsave(habitat_v3, file=here::here("outputs/", "habitat_stack_2019.png"),
       height = 16, width = 24, unit = "cm" )


###################################################################


############# Using LATITUDINAL data

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(tidytext)

# Load data

habitat <- read.csv(here::here("data", "using biomass-maxN", "PQs_latitudinal.csv"))

habitat_toplot <- habitat %>% 
  pivot_longer(cols= 6:13, names_to="Group") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(PQ, Transect, Group, Site, Area) %>%
  summarize(total=sum(value)) %>% 
  group_by(PQ, Transect, Group, Site, Area) %>% 
  mutate(percent_cover=(total*100/25)) %>% #25 points per PQ
  select(-total)


####(FIGURE X-1) ----
habitat_summary <-habitat_toplot%>%
  group_by(Site, Area, Group)%>%
  summarise(Percent_cover_site=mean(percent_cover,na.rm=T),
            stdev=sd(percent_cover))

habitat_plot <- ggplot()+
  geom_jitter(data=habitat_toplot, aes(x= reorder_within(Group, -percent_cover, Site, fun=mean), #note here that I included -Percent_cover to reorder descending
                                       y = percent_cover, colour=Group), width=0.15, alpha=0.4)+
  facet_wrap(.~Site, scales="free", nrow=3)+
  geom_errorbar(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site, ymin=Percent_cover_site-stdev,ymax=Percent_cover_site+stdev),fill = "black",width=.15)+
  geom_point(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site), size=4.5)+
  geom_point(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site, colour=Group), size=4)+
  theme_bw()+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

ggsave(habitat_plot, file=here::here("outputs/", "habitat_sites_2019.png"),
       height = 16, width = 24, unit = "cm" )

####(FIGURE X-2) ----
habitat_summary_habitat <-habitat_toplot%>%
  group_by(Area, Group)%>%
  summarise(Percent_cover_habitat=mean(percent_cover,na.rm=T),
            stdev=sd(percent_cover))

habitat_plot_hab <- ggplot()+
  geom_jitter(data=habitat_toplot, aes(x= reorder_within(Group, -percent_cover, Area, fun=mean), #note here that I included -Percent_cover to reorder descending
                                       y = percent_cover, colour=Group), width=0.15, alpha=0.4)+
  facet_wrap(.~Area, scales="free")+
  geom_errorbar(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Area, fun=mean), y = Percent_cover_habitat, ymin=Percent_cover_habitat-stdev,ymax=Percent_cover_habitat+stdev),fill = "black",width=.15)+
  geom_point(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Area, fun=mean), y = Percent_cover_habitat), size=4.5)+
  geom_point(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Area, fun=mean), y = Percent_cover_habitat, colour=Group), size=4)+
  theme_bw()+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

ggsave(habitat_plot_hab, file=here::here("outputs/", "habitat_latitudinal.png"),
       height = 16, width = 24, unit = "cm" )

####(FIGURE X-3) ----

habitat_v3<- ggplot(habitat_summary_habitat, aes(x=Area, y=Percent_cover_habitat, fill=Group)) +
  geom_bar(stat="identity", 
           position="stack",
           colour="black") + 
  theme_bw()+
geom_errorbar(aes(ymin=Percent_cover_habitat-stdev,ymax=Percent_cover_habitat+stdev), 
               size=0.5, width=.25, stat = "identity",
               position = "identity")

## Error bars not working

ggsave(habitat_v3, file=here::here("outputs/", "habitat_stack_latitudinal.png"),
       height = 16, width = 24, unit = "cm" )


###################### EXTRAS

habitat_plot_hab <- ggplot(habitat_toplot, aes(x=Group, y=percent_cover, fill=Group))+
  geom_boxplot()+
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0)) +
  facet_wrap(~ Habitat) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #panel.grid.minor  = element_blank(),
    #panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

habitat_plot_hab <- ggplot()+
  geom_jitter(data=habitat_toplot, aes(x= reorder_within(Group, -percent_cover, Habitat, fun=mean), #note here that I included -Percent_cover to reorder descending
                                       y = percent_cover, colour=Group), width=0.15, alpha=0.4)+
  facet_wrap(.~Habitat, scales="free")+
  geom_errorbar(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Habitat, fun=mean), y = Percent_cover_habitat, ymin=Percent_cover_habitat-stdev,ymax=Percent_cover_habitat+stdev),fill = "black",width=.15)+
  geom_point(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Habitat, fun=mean), y = Percent_cover_habitat), size=4.5)+
  geom_point(data=habitat_summary_habitat, aes(x= reorder_within(Group, -Percent_cover_habitat, Habitat, fun=mean), y = Percent_cover_habitat, colour=Group), size=4)+
  theme_bw()+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

ggsave(habitat_plot_hab, file=here::here("outputs/", "habitat_cross-shelf.png"),
       height = 16, width = 24, unit = "cm" )


####(FIGURE X-1) ----
habitat_summary <-habitat_toplot%>%
  group_by(Site, Habitat, Group)%>%
  summarise(Percent_cover_site=mean(percent_cover,na.rm=T),
            stdev=sd(percent_cover))

habitat_plot <- ggplot()+
  geom_jitter(data=habitat_toplot, aes(x= reorder_within(Group, -percent_cover, Site, fun=mean), #note here that I included -Percent_cover to reorder descending
                                       y = percent_cover, colour=Group), width=0.15, alpha=0.4)+
  facet_wrap(.~Site, scales="free", nrow=3)+
  geom_errorbar(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site, ymin=Percent_cover_site-stdev,ymax=Percent_cover_site+stdev),fill = "black",width=.15)+
  geom_point(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site), size=4.5)+
  geom_point(data=habitat_summary, aes(x= reorder_within(Group, -Percent_cover_site, Site, fun=mean), y = Percent_cover_site, colour=Group), size=4)+
  theme_bw()+
  theme(#axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom",
    axis.text.y = element_text(size=14),
    axis.text.x = element_blank())+
  xlab("")

ggsave(habitat_plot, file=here::here("outputs/", "habitat_sites.png"),
       height = 16, width = 24, unit = "cm" )