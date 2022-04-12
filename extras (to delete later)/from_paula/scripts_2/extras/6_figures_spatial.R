##############################################################
### Plot Functional Richness


## retrieve names of main input:
asb_fd_ind <- fd_ind_values
fd_details <- alpha_fd_indices$details


### Prepare data for plotting:

## get coordinates of species:
sp_faxes_coord <- fd_details$sp_faxes_coord

## set type, resolution and dimensions of file if to be saved:
device_file <- "png"
res_file <- 300
height_file <- 8 * c(2, 2, 3)
names(height_file) <- c("1", "3", "6")
width_file  <- 10 * c(2, 2, 3)
names(width_file) <- c("1", "3", "6")

## get number of dimensions in input:
nb_dim <- ncol(sp_faxes_coord)

## define arguments:
faxes               = NULL
faxes_nm            = NULL
range_faxes         = c(NA, NA)
plot_asb_nm <- c("Inshore", "Midshelf", "Offshore")
plot_sp_nm <- NULL



## define arguments values and prepare data for plotting:

# give faxes identity if faxes set to NULL:
if (is.null(faxes)) {
  faxes <- colnames(sp_faxes_coord)[1:min(c(4, nb_dim))]
}

# give faxes names if faxes set to NULL:
if (is.null(faxes_nm)) {
  faxes_nm <- faxes
}
names(faxes_nm) <- faxes

# get number of axes:
nb_faxes <- length(faxes)

# get combinations of axes on plot:
axes_plot <- utils::combn(faxes, 2)
plot_nb   <- ncol(axes_plot)

# set range of axes if c(NA, NA):
if (is.na(range_faxes[1]) && is.na(range_faxes[2])) {
  range_sp_coord  <- range(sp_faxes_coord)
  range_faxes <- range_sp_coord +
    c(-1, 1) * (range_sp_coord[2] - range_sp_coord[1]) * 0.1
}

# create a dataframe with species coordinates and option (vertices + label)...
# ... if required:
sp_faxes_coord_plot <- data.frame(sp_faxes_coord, label = "")

# if some species names to be plotted, adding a character variable to ...
# ... sp_faxes_coord:
if (! is.null(plot_sp_nm)) {
  sp_faxes_coord_plot[plot_sp_nm, "label"] <- plot_sp_nm
}

# get vertices of the convex hull of the species pool:
vert_pool <- fd_details$pool_vert_nm

# retrieve names and weights of species present in each assemblage:

# get names of assemblages:
pool <- "pool"
asb1 <- plot_asb_nm[1]
nm_asb <- asb1
asb2 <- plot_asb_nm[2]
nm_asb <- paste(nm_asb, asb2, sep = "_")
asb3 <- plot_asb_nm[3]
nm_asb <- paste(nm_asb, asb3, sep = "_")


sp_asb1 <- names(which(fd_details$asb_sp_occ[asb1, ] == 1))
sp_asb2 <- names(which(fd_details$asb_sp_occ[asb2, ] == 1))
sp_asb3 <- names(which(fd_details$asb_sp_occ[asb3, ] == 1))



### Plot FRic:

# list to store ggplot
panels_fric <- list()

# loop on combinations:
for (k in (1:plot_nb)) {
  
  # names of axes
  xy_k <- axes_plot[1:2, k]
  
  # get species coordinates along the 2 axes:
  sp_coord_xy <- as.matrix(sp_faxes_coord_plot[, xy_k])
  colnames(sp_coord_xy) <- c("x", "y")
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_coord_xy[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- fd_details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_coord_xy[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- fd_details$asb_vert_nm[[asb2]]
  asb_sp_coord2D_k[["asb3"]] <- sp_coord_xy[sp_asb3, ]
  vertices_nD_k[["asb3"]] <- fd_details$asb_vert_nm[[asb3]]
  
  # background = axes defined by range of values and names as specified:
  plot_k <- mFD::background.plot(range_faxes, faxes_nm = xy_k, color_bg = "grey95")
  
  # add species pool:
  plot_k <- mFD::pool.plot(ggplot_bg = plot_k,
                           sp_coord2D = sp_coord_xy,
                           vertices_nD = vert_pool,
                           plot_pool = TRUE,
                           color_ch = NA,
                           fill_ch = "white",
                           alpha_ch = 1,
                           shape_pool = 3,
                           size_pool = 0.7,
                           color_pool = "grey50",
                           fill_pool = NA,
                           shape_vert = 3,
                           size_vert = 0.7,
                           color_vert = "grey50",
                           fill_vert = NA)
  
  # plot 2D convex hulls and points for the 2 assemblages:
  plot_k <- mFD::fric.plot(ggplot_bg = plot_k,
                           asb_sp_coord2D = asb_sp_coord2D_k,
                           asb_vertices_nD = vertices_nD_k,
                           plot_sp = TRUE,
                           color_ch = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
                                        asb3 = "#DCE319FF"),
                           fill_ch = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
                                       asb3 = "#DCE319FF"),
                           alpha_ch = c(asb1 = 0.5, asb2 = 0.3, asb3 = 0.2),
                           shape_sp = c(asb1 = 15, asb2 = 17, asb3 = 16),
                           size_sp = c(asb1 = 1, asb2 = 1, asb3 = 1),
                           color_sp = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
                                        asb3 = "#DCE319FF"),
                           fill_sp = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
                                       asb3 = "#DCE319FF"),
                           shape_vert = c(asb1 = 15, asb2 = 17, asb3 = 16),
                           size_vert = c(asb1 = 1, asb2 = 1, asb3 = 1),
                           color_vert = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
                                          asb3 = "#DCE319FF"),
                           fill_vert = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
                                         asb3 = "#DCE319FF"))
  
  # add species names if needed:
  if (! is.null(plot_sp_nm)) {
    x <- NULL
    y <- NULL
    plot_k <- plot_k +
      ggrepel::geom_text_repel(data = sp_faxes_coord_plot,
                               ggplot2::aes_string(x = xy_k[1],
                                                   y = xy_k[2],
                                                   label = "label"),
                               size = size_sp_nm, colour= color_sp_nm,
                               fontface = fontface_sp_nm,
                               max.overlaps = Inf,
                               box.padding = grid::unit(2, 'lines'),
                               force = 5,
                               arrow = grid::arrow(length = grid::unit(0.02,
                                                                       'npc')),
                               segment.color = color_sp_nm)
  }
  
  # save plot in a list:
  panels_fric[[k]] <- plot_k
  
} # end of k



# create a caption to add to the patchwork object (legend):

#retrieve values to plot:
# top_fric <- c("Functional richness", asb1, "")
# values_fric <- c(round(asb_fd_ind[asb1, "fric"], 3), "")
# top_fric[3] <- asb2
# 
# values_fric[2] <- round(asb_fd_ind[asb2,"fric"], 3)
# top_fric[4] <- asb3
# values_fric[3] <- round(asb_fd_ind[asb3, "fric"], 3)

# customize position of texts in the plot
# spread_faxes <- (range_faxes[2] - range_faxes[1])
# hh <- c(1, 3, 4, 5, 6)
# vv <- 0.3

# plot window:
#x <- NULL
#y <- NULL
plot_caption <- ggplot2::ggplot(data.frame(x = range_faxes,                                            y = range_faxes),
                                ggplot2::aes(x = x, y = y)) +
  ggplot2::scale_x_continuous(limits = range_faxes, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = range_faxes, expand = c(0, 0)) +
  ggplot2::theme_void() + 
  ggplot2::theme(legend.position = "none") #+
#ggplot2::geom_rect(xmin = range_faxes[1], xmax = range_faxes[2],
#ymin = range_faxes[1], ymax = range_faxes[2],
#fill = "white", colour ="black")

# plot names of index and of assemblages:
# h   <- NULL
# v   <- NULL
# top <- NULL
# x <- NULL
# y <- NULL
#plot_caption <- plot_caption +
#   ggplot2::geom_text(data = data.frame(
#     h = range_faxes[1] + spread_faxes * 0.15 * hh[c(1,3:5)],
#      v = range_faxes[2] - spread_faxes * rep(0.2, 4),
#      top = top_fric),
#     ggplot2::aes(x = h, y = v, label = top),
#      size = 3, hjust = 0.5, fontface = "bold")

# plot FRic values:
# values_lab <- NULL
# data_caption <- data.frame(
#   h = range_faxes[1] + spread_faxes * 0.15 * hh[2:5],
#   v = range_faxes[2] - spread_faxes*rep(vv, 4),
#   values_lab = c("FRic", values_fric))
# plot_caption <- plot_caption +
#   ggplot2::geom_text(data = data_caption,
#                      ggplot2::aes(x = h, y = v, label = values_lab),
#                      size = 3, hjust = 0.5, fontface = "plain")
# vv <- vv + 0.1

# add text about dimensionality:
# nb <- NULL
# plot_caption <- plot_caption +
#   ggplot2::geom_text(data = data.frame(
#     h = range_faxes[1] + spread_faxes * 0.1,
#     v = range_faxes[2] - spread_faxes * vv,
#     nb = paste0("NB: Indices were computed in a ",
#                 nb_dim,"-dimensional space")),
#     ggplot2::aes(x = h, y = v, label = nb),
#     size = 3, hjust = 0, fontface = "italic")

# add legend (convex hull, asb species and pool species):

## plot legend:
# values_lab <- NULL

### retrieve arguments values:
# shape_sp = c(asb1 = 16, asb2 = 16, asb3 = 16)
# color_sp = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
#              asb3 = "#DCE319FF")
# fill_sp = c(asb1 = "#440154FF", asb2 = "#1F968BFF",
#             asb3 = "#DCE319FF")
# size_sp = c(asb1 = 1, asb2 = 1, asb3 = 1)
# alpha_ch = c(asb1 = 0.5, asb2 = 0.3, asb3 = 0.2)

### for 1st asb:
# plot_caption <- plot_caption +
#   ggplot2::geom_rect(xmin = range_faxes[1] + spread_faxes*0.10,
#                      xmax = range_faxes[1] + spread_faxes*0.15,
#                      ymin = range_faxes[2] - spread_faxes*0.51,
#                      ymax = range_faxes[2] - spread_faxes*0.55,
#                      fill = color_sp[["asb1"]], alpha = alpha_ch[["asb1"]]) + 
#   
#   ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
#                      y = range_faxes[2] - spread_faxes*0.525,
#                      label = paste0("convex hull of", sep = " ",
#                                     asb1),
#                      colour = color_sp[["asb1"]], size = 3) +
#   # 
#   ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
#                       y = range_faxes[2] - spread_faxes*0.58,
#                       size = size_sp[["asb1"]], shape = shape_sp[["asb1"]],
#                       color = color_sp[["asb1"]], fill = fill_sp[["asb1"]]) +
#  
#   
#     ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
#                      y = range_faxes[2] - spread_faxes*0.58,
#                      label = paste0("shape of species from", sep = " ",
#                                     asb1),
#                      colour = color_sp[["asb1"]], size = 3)

### for 2nd assemblage:

# plot_caption <- plot_caption +
#   ggplot2::geom_rect(xmin = range_faxes[1] + spread_faxes*0.10,
#                      xmax = range_faxes[1] + spread_faxes*0.15,
#                      ymin = range_faxes[2] - spread_faxes*0.64,
#                      ymax = range_faxes[2] - spread_faxes*0.68,
#                      fill = color_sp[["asb2"]], 
#                      alpha = alpha_ch[["asb2"]]) + 
#   
#   ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
#                      y = range_faxes[2] - spread_faxes*0.665,
#                      label = paste0("convex hull of", sep = " ", 
#                                     asb2),
#                      colour = color_sp[["asb2"]], size = 3) + 
#   
#   ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
#                       y = range_faxes[2] - spread_faxes*0.71,
#                       size = size_sp[["asb2"]], 
#                       shape = shape_sp[["asb2"]],
#                       color = color_sp[["asb2"]], 
#                       fill = fill_sp[["asb2"]]) + 
#   
#   ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
#                      y = range_faxes[2] - spread_faxes*0.71,
#                      label = paste0("shape of species from", sep = " ", 
#                                     asb2),
#                      colour = color_sp[["asb2"]], size = 3)

# for 3rd asb:
# plot_caption <- plot_caption +
#   ggplot2::geom_rect(xmin = range_faxes[1] + spread_faxes*0.10,
#                      xmax = range_faxes[1] + spread_faxes*0.15,
#                      ymin = range_faxes[2] - spread_faxes*0.77,
#                      ymax = range_faxes[2] - spread_faxes*0.81,
#                      fill = color_sp[["asb3"]], 
#                      alpha = alpha_ch[["asb3"]]) + 
#   
#   ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
#                      y = range_faxes[2] - spread_faxes*0.795,
#                      label = paste0("convex hull of", sep = " ", 
#                                     asb3),
#                      colour = color_sp[["asb3"]], size = 3) + 
#   
#   ggplot2::geom_point(x = range_faxes[1] + spread_faxes*0.125,
#                       y = range_faxes[2] - spread_faxes*0.84,
#                       size = size_sp[["asb3"]], 
#                       shape = shape_sp[["asb3"]],
#                       color = color_sp[["asb3"]], 
#                       fill = fill_sp[["asb3"]]) + 
#   
#   ggplot2::geom_text(x = range_faxes[1] + spread_faxes*0.35,
#                      y = range_faxes[2] - spread_faxes*0.84,
#                      label = paste0("shape of species from", sep = " ", 
#                                     asb3),
#                      colour = color_sp[["asb3"]], size = 3)

# create patchwork object:
patchwork_fric <- mFD::panels.to.patchwork(panels_fric,plot_caption) 

# title and caption
# tit_fric <- paste0("Functional Richness of ", asb1, ", ", asb2,
#                    " and ", asb3)

# create final patchwork object: 
patchwork_fric <- patchwork_fric 
#patchwork::plot_annotation(title = tit_fric)
#,caption = "made with mFD package")

# save plot:
file_fric <- paste0(nm_asb, "_", "FRic_", nb_dim, "D" , ".", device_file)
ggplot2::ggsave(filename = file_fric ,
                plot = patchwork_fric,
                device = device_file,
                scale = 1,
                height= height_file[as.character(plot_nb)],
                width = width_file[as.character(plot_nb)],
                units= "in",
                dpi = res_file)

# save for PC1 and PC2 only:
patchwork_fric_PC1_PC2 <- mFD::panels.to.patchwork(list(panels_fric[[1]]), plot_caption)

# title and caption
# tit_fric <- paste0("Functional Richness of ", asb1, ", ", asb2,
# " and ", asb3)

# create final patchwork object: 
patchwork_fric_PC1_PC2 <- patchwork_fric_PC1_PC2 +
  patchwork::plot_annotation()
# (title = tit_fric,caption = "made with mFD package")


plot(patchwork_fric_PC1_PC2)

## [PS] How to have the convex hull in the centre. Need help
# on this to make it "more pretty" 

#### CAP graph of FRic

# get scores

fric <- as.data.frame(scores(fric_dbrda_2, display = "sites")) 

rownames(fric) <- meta_data_3$Code

CAP1 <- fric$CAP1
CAP2 <- fric$CAP2

# add metadata
identical(as.character(meta_data_3$Code), rownames(fric)) # verify that data in same order
fric_habitat <- cbind(fric,select(meta_data_3, Habitat))


fric <- ggplot(fric_habitat, aes(x= CAP1, y= CAP2, color = Habitat)) + 
  stat_ellipse(aes(fill = Habitat), geom = "polygon", alpha = 0.2) +
  geom_point() +
  theme_classic() + 
  coord_cartesian(xlim=c(-10, 10), ylim=c(-5, 5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(color = "Habitat", fill = "Habitat") +
  ggtitle("Capscale Analysis - Functional Richness") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("CAP2") + 
  xlab("CAP1")

plot(fric)

#Graph CAP FOR tax beta diversity

# get scores

xy <- as.data.frame(scores(tax_diss_dbrda_transect, display = "sites")) 
CAP1 <- xy$CAP1
CAP2 <- xy$CAP2

# add metadata
identical(as.character(meta_data_for_functional$Code), rownames(xy)) # verify that data in same order
xy_habitat <- cbind(xy,select(meta_data_for_functional, Habitat))


tax_diss <- ggplot(xy_habitat, aes(x= CAP1, y= CAP2, color = Habitat)) + 
  stat_ellipse(aes(fill = Habitat), geom = "polygon", alpha = 0.2) +
  geom_point() +
  theme_classic() + 
  coord_cartesian(xlim=c(-5, 5), ylim=c(-5, 5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(color = "Habitat", fill = "Habitat") +
  ggtitle("Taxonomic beta diversity (Dissimilarity)") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("CAP2") + 
  xlab("CAP1")

plot(tax_diss)

#Graph CAP for tax beta diversity (turnover)

# get scores

turn <- as.data.frame(scores(tax_turn_dbrda_transect, display = "sites")) 
CAP1 <- turn$CAP1
CAP2 <- turn$CAP2

# add metadata
identical(as.character(meta_data_for_functional$Code), rownames(turn)) # verify that data in same order
turn_habitat <- cbind(turn,select(meta_data_for_functional, Habitat))


tax_turn <- ggplot(turn_habitat, aes(x= CAP1, y= CAP2, color = Habitat)) + 
  stat_ellipse(aes(fill = Habitat), geom = "polygon", alpha = 0.2) +
  geom_point() +
  theme_classic() + 
  coord_cartesian(xlim=c(-5, 5), ylim=c(-5, 5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(color = "Habitat", fill = "Habitat") +
  ggtitle("Taxonomic beta diversity (Turnover)") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("CAP2") + 
  xlab("CAP1")

tax_turn


#Graph CAP for functiona beta diversity

# get scores

func_diss <- as.data.frame(scores(func_diss_dbrda_transect, display = "sites")) 
CAP1 <- func_diss$CAP1
CAP2 <- func_diss$CAP2

# add metadata
identical(as.character(meta_data_for_functional$Code), rownames(func_diss)) # verify that data in same order
func_diss_habitat <- cbind(func_diss,select(meta_data_for_functional, Habitat))


func_diss <- ggplot(func_diss_habitat, aes(x= CAP1, y= CAP2, color = Habitat)) + 
  stat_ellipse(aes(fill = Habitat), geom = "polygon", alpha = 0.2) +
  geom_point() +
  theme_classic() + 
  coord_cartesian(xlim=c(-5, 5), ylim=c(-5, 5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(color = "Habitat", fill = "Habitat") +
  ggtitle("Functional beta diversity (Dissimilarity)") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("CAP2") + 
  xlab("CAP1")

plot(func_diss)

#### GRAPH

# get scores

func_turn <- as.data.frame(scores(func_turn_dbrda_transect, display = "sites")) 
CAP1 <- func_turn$CAP1
CAP2 <- func_turn$CAP2

# add metadata
identical(as.character(meta_data_for_functional$Code), rownames(func_turn)) # verify that data in same order
func_turn_habitat <- cbind(func_turn,select(meta_data_for_functional, Habitat))


func_turn <- ggplot(func_turn_habitat, aes(x= CAP1, y= CAP2, color = Habitat)) + 
  stat_ellipse(aes(fill = Habitat), geom = "polygon", alpha = 0.2) +
  geom_point() +
  theme_classic() + 
  coord_cartesian(xlim=c(-5, 5), ylim=c(-5, 5)) +
  geom_hline(yintercept = 0, linetype="dotted") + 
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(color = "Habitat", fill = "Habitat") +
  ggtitle("Functional beta diversity (Turnover)") +
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("CAP2") + 
  xlab("CAP1")

plot(func_turn)

##############################################################################
##               FIGURES
##############################################################################

# set general plotting parameters

ggplot2::theme_set(theme_classic())

# Color scale for habitat

colins    <-  "mediumseagreen"
colmid   <-  "lightsalmon1"
coloff  <-  "firebrick3"


#### Diversity indices ####

# Raw richness
gtot <- ggplot(fd_values_site, aes(habitat, sp_richn, fill=habitat)) +
  geom_boxplot(colour="black", width = 0.3) + # set boxplot colour by protection level
  scale_fill_manual(values = c(colins, colmid, coloff)) + # set legend items order
  scale_x_discrete(labels = c("Inshore", "Midshelf", "Offshore")) +
  labs(y = "Species richness", x = "") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") + # centralize title and hide legend
  coord_cartesian(ylim=c(0,30))

mytheme <-theme(panel.background=element_rect(fill="white"), panel.grid.minor = element_blank(), 
                axis.ticks = element_blank(), 
                panel.grid.major = element_blank(),axis.line = element_line(size = 1, colour = "black"),
                axis.text = element_text(size = (14), colour="black"), axis.title = element_text(size= (18)),
                legend.key = element_rect(fill = "white"))

print (gtot + mytheme)


