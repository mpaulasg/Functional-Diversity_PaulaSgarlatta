
#Loading data from both spatial/temporal

species_both <- read.csv(here::here("data", "species_both.csv")) %>% 
  mutate(type_data= if_else(data_1 == "temporal" & data_2 == "spatial", "both",
                            if_else(data_1 == "temporal" & data_2 == "no", "temporal","spatial"))) %>% 
  select(-data_1, -data_2)


# Add thermal aff to sp_faxes_coord
sp_3D_coord <- spatial_fd$details$sp_faxes_coord

sp_faxes_coord <- as.data.frame(sp_3D_coord) %>% 
  rownames_to_column("Species")

sp_faxes_coord <- inner_join(sp_faxes_coord, species_both, 
                             by="Species") 

# Change class of data type column: character to factor

sp_faxes_coord$type_data <- as.factor(sp_faxes_coord$type_data)
class(sp_faxes_coord$type_data) # ok

# Create a new assemblage*species df with assemblages being either spatial, temporal or both

asb_sp <- sp_faxes_coord[, c(1, 5)]

asb_sp_new <- asb_sp %>% 
  add_column(present = as.numeric(1)) %>% 
  pivot_wider(names_from = type_data, values_from = present)

asb_sp_new[is.na(asb_sp_new)] <- as.numeric(0)

asb_sp_new2 <- t(asb_sp_new)

colnames(asb_sp_new2) <- asb_sp_new2[1, ]

asb_sp_new2 <- asb_sp_new2[-1, ]


# create a dataframe that will contain the same values as asb_sp_new2 because (I don't knwo why) 
# I cannot convert character into numeric

asb_sp_new3 <- as.data.frame(matrix(nrow = nrow(asb_sp_new2), ncol = ncol(asb_sp_new2)))#(s)
colnames(asb_sp_new3) <- colnames(asb_sp_new2)
rownames(asb_sp_new3) <- rownames(asb_sp_new2)

for (i in (1:nrow(asb_sp_new3))) {
  for (j in (1:ncol(asb_sp_new3))) {
    asb_sp_new3[i, j] <- as.numeric(asb_sp_new2[i, j])
  }
}

asb_sp_new3 <- as.matrix(asb_sp_new3)


## settings ####

# vertices of all fe in 4D ----
pool_vert_nm <- spatial_fd$details$pool_vert_nm

# range of axes
range_faxes_coord <- range(sp_3D_coord[,1:3])
range_axes <- range_faxes_coord +
  c(-1, 1) * (range_faxes_coord[2] - range_faxes_coord[1]) * 0.1
spread_faxes <- range_axes[2] - range_axes[1]


## Compute FRic values #### 

# compute FRic for all habitats  ---
Fric <- alpha.fd.multidim(sp_faxes_coord = sp_3D_coord, 
                          asb_sp_w = asb_sp_new3,
                          ind_vect = c("fric"), 
                          scaling = TRUE, 
                          details_returned = TRUE)

# color code for data type

data_colors <- c(spatial = "lightsalmon1", temporal = "seagreen4", both = "#2C6BAA")


## plotting  ####

# list to store ggplot
ggplot_pc <- list()

# pairs of axes
pairs_axes <- list(c(1,2), c(1,3), c(2, 3))

for (z in 1:length(pairs_axes)) {
  
  # names of axes   
  xy <- pairs_axes[[z]]
  
  # background with axes range set + title
  ggplot_z <- background.plot(range_faxes = range_axes,
                              faxes_nm = paste0("PC", xy), 
                              color_bg = "grey95")
  
  
  # convex hull of global species pool
  ggplot_z <- pool.plot(ggplot_bg = ggplot_z,
                        sp_coord2D = sp_3D_coord[,xy],
                        vertices_nD = pool_vert_nm,
                        plot_pool = FALSE,
                        color_ch = "black", fill_ch = "white", alpha_ch = 1)
  
  # get names of assemblages:
  pool <- "pool"
  plot_asb_nm <- c("spatial", "temporal", "both")
  asb1 <- plot_asb_nm[1]
  nm_asb <- asb1
  asb2 <- plot_asb_nm[2]
  nm_asb <- paste(nm_asb, asb2, sep = "_")
  asb3 <- plot_asb_nm[3]
  nm_asb <- paste(nm_asb, asb3, sep = "_")
  
  
  sp_asb1 <- names(which(Fric$details$asb_sp_occ[asb1, ] == 1))
  sp_asb2 <- names(which(Fric$details$asb_sp_occ[asb2, ] == 1))
  sp_asb3 <- names(which(Fric$details$asb_sp_occ[asb3, ] == 1))
  
  # list with dataframes for plot:
  asb_sp_coord2D_k <- list()
  asb_sp_coord2D_k[["asb1"]] <- sp_3D_coord[sp_asb1, ]
  vertices_nD_k <- list()
  vertices_nD_k[["asb1"]] <- Fric$details$asb_vert_nm[[asb1]]
  asb_sp_coord2D_k[["asb2"]] <- sp_3D_coord[sp_asb2, ]
  vertices_nD_k[["asb2"]] <- Fric$details$asb_vert_nm[[asb2]]
  asb_sp_coord2D_k[["asb3"]] <- sp_3D_coord[sp_asb3, ]
  vertices_nD_k[["asb3"]] <- Fric$details$asb_vert_nm[[asb3]]
  
  # species present in space:
  sp_spatial <- sp_faxes_coord$Species[which(sp_faxes_coord$data_type == "spatial")]
  
  # species present in temporal:
  sp_temporal <- sp_faxes_coord$Species[which(sp_faxes_coord$data_type == "temporal")]
  
  # species present in both:
  sp_both <- sp_faxes_coord$Species[which(sp_faxes_coord$data_type == "both")]
  
  
  
  # vertices in spatial:
  vert_spatial <- Fric$details$asb_vert_nm$spatial
  
  # vertices in temporal:
  vert_temporal <- Fric$details$asb_vert_nm$temporal
  
  # vertices in both:
  vert_both <- Fric$details$asb_vert_nm$both
  
  
  # plot convex hull of assemblage but not species
  ggplot_z2 <-fric.plot(ggplot_bg = ggplot_z, 
                        asb_sp_coord2D = list(asb1 = sp_3D_coord[sp_spatial, xy], 
                                              asb2 = sp_3D_coord[sp_temporal, xy],
                                              asb3 = sp_3D_coord[sp_both, xy]),
                        asb_vertices_nD = list(asb1 = vert_spatial, 
                                               asb2 = vert_temporal,
                                               asb3=vert_both),
                        plot_sp = TRUE,
                        color_sp = data_colors,
                        fill_sp = c(asb1 = "white", asb2 = "white", asb3="white"),
                        size_sp = c(asb1 = 3, asb2 = 3, asb3=3),
                        shape_sp = c(asb1 = 16, asb2 = 16, asb=16),
                        color_vert = data_colors,
                        fill_vert = data_colors,
                        size_vert = c(asb1 = 3, asb2 = 3, asb3=3),
                        shape_vert = c(asb1 = 16, asb2 = 16, asb=16),
                        alpha_ch = c(asb1 = 0, asb2 = 0, asb3=0),
                        color_ch = c(asb1 = NA, asb2 = NA, asb3=NA),
                        fill_ch = c(asb1 = NA, asb2 = NA, asb3=NA))
  
  
  # ggplot stored in list
  ggplot_pc[[z]] <- ggplot_z2
  
  
}# end of z

## Compute Caption:

# plot white basic window:
plot_caption <- ggplot2::ggplot(data.frame(x = range_axes, 
                                           y = range_axes),
                                ggplot2::aes(x = x, y = y)) +
  ggplot2::scale_x_continuous(limits = range_axes, expand = c(0, 0)) +
  ggplot2::scale_y_continuous(limits = range_axes, expand = c(0, 0)) +
  ggplot2::theme_void() + ggplot2::theme(legend.position = "none") +
  ggplot2::geom_rect(xmin = range_axes[1], xmax = range_axes[2],
                     ymin = range_axes[1], ymax = range_axes[2],
                     fill = "white", colour ="white")


## merging all plots into a single figure and saving as png ####
figure <- panels.to.patchwork(ggplot_pc, plot_caption = plot_caption)

figure3 <- ggplot_pc [[1]]


ggsave(figure3, file=here::here("outputs",  "Figure3.jpeg"),
       height = 16, width = 30, unit = "cm" )