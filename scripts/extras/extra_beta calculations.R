

### Extra - beta calculations


## computing taxonomic and functional beta-diversity ####

# taxonomic dissimilarity = Jaccard index and its components
spatial_beta_taxo <- betapart::beta.pair(spatial_sp_occ, index.family = "jaccard")

# functional dissimilarity = Jaccard-like index and its components
spatial_beta_func <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_occ       = spatial_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
spatial_beta <- list (
  taxo_diss = spatial_beta_taxo$beta.jac,
  taxo_turn = spatial_beta_taxo$beta.jtu,
  func_diss = spatial_beta_func$pairasb_fbd_indices$jac_diss,
  func_turn = spatial_beta_func$pairasb_fbd_indices$jac_turn
)
spatial_beta$taxo_diss
# summary
cbind( min=lapply(spatial_beta, min), max=lapply(spatial_beta, max) )



## computing taxonomic and functional beta-diversity ####

# taxonomic dissimilarity = Jaccard index and its components ----
temporal_beta_taxo_nokelp <- betapart::beta.pair(nokelp_sp_occ, index.family = "jaccard")
temporal_beta_taxo_kelp <- betapart::beta.pair(kelp_sp_occ, index.family = "jaccard")

kelp_turnover <- temporal_beta_taxo_kelp$beta.jtu

write.csv(kelp_turnover, file=here::here("data", "raw_data", "kelp_turnover.csv"), 
          row.names = FALSE )

# functional beta no kelp sites ----
# functional dissimilarity = Jaccard-like index and its components
temporal_beta_func_nokelp <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_occ       = nokelp_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
temporal_beta_nokelp <- list (
  taxo_diss = temporal_beta_taxo_nokelp$beta.jac,
  taxo_turn = temporal_beta_taxo_nokelp$beta.jtu,
  func_diss = temporal_beta_func_nokelp$pairasb_fbd_indices$jac_diss,
  func_turn = temporal_beta_func_nokelp$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(temporal_beta_nokelp, min), max=lapply(temporal_beta_nokelp, max) )


# functional beta kelp sites ----
temporal_beta_func_kelp <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_occ       = kelp_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
temporal_beta_kelp <- list (
  taxo_diss = temporal_beta_taxo_kelp$beta.jac,
  taxo_turn = temporal_beta_taxo_kelp$beta.jtu,
  func_diss = temporal_beta_func_kelp$pairasb_fbd_indices$jac_diss,
  func_turn = temporal_beta_func_kelp$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(temporal_beta_kelp, min), max=lapply(temporal_beta_kelp, max) )


save(temporal_beta_nokelp, file=here::here("outputs/", "temporal_beta_nokelp.RData") )
save(temporal_beta_kelp, file=here::here("outputs/", "temporal_beta_kelp.RData") )




















save(spatial_beta, file=here::here("outputs/", "spatial_beta.RData") )