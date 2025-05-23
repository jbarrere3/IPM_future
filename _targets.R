#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "readxl", "cowplot",
                 "data.table", "factoextra", "terra", "ggmcmc", "R2jags", 
                 "betareg", "car", "scales", "MASS", "broom.mixed", "lme4", 
                 "modi", "ggridges", "purrr", "checkmate", "FD", "sf", 
                 "rnaturalearth", "rnaturalearthdata", "sinkr", "egg", "xtable", 
                 "spData", "ggnewscale")
for(i in 1:length(packages.in)){
  if(!(packages.in[i] %in% rownames(installed.packages()))){
    install.packages(packages.in[i])
  }
}  
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess", 
       future.globals.maxSize= 1048576000)
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 16)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load and format NFI and climate data ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # NFI data
  # -- File name
  tar_target(NFI_data_file, "data/NFI/FUNDIV_ba.csv", format = "file"), 
  # -- Read the file
  tar_target(NFI_data, fread(NFI_data_file)), 
  
  # Disturbance regime per NFI plot
  # -- File name
  tar_target(NFI_disturbance_file, "data/Disturbance/NFI_disturbance.csv", 
             format = "file"), 
  # -- Read the file
  tar_target(NFI_disturbance, fread(NFI_disturbance_file)), 
  
  # Climate files
  # -- File names
  tar_target(climate_files, list.files(
    "data/Climate", full.names = TRUE, recursive = TRUE), format = "file"), 
  # -- Extract climate for each NFI plot
  tar_target(NFI_climate, extract_climate_NFI(climate_files, NFI_disturbance)), 
  
  # Filter data based on species present in IPM and climate
  # -- select NFI plots
  tar_target(NFI_plots_selected, select_NFI_plots(
    NFI_data, NFI_climate, nclim = 10, nplots.per.clim = 1000)), 
  # -- subset all data
  tar_target(NFI_data_sub, subset(
    NFI_data, plotcode %in% NFI_plots_selected$plotcode)),
  tar_target(NFI_disturbance_sub, subset(
    NFI_disturbance, plotcode %in% NFI_plots_selected$plotcode)),
  tar_target(NFI_climate_sub, subset(
    NFI_climate, plotcode %in% NFI_plots_selected$plotcode)),
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load and format recruitment data ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Coefficients to model regional basal area
  # -- Raw file
  tar_target(coef_ba_reg_file, "data/Recruitment/coef_reg_ba.csv", format = "file"), 
  # -- Read the file
  tar_target(coef_ba_reg, fread(coef_ba_reg_file)), 
  
  # Climate pca (input to calculate regional basal area per plot)
  # -- Raw file
  tar_target(clim_pca_file, "data/Recruitment/clim_pca.rds", format = "file"), 
  # -- read rds file
  tar_target(clim_pca, readRDS(clim_pca_file)),
  
  # Dispersal kernel per species
  # -- Raw file
  tar_target(disp_kernel_file, "data/Recruitment/disp_kernel.csv", format = "file"), 
  # -- Read the file
  tar_target(disp_kernel, fread(disp_kernel_file)), 
  
  # New demographic parameters
  # -- Raw file
  tar_target(new_fit_list_file, "data/Recruitment/new_fit_list.rds", format = "file"), 
  # -- Load the rds file
  tar_target(new_fit_list, readRDS(new_fit_list_file)),
  
  # Forest cover per NFI plot
  # -- Raw file
  tar_target(NFI_forest_cover_file, 
             "data/Recruitment/FUNDIV_coord_landcover_1000.csv", format = "file"), 
  # -- Format raw file
  tar_target(NFI_forest_cover, get_forest_cover(NFI_forest_cover_file)),
  
  
  
  

  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Get distributions of storm and fire intensity ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # File containing disturbance intensity value
  tar_target(intensity_file, "data/Disturbance/jags_dominance.Rdata", 
             format = "file"), 
  
  # Extract intensity per disturbance as a quantile of posterior distribution
  tar_target(I_per_disturbance, get_I_per_disturbance(intensity_file, 0.85)),
  
  # Get the climate and disturbance dataframe for each plot and scenario
  tar_target(climate_dist_dflist, get_clim_dist_df(
    NFI_disturbance_sub, NFI_climate_sub, I_per_disturbance)),
  
  # Get the occurence of fire and storm per plot
  tar_target(dist_occurence, get_dist_occurence(climate_dist_dflist)),
  
  # Extend for the simulation the disturbances over several years
  tar_target(climate_dist_dflist_ext, extend_disturbance(climate_dist_dflist, 5)),

  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load and format traits data ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- Traits data file
  tar_target(traitsNFI_file, "data/Traits/traitsNFI.csv", format = "file"),
  tar_target(wood.density_file, "data/Traits/GlobalWoodDensityDatabase.xls", format = "file"),
  tar_target(shade.tolerance_file, "data/Traits/shade_tolerance_FrenchNFI.csv", format = "file"),
  tar_target(TRY_file, "data/Traits/TRY_data_request_21092.txt", format = "file"),
  # -- Compile traits from different databases
  tar_target(traits_compiled, compile_traits(
    wood.density_file, traitsNFI_file, shade.tolerance_file, TRY_file, sp.in.sim)),
  # -- Get all species included in the simulations for filtering
  tar_target(sp.in.sim, gsub("\\ ", "\\_", unique(species_list$species))),


  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare species objects and run simulations ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  # Extend the disturbance coefficients of matreex using traits
  tar_target(disturb_coef_ext, extend_disturb_coef(
    intensity_file, traits_compiled)),
   
  # Build the regional pool of each NFI plot
  tar_target(regional_pool, make_regional_pool(
    NFI_plots_selected, NFI_forest_cover, coef_ba_reg)),

  # Get the size distribution of each species in each plot
  tar_target(species_distrib, get_species_distrib(NFI_data_sub)),
  # Get the climatic margins
  tar_target(climate_margins, get_margins(NFI_climate_sub, NFI_plots_selected)),

  # Make species objects
  # -- create table listing all species to simulate
  tar_target(species_list, make_species_list(
    NFI_data_sub, NFI_plots_selected, regional_pool)),
  # -- Vector from 1 to number of species to make (useful for parallel computing)
  tar_target(ID.species, species_list$ID.species),
  # -- Make species via branching over ID.species
  tar_target(species_mu, make_species_mu_rds(
    species_list, climate_margins, disturb_coef_ext, new_fit_list, ID.species),
             pattern = map(ID.species), iteration = "vector", format = "file"),

  # Make simulations
  # -- Create a list of simulations to perform
  tar_target(simul_list, make_simul_list(NFI_plots_selected)),
  # -- Vector from 1 to number of simulations to make (for branching)
  tar_target(ID.simulation, c(1:dim(simul_list)[1])),
  # -- Make simulations with regional pool
  tar_target(simulations_pool, make_simulations(
    species_distrib, species_mu, species_list, climate_dist_dflist_ext,
    regional_pool, simul_list, disp_kernel, use_pool = TRUE, ID.simulation),
    pattern = map(ID.simulation), iteration = "list"),
  # -- Make simulations without regional pool
  tar_target(simulations_nopool, make_simulations(
    species_distrib, species_mu, species_list, climate_dist_dflist_ext,
    regional_pool, simul_list, disp_kernel, use_pool = FALSE, ID.simulation),
    pattern = map(ID.simulation), iteration = "list"),

  # Extract output of the simulations
  tar_target(sim_output_pool, bind_rows(simulations_pool, .id = NULL)),
  tar_target(sim_output_nopool, bind_rows(simulations_nopool, .id = NULL)),


  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Export plots ----
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # Plots for methods
  # -- Plot the position of species along functional and climatic space
  tar_target(fig_funclim_species, plot_funclim_species(
    NFI_data_sub, NFI_plots_selected, traits_compiled,
    "output/fig/methods/fig_funclim_species.jpg"), format = "file"),
  # -- Plot the map, and climate change on one plot
  tar_target(fig_map_clim_dist, plot_map_clim_dist(
    NFI_plots_selected, climate_dist_dflist,
    "output/fig/methods/fig_map_clim_dist.jpg"), format = "file"),


  # Plot for analyses
  # - Classify plots in succession stage
  tar_target(NFI_succession, classify_succession(NFI_data_sub, NFI_plots_selected)),
  # - Biogeo effect only for abundance as N
  tar_target(fig_biogeo_effect_N, plot_biogeo_effect_per.metric(
    sim_output_pool, NFI_succession, simul_list, NFI_data_sub, traits_compiled, "N",
    list(fig = "output/fig/analyses/fig_biogeo_effect_N.jpg", 
         resid = "output/fig/supplementary/fig_residuals_biogeo.jpg")), format = "file"),
  # - Local effect only for abundance as N
  tar_target(fig_pool_N, plot_pool_effect(
    regional_pool, sim_output_pool, sim_output_nopool, simul_list, NFI_data_sub, 
    traits_compiled, dist_occurence, "N", 
    list(fig = "output/fig/analyses/fig_pool_N.jpg", 
         resid = "output/fig/supplementary/fig_residuals_pool.jpg")), format = "file"),
  
  
  # Plots for supplementary material
  # - biogeo effect only for abundance as N but without the regional pool
  tar_target(fig_biogeo_effect_N_nopool, plot_biogeo_effect_per.metric(
    sim_output_nopool, NFI_succession, simul_list, NFI_data_sub, traits_compiled, "N",
    list(fig = "output/fig/supplementary/fig_biogeo_effect_N_nopool.jpg", 
         resid = "output/fig/supplementary/fig_residuals_biogeo_nopool.jpg")), format = "file"),
  # - plot the resulting size distribution per succession and climate
  tar_target(fig_distrib_succession, plot_succession_distrib(
    NFI_succession, "output/fig/supplementary/fig_distrib.jpg"), format = "file"),
  # - Plot the prediction of disturbance coefficients from functional traits
  tar_target(fig_disturb_coef, plot_disturb_coef(
    intensity_file, traits_compiled, "output/fig/supplementary/fig_coef.jpg"), 
    format = "file"),
  # - Plot the distribution of species richness in each climate
  tar_target(fig_richness_distrib, plot_richness_distrib(
    NFI_plots_selected, "output/fig/supplementary/fig_richness.jpg"), format = "file"),
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Export tables ----
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  tar_target(table_funclim_species, export_table_funclim_species(
    NFI_data_sub, NFI_plots_selected, traits_compiled, 
    "output/tables/table_species.tex"), format = "file")
)

