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
                 "rnaturalearth", "sinkr", "egg")
for(i in 1:length(packages.in)){
  if(!(packages.in[i] %in% rownames(installed.packages()))){
    install.packages(packages.in[i])
  }
}  
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 8)
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
  # -- Load and format traits data ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- Traits data file
  tar_target(traits_file, "data/Traits/traits.csv", format = "file"),
  tar_target(traitsNFI_file, "data/Traits/traitsNFI.csv", format = "file"),
  tar_target(wood.density_file, "data/Traits/GlobalWoodDensityDatabase.xls", format = "file"),
  tar_target(shade.tolerance_file, "data/Traits/shade_tolerance_FrenchNFI.csv", format = "file"),
  tar_target(TRY_file, "data/Traits/TRY_data_request_21092.txt", format = "file"),
  # -- Read traits data
  tar_target(traits, fread(traits_file)),
  tar_target(traitsNFI, fread(traitsNFI_file)),
  # -- Compile traits from different databases
  tar_target(traits_compiled, compile_traits(
    wood.density_file, traitsNFI_file, shade.tolerance_file, TRY_file, sp.in.sim)),
  # -- Get all species included in the simulations for filtering
  tar_target(sp.in.sim, gsub("\\ ", "\\_", unique(species_list$species))),
  
  
  
  
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
  # -- Prepare species objects and run simulations ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Extend the disturbance coefficients of matreex using traits
  tar_target(disturb_coef_ext, extend_disturb_coef(
    intensity_file, traits, traitsNFI)),
  
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
  # -- Make simulations
  tar_target(simulations, make_simulations(
    species_distrib, species_mu, species_list, climate_dist_dflist_ext,
    regional_pool, simul_list, disp_kernel, ID.simulation), 
    pattern = map(ID.simulation), iteration = "vector", format = "file"),

  # Extract output of the simulations
  tar_target(sim_output, get_simulations_output(simulations)),
  tar_target(sim_output_short, get_simulations_output_short(simulations)),
  

  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots ----
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # Plots for methods
  # -- Plot the position of species along functional and climatic space
  tar_target(fig_funclim_species, plot_funclim_species(
    NFI_data_sub, NFI_plots_selected, traits_compiled, 
    "output/fig/methods/fig_funclim_species.jpg"), format = "file"), 
  # -- Plot the map, and climate change on one plot
  tar_target(fig_map_clim_dist, plot_map_clim_dist(
    NFI_plots_selected, simul_list, climate_dist_dflist, 
    "output/fig/methods/fig_map_clim_dist.jpg"), format = "file"),
  
 
  # Plot for analyses
  # - Classify plots in succession stage
  tar_target(NFI_succession, classify_succession(NFI_data_sub, NFI_plots_selected)), 
  # - Biogeo effect only for abundance as N
  tar_target(fig_biogeo_effect_N, plot_biogeo_effect_per.metric(
    sim_output_short, NFI_succession, simul_list, traits_compiled, "N",
    "output/fig/analyses/fig_biogeo_effect_N.jpg"), format = "file"), 
  # - Local effect only for abundance as N
  tar_target(fig_local_N, plot_local_effect_permetric(
    regional_pool, sim_output_short, simul_list, NFI_succession, traits_compiled, 
    NFI_plots_selected, dist_occurence, "N", 
    "output/fig/analyses/fig_local_effect_N.jpg"), format = "file"), 
  
  # Plots for supplementary material
  # - plot the resulting size distribution per succession and climate
  tar_target(fig_distrib_succession, plot_succession_distrib(
    NFI_succession, "output/fig/supplementary/fig_distrib.jpg"), format = "file"), 
  # - Biogeo effect only for abundance as BA
  tar_target(fig_biogeo_effect_BA, plot_biogeo_effect_per.metric(
    sim_output_short, NFI_succession, simul_list, traits_compiled, "BA",
    "output/fig/supplementary/fig_biogeo_effect_BA.jpg"), format = "file"), 
  # - Local effect only for abundance as BA
  tar_target(fig_local_BA, plot_local_effect_permetric(
    regional_pool, sim_output_short, simul_list, NFI_succession, traits_compiled, 
    NFI_plots_selected, dist_occurence, "BA", 
    "output/fig/supplementary/fig_local_effect_BA.jpg"), format = "file")
)

