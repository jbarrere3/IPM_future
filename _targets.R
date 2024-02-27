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
                 "modi", "ggridges")
for(i in 1:length(packages.in)){
  if(!(packages.in[i] %in% rownames(installed.packages()))){
    install.packages(packages.in[i])
  }
}  
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
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
  
  # Filter data based on species present in IPM and disturbance regimes
  # -- select NFI plots
  tar_target(NFI_plots_selected, select_NFI_plots(
    NFI_data, NFI_disturbance, NFI_climate, nclim = 10, nplots.per.clim = 200)), 
  # -- subset all data
  tar_target(NFI_data_sub, subset(
    NFI_data, plotcode %in% NFI_plots_selected$plotcode)),
  tar_target(NFI_disturbance_sub, subset(
    NFI_disturbance, plotcode %in% NFI_plots_selected$plotcode)),
  tar_target(NFI_climate_sub, subset(
    NFI_climate, plotcode %in% NFI_plots_selected$plotcode)),
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load and format traits data ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- Traits data file
  tar_target(traits_file, "data/Traits/traits.csv", format = "file"),
  # -- Read traits data
  tar_target(traits, fread(traits_file)),
  # -- Get all species included in the simulations for filtering
  tar_target(sp.in.sim, gsub("\\ ", "\\_", unique(species_list$species))),
  # -- Get coordinates on the pca axis per species
  tar_target(pca1_per_species, get_pc1_per_species(traits, sp.in.sim)), 
  # -- Get demographic traits
  tar_target(traits_demo, get_traits_demo(sp.in.sim)),
  # -- Get coordinates on pca demographic axes per species
  tar_target(pca_demo_per_species, get_pc12_per_species(
    left_join(traits, traits_demo, by = "species"), sp.in.sim)), 
  
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Get distributions of storm and fire intensity ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # File containing disturbance intensity value
  tar_target(intensity_file, "data/Disturbance/jags_dominance.Rdata", 
             format = "file"), 
  
  # Fit model linking fire intensity to VPD
  tar_target(model_fire_vpd, fit_vpd_fire(intensity_file, NFI_data, NFI_climate)), 
  
  # Plot the relation between fire and vpd
  tar_target(fig_fire_vpd, plot_fire_vpd(
    model_fire_vpd, "output/fig/fig_fire_vpd.jpg"), format = "file"), 
  
  # Get distribution of storm intensity
  tar_target(Istorm_param, get_Istorm_param(intensity_file)),

  # Get the climate and disturbance dataframe for each plot and scenario
  tar_target(climate_dist_dflist, get_clim_dist_df(
    model_fire_vpd, Istorm_param, NFI_disturbance_sub, NFI_climate_sub)), 
  
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare species objects and run simulations ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Get the size distribution of each species in each plot
  tar_target(species_distrib, get_species_distrib(NFI_data_sub)),
  # Get the climatic margins
  tar_target(climate_margins, get_margins(NFI_climate_sub, NFI_plots_selected)), 
  
  # Make species objects
  # -- create table listing all species to simulate
  tar_target(species_list, make_species_list(NFI_data_sub, NFI_plots_selected)),
  # -- Vector from 1 to number of species to make (useful for parallel computing)
  tar_target(ID.species, species_list$ID.species),
  # -- Make species via branching over ID.species
  tar_target(species_mu, make_species_mu_rds(species_list, climate_margins, ID.species),
             pattern = map(ID.species), iteration = "vector", format = "file"), 
  
  # Make simulations
  # -- Create a list of simulations to perform
  tar_target(simul_list, make_simul_list(NFI_plots_selected)), 
  # -- Vector from 1 to number of simulations to make (for branching)
  tar_target(ID.simulation, c(1:dim(simul_list)[1])),
  # -- Make simulations 
  tar_target(simulations, make_simulations(
    species_distrib, species_mu, species_list, climate_dist_dflist,
    simul_list, ID.simulation), pattern = map(ID.simulation), 
    iteration = "vector", format = "file"), 
  
  # Extract output of the simulations
  tar_target(sim_output, get_simulations_output(simulations)), 
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Plots ---- 
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Plot of the simulations
  tar_target(fig_climate_effect, plot_climate_effect(
    sim_output, simul_list, pca1_per_species, "output/fig/analyses"), 
    format = "file"),
  tar_target(fig_temporalchange, plot_temporalchange(
    sim_output, simul_list, pca1_per_species, "output/fig/analyses_temporalchange"), 
    format = "file"),
  tar_target(fig_magnitudechange, plot_magnitudechange(
    sim_output, simul_list, pca_demo_per_species, "output/fig/analyses_magnitudechange"), 
    format = "file"), 
  
  # Plot characteristics of data along the gradient
  # -- Plot pca of recruitment traits
  tar_target(fig_pca_demo_traits, plot_traits_pca12(
    traits, traits_demo, sp.in.sim, "output/fig/description/pca_traits_demo.jpg"),
      format = "file"), 
  # -- Plot structure and composition along the climatic gradient
  tar_target(fig_str_compo, plot_str_compo_climate(
    NFI_climate_sub, NFI_data_sub, pca_demo_per_species, nclim = 10, 
    dir.out = "output/fig/description"), format = "file")
  
  
)

