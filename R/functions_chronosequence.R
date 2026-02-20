#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_chronosequence.R  
#' @description R script containing all functions relative to 
#               the chronosequence analysis for IPM validation
#' @author Julien Barrere, Georges Kunstler, Maxime Jaunatre
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to extract the climate class and functional composition of plots for chronosequence
#' @param NFI_plots_selected NFI plots selected for simulation
#' @param NFI_data Raw NFI data
#' @param climate_files climate raster files
#' @param traits_compiled traits information per species included in matreex
#' @param file_chronoseq file including path containing data fro chronosequence
get_sp_and_clim_chronoseq = function(NFI_plots_selected, NFI_data, climate_files, 
                                     traits_compiled, file_chronoseq){
  
  # Read template raster
  template_raster = rast(climate_files[1])
  
  # Restrict NFI plots to French plots and convert to sf
  data_plot = NFI_plots_selected %>%
    filter(plotcode %in% subset(NFI_data, country == "FG")$plotcode) %>%
    mutate(climate_num = as.numeric(gsub("\\clim", "", climate))) %>%
    dplyr::select(longitude, latitude, climate_num) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Convert to vector in terra
  data_plot_vect = vect(data_plot)
  
  # Extract mean value of climate class per raster cell
  dominant_raster <- rasterize(data_plot_vect, template_raster, 
                               field = "climate_num", 
                               fun = "mean")
  
  # Read chronosequence dataframe
  data_chronoseq = fread(file_chronoseq) %>%
    # convert species name
    mutate(species = gsub("\\ subsp.+", "", species), 
           species = gsub("\\ var.+", "", species), 
           species = gsub("\\ x\\ ", "\\ ", species), 
           species = gsub("\\ f\\.\\ .+", "", species), 
           species = gsub("\\ ", "\\_", species), 
           species = ifelse(species %in% c("Betula_pendula", "Betula_pubescens"), 
                            "Betula", species))
  
  # Plot level information of chronoseq dataframe
  data_chronoseq_plot = data_chronoseq %>%
    filter(status == "alive" & dbh >= 90) %>%
    # Binary variable indicating whether a species is in matreex
    mutate(isinmatreex = ifelse(species %in% fit_species, 1, 0)) %>%
    group_by(plotcode, longitude, latitude, age) %>%
    summarize(prop.isinmatreex = sum(ba.ha*isinmatreex, na.rm = TRUE)/sum(
      ba.ha, na.rm = TRUE)) %>% ungroup() %>%
    filter(prop.isinmatreex == 1) %>% dplyr::select(-prop.isinmatreex) %>%
    mutate(climate_num = extract(
      dominant_raster, .[, c("longitude", "latitude")])[, 2]) %>%
    group_by(plotcode) %>%
    mutate(climate_num = round(climate_num, digits = 0)) %>%
    ungroup() %>% drop_na() %>%
    mutate(climate = paste0("clim", climate_num))
  
  # Calculate functional diversity
  data.FD = data_chronoseq %>%
    # Remove plots with species not in matreex
    filter(plotcode %in% data_chronoseq_plot$plotcode & status == "alive" &
             dbh >= 90) %>%
    # Calculate abundance per species per plot
    group_by(plotcode, species) %>% 
    summarize(ba.ha.sp = sum(ba.ha, na.rm = TRUE)) %>% ungroup() %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Gather by functional axis
    gather(key = "axis", value = "trait_value", names(traits_compiled[[1]])) %>%
    # Calculate cwm for each plotcode
    group_by(plotcode, axis) %>%
    mutate(cwm = weighted.mean(trait_value, w = ba.ha.sp)) %>%
    ungroup() %>%
    # Calculate square distance of species to centroid along each axis
    mutate(zsq_per_axis = (trait_value - cwm)^2) %>%
    # Calculate the distance of each species to the centroid
    group_by(plotcode, ba.ha.sp, species) %>%
    summarize(z = sqrt(sum(zsq_per_axis))) %>% ungroup() %>%
    # Calculate functional dispersion
    group_by(plotcode) %>%
    summarize(FD = weighted.mean(z, w = ba.ha.sp))
  
  
  # Extract simulation output on forest composition
  data.sp = data_chronoseq %>%
    # Remove plots with species not in matreex
    filter(plotcode %in% data_chronoseq_plot$plotcode & status == "alive" &
             dbh >= 90) %>%
    # Calculate abundance per species per plot
    group_by(plotcode, species) %>% 
    summarize(ba.ha.sp = sum(ba.ha, na.rm = TRUE)) %>% ungroup() %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Calculate composition metrics
    group_by(plotcode) %>%
    mutate(p = ba.ha.sp/sum(ba.ha.sp), 
           plnp = p*log(p)) %>%
    summarise(H = -sum(plnp), 
              cwm_GrSurv = weighted.mean(GrSurv, w = ba.ha.sp), 
              cwm_ShadeDrought = weighted.mean(ShadeDrought, w = ba.ha.sp)) %>%
    # Add functional diversity
    left_join(data.FD, by = c("plotcode"))
  
  # Output 
  out = data_chronoseq_plot %>% left_join(data.sp, by = "plotcode") %>%
    # Add average climate per climate class
    left_join((NFI_plots_selected %>%
                 filter(plotcode %in% subset(NFI_data, country == "FG")$plotcode) %>%
                 group_by(climate) %>%
                 summarize(pca1 = mean(pca1, na.rm = TRUE), 
                           lambda = mean(lambda, na.rm = TRUE)) %>%
                 ungroup()), 
              by = "climate")
  
  # Return output
  return(out)
  
}


#' Build the regional pool of plots for chronosequence
#' @param plots_selected_chronoseq plot-level data for chronosequence selected for simulation
#' @param coef_ba_reg Coefficients to calculate the regional basal area
make_regional_pool_chronoseq = function(plots_selected_chronoseq, coef_ba_reg){
  
  # Initialize dataset for reg ba calculation
  data.in = expand.grid(plotcode = plots_selected_chronoseq$plotcode, 
                        species = coef_ba_reg$species) %>%
    # Add information on climate, parameters and land cover
    left_join((plots_selected_chronoseq %>% dplyr::select(plotcode, pca1, lambda)), 
              by = "plotcode") %>%
    left_join(coef_ba_reg, by = "species") %>%
    # Calculate regional basal area per species based on parameters
    mutate(ba_reg_th = a*exp(-(pca1 - b)^2/c), 
           ba_reg_th = ifelse(is.na(ba_reg_th), 0, ba_reg_th)) %>% 
    # Calculate probability of sampling for each species per plotcode
    group_by(plotcode) %>%
    mutate(prob = ba_reg_th/sum(ba_reg_th, na.rm = TRUE)) %>% ungroup() %>%
    # Correct regional basal area by forest cover
    mutate(ba_reg = ba_reg_th*0.6)
  
  # Initialize output
  out = data.frame(plotcode = character(0), species = character(0), 
                   ba_reg = character(0))
  
  # Loop on all plotcodes
  for(i in 1:length(unique(data.in$plotcode))){
    
    # Subset dataset
    # - right plotcode
    data.i = data.in %>% filter(plotcode == unique(data.in$plotcode)[i])
    
    # Sort species
    sp.i = unique(sample(x = data.i$species, size = 1 + rpois(1, unique(data.i$lambda)), 
                         replace = FALSE, prob = data.i$prob))
    
    # Complete the output dataset
    out = rbind(out, data.i %>% filter(species %in% sp.i) %>% 
                  dplyr::select(plotcode, species, ba_reg))
    
  }
  
  # Return output
  return(out)
  
}


#' Function to extract the distribution of each species in each plot for chronoseq
#' @param plots_selected_chronoseq plot-level info for chronoseq
#' @param file_chronoseq file containing tree-level info for chronosequence
get_species_distrib_chronoseq = function(plots_selected_chronoseq, file_chronoseq){
  
  # Initialize output list
  list.out = vector(mode = "list", length = length(unique(plots_selected_chronoseq$plotcode)))
  names(list.out) = unique(plots_selected_chronoseq$plotcode)
  
  # Read chronoseq file and restrict to plots of interest
  data_chronoseq = fread(file_chronoseq) %>%
    filter(plotcode %in% plots_selected_chronoseq$plotcode) %>%
    filter(status == "alive" & dbh >= 90) %>%
    # convert species name
    mutate(species = gsub("\\ subsp.+", "", species), 
           species = gsub("\\ var.+", "", species), 
           species = gsub("\\ x\\ ", "\\ ", species), 
           species = gsub("\\ f\\.\\ .+", "", species), 
           species = gsub("\\ ", "\\_", species), 
           species = ifelse(species %in% c("Betula_pendula", "Betula_pubescens"), 
                            "Betula", species))
  
  
  # Loop on all plotcodes
  for(i in 1:length(names(list.out))){
    
    # Identify species present in plotcode i
    species.i = unique(subset(data_chronoseq, plotcode == names(list.out)[i])$species)
    
    # Loop on all species present
    for(j in 1:length(species.i)){
      
      # Load parameters of species j
      eval(parse(text=paste0("fit.ij <- fit_", species.i[j])))
      
      # Get the delay for species i
      delay.i = as.numeric(fit.ij$info['delay'])
      
      # Initialize distribution of species j in plot i
      eval(parse(text = paste0("list.out[[i]]$", species.i[j], 
                               " = rep(0, (700 + delay.i))")))
      
      # Build the mesh of the species
      mesh.ij  = data.frame(
        mesh.id = c(1:700)+delay.i, 
        dbh.mesh = seq(from = 90, to = get_maxdbh(fit.ij) * 1.1, length.out = 700))
      
      # Subset tree-level NFI data
      data.ij = data_chronoseq %>%
        filter(plotcode == names(list.out)[i] & species == species.i[j]) %>%
        dplyr::select(dbh, ba.ha)
      
      # Loop on all trees present in the plot from species j
      for(k in 1:dim(data.ij)[1]){
        
        # Position of tree k in the mesh
        id.k = (mesh.ij %>% 
                  mutate(diff = abs(dbh.mesh - data.ij$dbh[k])) %>%
                  arrange(diff))$mesh.id[1]
        
        # Add basal area of tree k in the output list
        list.out[[i]][[j]][id.k] = list.out[[i]][[j]][id.k] + data.ij$ba.ha[k]
      }
      
    }
    
  }
  
  # Return output list
  return(list.out)
  
}


#' Prepare average climate data for each climate class
#' @param NFI_plots_selected Plot-level information for plots included in sim
#' @param NFI_data Raw NFI data
#' @param NFI_climate Plot-level climate data
get_meanclimate_chronoseq = function(NFI_plots_selected, NFI_data, NFI_climate){
  
  NFI_climate %>%
    filter(plotcode %in% subset(NFI_data, country == "FG")$plotcode) %>%
    dplyr::select(all_of(c("plotcode", paste0("hist_sgdd_", c(1991:2000)), 
                           paste0("hist_wai_", c(1991:2000))))) %>%
    pivot_longer(names_to = "variable", values_to = "value", colnames(.)[2:dim(.)[2]]) %>%
    mutate(variable = gsub("hist\\_", "", variable)) %>%
    separate(col = "variable", into = c("var", "year"), sep = "_") %>%
    left_join(NFI_plots_selected %>% dplyr::select(plotcode, climate), 
              by = "plotcode") %>%
    group_by(climate, var) %>%
    summarize(mean = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(names_from = "var", values_from = "mean") %>%
    mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
           waib = 1/(1 + wai)) %>%
    drop_na()
  
}

#' Make mu simulations for chronosquence
#' @param plots_selected_chronoseq Plot level info on plots selected for chronoseq
#' @param species_distrib_chronoseq List of species distribution (per plot)
#' @param species_list df with information on all species object
#' @param species_mu vector containing all species mu rds files created
#' @param meanclimate_chronoseq Mean climate in each climate category
#' @param disp_kernel dispersal kernel for each species
#' @param regional_pool_chronoseq Regional pool per plot (chronosequence version)
#' @param use_pool boolean to indicate whether the regional pool should be used
#' @param ID.simulation ID in simul_dist of the simulation to perform
make_simulations_chronoseq = function(
    plots_selected_chronoseq, species_distrib_chronoseq, species_list, species_mu,
    meanclimate_chronoseq, disp_kernel, regional_pool_chronoseq, use_pool, ID.simulation){
  
  # Print simulation ID
  print(ID.simulation)
  
  # Identify the plotcode, climate, ssp, etc
  plot.in = plots_selected_chronoseq$plotcode[ID.simulation]
  clim.in = plots_selected_chronoseq$climate[ID.simulation]
  
  # Extract climate and distributions
  climate.in = as.numeric(meanclimate_chronoseq %>% 
                            filter(climate == clim.in) %>% 
                            dplyr::select(-climate))
  names(climate.in) = colnames(meanclimate_chronoseq)[-1]
  distributions.in = species_distrib_chronoseq[[as.character(plot.in)]]
  
  # Final list of species based on regional pool and distribution.in
  sp.reg = filter(regional_pool_chronoseq, plotcode == plot.in)$species
  sp.tot = unique(c(names(distributions.in), sp.reg))
  # Compile to make a regional pool
  reg.pool.in = (data.frame(species = sp.tot) %>%
                   left_join(filter(regional_pool_chronoseq, plotcode == plot.in), by = "species") %>%
                   mutate(ba_reg = ifelse(is.na(ba_reg), 0, ba_reg)))$ba_reg
  names(reg.pool.in) = sp.tot
  # Vector of migration rates
  mig.rate.in = 1 - left_join(data.frame(species = sp.tot), disp_kernel, 
                              by = "species")$p30
  names(mig.rate.in) = sp.tot
  
  # Initialize list of species
  list.sp = vector(mode = "list", length = length(sp.tot))
  names(list.sp) = paste0("mu_", sp.tot)
  
  # Check that the species object were created for this climate
  if(all(sp.tot %in% subset(species_list, climate == clim.in)$species)){
    # Loop on all species of the list
    for(i in 1:length(names(list.sp))){
      
      # Get the index of mu from species list
      id.mu = subset(species_list, species == sp.tot[i] & climate == clim.in)$ID.species
      
      # Read the species
      list.sp[[i]] = readRDS(species_mu[id.mu])
      
      # Convert the mesh from basal area per ha to number of tree per ha
      # -- mesh in basal area (m2)
      mesh_ba.i = pi*(list.sp[[i]]$IPM$mesh/2000)^2
      # -- If species is present in the data (and not only in reg pool)
      if(sp.tot[i] %in% names(distributions.in)){
        # -- divide size distribution in basal area per m2 by mesh in m2
        distrib.i = distributions.in[[sp.tot[i]]]/mesh_ba.i
        # -- Set na value to 0 (delay)
        distrib.i[which(is.na(distrib.i))] = 0
        # -- If species only in regional pool, set distribution to 0
      } else {
        distrib.i = 0*mesh_ba.i
      }
      
      # Use the new distribution in N/ha to initialize the population 
      list.sp[[i]]$init_pop <- def_init_k(distrib.i)
      
    }
    
    # Make forest
    # - With the regional pool
    if(use_pool){
      forest.in = forest(species = list.sp, harv_rules = c(
        Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1), 
        regional_abundance = reg.pool.in, migration_rate = mig.rate.in)
    }
    # - Without the regional pool
    if(!use_pool){
      forest.in = forest(species = list.sp[paste0("mu_", names(distributions.in))], 
                         harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1), 
                         migration_rate = mig.rate.in[names(distributions.in)])
    }
    
    
    # Different code depending on whether the scenario includes disturbances or not
    try(sim.in <- sim_deter_forest(
      forest.in, tlim = 170, climate = climate.in, 
      equil_dist = 170,  equil_time = 170, 
      verbose = TRUE, correction = "cut"), silent = TRUE)
    
    
    
    # If the simulation failed, return an empty object
    if(!exists("sim.in")) sim.in = list()
    # Otherwise, extract lag and add it to the distribution
    else{
      # Extract the mesh for the delay at equilibrium
      sim.in = sim.in %>%
        filter(!equil & var %in% c("BAsp", "N")) %>%
        mutate(ID.simulation = ID.simulation) %>% 
        dplyr::select(ID.simulation, species, time, var, value) %>%
        spread(key = "var", value = "value") %>%
        rename(BA = BAsp)
    } 
    
  }else{
    sim.in = list()
  }
  
  # Return output list
  return(sim.in)
}

make_simulations_chronoseq_lag= function(
    plots_selected_chronoseq, species_distrib_chronoseq, species_list, species_mu,
    meanclimate_chronoseq, disp_kernel, regional_pool_chronoseq, use_pool, ID.simulation){
  
  # Print simulation ID
  print(ID.simulation)
  
  # Identify the plotcode, climate, ssp, etc
  plot.in = plots_selected_chronoseq$plotcode[ID.simulation]
  clim.in = plots_selected_chronoseq$climate[ID.simulation]
  
  # Extract climate and distributions
  climate.in = as.numeric(meanclimate_chronoseq %>% 
                            filter(climate == clim.in) %>% 
                            dplyr::select(-climate))
  names(climate.in) = colnames(meanclimate_chronoseq)[-1]
  distributions.in = species_distrib_chronoseq[[as.character(plot.in)]]
  
  # Final list of species based on regional pool and distribution.in
  sp.reg = filter(regional_pool_chronoseq, plotcode == plot.in)$species
  sp.tot = unique(c(names(distributions.in), sp.reg))
  # Compile to make a regional pool
  reg.pool.in = (data.frame(species = sp.tot) %>%
                   left_join(filter(regional_pool_chronoseq, plotcode == plot.in), by = "species") %>%
                   mutate(ba_reg = ifelse(is.na(ba_reg), 0, ba_reg)))$ba_reg
  names(reg.pool.in) = sp.tot
  # Vector of migration rates
  mig.rate.in = 1 - left_join(data.frame(species = sp.tot), disp_kernel, 
                              by = "species")$p30
  names(mig.rate.in) = sp.tot
  
  # Initialize list of species
  list.sp = vector(mode = "list", length = length(sp.tot))
  names(list.sp) = paste0("mu_", sp.tot)
  
  # Check that the species object were created for this climate
  if(all(sp.tot %in% subset(species_list, climate == clim.in)$species)){
    # Loop on all species of the list
    for(i in 1:length(names(list.sp))){
      
      # Get the index of mu from species list
      id.mu = subset(species_list, species == sp.tot[i] & climate == clim.in)$ID.species
      
      # Read the species
      list.sp[[i]] = readRDS(species_mu[id.mu])
      
      # Convert the mesh from basal area per ha to number of tree per ha
      # -- mesh in basal area (m2)
      mesh_ba.i = pi*(list.sp[[i]]$IPM$mesh/2000)^2
      # -- If species is present in the data (and not only in reg pool)
      if(sp.tot[i] %in% names(distributions.in)){
        # -- divide size distribution in basal area per m2 by mesh in m2
        distrib.i = distributions.in[[sp.tot[i]]]/mesh_ba.i
        # -- Set na value to 0 (delay)
        distrib.i[which(is.na(distrib.i))] = 0
        # -- If species only in regional pool, set distribution to 0
      } else {
        distrib.i = 0*mesh_ba.i
      }
      
      # Use the new distribution in N/ha to initialize the population 
      list.sp[[i]]$init_pop <- def_init_k(distrib.i)
      
    }
    
    # Make forest
    # - With the regional pool
    if(use_pool){
      forest.in = forest(species = list.sp, harv_rules = c(
        Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1), 
        regional_abundance = reg.pool.in, migration_rate = mig.rate.in)
    }
    # - Without the regional pool
    if(!use_pool){
      forest.in = forest(species = list.sp[paste0("mu_", names(distributions.in))], 
                         harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1), 
                         migration_rate = mig.rate.in[names(distributions.in)])
    }
    
    # Different code depending on whether the scenario includes disturbances or not
    try({
      sim.in.50 <- sim_deter_forest(
        forest.in, tlim = 50, climate = climate.in,
        equil_dist = 51,  equil_time = 51,
        verbose = TRUE, correction = "cut");
      for(i in 1:length(names(list.sp))){
        
        # Get the index of mu from species list
        id.mu = subset(species_list, species == sp.tot[i] & climate == clim.in)$ID.species
        
        # Read the species
        list.sp[[i]] = readRDS(species_mu[id.mu])
        
        # Convert the mesh from basal area per ha to number of tree per ha
        # -- mesh in basal area (m2)
        mesh_ba.i = pi*(list.sp[[i]]$IPM$mesh/2000)^2
        # -- If species is present in the data (and not only in reg pool)
        if(sp.tot[i] %in% names(distributions.in)){
          # -- divide size distribution in basal area per m2 by mesh in m2
          distrib.i = distributions.in[[sp.tot[i]]]/mesh_ba.i
          # -- Set na value to 0 (delay)
          distrib.i[which(is.na(distrib.i))] = 0
          # -- If species only in regional pool, set distribution to 0
        } else {
          distrib.i = 0*mesh_ba.i
        }
        
        # Use the new distribution in N/ha to initialize the population
        distrib_t50 <- sim.in.50 %>%
          dplyr::filter(var == "n", time == 50, species == sp.tot[i]) %>%
          mutate(n2 =distrib.i) %>% mutate(n2 = if_else(size >0, n2, value) )  %>%
          pull(n2)
        
        list.sp[[i]]$init_pop  <- def_init_k(distrib_t50)
      };
      forest.in2 = forest(species = list.sp, harv_rules = c(
        Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
        regional_abundance = reg.pool.in, migration_rate = mig.rate.in);
      sim.in <- sim_deter_forest(
        forest.in2, tlim = 170, climate = climate.in,
        equil_dist = 170,  equil_time = 170,
        verbose = TRUE, correction = "cut")
    }, silent = TRUE)
    
    # If the simulation failed, return an empty object
    if(!exists("sim.in")) sim.in = list()
    # Otherwise, extract lag and add it to the distribution
    else{
      # Extract the mesh for the delay at equilibrium
      sim.in = sim.in %>%
        filter(!equil & var %in% c("BAsp", "N")) %>%
        mutate(ID.simulation = ID.simulation) %>% 
        dplyr::select(ID.simulation, species, time, var, value) %>%
        spread(key = "var", value = "value") %>%
        rename(BA = BAsp)
    } 
    
  }else{
    sim.in = list()
  }
  
  # Return output list
  return(sim.in)
}


#' Function to plot the change in species composition in simulation vs chronoseq
#' @param plots_selected_chronoseq NFI plots selected for the chronosequence analysis
#' @param sim_output_chronoseq output of the simulations for chronosequence
#' @param sp_and_clim_chronoseq Forest composition in Frnech NFI plots
#' @param traits_compiled Functional traits value per species
#' @param file.out Name of the file to save, including path
plot_chronosequence = function(plots_selected_chronoseq, sim_output_chronoseq, 
                               sp_and_clim_chronoseq, traits_compiled, 
                               file.out, metrics = c("div", "traits")){
  metrics <- match.arg(metrics)
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Calculate functional diversity
  data.FD = sim_output_chronoseq %>%
    # Add plotcode and climate
    left_join(data.frame(plotcode = plots_selected_chronoseq$plotcode, 
                         ID.simulation = c(1:dim(plots_selected_chronoseq)[1])), 
              by = "ID.simulation") %>%
    left_join(plots_selected_chronoseq %>% dplyr::select(plotcode, climate), 
              by = "plotcode") %>%
    # Convert time in age (all plots have 20 years-old)
    mutate(age = time + 20) %>%
    # Select columns of interest
    dplyr::select(plotcode, climate, age, species, BA) %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") 
  
  # Calculate species diversity and cwm
  data.sp = data.FD %>%
    # Calculate composition metrics
    group_by(plotcode, climate, age) %>%
    mutate(p = BA/sum(BA), 
           plnp = p*log(p)) %>%
    summarise(H = -sum(plnp), 
              cwm_GrSurv = weighted.mean(GrSurv, w = BA), 
              cwm_ShadeDrought = weighted.mean(ShadeDrought, w = BA))
  
  # Finish calculation of functional diversity
  data.FD = data.FD %>%
    # Gather by functional axis
    gather(key = "axis", value = "trait_value", names(traits_compiled[[1]])) %>%
    # Calculate cwm for each time step and each simulation
    group_by(plotcode, climate, age, axis) %>%
    mutate(cwm = weighted.mean(trait_value, w = BA)) %>%
    ungroup() %>%
    # Calculate square distance of species to centroid along each axis
    mutate(zsq_per_axis = (trait_value - cwm)^2) %>%
    # Calculate the distance of each species to the centroid
    group_by(plotcode, climate, age, species, BA) %>%
    summarize(z = sqrt(sum(zsq_per_axis))) %>% ungroup() %>%
    # Calculate functional dispersion
    group_by(plotcode, climate, age) %>%
    summarize(FD = weighted.mean(z, w = BA))
  
  # Join all data together And average per decade and climate
  data.out = data.sp %>% 
    # Keep only one simulation every 20 years
    filter(age %in% unique(round(c(10:200)/20, digits = 0)*20)) %>%
    # Remove NA in species diversity
    mutate(H = ifelse(is.na(H), 0, H)) %>%
    # Add functional diversity
    left_join(data.FD, by = c("plotcode", "climate", "age")) %>%
    # Add the source of data (here simulation)
    mutate(source = "simulation") %>%
    # Add chronosequence data
    rbind(sp_and_clim_chronoseq %>%
            mutate(age = round(age/20, digits = 0)*20) %>%
            mutate(source = "chronosequence") %>%
            dplyr::select(plotcode, climate, age, H, cwm_GrSurv, cwm_ShadeDrought, 
                          FD, source)) %>%
    # Average composition
    pivot_longer(values_to = "trait", names_to = "comp_var", cols = all_of(c(
      "H", "FD", "cwm_GrSurv", "cwm_ShadeDrought"))) %>%
    group_by(age, climate, comp_var, source) %>%
    summarize(mean = mean(trait, na.rm = TRUE), 
              lwr = quantile(trait, 0.05, na.rm = TRUE), 
              upr = quantile(trait, 0.95, na.rm = TRUE), 
              se = sd(trait, na.rm = TRUE),
              n = n()) %>%
    ungroup() # %>% filter(n > 9)
  
  if(metrics == "div"){
    data.out <- filter(data.out, comp_var %in% c("H", "FD"))
  } else {
    data.out <- filter(data.out, comp_var %in% c("cwm_GrSurv", "cwm_ShadeDrought"))
  }
  
  # Plot chronosequence vs simulations
  plot.out = data.out %>%
    filter(age < 200) %>%
    mutate(climate = factor(climate, levels = paste0("clim", c(1:10))), 
           age = ifelse(source == "simulation", age+3, age-3)) %>%
    ggplot(aes(x = age, color = source, fill = source)) + 
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) + 
    geom_point(aes(y = mean), shape = 21) +
    scale_color_manual(values = c("#780000", "#003049")) + 
    scale_fill_manual(values = c("#C1121F", "#669BBC")) +
    facet_grid(comp_var ~ climate, scales = "free_y") + 
    theme_bw()
  
  # Save the plot
  ggsave(file.out, plot.out, width = 38, height = 16 , units = "cm", 
         dpi = 600, bg = "white")
  
  # Return file saved
  return(file.out)
}


#' @param species_distrib_chronoseq Species distribution in each plot for 
#' chronoseq simulations
#' @param plot_per_climate list of num vector with ID of plot in each climate
#' @param climate names of the different climate
#' @import purrr
getmeandistrib_chronoseq <- function(species_distrib_chronoseq, 
                                     plot_per_climate, 
                                     mean_plot_ID.simulation_chronoseq, 
                                     species_list){
  
  res <- vector("list", length(plot_per_climate))
  names(res) <- names(plot_per_climate)
  
  for(c in names(plot_per_climate)){
    clim_distrib <- species_distrib_chronoseq[as.character(plot_per_climate[[c]])]
    sp <- map(clim_distrib, names) |> flatten_chr() |> unique()
    sp <- sp[sp %in% with(species_list, species[climate == c])]
    sum_distrib <- map(sp, function(s) map(clim_distrib, pluck(s)) |> compact() |> reduce(`+`))
    names(sum_distrib) <- sp
    mean_distrib <- map(sum_distrib, ~ .x / length(clim_distrib))
    res[[c]] <- mean_distrib
  }
  
  names(res) <- as.character(unname(
    mean_plot_ID.simulation_chronoseq[names(res)]
  ))
  
  return(res)
}



