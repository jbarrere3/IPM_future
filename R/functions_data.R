#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Manage NFI, Climate and Disturbance data  ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Extract climate data for each NFI plot
#' @param climate_files character containing the name of all climate files
#' @param NFI_coord dataframe containing NFI coordinates
extract_climate_NFI = function(climate_files, NFI_coord){
  
  # Initialize output
  out = NFI_coord %>%
    dplyr::select(plotcode, longitude, latitude)
  
  # Loop on all climate files
  for(i in 1:length(climate_files)){
    
    # File i
    file.i = climate_files[i]
    
    # Get scenario i
    scenario.i = gsub("\\/.+", "", gsub("data\\/Climate\\/", "", file.i))
    
    # Get variable i
    var.i = gsub("\\_.+", "", gsub(".+\\/", "", file.i))
    
    # Get year i
    year.i = gsub("\\..+", "", gsub(".+\\_", "", file.i))
    
    # Read raster for file i
    rast.i = rast(file.i)
    
    # Extract data from file i
    out$value = extract(rast.i, out[, c("longitude", "latitude")])[, 2]
    colnames(out)[dim(out)[2]] = paste(scenario.i, var.i, year.i, sep = "_")
  }
  
  # Return output
  return(out)
  
}

#' Function to get parameter values per iteration for each species and disturbance
#' @param rdata.file rdata contianing the model outputs
get_I_from_rdata <- function(rdata.file){
  
  # Load rdata
  load(rdata.file)
  
  # Loop on all disturbances
  for(i in 1:length(names(jags.list))){
    
    # Format the table containing parameter value per species, disturbance and iteration
    param.table.i <- ggs(as.mcmc(jags.list[[i]])) %>%
      # Add disturbance
      mutate(disturbance = names(jags.list)[i]) %>%
      # Extract information on parameter, species and country
      mutate(Param = gsub("\\[.+", "", Parameter), 
             plot = as.integer(ifelse(Param == "I", gsub(".+\\[", "", gsub("\\]", "", Parameter)), NA_integer_))) %>%
      # Remove the estimation of intensity and deviance
      filter(Param == "I") %>%
      dplyr::select(Iteration, Chain, disturbance, plot, I = value)
    
    # Store table in the final table
    if(i == 1) param.table <- param.table.i
    else param.table <- rbind(param.table, param.table.i)
  }
  
  # return the list
  return(param.table)
  
}

#' Function to fit a model linking fire intensity to VPD
#' @param intensity_file file containing fire intensity values
#' @param NFI_data dataframe containin NFI data
#' @param NFI_climate data frame containing climate data for NFI
fit_vpd_fire = function(intensity_file, NFI_data, NFI_climate){
  
  # Load sensitivity file
  load(intensity_file)
  
  # Get intensity
  I_data = get_I_from_rdata(intensity_file)
  
  # Format dataset
  I.per.plot = I_data %>%
    filter(disturbance == "fire") %>%
    left_join(corresp.tables$fire$plotcode.table, by = "plot") %>%
    mutate(I.logit = log(I/(1 - I))) %>%
    group_by(plotcode) %>%
    summarize(I.mean = mean(I, na.rm = TRUE), 
              I.lwr = quantile(I, 0.025, na.rm = TRUE), 
              I.upr = quantile(I, 0.975, na.rm = TRUE), 
              I.logit.mean = mean(I.logit, na.rm = TRUE), 
              I.logit.var = var(I.logit, na.rm = TRUE))
  
  
  # Plot data 
  data.plot = NFI_data %>%
    filter(plotcode %in% I.per.plot$plotcode) %>%
    dplyr::select(plotcode, country, longitude, latitude, year1, year2) %>%
    distinct()
  
  # Extract vpd for each plot
  colnames.vpd = colnames(NFI_climate)[grep("vpd", colnames(NFI_climate))]
  data.vpd = NFI_climate %>%
    filter(plotcode %in% data.plot$plotcode) %>%
    dplyr::select(plotcode, colnames.vpd) %>%
    gather(key = "variable", value = "vpd", colnames.vpd) %>%
    separate(col = "variable", into = c("scen", "var", "year"), sep = "\\_") %>%
    filter(scen == "hist") %>%
    left_join((data.plot %>% dplyr::select(plotcode, year1, year2)), 
              by = "plotcode") %>%
    mutate(year = as.numeric(year)) %>%
    filter(year >= year1 & year <= year2) %>%
    group_by(plotcode) %>%
    summarize(vpd.mean = mean(vpd, na.rm = TRUE))
  
  # Add vpd to plot data
  data.plot = data.plot %>% left_join(data.vpd, by = "plotcode")
  
  # Final dataset
  data = data.plot %>% 
    left_join(I.per.plot, by = "plotcode") %>% 
    mutate(weight = 1/I.logit.var) %>%
    drop_na()
  
  # Fit model
  model = betareg(I.mean ~ vpd.mean, weights = weight, data = data, link = "logit")
  
  # Output list : model and the data
  out = list(model = model, data = data)
  
  # Return output list
  return(out)
  
}


#' Function to extract beta parameters of storm intensity distribution
#' @param intensity_file rdata file containing posterior of intensity
get_Istorm_param = function(intensity_file){
  
  # Load sensitivity file
  load(intensity_file)
  
  # Get intensity
  I_data = get_I_from_rdata(intensity_file)
  
  # Format dataset
  I.per.plot = I_data %>%
    filter(disturbance == "storm") %>%
    left_join(corresp.tables$storm$plotcode.table, by = "plot") %>%
    group_by(plotcode) %>%
    summarize(I.mean = mean(I, na.rm = TRUE))
  
  # Extract beta parameters for storm
  out = fitdistr(I.per.plot$I.mean, densfun = "beta", 
                 start = list(shape1 = 1, shape2 = 2))[[1]]
  
  # Return output
  return(out)
  
}



#' Function to prepare disturbance and climate dataframe for each plot
#' @param model_fire_vpd model and data linking disturbance intensity to vpd
#' @param Istorm_param parameters of storm intensity beta distribution
#' @param NFI_disturbance data of disturbance probability per plot
#' @param NFI_climate data of climate per plot
get_clim_dist_df = function(model_fire_vpd, Istorm_param, NFI_disturbance, 
                            NFI_climate){
  
  # Initialize output list
  list.out = vector(mode='list', length=dim(NFI_disturbance)[1])
  names(list.out) = NFI_disturbance$plotcode
  
  # Table matching years and "time"
  year.table = data.frame(name = colnames(NFI_climate)) %>%
    filter(!(name %in% c("plotcode", "longitude", "latitude"))) %>%
    mutate(year = as.numeric(gsub(".+\\_", "", name))) %>%
    dplyr::select(year) %>%
    distinct() %>%
    arrange(year) %>%
    mutate(t = c(1:dim(.)[1]))
  
  # scenarios for which to generate climate files
  ssp.in = c("ssp126", "ssp370", "ssp585")
  
  # Inside the loop
  for(i in 1:length(names(list.out))){
    
    # Printer
    if(floor(i/100) == i/100) print(paste0("plot ", i, "/", length(names(list.out))))
    
    # Prepare climate dataframe for plot i
    data.climate.i = NFI_climate %>%
      filter(plotcode == names(list.out)[i]) %>%
      dplyr::select(colnames(.)[grep("_", colnames(.))]) %>%
      gather(key = "variable", value = "value", colnames(.)) %>%
      separate(col = "variable", into = c("scenario", "var", "year"), sep = "_") %>%
      mutate(year = as.numeric(year)) %>%
      left_join(year.table, by = "year") %>%
      dplyr::select(t, scenario, var, value) %>%
      spread(key = "var", value = "value") %>%
      mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
             waib = 1/(1 + wai)) %>%
      rename(vpd.mean = vpd) %>%
      mutate(Ifire.mu = predict(model_fire_vpd$model, newdata = ., type = "response"), 
             Ifire.phi = predict(model_fire_vpd$model, newdata = ., type = "precision")) %>%
      mutate(afire = Ifire.mu*Ifire.phi, bfire = (1 - Ifire.mu)*Ifire.phi) 
    
    # Prepare disturbance dataframe for plot i
    data.disturbance.i = NFI_disturbance %>%
      filter(plotcode == names(list.out)[i]) %>%
      dplyr::select(colnames(.)[grep("_", colnames(.))]) %>%
      gather(key = "variable", value = "value", colnames(.)) %>%
      separate(col = "variable", into = c("disturbance", "year"), sep = "_") %>%
      mutate(year = as.numeric(year)) %>%
      left_join(year.table, by = "year") %>%
      dplyr::select(t, disturbance, value) %>%
      spread(key = "disturbance", value = "value") %>%
      drop_na()
    
    # Loop on all scenario
    for(j in 1:length(ssp.in)){
      # Initialize the list for this scenario
      eval(parse(text = paste0("list.out[[i]]$", ssp.in[j], " = list()")))
      
      # subset climate table and disturbance table to the right scenario
      data.climate.ij = data.climate.i %>% filter(scenario %in% c("hist", ssp.in[j]))
      
      # Initialize the disturbance table for scenario j
      disturbance.ij = data.frame(t = year.table$t, IsSurv = FALSE) %>%
        mutate(storm = rbinom(dim(.)[1], 1, data.disturbance.i$storm), 
               fire = rbinom(dim(.)[1], 1, data.disturbance.i$fire)) %>%
        mutate(type = case_when(fire == 1 & storm == 1 ~ sample(c("fire", "storm"), 1), 
                                fire == 1 & storm == 0 ~ "fire", 
                                fire == 0 & storm == 1 ~ "storm", 
                                TRUE ~ "none"), 
               intensity = NA_real_) %>%
        filter(type != "none") %>%
        dplyr::select(type, IsSurv, t, intensity)
      
      # Fill the intensity column if there were disturbances
      if(dim(disturbance.ij)[1] > 0){
        # Loop on all disturbance events
        for(k in 1:dim(disturbance.ij)[1]){
          # If it is a fire disturbance
          if(disturbance.ij$type[k] == "fire"){
            # Get fire intensity parameters 
            param.k = (data.climate.ij %>% 
                         filter(t == disturbance.ij$t[k]))[, c("afire", "bfire")]
            # Attribute an intensity value based on the parameters
            disturbance.ij$intensity[k] = rbeta(1, shape1 = param.k$afire, 
                                                shape2 = param.k$bfire)
          }
          # If it is a storm disturbance
          if(disturbance.ij$type[k] == "storm"){
            # Attribute an intensity value based on the parameters
            disturbance.ij$intensity[k] = rbeta(1, shape1 = Istorm_param[1], 
                                                shape2 = Istorm_param[2])
          }
        }
      }
      
      # Finish formatting the climate table for this scenario
      data.climate.ij = data.climate.ij  %>% 
        dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, t)
      
      # Add to the output list
      # -- climate dataframe
      eval(parse(text = paste0(
        "list.out[[i]]$", ssp.in[j], "$climate = data.climate.ij")))
      # -- disturbance dataframe
      eval(parse(text = paste0(
        "list.out[[i]]$", ssp.in[j], "$disturbance = disturbance.ij")))
    }
    
  }
  
  # Return output list 
  return(list.out)
}

#' Function to select plots based on whether species can be simulated or not
#' @param NFI_data data of individual trees in NFI
#' @param NFI_disturbance data of disturbance probability per plot
#' @param NFI_climate data of climate per plot
#' @param nclim Number of climate in which to divide data
#' @param nplots.per.clim Maximum number of plots to select per climate
select_NFI_plots = function(NFI_data, NFI_disturbance, NFI_climate, 
                            nclim = 10, nplots.per.clim = 100){
  
  # Get mean climate per plot
  meanclim = NFI_climate %>% dplyr::select(plotcode)
  meanclim$sgdd = rowMeans(NFI_climate %>% dplyr::select(grep("sgdd", colnames(.))))  
  meanclim$wai = rowMeans(NFI_climate %>% dplyr::select(grep("wai", colnames(.))))  
  
  # Make a pca
  pca <- prcomp((meanclim %>% dplyr::select(-plotcode)), 
                center = T, scale = T)
  
  # Add pca coordinate to the dataframe
  meanclim$pca1 = get_pca_ind(pca)[[1]][, 1]
  
  # Group plots per climate
  meanclim = meanclim %>%
    mutate(quartile = ntile(pca1, n = nclim)) %>%
    mutate(climate = paste0("clim", quartile)) %>%
    dplyr::select(plotcode, sgdd, wai, pca1, climate)
  
  # Identify the sum of disturbances per scenario
  # -- Initialize with plotcode and climate
  disturbance.sum = NFI_disturbance %>% 
    dplyr::select(plotcode) %>%
    left_join(meanclim, by = "plotcode")
  # -- Add sum of fire probability
  disturbance.sum$fire = rowSums(
    NFI_disturbance %>% dplyr::select(colnames(.)[grep("fire", colnames(.))]))
  # -- Add sum of storm probability
  disturbance.sum$storm = rowSums(
    NFI_disturbance %>% dplyr::select(colnames(.)[grep("storm", colnames(.))]))
  # -- Add binary to indicate occurence of fire and storm per plot
  disturbance.sum = disturbance.sum %>%
    mutate(fire.bin = ifelse(fire > 0, 1, 0), storm.bin = ifelse(storm > 0, 1, 0))
  
  # Format NFI data to add compatibility in fire and storm sensitivity
  # -- Species which are included in the IPM
  species.ipm = fit_species[which(!(fit_species %in% c("Carpinus_betulus", "Quercus_ilex")))]
  # -- Identify species that are storm compatible and within IPM
  species.storm = (disturb_coef %>% 
                     filter(disturbance == "storm" & species %in% species.ipm))$species
  # -- Identify species that are storm compatible and within IPM
  species.fire = (disturb_coef %>% 
                    filter(disturbance == "fire" & species %in% species.ipm))$species
  # -- Format NFI data to get the proportion of basal area that can be simulated
  simulable.data = NFI_data %>%
    mutate(ba_ha = ba_ha/1000000, species = gsub("\\ ", "\\_", species)) %>%
    mutate(species = gsub("\\_sp", "", species)) %>%
    left_join((disturbance.sum %>% dplyr::select(plotcode, fire.bin, storm.bin)), 
              by = "plotcode") %>%
    mutate(fire.sp = ifelse(species %in% species.fire, 1, 0), 
           storm.sp = ifelse(species %in% species.storm, 1, 0), 
           ipm.sp = ifelse(species %in% species.ipm, 1, 0), 
           fire.ok = ifelse(fire.bin == 1 & fire.sp == 0, 0, 1), 
           storm.ok = ifelse(storm.bin == 1 & storm.sp == 0, 0, 1),
           simulable.sp = ipm.sp*fire.ok*storm.ok, 
           ba_ha.simulable = ba_ha*simulable.sp)
  
  # Final dataset to export
  data.out = meanclim %>%
    left_join((disturbance.sum %>% dplyr::select(plotcode, fire.bin, storm.bin)), 
              by = "plotcode") %>%
    left_join((simulable.data %>%
                 group_by(plotcode) %>%
                 summarize(prop.simul = sum(ba_ha.simulable)/sum(ba_ha))), 
              by = "plotcode") %>%
    mutate(selected = 0)
  
  # Loop on all climates
  for(i in 1:nclim){
    
    # Identify plotcode for climate i that are 100% simulable
    plot.simulable.i = (data.out %>%
                          filter(climate == paste0("clim", i) & prop.simul == 1))$plotcode
    
    # If more than the max number of plot per climate, sample some of them
    if(length(plot.simulable.i) > nplots.per.clim){
      plot.simulable.i = sample(plot.simulable.i, nplots.per.clim, replace = FALSE)
    }
    
    # Complete final dataset
    data.out = data.out %>% 
      mutate(selected = ifelse(plotcode %in% plot.simulable.i, 1, selected))
  }
  
  # Filter data.out
  data.out = data.out %>%
    filter(selected == 1) %>%
    dplyr::select(-prop.simul, -selected)
  
  # Return output
  return(data.out)
}


#' Extract the climatic margins for each climate
#' @param NFI_climate_sub Climate data for NFI plots, with only selected plots
#' @param NFI_plots_selected data frame with information on the plots selected
get_margins = function(NFI_climate_sub, NFI_plots_selected){
  
  # Identify colnames of sgdd and of wai
  sgdd.col = colnames(NFI_climate_sub)[grep("sgdd", colnames(NFI_climate_sub))]
  wai.col = colnames(NFI_climate_sub)[grep("wai", colnames(NFI_climate_sub))]
  
  # Identify all the different climates
  climate.in = unique(NFI_plots_selected$climate)
  
  # Initialize the output list
  list.out = vector(mode = "list", length = length(climate.in))
  names(list.out) = climate.in
  
  # Loop on all climates
  for(i in 1:length(climate.in)){
    
    # Extract margins
    margins.i = NFI_climate_sub %>%
      left_join((NFI_plots_selected %>% dplyr::select(plotcode, climate)), 
                by = "plotcode") %>%
      filter(climate == climate.in[i]) %>%
      dplyr::select(sgdd.col, wai.col) %>%
      gather(key = "variable", value = "value", sgdd.col, wai.col) %>%
      separate(col = "variable", into = c("scen", "var", "year"), sep = "\\_") %>%
      group_by(var) %>%
      summarize(min = min(value, na.rm = TRUE), 
                mean = mean(value, na.rm = TRUE), 
                max = max(value, na.rm = TRUE)) %>%
      gather(key = "stat", value = "value", "min", "mean", "max") %>%
      mutate(N = case_when(
        var == "sgdd" & stat == "min" ~ 3, 
        var == "sgdd" & stat == "max" ~ 1, 
        var == "wai" & stat == "min" ~ 1, 
        var == "wai" & stat == "max" ~ 3, 
        stat == "mean" ~ 2)) %>%
      dplyr::select(-stat) %>%
      spread(key = "var", value = "value") %>%
      mutate(PC1 = 0, PC2 = 0, sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
             waib = 1/(1 + wai), SDM = 0) %>%
      dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
    
    # Add to the output list
    list.out[[i]] = margins.i
    
  }
  
  # Return output list
  return(list.out)
  
}


#' Prepare the list of species objects to make
#' @param NFI_data_sub Tree-level data for NFI plots, with only selected plots
#' @param NFI_plots_selected data frame with information on the plots selected
make_species_list = function(NFI_data_sub, NFI_plots_selected){
  
  NFI_data_sub %>%
    left_join((NFI_plots_selected %>% dplyr::select(plotcode, climate)), 
              by = "plotcode") %>%
    dplyr::select(species, climate) %>%
    distinct() %>%
    arrange(climate, species) %>%
    mutate(species = gsub("\\ sp", "", species)) %>%
    mutate(ID.species = c(1:dim(.)[1])) %>%
    dplyr::select(ID.species, species, climate)
  
}


#' Function to make a species object, save it as rds and return filename
#' @param species_list table containing all species to create
#' @param climate_margins climatic margins for each climate (list)
#' @param ID.species.in ID of the species to make in species_list
make_species_mu_rds = function(species_list, climate_margins, ID.species.in){
  
  # Identify the species and case study for iteration i
  species.in = gsub("\\ ", "\\_", species_list$species[ID.species.in])
  climate.in = species_list$climate[ID.species.in]
  
  # Load demographic parameter of the species 
  eval(parse(text=paste0("fit.in <- fit_", species.in)))
  
  # From clim.in, build the margins file for iteration i
  margins.in = climate_margins[[climate.in]]
  
  # Make the mu matrix
  mu.in <- make_mu_gr(
    species = species.in, fit = fit.in, climate = margins.in, 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit.in) * 1.1),
    verbose = TRUE, stepMu = 0.001)
  
  # Create species object from random distribution
  sp.in = species(IPM = mu.in, init_pop = def_initBA(20),
                  harvest_fun = def_harv)
  # Update disturbance function
  sp.in$disturb_fun = disturb_fun
  # Add disturbance coefficients
  sp.in$disturb_coef  <- filter(
    matreex::disturb_coef, species == species.in)
  
  # Name of the file to save
  file.in = paste0("rds/", climate.in, "/species_mu/", species.in, ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sp.in, file.in)
  
  # Return output list
  return(file.in)
  
}



#' Function to extract the distribution of each species in each plot
#' @param NFI_data_sub tree-level NFI data for the plots selected
get_species_distrib = function(NFI_data_sub){
  
  # Initialize output list
  list.out = vector(mode = "list", length = length(unique(NFI_data_sub$plotcode)))
  names(list.out) = unique(NFI_data_sub$plotcode)
  
  # Correct species and basal area in NFI data 
  NFI_data.in = NFI_data_sub %>%
    mutate(ba_ha = ba_ha/1000000) %>%
    mutate(species = gsub("\\ ", "\\_", species), 
           species = gsub("\\_sp", "", species))
  
  # Loop on all plotcodes
  for(i in 1:length(names(list.out))){
    
    # Identify species present in plotcode i
    species.i = unique(subset(NFI_data.in, plotcode == names(list.out)[i])$species)
    
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
      data.ij = NFI_data.in %>%
        filter(plotcode == names(list.out)[i] & species == species.i[j]) %>%
        dplyr::select(dbh, ba_ha)
      
      # Loop on all trees present in the plot from species j
      for(k in 1:dim(data.ij)[1]){
        
        # Position of tree k in the mesh
        id.k = (mesh.ij %>% 
                  mutate(diff = abs(dbh.mesh - data.ij$dbh[k])) %>%
                  arrange(diff))$mesh.id[1]
        
        # Add basal area of tree k in the output list
        list.out[[i]][[j]][id.k] = list.out[[i]][[j]][id.k] + data.ij$ba_ha[k]
      }
      
    }
    
  }
  
  # Return output list
  return(list.out)
  
}


#' Make a list of simulations to make
#' @param NFI_plots_selected Information on the NFI plots selected
#' @param ssp.in ssp scenarios to simulate
make_simul_list = function(
  NFI_plots_selected, ssp.in = c("ssp126", "ssp370", "ssp585")){
  
  expand.grid(plotcode = NFI_plots_selected$plotcode, 
              ssp = ssp.in, 
              dist = c("dist", "nodist")) %>%
    left_join((NFI_plots_selected %>% dplyr::select(plotcode, climate, pca1)), 
              by = "plotcode") %>%
    mutate(ID.simulation = c(1:dim(.)[1])) %>%
    dplyr::select(ID.simulation, plotcode, pca1, ssp, dist, climate)
  
}


#' Make mu simulations with dist and changing climate
#' @param species_distrib List of species distribution (per plot)
#' @param species_list df with information on all species object
#' @param species_mu vector containing all species mu rds files created
#' @param climate_dist_dflist list of climate and disturbance df
#' @param simul_list df listing simulations to perform
#' @param ID.simulation ID in simul_dist of the simulation to perform
make_simulations = function(
  species_distrib, species_mu, species_list, climate_dist_dflist,
  simul_list, ID.simulation){
  
  # Identify the plotcode, climate, ssp, etc
  plot.in = simul_list$plotcode[ID.simulation]
  ssp.in = simul_list$ssp[ID.simulation]
  clim.in = simul_list$climate[ID.simulation]
  dist.in = simul_list$dist[ID.simulation]
  
  
  # Extract climate, disturbance and distributions
  climate.in = climate_dist_dflist[[plot.in]][[ssp.in]]$climate
  disturbance.in = climate_dist_dflist[[plot.in]][[ssp.in]]$disturbance
  distributions.in = species_distrib[[plot.in]]
  
  # Initialize list of species
  list.sp = vector(mode = "list", length = length(names(distributions.in)))
  names(list.sp) = paste0("mu_", names(distributions.in))
  
  # Correct species name in species list
  species_list = species_list %>% mutate(species = gsub("\\ ", "\\_", species))
  
  # Loop on all species of the list
  for(i in 1:length(names(list.sp))){
    
    # Name of species i
    spname.i = names(distributions.in)[i]
    
    # Get the index of mu from species list
    id.mu = subset(species_list, species == spname.i & climate == clim.in)$ID.species
    
    # Read the species
    list.sp[[i]] = readRDS(species_mu[id.mu])
    
    # Change initial distribution
    list.sp[[i]]$init_pop <- def_init_k(distributions.in[[i]]*0.03)
    
  }
  
  # Make forest
  forest.in = forest(species = list.sp, harv_rules = c(
    Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1))
  
  # Different code depending on whether the scenario includes disturbances or not
  if(dist.in == "dist"){
    # Run simulation with a changing climate and disturbance
    sim.in = sim_deter_forest(
      forest.in, tlim = max(climate.in$t), climate = climate.in, 
      equil_dist = max(climate.in$t),  equil_time = max(climate.in$t), 
      verbose = TRUE, correction = "cut", disturbance = disturbance.in)
  } else {
    # Run simulation with a changing climate but no disturbance
    sim.in = sim_deter_forest(
      forest.in, tlim = max(climate.in$t), climate = climate.in, 
      equil_dist = max(climate.in$t),  equil_time = max(climate.in$t), 
      verbose = TRUE, correction = "cut")
  }
  
  # Name of the file to save
  file.in = paste0("rds/", clim.in, "/simulations/simul", ID.simulation, ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sim.in, file.in)
  
  # Return output list
  return(file.in)
}


#' Disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Not used in this case.
#' \describe{
#' \item{qmd}{Forest Quadratic Mean Diameter}
#' }
#' @author Maxime Jeaunatre
#'
disturb_fun <- function(x, species, disturb = NULL, ...){
  
  dots <- list(...)
  qmd <- dots$qmd 
  size <- species$IPM$mesh
  coef <- species$disturb_coef
  if(any(disturb$type %in% coef$disturbance)){
    coef <- subset(coef, disturbance == disturb$type)
  } else {
    stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                 sp_name(species), disturb$type))
  }
  
  # edits for delay
  size[size == 0] <- min(size[size !=0])
  
  logratio <-  log(size / qmd)
  dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
  logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
  Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled + 
                    coef$b * disturb$intensity ^(coef$c * dbh.scaled))
  
  return(x* Pkill) # always return the mortality distribution
}



#' Function to extract for each timestep prodictivity, mean dbh and diversity
#' @param simulations simulation data
get_simulations_output = function(simulations){
  
  
  # Loop on all simulations
  for(ID.simulation in 1:length(simulations)){
    
    print(ID.simulation)
    
    # Read simulation
    sim.in = readRDS(simulations[ID.simulation])
    
    # Initialize dataset for productivity
    data.prod.in = sim.in %>%
      filter(!equil & var == "BAsp") %>%
      dplyr::select(species, time, BA = value) %>%
      rbind((sim.in %>%
               filter(!equil & var == "BAsp") %>%
               group_by(time) %>% summarize(BA = sum(value)) %>% 
               mutate(species = "all") %>% ungroup() %>%
               dplyr::select(species, time, BA))) %>%
      mutate(prod = 0, ID = c(1:dim(.)[1]))
    
    # Loop on all time step
    for(t in 2:max(data.prod.in$time)){
      # Loop on all species
      for(s in unique(data.prod.in$species)){
        id.t = subset(data.prod.in, species == s & time == t)$ID
        id.t1 = subset(data.prod.in, species == s & time == (t-1))$ID
        if(!is.na(data.prod.in$BA[id.t])){
          if(data.prod.in$BA[id.t] > data.prod.in$BA[id.t1]){
            data.prod.in$prod[id.t] = data.prod.in$BA[id.t] - data.prod.in$BA[id.t1]
          }
        }else{data.prod.in$prod[id.t] = NA_real_}
      }
    }
    
    # Data with structural information
    # data.dbh.in = sim.in %>%
    #   filter(var == "n" & !equil) %>%
    #   dplyr::select(size, time, species, value) %>%
    #   rbind((sim.in %>%
    #            filter(var == "n" & !equil) %>%
    #            group_by(size, time) %>% summarize(value = mean(value)) %>%
    #            ungroup() %>% mutate(species = "all") %>%
    #            dplyr::select(size, time, species, value))) %>%
    #   group_by(size, time, species) %>%
    #   summarize(ntot = sum(value)) %>%
    #   ungroup() %>% group_by(time, species) %>%
    #   filter(size > 0) %>%
    #   mutate(ntot_size = ntot*size) %>%
    #   summarize(dbh.mean = weighted.mean(size, w = ntot))
    
    # Diversity along time
    data.div.in = sim.in %>%
      filter(var == "BAsp" & !equil) %>%
      group_by(time) %>%
      mutate(p = value/sum(.$value), 
             plnp = p*log(p)) %>%
      summarise(H = -sum(plnp)) 
    
    # Final dataset
    data.out = data.prod.in %>%
      #left_join(data.dbh.in, by = c("time", "species")) %>%
      left_join(data.div.in, by = "time") %>%
      mutate(ID.simulation = ID.simulation) %>%
      dplyr::select(ID.simulation, time, species, prod, H)
    
    # Addd to final dataframe
    if(ID.simulation == 1) out = data.out
    else out = rbind(out, data.out)
    
  }
  
  # Return output
  return(out)
  
  
}

