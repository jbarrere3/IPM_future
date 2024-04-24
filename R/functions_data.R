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


#' Function to get parameter values per iteration for each species and disturbance
#' @param rdata.file rdata contianing the model outputs
get_param_from_rdata <- function(rdata.file){
  
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
             sp = as.integer(ifelse(Param == "a0", gsub(".+\\[", "", gsub("\\,.+", "", Parameter)), 
                                    gsub(".+\\[", "", gsub("\\]", "", Parameter)))), 
             co = as.integer(ifelse(Param == "a0", gsub(".+\\,", "", gsub("\\]", "", Parameter)), NA_integer_))) %>%
      # Remove the estimation of intensity and deviance
      filter(Param != "I") %>%
      filter(Param != "deviance") %>%
      # Add name of the country and species, and weight of each species per country
      left_join(corresp.tables[[i]]$country.table, by = "co") %>%
      left_join(corresp.tables[[i]]$species.table, by = "sp") %>%
      left_join(weight.tables[[i]], by = c("species", "country")) %>%
      # No weight for the parameters that do not rely on the country
      mutate(weight = ifelse(Param == "a0", weight, 1)) %>%
      mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
      # Summarize Parameter value per country (only apply to a0)
      group_by(disturbance, Iteration, Chain, Param, species) %>%
      summarize(val = sum(value*weight, na.rm = TRUE)/sum(weight, na.rm = TRUE)) %>%
      # Format to get one column per parameter
      spread(key = Param, value = val) %>%
      # Set a1 to 0 (dominance effect) if disturbance is not storm or snow
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, 0)) %>%
      # Add parameters to scale dbh and logratio
      mutate(dbh.intercept = scale.tables[[i]]$dbh.intercept, 
             dbh.slope = scale.tables[[i]]$dbh.slope, 
             logratio.intercept = scale.tables[[i]]$logratio.intercept, 
             logratio.slope = scale.tables[[i]]$logratio.slope)
    
    # Store table in the final table
    if(i == 1) param.table <- param.table.i
    else param.table <- rbind(param.table, param.table.i)
  }
  
  # return the list
  return(param.table)
  
}

#' #' Function to fit a model linking fire intensity to VPD
#' #' @param intensity_file file containing fire intensity values
#' #' @param NFI_data dataframe containin NFI data
#' #' @param NFI_climate data frame containing climate data for NFI
#' fit_vpd_fire = function(intensity_file, NFI_data, NFI_climate){
#'   
#'   # Load sensitivity file
#'   load(intensity_file)
#'   
#'   # Get intensity
#'   I_data = get_I_from_rdata(intensity_file)
#'   
#'   # Format dataset
#'   I.per.plot = I_data %>%
#'     filter(disturbance == "fire") %>%
#'     left_join(corresp.tables$fire$plotcode.table, by = "plot") %>%
#'     mutate(I.logit = log(I/(1 - I))) %>%
#'     group_by(plotcode) %>%
#'     summarize(I.mean = mean(I, na.rm = TRUE), 
#'               I.lwr = quantile(I, 0.025, na.rm = TRUE), 
#'               I.upr = quantile(I, 0.975, na.rm = TRUE), 
#'               I.logit.mean = mean(I.logit, na.rm = TRUE), 
#'               I.logit.var = var(I.logit, na.rm = TRUE))
#'   
#'   
#'   # Plot data 
#'   data.plot = NFI_data %>%
#'     filter(plotcode %in% I.per.plot$plotcode) %>%
#'     dplyr::select(plotcode, country, longitude, latitude, year1, year2) %>%
#'     distinct()
#'   
#'   # Extract vpd for each plot
#'   colnames.vpd = colnames(NFI_climate)[grep("vpd", colnames(NFI_climate))]
#'   data.vpd = NFI_climate %>%
#'     filter(plotcode %in% data.plot$plotcode) %>%
#'     dplyr::select(plotcode, colnames.vpd) %>%
#'     gather(key = "variable", value = "vpd", colnames.vpd) %>%
#'     separate(col = "variable", into = c("scen", "var", "year"), sep = "\\_") %>%
#'     filter(scen == "hist") %>%
#'     left_join((data.plot %>% dplyr::select(plotcode, year1, year2)), 
#'               by = "plotcode") %>%
#'     mutate(year = as.numeric(year)) %>%
#'     filter(year >= year1 & year <= year2) %>%
#'     group_by(plotcode) %>%
#'     summarize(vpd.mean = mean(vpd, na.rm = TRUE))
#'   
#'   # Add vpd to plot data
#'   data.plot = data.plot %>% left_join(data.vpd, by = "plotcode")
#'   
#'   # Final dataset
#'   data = data.plot %>% 
#'     left_join(I.per.plot, by = "plotcode") %>% 
#'     mutate(weight = 1/I.logit.var) %>%
#'     drop_na()
#'   
#'   # Fit model
#'   model = betareg(I.mean ~ vpd.mean, weights = weight, data = data, link = "logit")
#'   
#'   # Output list : model and the data
#'   out = list(model = model, data = data)
#'   
#'   # Return output list
#'   return(out)
#'   
#' }

#' Function to extract reference intensity per disturbance
#' @param intensity_file rdata file containing posterior of intensity
#' @param quantile.ref reference quantile to take from posterior I distribution
get_I_per_disturbance = function(intensity_file, quantile.ref){
  
  # Load sensitivity file
  load(intensity_file)
  
  # Get intensity per disturbance
  out = get_I_from_rdata(intensity_file) %>% 
    filter(disturbance %in% c("storm", "fire")) %>% 
    group_by(disturbance) %>% 
    summarize(quantile = quantile(I, quantile.ref))
  
  # Return output
  return(out)
  
}



#' #' Function to extract beta parameters of storm intensity distribution
#' #' @param intensity_file rdata file containing posterior of intensity
#' get_Istorm_param = function(intensity_file){
#'   
#'   # Load sensitivity file
#'   load(intensity_file)
#'   
#'   # Get intensity
#'   I_data = get_I_from_rdata(intensity_file)
#'   
#'   # Format dataset
#'   I.per.plot = I_data %>%
#'     filter(disturbance == "storm") %>%
#'     left_join(corresp.tables$storm$plotcode.table, by = "plot") %>%
#'     group_by(plotcode) %>%
#'     summarize(I.mean = mean(I, na.rm = TRUE))
#'   
#'   # Extract beta parameters for storm
#'   out = fitdistr(I.per.plot$I.mean, densfun = "beta", 
#'                  start = list(shape1 = 1, shape2 = 2))[[1]]
#'   
#'   # Return output
#'   return(out)
#'   
#' }


#' Function to prepare disturbance and climate dataframe for each plot
#' @param NFI_disturbance data of disturbance probability per plot
#' @param NFI_climate data of climate per plot
#' @param I_per_disturbance dataframe listing the intensity per disturbance agent
get_clim_dist_df = function(NFI_disturbance, NFI_climate, I_per_disturbance){
  
  # Initialize output list
  list.out = vector(mode='list', length=dim(NFI_disturbance)[1])
  names(list.out) = NFI_disturbance$plotcode
  
  # Table matching years and "time"
  year.table = data.frame(name = colnames(NFI_disturbance)) %>%
    filter(!(name %in% c("plotcode", "longitude", "latitude"))) %>%
    mutate(year = as.numeric(gsub(".+\\_", "", gsub("\\_mpi.+", "", name)))) %>%
    dplyr::select(year) %>%
    distinct() %>%
    arrange(year) %>%
    mutate(t = c(1:dim(.)[1]))
  
  # scenarios for which to generate climate files
  ssp.in = c("ssp126", "ssp585")
  rcp.in = c("rcp26", "rcp85")
  
  # Reference intensity for storm and fire
  Ifire = as.numeric(subset(I_per_disturbance, disturbance == "fire")$quantile)
  Istorm = as.numeric(subset(I_per_disturbance, disturbance == "storm")$quantile)
  
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
      filter(year %in% year.table$year & scenario %in% ssp.in & var != "vpd") %>%
      left_join(year.table, by = "year") %>%
      dplyr::select(t, scenario, var, value) %>%
      spread(key = "var", value = "value") %>%
      mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
             waib = 1/(1 + wai))
    
    # Prepare disturbance dataframe for plot i
    data.disturbance.i = NFI_disturbance %>%
      filter(plotcode == names(list.out)[i]) %>%
      dplyr::select(colnames(.)[grep("_", colnames(.))]) %>%
      gather(key = "variable", value = "value", colnames(.)) %>%
      separate(col = "variable", into = c("disturbance", "year", "model", "rcp"), sep = "_") %>%
      mutate(year = as.numeric(year)) %>%
      left_join(year.table, by = "year") %>%
      dplyr::select(t, rcp, disturbance, value)
    # Manage rcp for storm by duplicating observations 
    data.disturbance.i = data.disturbance.i %>%
      mutate(rcp = ifelse(is.na(rcp), rcp.in[1], rcp)) %>%
      rbind((data.disturbance.i %>% filter(disturbance == "storm") %>% 
               mutate(rcp = rcp.in[2]))) %>%
      # Finish formatting
      spread(key = "disturbance", value = "value") %>%
      # Convert rcp in ssp
      left_join(data.frame(scenario = ssp.in, rcp = rcp.in), 
                by = "rcp") %>% 
      dplyr::select(t, scenario, fire, storm) %>%
      drop_na()
    
    # Loop on all scenario
    for(j in 1:length(ssp.in)){
      # Initialize the list for this scenario
      eval(parse(text = paste0("list.out[[i]]$", ssp.in[j], " = list()")))
      
      # subset climate table and disturbance table to the right scenario
      data.climate.ij = data.climate.i %>% filter(scenario %in% c("hist", ssp.in[j]))
      data.disturbance.ij = data.disturbance.i %>% filter(scenario %in% c("hist", ssp.in[j]))
      
      # Initialize the disturbance table for scenario j
      disturbance.ij = data.frame(t = year.table$t, IsSurv = FALSE) %>%
        mutate(storm = rbinom(dim(.)[1], 1, data.disturbance.ij$storm), 
               fire = rbinom(dim(.)[1], 1, data.disturbance.ij$fire)) %>%
        mutate(type = case_when(fire == 1 & storm == 1 ~ sample(c("fire", "storm"), 1), 
                                fire == 1 & storm == 0 ~ "fire", 
                                fire == 0 & storm == 1 ~ "storm", 
                                TRUE ~ "none")) %>%
        filter(type != "none") %>%
        mutate(intensity = ifelse(type == "fire", Ifire, Istorm)) %>%
        dplyr::select(type, IsSurv, t, intensity)
      
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


#' Extend the climate list to apply disturbances over several years
#' @param climate_dist_dflist list of distrurbance and climate per plot simulated
#' @param dist.duration Over how many years should we apply the disturbance
extend_disturbance = function(climate_dist_dflist, dist.duration){
  
  # Vector of years simulated
  t.vec = climate_dist_dflist[[1]][[1]]$climate$t
  
  # Initialize output
  out = climate_dist_dflist
  
  # Loop on all plotcodes
  for(i in 1:length(names(climate_dist_dflist))){
    # Loop on all ssp scenarios
    for(j in 1:length(names(climate_dist_dflist[[i]]))){
      # Extract disturbance table ij
      dist.ij = climate_dist_dflist[[i]][[j]]$disturbance
      # Extend disturbance only if disturbance occured
      if(dim(dist.ij)[1] > 0){
        # Initialize the new disturbance dataframe
        dist.ij.new = dist.ij
        # Loop on all disturbances
        for(k in 1:dim(dist.ij)[1]){
          # Initialize extra years for disturbance k
          dist.ijk = data.frame(type = dist.ij$type[k], IsSurv = FALSE, 
                                t = c((dist.ij$t[k]+1):(dist.ij$t[k] + dist.duration - 1)), 
                                intensity = dist.ij$intensity[k])
          # If there are other disturbances after, remove overlapping disturbances
          if(k < dim(dist.ij)[1]) dist.ijk = filter(dist.ijk, t < dist.ij$t[k+1])
          # If it is the last disturbance, remove occurrences beyond maximum time
          if(k == dim(dist.ij)[1]) dist.ijk = filter(dist.ijk, t <= max(t.vec))
          # Add extra years to the new disturbance dataframe
          dist.ij.new = rbind(dist.ij.new, dist.ijk)
        }
        # Order the new disturbance dataframe
        dist.ij.new = dist.ij.new %>% arrange(t)
        # Replace by the new disturbance dataset
        out[[i]][[j]]$disturbance = dist.ij.new
      }
    }
  }
  
  # Return the new list 
  return(out)
}


#' Function to select plots based on whether species can be simulated or not
#' @param NFI_data data of individual trees in NFI
#' @param NFI_climate data of climate per plot
#' @param nclim Number of climate in which to divide data
#' @param nplots.per.clim Maximum number of plots to select per climate
select_NFI_plots = function(NFI_data, NFI_climate, 
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
  
  # Calculate the distribution of richness in each climate
  # - Calculate richness per plotcode and join climate
  distr.richness_per_clim = NFI_data %>%
    dplyr::select(plotcode, species) %>%
    distinct() %>%
    group_by(plotcode) %>% 
    summarize(richness = n()) %>%
    left_join(meanclim, by = "plotcode")
  # - Initialize data containing poisson parameter of richness distribution per clim
  data.distr = data.frame(climate = paste0("clim", c(1:nclim)), 
                          lambda = NA_real_)
  # - Loop on all climates to extract distribution
  for(i in 1:dim(data.distr)[1]) data.distr$lambda[i] = fitdistr(
    subset(distr.richness_per_clim, climate == data.distr$climate[i])$richness - 1, 
    densfun = "Poisson")$estimate
  
  
  # Format NFI data to add compatibility to IPM species: 
  # - Species that we can simulate
  sp.simulable = fit_species[-grep("Juniperus", fit_species)]
  # - Determine which plots can be simulated
  simulable.data = NFI_data %>%
    mutate(ba_ha = ba_ha/1000000, species = gsub("\\ ", "\\_", species)) %>%
    mutate(species = gsub("\\_sp", "", species)) %>%
    mutate(species = ifelse(species == "Betula_sp", "Betula", species)) %>%
    mutate(ipm.sp = ifelse(species %in% sp.simulable, 1, 0), 
           ba_ha.simulable = ba_ha*ipm.sp)
  
  # Final dataset to export
  data.out = meanclim %>%
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
    # Only keep the plots selected
    filter(selected == 1) %>%
    # Add plot coordinates
    left_join(NFI_climate[, c("plotcode", "latitude", "longitude")], 
              by = "plotcode") %>%
    # Add distribution of richness
    left_join(data.distr, by = "climate") %>%
    # Final selection of columns
    dplyr::select(plotcode, longitude, latitude, sgdd, pca1, climate, lambda)
  
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


#' get forest cover per plotcode from land cover file
#' @param land_cover_file file containing data on land cover per plotcode
get_forest_cover = function(land_cover_file){
  
  fread(land_cover_file) %>%
    mutate(forest_cov_1km = frac_82 + frac_83) %>%
    dplyr::select(plotcode, forest_cov_1km)
  
}


#' Function to identify the occurence of disturbances in each plotcode
#' @param climate_dist_dflist list of climate and disturbance df for simulations
get_dist_occurence = function(climate_dist_dflist){
  
  # Initialize the output dataframe
  out = data.frame(plotcode = names(climate_dist_dflist), 
                   fire.bin = 0, storm.bin = 0)
  
  # Loop on all plotcodes
  for(i in 1:dim(out)[1]){
    # Loop on all ssp scenarios
    for(j in 1:length(names(climate_dist_dflist[[i]]))){
      # Check the occurence of fire or of storm
      if("storm" %in% climate_dist_dflist[[i]][[j]]$disturbance$type) out$storm.bin[i] = 1
      if("fire" %in% climate_dist_dflist[[i]][[j]]$disturbance$type) out$fire.bin[i] = 1
    }
  }
  
  # Return the output
  return(out)
  
}

#' Build the regional pool of the plots selected
#' @param NFI_plots_selected dataframe listing the NFI plots selected for simulations
#' @param NFI_forest_cover Proportion of forest cover in a 1km radius around each plot
#' @param coef_ba_reg Coefficients to calculate the regional basal area
make_regional_pool = function(NFI_plots_selected, NFI_forest_cover, coef_ba_reg){
  
  # Initialize dataset for reg ba calculation
  data.in = expand.grid(plotcode = NFI_plots_selected$plotcode, 
                        species = coef_ba_reg$species) %>%
    # Add information on climate, parameters and land cover
    left_join(NFI_plots_selected, by = "plotcode") %>%
    left_join(coef_ba_reg, by = "species") %>%
    left_join(NFI_forest_cover, by = "plotcode") %>%
    # In case forest cover is missing for some plots, replace by its average
    mutate(forest_cov_1km = ifelse(is.na(forest_cov_1km), 0.6, forest_cov_1km)) %>%
    # Calculate regional basal area per species based on parameters
    mutate(ba_reg_th = a*exp(-(pca1 - b)^2/c), 
           ba_reg_th = ifelse(is.na(ba_reg_th), 0, ba_reg_th)) %>% 
    # Calculate probability of sampling for each species per plotcode
    group_by(plotcode) %>%
    mutate(prob = ba_reg_th/sum(ba_reg_th, na.rm = TRUE)) %>% ungroup() %>%
    # Correct regional basal area by forest cover
    mutate(ba_reg = ba_reg_th*forest_cov_1km)
  
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


#' Calculate dqm, classify plots in succession stage and fit distributions
#' @param NFI_data_sub tree level data for the plots simulated
#' @param NFI_plots_selected Plot level data for the plots simulated
classify_succession = function(NFI_data_sub, NFI_plots_selected){
  
  # Compile data to make classes of quadratic diameter per climate
  data_dqm = NFI_data_sub %>%
    mutate(dbh1.weight = (dbh^2)*Nha) %>%
    group_by(plotcode) %>%
    summarize(dqm = sqrt(sum(dbh1.weight, na.rm = TRUE)/sum(Nha, na.rm = TRUE))) %>%
    left_join(NFI_plots_selected[, c("plotcode", "climate")], by = "plotcode") %>%
    ungroup() %>% group_by(climate) %>%
    mutate(dqm033 = quantile(dqm, 0.33), dqm066 = quantile(dqm, 0.66)) %>%
    ungroup() %>%
    mutate(dqm_class = case_when(dqm < dqm033 ~ "succession1", 
                                 dqm > dqm066 ~ "succession3", 
                                 TRUE ~ "succession2")) %>%
    dplyr::select(plotcode, climate, dqm, dqm_class)
  
  
  # Initialize output dataframe
  out = data.frame(plotcode = unique(NFI_data_sub$plotcode), 
                   shape = NA_real_, scale = NA_real_)
  
  
  # Loop on all plot codes
  for(i in 1:dim(out)[1]){
    print(i)
    if(length(which(NFI_data_sub$plotcode == out$plotcode[i])) > 1){
      
      # - Vector if number of trees per ha
      Nha.i = round(subset(NFI_data_sub, plotcode == out$plotcode[i])$Nha, 
                    digits = 0)
      # - Vector of dbh
      dbh.i = round(subset(NFI_data_sub, plotcode == out$plotcode[i])$dbh, 
                    digits = 0)
      # - Initialize final vector of dbh
      vec.i = c()
      # - Loop on all trees
      for(j in 1:length(dbh.i)) vec.i = c(vec.i, rep(dbh.i[j], Nha.i[j]))
      # - Size distribution of plot i
      distr.i = 
        # - Add to final output
        try(out$shape[i] <- as.numeric(fitdistr(
          vec.i, densfun = "weibull", start = list(shape = 1, scale = 1))[[1]][1]), 
          silent = TRUE)
      try(out$scale[i] <- as.numeric(fitdistr(
        vec.i, densfun = "weibull", start = list(shape = 1, scale = 1))[[1]][2]), 
        silent = TRUE)
      
    }
  }
  
  # Add quadratic diameter and succession stage
  out = out %>%
    left_join(data_dqm, by = "plotcode") %>%
    dplyr::select(plotcode, climate, dqm_class, dqm, shape, scale)
  
  # Return output
  return(out)
  
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Prepare and run simulations  ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function to extend disturbance coefficients to all species
#' @param intensity_file posterior of disturbance coefficients from Barrere et al. 2023
#' @param traits dataframe containing wood density data per species
#' @param traitsNFI dataframe containing traits calculated from NFI data (bark thickness)
extend_disturb_coef = function(intensity_file, traits, traitsNFI){
  
  # Get the parameters of the global models
  param = get_param_from_rdata(intensity_file) 
  
  # Extract the parameters for Other broadleaf and other conifer
  param.other = param %>%
    filter(disturbance %in% c("storm", "fire")) %>%
    filter(species %in% c("Other broadleaf", "Other conifer")) %>%
    group_by(disturbance, species) %>%
    summarize(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c),
              dbh.intercept = mean(dbh.intercept),
              dbh.slope = mean(dbh.slope),
              logratio.intercept = mean(logratio.intercept),
              logratio.slope = mean(logratio.slope))
  
  
  
  # Fit a model with coefficients as a function of bark thickness
  # - Prepare data
  data.bt.fit = param %>%
    filter(disturbance == "fire") %>%
    ungroup() %>%
    group_by(species)  %>%
    summarise(a0 = mean(a0), b = mean(b), c = mean(c)) %>%
    left_join(dplyr::select(traitsNFI, "species", "BT" = "bark.thickness_mm"), 
              by = "species") %>%
    drop_na()
  # - Fit model
  mod.fire = lm(cbind(a0, b, c) ~ BT, data = data.bt.fit)
  # - Scale parameters
  scale.fire = param %>% ungroup() %>%
    filter(disturbance == "fire") %>%
    dplyr::select(logratio.intercept, logratio.slope, dbh.intercept, dbh.slope) %>%
    distinct()
  
  # Fit a model with coefficients as a function of bark thickness
  # - Prepare data
  data.wd.fit.storm = param %>%
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    filter(disturbance == "storm") %>%
    ungroup() %>%
    group_by(species)  %>%
    summarise(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c)) %>%
    left_join(dplyr::select(traits, "species", "WD" = "wood.density"), 
              by = "species") %>%
    drop_na()
  # - Fit model
  mod.storm = lm(cbind(a0, a1, b, c) ~ WD, data = data.wd.fit.storm)
  # - Scale parameters
  scale.storm = param %>% ungroup() %>%
    filter(disturbance == "storm") %>%
    dplyr::select(logratio.intercept, logratio.slope, dbh.intercept, dbh.slope) %>%
    distinct()
  
  # Slightly change disturb_coef to include the right name for Betula
  disturb_coef.2 = matreex::disturb_coef %>%
    filter(!(species %in% c("Betula_pubescens", "Betula_sp"))) %>%
    mutate(species = ifelse(species == "Betula_pendula", "Betula", species))
  # Identify species for which we already have fire or storm parameters
  sp.fire = unique(subset(disturb_coef.2, disturbance == "fire")$species)
  sp.storm = unique(subset(disturb_coef.2, disturbance == "storm")$species)
  # Identify species not in storm and fire
  sp.nostorm = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.storm)))]
  sp.nofire = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.fire)))]
  # Predict parameters based on traits
  out = disturb_coef.2 %>%
    # Add missing storm
    rbind((data.frame(disturbance = "storm", species = sp.nostorm) %>%
             left_join(dplyr::select(traits, "species", "WD" = "wood.density"), 
                       by = "species") %>%
             cbind(predict(mod.storm, newdata = .)) %>%
             drop_na() %>%
             mutate(dbh.intercept = scale.storm$dbh.intercept, 
                    dbh.slope = scale.storm$dbh.slope, 
                    logratio.intercept = scale.storm$logratio.intercept, 
                    logratio.slope = scale.storm$logratio.slope) %>%
             dplyr::select(-WD))) %>%
    # Add missing fire
    rbind((data.frame(disturbance = "fire", species = sp.nofire, a1 = 0) %>%
             left_join((traitsNFI %>% mutate(species = gsub("\\ ", "\\_", species)) %>%
                          dplyr::select("species", "BT" = "bark.thickness_mm")), 
                       by = "species") %>%
             cbind(predict(mod.fire, newdata = .)) %>%
             drop_na() %>%
             mutate(dbh.intercept = scale.storm$dbh.intercept, 
                    dbh.slope = scale.storm$dbh.slope, 
                    logratio.intercept = scale.storm$logratio.intercept, 
                    logratio.slope = scale.storm$logratio.slope) %>%
             dplyr::select(disturbance, species, a0, a1, b, c, dbh.intercept, 
                           dbh.slope, logratio.intercept, logratio.slope)))
  # Update species for which we don't have parameters
  sp.fire = unique(subset(out, disturbance == "fire")$species)
  sp.storm = unique(subset(out, disturbance == "storm")$species)
  sp.nostorm = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.storm)))]
  sp.nofire = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.fire)))]
  
  # Add to the final dataset using the parameters of other broadleaf / connifer
  out = out %>%
    rbind(rbind(data.frame(disturbance = "storm", sp = sp.nostorm), 
                data.frame(disturbance = "fire", sp = sp.nofire)) %>%
            mutate(genus = gsub("\\_.+", "", sp)) %>%
            mutate(species = ifelse(genus %in% c(
              "Acer", "Alnus", "Betula", "Carpinus", "Fagus", "Fraxinus", 
              "Populus", "Prunus", "Quercus", "Salix"), 
              "Other broadleaf", "Other conifer")) %>%
            left_join(param.other, by = c("species", "disturbance")) %>%
            dplyr::select(disturbance, species = sp, a0, a1, b, c, dbh.intercept, dbh.slope, 
                          logratio.intercept, logratio.slope))
  
  # Return output
  return(out)
  
}





#' Prepare the list of species objects to make
#' @param NFI_data_sub Tree-level data for NFI plots, with only selected plots
#' @param NFI_plots_selected data frame with information on the plots selected
#' @param regional_pool dataframe listing the regional pool per plot
make_species_list = function(NFI_data_sub, NFI_plots_selected, regional_pool){
  
  rbind(
    # Species per climate in the regional pool
    (regional_pool %>% 
             left_join(NFI_plots_selected[, c("plotcode", "climate")], 
                       by = "plotcode") %>%
             dplyr::select(species, climate)), 
    # Species per climate in the NFI data
    (NFI_data_sub %>%
       left_join((NFI_plots_selected %>% dplyr::select(plotcode, climate)), 
                 by = "plotcode") %>%
       mutate(species = gsub("\\ ", "\\_", species), 
              species = gsub("Betula\\_sp", "Betula", species)) %>%
       dplyr::select(species, climate))) %>%
    # Remove duplicates
    distinct() %>%
    arrange(climate, species) %>%
    mutate(species = gsub("\\ sp", "", species)) %>%
    mutate(ID.species = c(1:dim(.)[1])) %>%
    dplyr::select(ID.species, species, climate)
  
}


#' Function to make a species object, save it as rds and return filename
#' @param species_list table containing all species to create
#' @param climate_margins climatic margins for each climate (list)
#' @param disturb_coef_ext Extension of the disturb_coef data of matreex for all species
#' @param new_fit_list New demographic parameters
#' @param ID.species.in ID of the species to make in species_list
make_species_mu_rds = function(species_list, climate_margins, disturb_coef_ext, 
                               new_fit_list, ID.species.in){
  
  # Identify the species and case study for iteration i
  species.in = gsub("\\ ", "\\_", species_list$species[ID.species.in])
  climate.in = species_list$climate[ID.species.in]
  
  # Load demographic parameter of the species 
  fit.in = new_fit_list[[species.in]]
  
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
  sp.in$disturb_coef  <- filter(disturb_coef_ext, species == species.in)
  
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
  NFI_plots_selected, ssp.in = c("ssp126", "ssp585")){
  
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
  regional_pool, simul_list, disp_kernel, ID.simulation){
  
  # Print simulation ID
  print(ID.simulation)
  
  # Identify the plotcode, climate, ssp, etc
  plot.in = simul_list$plotcode[ID.simulation]
  ssp.in = simul_list$ssp[ID.simulation]
  clim.in = simul_list$climate[ID.simulation]
  dist.in = simul_list$dist[ID.simulation]
  
  
  # Extract climate, disturbance and distributions
  climate.in = climate_dist_dflist[[plot.in]][[ssp.in]]$climate
  disturbance.in = climate_dist_dflist[[plot.in]][[ssp.in]]$disturbance
  distributions.in = species_distrib[[plot.in]]
  
  # Final list of species based on regional pool and distribution.in
  sp.reg = filter(regional_pool, plotcode == plot.in)$species
  sp.tot = unique(c(names(distributions.in), sp.reg))
  # Compile to make a regional pool
  reg.pool.in = (data.frame(species = sp.tot) %>%
                   left_join(filter(regional_pool, plotcode == plot.in), by = "species") %>%
                   mutate(ba_reg = ifelse(is.na(ba_reg), 0, ba_reg)))$ba_reg
  names(reg.pool.in) = sp.tot
  # Vector of migration rates
  mig.rate.in = 1 - left_join(data.frame(species = sp.tot), disp_kernel, 
                              by = "species")$p30
  names(mig.rate.in) = sp.tot
  
  # Initialize list of species
  list.sp = vector(mode = "list", length = length(sp.tot))
  names(list.sp) = paste0("mu_", sp.tot)
  
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
  forest.in = forest(species = list.sp, harv_rules = c(
    Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1), 
    regional_abundance = reg.pool.in, migration_rate = mig.rate.in)
  
  # Different code depending on whether the scenario includes disturbances or not
  if(dist.in == "dist"){
    # Run simulation with a changing climate and disturbance
    try(sim.in <- sim_deter_forest(
      forest.in, tlim = max(climate.in$t), climate = climate.in, 
      equil_dist = max(climate.in$t),  equil_time = max(climate.in$t), 
      verbose = TRUE, correction = "cut", disturbance = disturbance.in), 
      silent = TRUE)
  } else {
    # Run simulation with a changing climate but no disturbance
    try(sim.in <- sim_deter_forest(
      forest.in, tlim = max(climate.in$t), climate = climate.in, 
      equil_dist = max(climate.in$t),  equil_time = max(climate.in$t), 
      verbose = TRUE, correction = "cut"), silent = TRUE)
  }
  
  # If the simulation failed, return an empty object
  if(!exists("sim.in")) sim.in = list()
  
  # Name of the file to save
  file.in = paste0("rds/", clim.in, "/simulations/simul", ID.simulation, ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(file.in)
  saveRDS(sim.in, file.in)
  
  # Return output list
  return(file.in)
}





#' Function to extract for each timestep prodictivity, mean dbh and diversity
#' @param simulations simulation data
get_simulations_output = function(simulations){
  
  
  # Loop on all simulations
  for(ID.simulation in 1:length(simulations)){
    
    print(ID.simulation)
    
    # Read simulation
    sim.in = readRDS(simulations[ID.simulation])
    
    # Check that the simulation didn't fail
    if(length(sim.in) > 0){
      
      # Check that there was no exponential increase in basal area
      if(!any(is.na(sim.in$value))){
        
        
        
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
        data.dbh.in = sim.in %>%
          filter(var == "n" & !equil) %>%
          dplyr::select(size, time, species, value) %>%
          rbind((sim.in %>%
                   filter(var == "n" & !equil) %>%
                   group_by(size, time) %>% summarize(value = mean(value)) %>%
                   ungroup() %>% mutate(species = "all") %>%
                   dplyr::select(size, time, species, value))) %>%
          group_by(size, time, species) %>%
          summarize(ntot = sum(value)) %>%
          ungroup() %>% group_by(time, species) %>%
          filter(size > 0) %>%
          mutate(ntot_size = ntot*size) %>%
          summarize(dbh.mean = weighted.mean(size, w = ntot), 
                    dbh.var = weighted.var(size, w = ntot))
        
        # Diversity along time
        data.div.in = sim.in %>%
          filter(var == "BAsp" & !equil) %>%
          group_by(time) %>%
          mutate(p = value/sum(value), 
                 plnp = p*log(p)) %>%
          summarise(H = -sum(plnp)) 
        
        # Final dataset
        data.out = data.prod.in %>%
          left_join(data.dbh.in, by = c("time", "species")) %>%
          left_join(data.div.in, by = "time") %>%
          mutate(ID.simulation = ID.simulation) %>%
          dplyr::select(ID.simulation, time, species, prod, H, BA, 
                        dbh.mean, dbh.var)
        
        # Addd to final dataframe
        if(ID.simulation == 1) out = data.out
        else out = rbind(out, data.out)
        
        
      }
    }
    
  }
  
  # Return output
  return(out)
  
}


#' Faster function to extract for each timestep species composition
#' @param simulations simulation data
get_simulations_output_short = function(simulations){
  
  
  # Loop on all simulations
  for(ID.simulation in 1:length(simulations)){
    
    print(ID.simulation)
    
    # Read simulation
    sim.in = readRDS(simulations[ID.simulation])
    
    # Check that the simulation didn't fail
    if(length(sim.in) > 0){
      
      # Check that there was no exponential increase in basal area
      if(!any(is.na(sim.in$value))){
        
        
        # Extract mean basal area and number of individuals per time and species
        data.out = sim.in %>%
          filter(!equil & var %in% c("BAsp", "N")) %>%
          mutate(ID.simulation = ID.simulation) %>% 
          dplyr::select(ID.simulation, species, time, var, value) %>%
          spread(key = "var", value = "value") %>%
          rename(BA = BAsp)
        
        # Addd to final dataframe
        if(ID.simulation == 1) out = data.out
        else out = rbind(out, data.out)
        
        
      }
    }
    
  }
  
  # Return output
  return(out)
  
}







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## functions to manage traits data  -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function to make imputation of functional traits from a trait df
#' @param traits.in input trait dataset
input_traits = function(traits.in){
  
  
  traits.matrix = as.matrix(traits.in[, -1])
  rownames(traits.matrix) = traits.in[, 1]
  
  # Try imputation by DINEOF
  test.eof = eofRecon(eof(traits.matrix, recursive = TRUE))
  
  # Output dataframe
  matrix.out = traits.matrix
  matrix.out[which(is.na(traits.matrix))] = test.eof[which(is.na(traits.matrix))]
  out = traits.in[, 1] %>% cbind(matrix.out)
  colnames(out) = colnames(traits.in)
  out = as.data.frame(out)
  for(i in 2:dim(out)[2]) out[, i] = as.numeric(out[, i])
  
  # Return output
  return(out)
}

#' Function to compile functional traits data and make imputation and PCA
#' @param wood.density_file WD data from Chave et al. 
#' @param traitsNFI_file traits calculated from NFI data (see Barrere et al. 2023)
#' @param shade.tolerance_file ST from Niinemets & Valladares 2006
#' @param TRY_file TRY request used in Barrere et al. 2023
#' @param sp.in.sim species included in the simulations
compile_traits = function(wood.density_file, traitsNFI_file, shade.tolerance_file, 
                          TRY_file, sp.in.sim){
  
  # Initilize output list
  list.out = list()
  
  # Wood density (Chave et al. 2008 + Dryad to quote)
  # - Read file
  wood.density <- read_xls(wood.density_file, sheet = "Data")
  # - CHange colnames
  colnames(wood.density) <- c("n", "family", "species", "wood.density", "region", "reference")
  # - Extract data for all simulated species
  wood.density = wood.density %>%
    mutate(species = gsub("\\ ", "\\_", species), 
           species = ifelse(species %in% c("Betula_pubescens", "Betula_pendula"), 
                            "Betula", species)) %>%
    filter(species %in% sp.in.sim) %>%
    group_by(species) %>%
    summarise(wood.density = mean(wood.density, na.rm = TRUE))
  
  # Shade tolerance (Niimenets)
  shade.tolerance <- fread(shade.tolerance_file) %>%
    mutate(shade.tolerance = as.numeric(gsub("\\,", "\\.", shade.tolerance.mean))) %>%
    dplyr::select(species, shade.tolerance) %>%
    mutate(species = gsub("\\ ", "\\_", species), 
           species = ifelse(species %in% c("Betula_pubescens", "Betula_pendula"), 
                            "Betula", species)) %>%
    filter(species %in% sp.in.sim)  %>%
    drop_na()
  
  # NFI traits
  NFI_traits = fread(traitsNFI_file) %>%
    mutate(species = gsub("\\ ", "\\_", species), 
           species = ifelse(species %in% c("Betula_pubescens", "Betula_pendula"), 
                            "Betula", species)) %>%
    filter(species %in% sp.in.sim) %>%
    rename(bark.thickness = bark.thickness_mm)
  
  # Leaf thickness from TRY
  traits.TRY <- fread(TRY_file) %>%
    filter(TraitID == 46) %>%
    filter(!is.na(StdValue)) %>%
    mutate(species = gsub("\\ ", "\\_", AccSpeciesName), 
           species = ifelse(species %in% c("Betula_pubescens", "Betula_pendula"), 
                            "Betula", species)) %>%
    filter(species %in% sp.in.sim)  %>%
    group_by(species) %>%
    summarize(leaf.thickness = mean(StdValue, na.rm = TRUE), 
              n = n()) %>%
    ungroup() %>% filter(n >= 5) %>%
    dplyr::select(-n) 
  
  
  # Add to the output list raw traits databases (with NA)
  # - Create element in the list
  list.out$traits_raw = list()
  # - traits related to the growth-survival trade-off
  list.out$traits_raw$GrSurv = data.frame(species = sp.in.sim) %>%
    left_join(wood.density, by = "species") %>%
    left_join((NFI_traits %>% dplyr::select(-bark.thickness)), 
              by = "species")
  # - traits related to the shadetolerance drought tolerance trade-off
  list.out$traits_raw$ShadeDrought = data.frame(species = sp.in.sim) %>%
    left_join(shade.tolerance, by = "species") %>%
    left_join((NFI_traits %>% dplyr::select(species, bark.thickness)), 
              by = "species") %>%
    left_join(traits.TRY, by = "species") %>%
    left_join((climate_species %>% filter(N == 2) %>%
                 dplyr::select(species = sp, optimum = PC1)), by = "species")
  
  # Add to the output list inputed traits database
  # - Crerate element in the list
  list.out$traits_imputed = list()
  # - Traits related to the growth survival trade-off
  list.out$traits_imputed$GrSurv = input_traits(list.out$traits_raw$GrSurv)
  # - Traits related to the shade - drought tolerance trade-off
  list.out$traits_imputed$ShadeDrought = input_traits(list.out$traits_raw$ShadeDrought)
  
  # Make PCAs based on imputed traits
  # - Crerate element in the list
  list.out$pca = list()
  # - Traits related to the growth survival trade-off
  list.out$pca$GrSurv = prcomp((list.out$traits_imputed$GrSurv %>% 
                                  dplyr::select(-species)), 
                               center = T, scale = T)
  # - Traits related to the shade - drought tolerance trade-off
  list.out$pca$ShadeDrought = prcomp((list.out$traits_imputed$ShadeDrought %>% 
                                        dplyr::select(-species)), 
                                     center = T, scale = T)
  
  # Extract traits coordinates for each pca
  # - Crerate element in the list
  list.out$traits_pca_coord = list()
  # - Traits related to the growth survival trade-off
  list.out$traits_pca_coord$GrSurv = data.frame(
    trait = as.character(rownames(get_pca_var(list.out$pca$GrSurv)[[1]])), 
    pca1 = as.numeric(get_pca_var(list.out$pca$GrSurv)[[1]][, 1]))
  # - Traits related to the shade - drought tolerance trade-off
  list.out$traits_pca_coord$ShadeDrought = data.frame(
    trait = as.character(rownames(get_pca_var(list.out$pca$ShadeDrought)[[1]])), 
    pca1 = as.numeric(get_pca_var(list.out$pca$ShadeDrought)[[1]][, 1]))
  
  # Make a dataframe with the names of each axis
  list.out$title_axes = data.frame(
    axis = c("GrSurv", "ShadeDrought"), title = c(
      paste0("GrSurv Axis (", round(summary(
        list.out$pca$GrSurv)$importance[2, 1]*100, digits = 2), "%)",
        "\nHigh growth <--> High survival"), 
      paste0("ShadeDrought Axis (", round(summary(
        list.out$pca$ShadeDrought)$importance[2, 1]*100, digits = 2), "%)",
        "\nHigh shade tol. <--> High drought tol."))
  )
  
  # Data frame with the coordinates of each species along each axis
  list.out$species_coord = data.frame(species = sp.in.sim) %>%
    # Add coordinates on the growth survival axis
    left_join(data.frame(
      species = as.character(rownames(get_pca_ind(list.out$pca$GrSurv)[[1]])), 
      GrSurv = as.numeric(get_pca_ind(list.out$pca$GrSurv)[[1]][, 1])), 
      by = "species") %>%
    # Add coordinates on the drought-shade axis
    left_join(data.frame(
      species = as.character(rownames(get_pca_ind(list.out$pca$ShadeDrought)[[1]])), 
      ShadeDrought = as.numeric(get_pca_ind(list.out$pca$ShadeDrought)[[1]][, 1])), 
      by = "species")
  
  # Return the output list
  return(list.out)
  
  
}



#' Function to extract demographic traits from the model parameters
#' @param sp.in.sim Species included in the simulations
#' @param new_fit_list New list of parameters (based on the new recruitment fit)
get_traits_demo = function(sp.in.sim, new_fit_list, coef_ba_reg){
  
  # Load climatic optimum for each species
  data("climate_species")
  
  # Calculate the mean dbh, competition (hetero and conspecific)
  dbh.ref = 250
  BASP.ref = 20
  BANONSP.ref = 10
  BATOT.ref = BASP.ref + BANONSP.ref
  
  # Provide mean value of competition across all species
  data_species = data.frame(species = sp.in.sim, 
                            BATOTSP = BASP.ref, 
                            BATOTNonSP = BANONSP.ref, 
                            logBATOTSP = log(BASP.ref),
                            BATOTcomp = BATOT.ref,
                            size = dbh.ref,
                            logsize = log(dbh.ref),
                            intercept = 1) %>%
    # Add climatic data
    left_join((climate_species %>%
                 filter(N == 2) %>%
                 dplyr::select(species = sp, wai, wai2, 
                               waib, sgdd, sgdd2, sgddb, PC1)), 
              by = "species") %>%
    # Add regional basal area and dispersion kernel
    left_join(coef_ba_reg, by = "species") %>%
    mutate(BA.reg = a*exp(-(PC1 - b)^2/c)) %>%
    left_join(disp_kernel, by = "species") %>%
    # Calculate fecundity
    mutate(Fec = log(p30*BATOTSP + (1 - p30)*BA.reg*0.9)) %>%
    dplyr::select(-p30, -a, -b, -c, -PC1) %>%
    # Add interactions
    mutate(`logsize:sgdd` = logsize*sgdd, `logsize:wai` = logsize*wai, 
           `size:sgdd` = size*sgdd, `size:wai` = size*wai, 
           `BATOTcomp:sgdd` = BATOTcomp*sgdd, `BATOTcomp:wai` = BATOTcomp*wai, 
           `logsize:sgddb` = logsize*sgddb, `logsize:waib` = logsize*waib, 
           `size:sgddb` = size*sgddb, `size:waib` = size*waib, 
           `BATOTcomp:sgddb` = BATOTcomp*sgddb, `BATOTcomp:waib` = BATOTcomp*waib)
  
  # Initialize the output dataset
  traits_demo = data.frame(species = sp.in.sim, recruitment = NA_real_, 
                           delay = NA_real_, growth = NA_real_, survival = NA_real_)
  
  # Loop on all species to gather traits
  for(i in 1:dim(traits_demo)[1]){
    
    # Species i
    sp.i = traits_demo$species[i]
    
    # Give the delay from fit list
    traits_demo$delay[i] = as.numeric(new_fit_list[[i]]$info["delay"])
    
    # Demographic parameters for species i
    param.rec.i = new_fit_list[[i]]$rec$params_m
    names(param.rec.i)[which(names(param.rec.i) == "logBATOTSP")] = "Fec"
    param.gr.i = new_fit_list[[i]]$gr$params_m
    param.sv.i = new_fit_list[[i]]$sv$params_m
    
    # Associated vector of variables
    val.rec.i = as.vector(subset(data_species, species == sp.i)[, names(param.rec.i)])
    val.gr.i = as.vector(subset(data_species, species == sp.i)[, names(param.gr.i)])
    val.sv.i = as.vector(subset(data_species, species == sp.i)[, names(param.sv.i)])
    
    # Get the demographic parameters
    traits_demo$recruitment[i] = exp(sum(param.rec.i*val.rec.i))
    traits_demo$growth[i] = exp(sum(param.gr.i*val.gr.i))
    traits_demo$survival[i] = 1 - plogis(sum(param.sv.i*val.sv.i))
    
  }
  
  # Return final dataframe
  return(traits_demo)
  
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## functions for matreex -------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


#' species are initiated with their own functions and only the regional abundance
#' is set in the forest object as well as the migration rate

#' Constructor for forest class
#'
#' Only used in the matreex package
#'
#' @param species List of species created with matreex package.
#' @param harv_rules Vector for harvest rules at the scale of the forest.
#' \describe{
#'   \item{Pmax}{maximum proportion of BAcut / BA}
#'   \item{dBAmin}{the minimum BA to perform cut}
#'   \item{freq}{Frequence at which the harvest will be executed.}
#' }
#' @param regional_abundance list of vector with size distribution for each species.
#' This is not a direct species relative abundance, but I don't know how to implement this...help ,
#' @param migration_rate numeric vector with a migration rate in percentage between 1 et and 0.
#'
#' @importFrom purrr map_chr
#'
#' @keywords internal
#' @export
new_forest <- function(species = list(),
                       harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
                       regional_abundance = NULL,
                       migration_rate = NULL
){
  
  sp <- map_chr(species, sp_name)
  names(species) <- sp
  if(!is_null(regional_abundance)){
    names(regional_abundance) <- sp
    names(migration_rate) <- sp
  }
  forest <- list(
    species = species, harv_rules = harv_rules,
    info = list(species = sp,
                clim_lab = map_chr(species, climatic)),
    regional_abundance = regional_abundance,
    migration_rate = migration_rate
  )
  
  if(!is_null(regional_abundance)){
    class(forest) <- c("forest", "reg_forest")
  } else {
    class(forest) <- "forest"
  }
  
  return(forest)
}

#' validator for forest class.
#'
#' @param x forest class object
#'
#' @import checkmate
#'
#' @noRd
validate_forest <- function(x){
  
  regional <- inherits(x, "reg_forest")
  values <- unclass(x)
  names <- attr(x, "names")
  
  #map(values$species, validate_species)
  # TODO check forest harv rules
  
  clim_lab <- values$info$clim_lab
  if(length(unique(clim_lab)) > 1){
    clim_ipm <- clim_lab[clim_lab != "mu_gr"]
    if(length(clim_ipm) > 1){ # D & F
      stop(paste0("Some ipm species are not defined with the same climatic name.",
                  "Check it with : map_chr(species, climatic)"))
    }
  }
  
  # check the regional pool settings
  if(regional){
    
    assertNumeric(values$migration_rate, len = length(values$species),
                  lower = 0, upper = 1)
    if(all(values$migration_rate == 0)){
      warning("All migration rate are 0, the regional pool of this forest is deleted")
      x$regional_abundance <- NULL
      x$migration_rate <- NULL
      class(x) <- "forest"
      
      return(invisible(x))
    }
    
    assertSubset(names(values$migration_rate), names(values$species))
    # length_meshs <- map_dbl(values$species, ~ length(.x$IPM$mesh))
    
    # assertList(values$regional_abundance, types = "numeric",
    # len = length(values$species))
    # if(any(lengths(values$regional_abundance) != length_meshs)){
    # stop("regional abundance numeric vector should be the length of the species mesh.")
    # }
    
    assertSubset(names(values$regional_abundance), names(values$species))
  }
  
  invisible(x)
}

#' Create a new forest for simulation
#'
#' A forest is a group of one of multiple species to silumate along time using
#' the IPM defined for each species and harvest rules.
#'
#' @inheritParams new_forest
#'
#' @export
forest <- function(species = list(),
                   harv_rules = c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1),
                   regional_abundance = NULL,
                   migration_rate = NULL
){
  
  res <- validate_forest(new_forest(
    species = species,
    harv_rules = harv_rules,
    regional_abundance = regional_abundance,
    migration_rate = migration_rate
  ))
  
  return(res)
}

