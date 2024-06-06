#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot for methods ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to plot the functional and climate position of each species
#' @param NFI_data_sub Subset of NFI data for the plots selected
#' @param NFI_plots_selected NFI plots that were selected for the analyses
#' @param traits_compiled list compiling traits data, pca and functional position
#' @param file.out Name of the file to save, including path
plot_funclim_species = function(NFI_data_sub, NFI_plots_selected, 
                                traits_compiled, file.out){
  
  # Create output directory if necessary
  create_dir_if_needed(file.out)
  
  ## -- Step 1 - Species plot
  
  # Compile species level data
  data.species = NFI_data_sub %>%
    dplyr::select(plotcode, species) %>%
    distinct() %>%
    left_join(NFI_plots_selected[, c("plotcode", "pca1")], by = "plotcode") %>%
    group_by(species) %>%
    summarize(Cold = quantile(pca1, 0.975, na.rm = TRUE), 
              Optimum = mean(pca1, na.rm = TRUE), 
              Hot = quantile(pca1, 0.025, na.rm = TRUE)) %>%
    left_join((traits_compiled$species_coord %>%
                 mutate(species = gsub("\\_", "\\ ", species)) %>%
                 mutate(species = ifelse(species == "Betula", "Betula sp", species))), 
              by = "species") %>%
    gather(key = "axis", value = "FunValue", "GrSurv", "ShadeDrought") %>%
    left_join(traits_compiled$title_axes, by = "axis")
  
  # Compile climate data
  data.clim = NFI_plots_selected %>%
    group_by(climate) %>%
    summarize(Hot = min(pca1, na.rm = TRUE), 
              Cold = max(pca1, na.rm = TRUE)) %>%
    merge(traits_compiled$title_axes) %>%
    left_join((traits_compiled$species_coord %>%
                 gather(key = "axis", value = "value", 
                        "GrSurv", "ShadeDrought") %>%
                 group_by(axis) %>%
                 summarize(min = min(value), max = max(value))), 
              by = "axis") %>%
    mutate(climate = factor(climate, levels = paste0(
      "clim", c(1:length(unique(.$climate))))), 
      FunValue = NA_real_) %>%
    mutate(Cold = ifelse(Cold > max, max, Cold))
  
  
  ## -- Step 2 - Functional traits plot
  
  # Initialize output list
  plotlist.out = list()
  
  # Loop on the two variables
  for(i in 1:length(names(traits_compiled$traits_raw))){
    
    # Identify the variable name i
    name.i = names(traits_compiled$traits_raw)[i]
    
    # Table i of pca coordinates 
    coord.i = traits_compiled$traits_pca_coord[[i]]
    
    # Species data subsetted for axis i
    df.i = subset(data.species, axis == name.i)
    
    # Make trait plot for variable i
    plot.trait.i = coord.i  %>%
      mutate(pos = c(1:dim(.)[1]), 
             pca1 = pca1*(min(abs(range(df.i$FunValue)))/max(abs(coord.i$pca1)))) %>%
      ggplot() + 
      geom_segment(aes(x = pos, xend = pos, y = 0, yend = pca1), 
                   arrow = arrow(length = unit(0.2, "cm"), type = "closed")) + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      scale_x_continuous(breaks = c(1:dim(coord.i)[1]), 
                         label = coord.i$trait, 
                         limits = c(0.5, dim(coord.i)[1] + 0.5)) + 
      ylab(traits_compiled$title_axes$title[i]) + xlab("") +
      ylim(range(df.i$FunValue)) +
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    # Make species plot i
    plot.species.i = df.i %>%
      ggplot(aes(y = FunValue)) + 
      geom_rect(data = (data.clim %>% filter(axis == name.i)), 
                aes(fill = climate, xmin = Hot, xmax = Cold, ymin = min, ymax = max), 
                color = "#495057", inherit.aes = TRUE, alpha = 0.5) + 
      geom_segment(aes(x = Hot, xend = Cold, yend = FunValue)) + 
      geom_point(aes(x = Optimum), shape = 21, color = "black", fill = "grey") + 
      geom_point(aes(x = Hot), shape = 21, color = "black", fill = "red") + 
      geom_point(aes(x = Cold), shape = 21, color = "black", fill = "blue") + 
      scale_fill_manual(values = colorRampPalette(c("#FFBF69", "#C5D86D", "#AEB8FE"))(10))+ 
      geom_text(aes(x = Optimum, y = FunValue - 0.1, label = species), size = 2) + 
      xlab("Position along the climate axis\n(hot-dry to cold-wet)") + 
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), 
            axis.title.y = element_blank(), 
            legend.position = "none")
    
    # Assemble in one single plot
    plot.i = ggarrange(plot.trait.i, plot.species.i, nrow = 1, widths = c(0.3, 1))
    
    
    
    # Add to the output list
    eval(parse(text = paste0("plotlist.out$", name.i, " = plot.i")))
    
  }
  
  # Make the final plot
  plot.out = plot_grid(plotlist = plotlist.out, nrow = 1, scale = 0.95, 
                       labels = c("(a)", "(b)"), align = "h")
  
  # Save the plot
  ggsave(file.out, plot.out, width = 25, height = 10 , units = "cm", 
         dpi = 600, bg = "white")
  
  
  # Return file saved
  return(file.out)
  
}






#' Plot the change in disturbance frequency along the climatic gradient
#' @param NFI_plots_selected df with information on the NFI plots selected
#' @param climate_dist_dflist climate and disturbance regimes per plot simulated
#' @param file.out Name of the file to save, including path
plot_map_clim_dist = function(NFI_plots_selected, climate_dist_dflist, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Vector of color for plotting
  color.vec = colorRampPalette(c("#FFBF69", "#C5D86D", "#AEB8FE"))(
    length(unique(NFI_plots_selected$climate)))
  names(color.vec) = paste0("clim", c(1:length(color.vec)))
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%
  ## -- First step : Map 
  
  # Convert the NFI dataframe in sf format
  data_sf = NFI_plots_selected  %>%
    mutate(climate = factor(climate, levels = names(color.vec))) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Make the map
  plot.map = ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(fill = "#343A40", color = "gray", show.legend = F, size = 0.2) + 
    geom_sf(data = data_sf, aes(color = climate), size = 1, shape = 20) +
    scale_color_manual(values = color.vec) +
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(), 
          legend.title = element_blank(), 
          legend.key = element_blank()) + 
    guides(color = guide_legend(override.aes = list(size=5, alpha = 0.85)))
  
  
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%
  ## -- Second step : climate and disturbance 
  
  
  # Identify the duration of simulations
  t.sim = dim(climate_dist_dflist[[1]][[1]]$climate)[1]
  
  # Periods simulated and plotted
  years.sim = c(1991:2100)
  periods = c("2020 - 2060", "2060 - 2100")
  
  # Initialize the output list
  list.out = vector(mode = "list", length = length(periods)) 
  names(list.out) = periods
  
  # Loop on all periods
  for(j in 1:length(list.out)){
    
    # Range of the period
    years.range.j = as.numeric(strsplit(periods[j], " - ")[[1]])
    
    # Years in period j
    years.j = c(years.range.j[1]:years.range.j[2])
    
    # Convert in a period of time based on the period simulated
    time.j = which(years.sim %in% years.j)
    
    # Initialize the data set for period j
    list.out[[j]] = expand.grid(plotcode = names(climate_dist_dflist), 
                                ssp = names(climate_dist_dflist[[1]]),
                                freq.storm = 0, freq.fire = 0, sgdd = NA, wai = NA)
    
    # Loop on all plotcode - ssp combination
    for(i in 1:dim(list.out[[j]])[1]){
      
      # Disturbance dataframe corresponding to simulation i
      dist.df.ij = climate_dist_dflist[[
        list.out[[j]]$plotcode[i]]][[list.out[[j]]$ssp[i]]]$disturbance %>%
        filter(t %in% time.j)
      
      # Count the occurrences of storm 
      if("storm" %in% dist.df.ij$type) list.out[[j]]$freq.storm[i] = length(
        which(dist.df.ij$type == "storm"))/length(time.j)
      
      # Count the occurrences of fire 
      if("fire" %in% dist.df.ij$type) list.out[[j]]$freq.fire[i] = length(
        which(dist.df.ij$type == "fire"))/length(time.j)
      
      # Average climate for plotcode i and period j
      # - sgdd
      list.out[[j]]$sgdd[i] = mean(climate_dist_dflist[[
        list.out[[j]]$plotcode[i]]][[list.out[[j]]$ssp[i]]]$climate$sgdd[time.j])
      # - wai
      list.out[[j]]$wai[i] = mean(climate_dist_dflist[[
        list.out[[j]]$plotcode[i]]][[list.out[[j]]$ssp[i]]]$climate$wai[time.j])
      
    }
    
  }
  
  # Combine the data into a dataframe
  data.out = bind_rows(list.out, .id = "period") %>%
    left_join(NFI_plots_selected[, c("plotcode", "climate", "pca1")], 
              by = "plotcode") %>%
    # Calculate average climate and disturbance per ssp, climate and period
    group_by(climate, period, ssp) %>%
    summarize(Hot = min(pca1, na.rm = TRUE), 
              Cold = max(pca1, na.rm = TRUE), 
              pca1.mean = mean(pca1, na.rm = TRUE), 
              sgdd.mean = mean(sgdd, na.rm = TRUE), 
              sgdd.sd = sd(sgdd, na.rm = TRUE), 
              wai.mean = mean(wai, na.rm = TRUE), 
              wai.sd = sd(wai, na.rm = TRUE), 
              freq.fire.mean = mean(freq.fire, na.rm = TRUE), 
              freq.fire.sd = sd(freq.fire, na.rm = TRUE)/sqrt(n()), 
              freq.storm.mean = mean(freq.storm, na.rm = TRUE), 
              freq.storm.sd = sd(freq.storm, na.rm = TRUE)/sqrt(n())) %>%
    # Calculate upper and lower boundaries for frequency and climate
    mutate(sgdd_upr = sgdd.mean + sgdd.sd, sgdd_lwr = sgdd.mean - sgdd.sd, 
           wai_upr = wai.mean + wai.sd, wai_lwr = wai.mean - wai.sd, 
           freq.fire_upr = freq.fire.mean + freq.fire.sd, 
           freq.storm_upr = freq.storm.mean + freq.storm.sd, 
           freq_upr = max(freq.fire_upr, freq.storm_upr)) %>%
    ungroup() %>% 
    mutate(sgdd.upr = max(sgdd_upr), sgdd.lwr = min(sgdd_lwr), 
           wai.upr = max(wai_upr), wai.lwr = min(wai_lwr), 
           freq.lwr = 0, freq.upr = max(freq_upr)) %>%
    ungroup() %>% dplyr::select(-sgdd_upr, -sgdd_lwr, -wai_upr, -wai_lwr, -freq_upr, 
                                - freq.fire_upr, - freq.storm_upr) %>%
    # Reduce the size of rectangles for the most extreme climates
    mutate(Hot = ifelse(climate == "clim1", pca1.mean, Hot), 
           Cold = ifelse(climate == "clim10", pca1.mean, Cold))
  
  
  # Plot the change in disturbance frequency
  plot.dist = data.out %>%
    gather(key = "metric", value = "value", "freq.fire.mean", "freq.storm.mean", 
           "freq.fire.sd", "freq.storm.sd") %>%
    mutate(metric = gsub("freq\\.", "", metric)) %>%
    separate(col = "metric", into = c("disturbance", "variable"), "\\.") %>%
    spread(key = "variable", value = "value") %>%
    ggplot(aes(x = pca1.mean, y = mean, color = disturbance, linetype = ssp)) + 
    geom_rect(data = mutate(data.out, disturbance = NA_character_, mean = NA_real_), 
              aes(fill = climate, xmin = Hot, xmax = Cold, ymin = freq.lwr, 
                  ymax = freq.upr), color = NA, inherit.aes = TRUE, 
              alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = color.vec) + 
    geom_ribbon(aes(ymin = mean-sd, ymax = mean+sd, alpha = ssp), fill = "grey") + 
    geom_line() + 
    facet_wrap(~ period) + 
    scale_color_manual(values = c('fire' = "red", 'storm' = "blue")) + 
    scale_linetype_manual(values = c('ssp126' = "dashed", 'ssp585' = "solid")) +
    scale_alpha_manual(values = c('ssp126' = 0.35, 'ssp585' = 0.65)) +
    xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
    ylab("Disturbance frequency") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank())
  
  # Plot the change in sgdd
  plot.sgdd = data.out %>%
    ggplot(aes(x = pca1.mean, y = sgdd.mean, ymin = sgdd.mean - sgdd.sd, 
               ymax = sgdd.mean + sgdd.sd, linetype = ssp, alpha = ssp)) + 
    geom_rect(aes(fill = climate, xmin = Hot, xmax = Cold, ymin = sgdd.lwr, 
                  ymax = sgdd.upr), color = NA, alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = color.vec) + 
    geom_ribbon() +
    geom_line() + 
    facet_wrap(~ period) + 
    scale_alpha_manual(values = c('ssp126' = 0.2, 'ssp585' = 0.4)) + 
    scale_linetype_manual(values = c('ssp126' = "dashed", 'ssp585' = "solid")) +
    xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
    ylab("Average sum of \ndegree days (sgdd)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank())
  
  # Plot the change in wai
  plot.wai = data.out %>%
    ggplot(aes(x = pca1.mean, y = wai.mean, ymin = wai.mean - wai.sd, 
               ymax = wai.mean + wai.sd, linetype = ssp, alpha = ssp)) + 
    geom_rect(aes(fill = climate, xmin = Hot, xmax = Cold, ymin = wai.lwr, 
                  ymax = wai.upr), color = NA, alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = color.vec) + 
    geom_ribbon() +
    geom_line() + 
    facet_wrap(~ period) + 
    scale_alpha_manual(values = c('ssp126' = 0.2, 'ssp585' = 0.4)) + 
    scale_linetype_manual(values = c('ssp126' = "dashed", 'ssp585' = "solid")) +
    xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
    ylab("Average water\navailability index (wai)") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank())
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%
  ## -- Step 3 : assemble plots
  
  # Make the final plot
  plot.out = plot_grid(
    plot.map, plot_grid(plot.sgdd, plot.wai, plot.dist, ncol = 1, scale = 0.9,
                        labels = c("(b)", "(c)", "(d)"), align = "v"), 
    nrow = 1, labels = c("(a)", ""))
  
  
  # Save the plot
  ggsave(file.out, plot.out, width = 30, height = 18 , units = "cm", 
         dpi = 600, bg = "white")
  
  
  # Return file saved
  return(file.out)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot for statistical analyses ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#' Function to plot the effect of regional pool on compositional shifts
#' @param regional_pool regional pool per plotcode
#' @param sim_output_pool data frame containing output of the simulations with regional pool
#' @param sim_output_nopool data frame containing output of the simulations without regional pool
#' @param simul_list data frame containing climate and plotcode per simulation
#' @param NFI_data_sub subset of NFI data with the plots selected only
#' @param traits_compiled List with traits data, pca and functional axes
#' @param dist_occurence Occurrence of disturbance in each plot simulated
#' @param metric.ref Reference metric to use to quantify abundance ("BA" or "N")
#' @param file.out Name of the file to save, including path
plot_pool_effect = function(
    regional_pool, sim_output_pool, sim_output_nopool, simul_list, NFI_data_sub, 
    traits_compiled, dist_occurence, metric.ref, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Merge simulation output with or without regional pool
  sim_output = bind_rows(list(with_pool = sim_output_pool, 
                              without_pool = sim_output_nopool), .id = "dispersal")
  
  # Calculate functional diversity
  data.FD = sim_output %>%
    gather(key = "metric", value = "weight", "N", "BA") %>%
    # Only focus on the second half of the simulations
    filter(time >= 70) %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Gather by functional axis
    gather(key = "axis", value = "trait_value", names(traits_compiled[[1]])) %>%
    # Calculate cwm for each time step and each simulation
    group_by(dispersal, ID.simulation, time, metric, axis) %>%
    mutate(cwm = weighted.mean(trait_value, w = weight)) %>%
    ungroup() %>%
    # Calculate square distance of species to centroid along each axis
    mutate(zsq_per_axis = (trait_value - cwm)^2) %>%
    # Calculate the distance of each species to the centroid
    group_by(dispersal, ID.simulation, time, species, metric, weight) %>%
    summarize(z = sqrt(sum(zsq_per_axis))) %>% ungroup() %>%
    # Calculate functional dispersion
    group_by(dispersal, ID.simulation, metric, time) %>%
    summarize(FD = weighted.mean(z, w = weight))
  
  
  # Extract simulation output on forest composition
  data.sp = sim_output %>%
    gather(key = "metric", value = "weight", "N", "BA") %>%
    # Only focus on the second half of the simulations
    filter(time >= 70) %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Calculate composition metrics
    group_by(dispersal, ID.simulation, metric, time) %>%
    mutate(p = weight/sum(weight), 
           plnp = p*log(p)) %>%
    summarise(H = -sum(plnp), 
              traits.pca1.mean = weighted.mean(GrSurv, w = weight), 
              traits.pca2.mean = weighted.mean(ShadeDrought, w = weight)) %>%
    # Add functional diversity
    left_join((data.FD %>% mutate(ID.simulation = as.integer(ID.simulation),
                                  time = as.integer(time))),
              by = c("dispersal", "ID.simulation", "time", "metric")) %>%
    # Average metrics
    ungroup() %>% group_by(dispersal, ID.simulation, metric) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              FD = mean(FD, na.rm = TRUE),
              cwm1 = mean(traits.pca1.mean, na.rm = TRUE), 
              cwm2 = mean(traits.pca2.mean, na.rm = TRUE)) %>%
    drop_na() 
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("H", "FD", "cwm1", "cwm2"), 
    title = c("Species_diversity", "Functional_diversity", 
              traits_compiled$title_axes$axis), 
    label = c("Species diversity\n(Shannon index)", "Functional diversity",
              traits_compiled$title_axes$title)) %>%
    mutate(label = gsub("\\(.+\\)", "", label))
  
  # Join all data together
  data = data.sp %>%
    left_join(simul_list, by = "ID.simulation") %>%
    dplyr::select(dispersal, ID.simulation, plotcode, metric, climate, ssp, 
                  dist, pca1, H, cwm1, cwm2, FD) %>%
    gather(key = "variable", value = "value", data.var$var) %>%
    drop_na()
  
  # Calculate the change relative to the reference scenario (ssp126 and no dist)
  data = data %>%
    left_join((data %>% filter(dist == "nodist" & ssp == "ssp126") %>% ungroup() %>%
                 dplyr::select(dispersal, plotcode, variable, metric, value.ref = value)), 
              by = c("dispersal", "plotcode", "variable", "metric")) %>%
    mutate(var.change = value-value.ref)  %>%
    ungroup() %>% 
    filter(metric == metric.ref) %>%
    dplyr::select(plotcode, climate, ssp, dist, 
                  var.change, dispersal, variable) %>%
    filter(!is.infinite(var.change)) %>%
    filter(!(ssp == "ssp126" & dist == "nodist")) %>%
    mutate(scenario = paste0(ssp, "_", dist)) %>%
    dplyr::select(- ssp, -dist) %>%
    spread(key = "scenario", value = "var.change") %>%
    rename(CC_nodist = ssp585_nodist, CC_dist = ssp585_dist)
  
  # Calculate local richness and number of new species from the pool
  # -- Get species presence in regional pool per plotcode
  sp.pool = regional_pool[, c("plotcode", "species")] %>%
    mutate(pool.pres = 1) %>%
    spread(key = "species", value = "pool.pres") %>%
    gather(key = "species", value = "pool.pres", unique(regional_pool$species)) %>%
    mutate(pool.pres = ifelse(is.na(pool.pres), 0, pool.pres))
  # -- Get species presence in NFI plots
  sp.local = NFI_data_sub[, c("plotcode", "species")] %>%
    mutate(nfi.pres = 1) %>%
    distinct() %>%
    spread(key = "species", value = "nfi.pres") %>%
    gather(key = "species", value = "nfi.pres", unique(NFI_data_sub$species)) %>%
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    mutate(species = ifelse(species == "Betula_sp", "Betula", species)) %>%
    mutate(nfi.pres = ifelse(is.na(nfi.pres), 0, nfi.pres))
  # -- Plot level information on richness and new species
  sp.plot_data = sp.pool %>% 
    left_join(sp.local, by = c("plotcode", "species")) %>%
    mutate(new.from.pool = ifelse((pool.pres == 1 & nfi.pres == 0), 1, 0)) %>%
    group_by(plotcode) %>%
    summarize(R.init = sum(nfi.pres), R.new_from_pool = sum(new.from.pool)) 
  
  # -- Format data to plot the difference between with and without pool
  data.diffpool = data %>%
    dplyr::select(-CC_nodist) %>%
    spread(key = "dispersal", value = "CC_dist") %>%
    mutate(diff_pool = with_pool - without_pool) %>%
    left_join(sp.plot_data, by = "plotcode") %>%
    mutate(R.init2 = log(R.init)) %>%
    left_join((dist_occurence %>% mutate(dist = ifelse(
      storm.bin == 1 | fire.bin == 1, "disturbed", "undisturbed")) %>%
        dplyr::select(-fire.bin, -storm.bin)), by = "plotcode") %>%
    mutate(dist = factor(dist, levels = c("disturbed", "undisturbed")))
  # Make stat analysis
  # -- Initialize stat vector
  var.vec = unique(data.diffpool$variable)
  # -- List of models
  mod.list = vector(mode = "list", length = length(var.vec))
  names(mod.list) = var.vec
  # -- List of newdata
  newdata.list = vector(mode = "list", length = length(var.vec))
  names(newdata.list) = var.vec
  # -- Loop on all variables to test
  for(i in 1:length(var.vec)){
    # Data to fit model i
    data.i = subset(data.diffpool, variable == var.vec[i])
    # Fit model i
    mod.list[[i]] = lm(diff_pool ~ R.init*dist + R.init2*dist, data = data.i)
    # Initialize newdata for model i
    newdata.list[[i]] = expand.grid(dist = unique(data.i$dist), 
                                    R.init = seq(from = 1, to = 5, length.out = 100)) %>% 
      mutate(diff_pool = 0, R.init2 = log(R.init)) %>% arrange(dist, R.init)
    # Get predictions
    form.i = formula("diff_pool ~ R.init*dist + R.init2*dist")
    # -- Generate matrix
    X <- model.matrix(form.i, newdata.list[[i]])
    # -- Extract fixed effects
    beta = coefficients(mod.list[[i]])
    # -- Extract variance covariance matrix
    vc = vcov(mod.list[[i]])
    # -- Prediction
    pred <- X %*% beta
    # -- Standard error of prediction
    pred.se <- sqrt(diag(X %*% vc %*% t(X))) 
    # -- Criteria to calculate confidence interval
    crit <- -qnorm(0.05/2)
    # -- Calculate confidence interval
    lwr <- pred-crit*pred.se 
    upr <- pred+crit*pred.se
    # -- Finish formatting
    newdata.list[[i]] = newdata.list[[i]] %>%
      dplyr::select(-diff_pool) %>%
      # Add predictions not in logit format
      mutate(fit = pred, lwr = lwr, upr = upr)
  }
  # Format predictions
  data.diffpool.fit = bind_rows(newdata.list, .id = "var") %>%
    left_join(data.var, by = "var")  %>%
    mutate(label = factor(label, levels = data.var$label))
  
  # Plot raw data along with predictions
  plot.diffpool = data.diffpool %>%
    rename(var = variable) %>%
    group_by(var, dist, R.init) %>%
    summarize(mean = mean(diff_pool, na.rm = TRUE), 
              se = sd(diff_pool, na.rm = TRUE)/sqrt(n()), 
              n = n()) %>%
    filter(n >= 10) %>%
    group_by(var, R.init) %>%
    mutate(mean.max = max(mean), lwr.data = min(mean - se), upr.data = max(mean+se)) %>%
    # Add lower and upper value of the predictions
    left_join((data.diffpool.fit %>% mutate(R.init = round(R.init*10, 0)/10) %>% 
                 filter(R.init %in% c(1:5)) %>% group_by(var, R.init) %>% 
                 summarize(lwr.fit = min(lwr), upr.fit = max(upr))), 
              by = c("var", "R.init")) %>%
    # Calculate the range per variable plotted
    ungroup() %>% group_by(var) %>%
    mutate(min = min(min(lwr.data), min(lwr.fit)), 
           max = max(max(upr.data), max(upr.fit)), 
           range = max - min) %>%
    ungroup() %>% dplyr::select(-max, -min) %>%
    # Add significance
    left_join((data.diffpool.fit %>% mutate(R.init = round(R.init*10, 0)/10) %>% 
                 filter(R.init %in% c(1:5)) %>% group_by(var, R.init, dist) %>% 
                 summarize(lwr = mean(lwr), upr = mean(upr)) %>%
                 mutate(signif = ifelse(lwr > 0 | upr < 0, "*", "")) %>% 
                 dplyr::select(-lwr, -upr)), by = c("var", "R.init", "dist")) %>%
    # Attribute label position depending on which is upper and lower
    group_by(var, R.init) %>%
    mutate(label.pos = ifelse(mean == mean.max, max(upr.fit, upr.data) + 0.15*range, 
                              min(lwr.fit, lwr.data) - 0.15*range), 
           label.n = ifelse(mean == mean.max, paste0(" \n", signif, "\n (", n, ") "), 
                            paste0(" (", n, ") \n", signif, "\n "))) %>%
    left_join(data.var, by = "var")  %>%
    mutate(label = factor(label, levels = data.var$label)) %>%
    ggplot(aes(x = R.init, fill = dist, color = dist)) + 
    geom_ribbon(data = data.diffpool.fit, inherit.aes = TRUE, 
                aes(ymin = lwr, ymax = upr), alpha = 0.3, color = NA) + 
    geom_line(data = data.diffpool.fit, inherit.aes = TRUE, aes(y = fit)) +
    geom_errorbar(aes(y = mean, ymin = mean - se, ymax = mean + se), 
                  position = position_dodge(0.1), color = "black", width = 0.1) + 
    geom_point(aes(y = mean), position = position_dodge(0.1), color = "black", shape = 21) + 
    geom_text(aes(y = label.pos, label = label.n), size = 2.5, lineheight = 0.8) +
    facet_wrap( ~ label, scale = "free", nrow = 1) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    scale_color_manual(values = c(`disturbed` = "#2A6F97", `undisturbed` = "#61A5C2")) + 
    scale_fill_manual(values = c(`disturbed` = "#2A6F97", `undisturbed` = "#61A5C2")) + 
    xlab("Initial species richness") + 
    ylab("Effect of the regional pool\non composition shift") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          legend.title = element_blank(), 
          strip.background = element_blank())
  
  
  # Save the plot
  ggsave(file.out, plot.diffpool, width = 30, 
         height = 7 , units = "cm", dpi = 600, bg = "white")
  
  # Return the file saved 
  return(file.out)
}






#' Plot the effect of climate and disturbances on temporal change in different variables
#' @param sim_output_short formatted outputs of the simulations
#' @param NFI_succession Succession stage per climate and NFI plot
#' @param simul_list Informations on the simulations made
#' @param traits_compiled List of traits information
#' @param metric.ref "N" or "BA" : which metric should be used as abundance
#' @param file.out Name of the file to save, including path 
plot_biogeo_effect_per.metric = function(
    sim_output_short, NFI_succession, simul_list, NFI_data_sub, traits_compiled, 
    metric.ref, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Calculate functional diversity
  data.FD = sim_output_short %>%
    gather(key = "metric", value = "weight", "N", "BA") %>%
    # Only focus on the second half of the simulations
    filter(time >= 70) %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Gather by functional axis
    gather(key = "axis", value = "trait_value", names(traits_compiled[[1]])) %>%
    # Calculate cwm for each time step and each simulation
    group_by(ID.simulation, time, metric, axis) %>%
    mutate(cwm = weighted.mean(trait_value, w = weight)) %>%
    ungroup() %>%
    # Calculate square distance of species to centroid along each axis
    mutate(zsq_per_axis = (trait_value - cwm)^2) %>%
    # Calculate the distance of each species to the centroid
    group_by(ID.simulation, time, species, metric, weight) %>%
    summarize(z = sqrt(sum(zsq_per_axis))) %>% ungroup() %>%
    # Calculate functional dispersion
    group_by(ID.simulation, metric, time) %>%
    summarize(FD = weighted.mean(z, w = weight))
  
  
  # Extract simulation output on forest composition
  data.sp = sim_output_short %>%
    gather(key = "metric", value = "weight", "N", "BA") %>%
    # Only focus on the second half of the simulations
    filter(time >= 70) %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Calculate composition metrics
    group_by(ID.simulation, metric, time) %>%
    mutate(p = weight/sum(weight), 
           plnp = p*log(p)) %>%
    summarise(H = -sum(plnp), 
              traits.pca1.mean = weighted.mean(GrSurv, w = weight), 
              traits.pca2.mean = weighted.mean(ShadeDrought, w = weight)) %>%
    # Add functional diversity
    left_join((data.FD %>% mutate(ID.simulation = as.integer(ID.simulation),
                                  time = as.integer(time))),
              by = c("ID.simulation", "time", "metric")) %>%
    # Average metrics
    ungroup() %>% group_by(ID.simulation, metric) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              FD = mean(FD, na.rm = TRUE),
              cwm1 = mean(traits.pca1.mean, na.rm = TRUE), 
              cwm2 = mean(traits.pca2.mean, na.rm = TRUE)) %>%
    drop_na() 
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("H", "FD", "cwm1", "cwm2"), 
    title = c("Species_diversity", "Functional_diversity", 
              traits_compiled$title_axes$axis), 
    label = c("Species diversity\n(Shannon index)", "Functional diversity",
              traits_compiled$title_axes$title)
  )
  
  # Join all data together
  data = data.sp %>%
    left_join(simul_list, by = "ID.simulation") %>%
    left_join(NFI_succession[, c("plotcode", "dqm_class")], by = "plotcode") %>%
    dplyr::select(ID.simulation, plotcode, metric, climate, dqm_class, ssp, 
                  dist, pca1, H, cwm1, cwm2, FD) %>%
    gather(key = "variable", value = "value", data.var$var) %>%
    drop_na()
  
  
  # Calculate the change relative to the reference scenario (ssp126 and no dist)
  data = data %>%
    left_join((data %>% filter(dist == "nodist" & ssp == "ssp126") %>%
                 dplyr::select(plotcode, variable, metric, value.ref = value)), 
              by = c("plotcode", "variable", "metric")) %>%
    # mutate(var.change = 100*(value-value.ref)/value.ref, 
    mutate(var.change = value-value.ref, 
           scenario = case_when(
             dist == "nodist" & ssp == "ssp126" ~ "reference", 
             dist == "nodist" & ssp == "ssp585" ~ "Climate change only", 
             dist == "dist" & ssp == "ssp585" ~ "Disturbance and climate change"
           )) %>%
    as.data.frame()
  
  # Prepare data for mapping
  data.map = data %>%
    filter(metric == metric.ref & scenario == "Disturbance and climate change") %>%
    ungroup() %>% dplyr::select(plotcode, variable, var.change) %>%
    left_join((NFI_data_sub[, c("plotcode", "longitude", "latitude")] %>% distinct), 
              by = "plotcode") %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Loop on all succession stage
  for(i in 1:length(unique(data$dqm_class))){
    
    # Succession stage i
    succ.i = unique(data$dqm_class)[i]
    
    # Loop on all variables for which to make analyses
    for(j in 1:dim(data.var)[1]){
      
      # Loop on all abundance metrics
      for(m in 1:2){
        
        # Metric m
        metric.m = unique(data$metric)[m]
        
        # Restrict the dataset to the variable j
        data.ijm = data %>%
          ungroup() %>% 
          filter(variable == data.var$var[j] & dqm_class == succ.i & 
                   metric == metric.m) %>%
          mutate(pca1sq = pca1^2, 
                 ssp = factor(ssp, levels = c("ssp126", "ssp585"))) %>%
          dplyr::select(plotcode, dqm_class, climate, ssp, dist, pca1, pca1sq, 
                        var.change, scenario) %>%
          filter(!is.infinite(var.change)) %>%
          filter(!(scenario %in% c("reference", "Disturbance only"))) %>%
          filter(var.change < quantile(.$var.change, 0.99, na.rm = TRUE) & 
                   var.change > quantile(.$var.change, 0.01, na.rm = TRUE)) 
        
        # Fit a first model
        mod = lmer(var.change ~ pca1*scenario + pca1sq*scenario + (1|plotcode), 
                   data = data.ijm)
        # Initialize boolean to stop model selection and counter of interactions removed
        selection.over = FALSE; k=0
        # Start while loop
        while(!selection.over){
          print(k)
          # Simple and double interaction terms
          terms.all.in = rownames(Anova(mod))[grep(":", rownames(Anova(mod)))]
          terms.double.in = rownames(Anova(mod))[grep(":.+:", rownames(Anova(mod)))]
          terms.simple.in = setdiff(terms.all.in, terms.double.in)
          terms.noint.in = setdiff(rownames(Anova(mod)), terms.all.in)
          # Initialize formula of the model
          form = paste(gsub("\\:", "\\*", c(terms.noint.in, terms.all.in)), collapse = " + ")
          # First configuration : there are double interactions left
          if(length(terms.double.in) > 0){
            # First sub-configuration : all double interactinos are significant
            if(!any(Anova(mod)[terms.double.in, 3] > 0.05)){
              # We stop model selection, no interactions can be removed
              selection.over = TRUE
            }else{
              # Otherwise, we remove the least significant double interaction 
              # -- increase counter
              k = k+1
              # -- remove from double iteractions the least significant
              terms.double.in = terms.double.in[
                -which(Anova(mod)[terms.double.in, 3] == max(Anova(mod)[terms.double.in, 3]))]
              # -- update vector containing all terms
              terms.all.in = c(terms.simple.in, terms.double.in)
              # -- New formula
              form = paste(gsub("\\:", "\\*", c(terms.noint.in, terms.all.in)), collapse = " + ")
              # -- Fit model with the new formula 
              eval(parse(text = paste0(
                "mod = lmer(var.change ~ ", form, " + (1|plotcode), data = data.ijm)")))
            }
          }else{
            # Second configuration : no more double interactions
            # First sub-configuration : all simple interactinos are significant
            if(!any(Anova(mod)[terms.simple.in, 3] > 0.05)){
              # We stop model selection, no interactions can be removed
              selection.over = TRUE
            }else{
              # Otherwise, we remove the least significant double interaction 
              # -- increase counter
              k = k+1
              # -- remove from double iteractions the least significant
              terms.simple.in = terms.simple.in[
                -which(Anova(mod)[terms.simple.in, 3] == max(Anova(mod)[terms.simple.in, 3]))]
              # -- update vector containing all terms
              terms.all.in = c(terms.simple.in, terms.double.in)
              # -- Only make new model if there are terms left in the formula
              if(length(terms.all.in) > 0){
                # -- New formula
                form = paste(gsub("\\:", "\\*", c(terms.noint.in, terms.all.in)), collapse = " + ")
                # -- Fit model with the new formula 
                eval(parse(text = paste0(
                  "mod = lmer(var.change ~ ", form, " + (1|plotcode), data = data.ijm)")))
              } else { 
                # Otherwise, stop model selection
                selection.over = TRUE
              }
            }
          }
        }
        
        
        # Make predictions based on fixed effects
        # -- Extract fixed effects
        beta = fixef(mod)
        # -- Extract variance vocariance matrix
        v = vcov(mod)
        # -- Initialize data for predictions
        newdata <- expand.grid(
          pca1 = seq(from = quantile(data.ijm$pca1, 0.01), 
                     to = quantile(data.ijm$pca1, 0.99), length.out = 100), 
          scenario = unique(data.ijm$scenario)[order(unique(data.ijm$scenario))]) %>%
          mutate(pca1sq = pca1^2, var.change = 0)
        # -- Same formula without random plot
        form <- formula(paste0("var.change ~ ", form))
        # -- Generate matrix
        X <- model.matrix(form, newdata)
        # -- Prediction
        pred <- X %*% beta
        # -- Standard error of prediction
        pred.se <- sqrt(diag(X %*% v %*% t(X))) 
        # -- Criteria to calculate confidence interval
        crit <- -qnorm(0.05/2)
        # -- Calculate confidence interval
        lwr <- pred-crit*pred.se 
        upr <- pred+crit*pred.se
        # -- Add to the prediction dataset
        newdata = newdata %>% mutate(fit = pred, lwr = lwr, upr = upr)
        # -- Add succession stage, metric  and variable
        newdata$dqm_class = succ.i
        newdata$metric = metric.m
        newdata$var = data.var$var[j]
        # Prepare points data for plotting
        data.points.ijm = data.ijm %>%
          group_by(climate, dqm_class, scenario) %>%
          summarize(pca1 = mean(pca1), 
                    fit = mean(var.change, na.rm = TRUE), 
                    lwr = quantile(var.change, 0.05, na.rm = TRUE), 
                    upr = quantile(var.change, 0.95, na.rm = TRUE), 
                    sd = sd(var.change, na.rm = TRUE)/sqrt(n())) %>%
          mutate(lwr = fit - sd, upr = fit + sd) %>%
          mutate(x.pos = case_when(
            scenario == "Disturbance only" ~ pca1 - 0.02*diff(range(.$pca1)), 
            scenario == "Disturbance and climate change" ~ pca1 + 0.02*diff(range(.$pca1)), 
            TRUE ~ pca1), 
            var = data.var$var[j]) %>%
          left_join(data.var, by = "var") %>%
          mutate(s = factor(dqm_class, levels = paste0(
            c("early", "late"), "-succession")), metric = metric.m) 
        
        # Add to final dataset
        if(j == 1 & i == 1 & m == 1){
          data.fit.i = newdata
          data.points.i = data.points.ijm
        }else{
          data.fit.i = rbind(data.fit.i, newdata)
          data.points.i = rbind(data.points.i, data.points.ijm)
        }
        
        
      }
      
    }
    
    
    
    
  }
  
  # Initialize the output plotlist
  plotlist.out = vector(mode = "list", length = dim(data.var)[1])
  names(plotlist.out) = data.var$var
  
  # Loop on all variables to make the plot
  for(k in 1:length(plotlist.out)){
    
    # Get the ylab
    if(data.var$var[k] == "FD") ylab.k = paste0(
      "Change in \n", gsub("Ax.+\\)", "Axis", data.var$label[k]))
    else ylab.k = paste0("Change in ", gsub("Ax.+\\)", "Axis", data.var$label[k]))
    
    # Subset data for variable k
    data.k = subset(data.map, variable == data.var$var[k])
    
    # Vector of points for color
    quant.0 = -min(data.k$var.change, na.rm = TRUE)/diff(range(data.k$var.change, na.rm = TRUE))
    point.25 = as.numeric(quantile(filter(data.k, var.change < 0)$var.change, 0.5, na.rm = TRUE))
    point.75 = as.numeric(quantile(filter(data.k, var.change > 0)$var.change, 0.5, na.rm = TRUE))
    quant.25 = point.25-min(data.k$var.change, na.rm = TRUE)/diff(range(data.k$var.change, na.rm = TRUE))
    quant.75 = point.75-min(data.k$var.change, na.rm = TRUE)/diff(range(data.k$var.change, na.rm = TRUE))
    vec.quant.k = c(0, quant.25, quant.0, quant.75, 1)
    
    # Map plot
    plot.map.k = ne_countries(scale = "medium", returnclass = "sf") %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(fill = "#343A40", color = "gray", show.legend = F, size = 0.2) + 
      geom_sf(data = subset(data.map, variable == data.var$var[k]), 
              aes(color = var.change), size = 1, shape = 20) +
      scale_color_gradientn(
        colors = c('#1D3461', '#1368AA', 'white', '#F29479', '#CB1B16'),
        values = vec.quant.k) +
      coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) + 
      # guides(color = guide_legend(override.aes = list(size = 0.5))) +
      theme(panel.background = element_rect(color = 'black', fill = 'white'), 
            panel.grid = element_blank(), 
            legend.title = element_blank(), 
            legend.key = element_blank(), 
            legend.position = c(0.25, 0.9), 
            legend.key.width = unit(0.4, "cm"), 
            legend.key.height = unit(0.15, "cm"), 
            legend.direction = "horizontal", 
            legend.text = element_text(size = 6)) 
    
    
    # Baseline plot
    plot.base.k = data.fit.i %>%
      filter(var == data.var$var[k]) %>%
      filter(dqm_class != "succession2") %>%
      filter(metric == metric.ref) %>%
      mutate(s = factor(dqm_class, levels = paste0(
        c("early", "late"), "-succession"))) %>%
      ggplot(aes(x = pca1, y = fit, group = scenario, 
                 color = scenario, fill = scenario, ymin = lwr, ymax = upr)) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_errorbar(data = (subset(data.points.i, var == data.var$var[k]) %>%
                              filter(metric == metric.ref)), color = "black", 
                    aes(x = x.pos), inherit.aes = TRUE, alpha = 0.5, width = 0.1) +
      geom_point(data = (subset(data.points.i, var == data.var$var[k]) %>%
                           filter(metric == metric.ref)), color = "black",
                 aes(x = x.pos), inherit.aes = TRUE, shape = 21) + 
      geom_line() + 
      geom_ribbon(alpha = 0.3, color = NA) + 
      xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
      ylab(ylab.k) +
      scale_color_manual(values = c('Disturbance only' = "#001219", 
                                    'Climate change only' = "#CA6702", 
                                    'Disturbance and climate change' = "#9B2226")) +
      scale_fill_manual(values = c('Disturbance only' = "#005F73", 
                                   'Climate change only' = "#EE9B00", 
                                   'Disturbance and climate change' = "#AE2012")) +
      facet_wrap( ~ s, ncol = 1) + 
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.key = element_blank(), 
            strip.background = element_blank(), 
            strip.text = element_text(size = 12))
    
    # Add different legend depending on the position
    if(k == 1){
      # Extract legend first
      plot.legend = get_legend(plot.base.k + theme(
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 16)))
    }
    
    # Final plot 
    plotlist.out[[k]] = plot_grid(plot.base.k + theme(legend.position = "none"), 
                                  plot.map.k, ncol = 1, align = "v")
  }
  # Plot predictions of the model
  plot.out = plot_grid(
    plot.legend, 
    plot_grid(plotlist = plotlist.out, nrow = 1, align = "hv", scale = 0.95), 
    ncol = 1, rel_heights = c(0.03, 1)
  ) 
  
  # -- Save the plots
  ggsave(file.out, plot.out, width = 35, height = 20 , units = "cm", dpi = 600, bg = "white")
  
  # Return the files saved
  return(file.out)
  
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot for supplementary material ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to plot the theoretical distribution per climate and dqm class
#' @param NFI_succession Succession stage per NFI plot and climate
#' @param file.out Name of the file to save, including path
plot_succession_distrib = function(NFI_succession, file.out){
  
  # Create output directory if needed 
  create_dir_if_needed(file.out)
  
  # Number of iterations for plotting 
  n.iter = 2000
  
  # Average parameters per climate and succession stage
  param.out = NFI_succession %>%
    group_by(climate, dqm_class) %>%
    summarize(shape = mean(shape, na.rm = TRUE), 
              scale = mean(scale, na.rm = TRUE))
  
  
  
  # Plot distribution per climate and per succession stage
  # - Build data set
  data.plot = expand.grid(climate = unique(param.out$climate), 
                          dqm_class = unique(param.out$dqm_class), 
                          iteration = c(1:n.iter), 
                          dbh.simulated = NA_real_) 
  # - Loop on all climate succession combination
  for(j in 1:dim(param.out)[1]){
    # Identify the IDs in the plot data
    id.i = which(data.plot$climate == param.out$climate[j] & 
                   data.plot$dqm_class == param.out$dqm_class[j])
    # Fill the plot data with simulated dbh
    data.plot$dbh.simulated[id.i] = rweibull(n.iter, scale = param.out$scale[j], 
                                             shape = param.out$shape[j])
  }
  
  # - Plot the distributions
  plot.out = data.plot %>%
    mutate(climate = factor(climate, levels = paste0(
      "clim", c(1:length(unique(data.plot$climate)))))) %>%
    ggplot(aes(x = dbh.simulated)) + 
    geom_histogram(aes(fill = climate), color = "black", show.legend = FALSE) +
    scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
    facet_grid(dqm_class ~ climate) + 
    theme_bw()
  
  # Export the plot
  ggsave(file.out, plot.out, width = 25, height = 10, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return the file saved
  return(file.out)
  
}






