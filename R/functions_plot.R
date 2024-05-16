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
#' @param simul_list dataframe listing the properties of each simulation
#' @param climate_dist_dflist climate and disturbance regimes per plot simulated
#' @param file.out Name of the file to save, including path
plot_map_clim_dist = function(NFI_plots_selected, simul_list, 
                              climate_dist_dflist, file.out){
  
  #• Create output directory if needed
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
  
  # Initialize the output dataframe
  out = simul_list %>%
    mutate(freq.fire_half1 = 0, freq.storm_half1 = 0, 
           freq.fire_half2 = 0, freq.storm_half2 = 0, 
           wai_change = 0, sgdd_change = 0) %>%
    filter(dist == "dist")
  
  # Loop on all plotcodes
  for(i in 1:dim(out)[1]){
    
    # Disturbance dataframe corresponding to simulation i
    dist.df.i = climate_dist_dflist[[out$plotcode[i]]][[out$ssp[i]]]$disturbance
    dist.df.i_1 = subset(dist.df.i, t < t.sim/2)
    dist.df.i_2 = subset(dist.df.i, t > t.sim/2)
    
    # Count the occurrences of storm the first half of simulations
    if("storm" %in% dist.df.i_1$type){
      out$freq.storm_half1[i] = length(which(dist.df.i_1$type == "storm"))/(0.5*t.sim)
    } 
    # Count the occurrences of fire the first half of simulations
    if("fire" %in% dist.df.i_1$type){
      out$freq.fire_half1[i] = length(which(dist.df.i_1$type == "fire"))/(0.5*t.sim)
    } 
    # Count the occurrences of storm the second half of simulations
    if("storm" %in% dist.df.i_2$type){
      out$freq.storm_half2[i] = length(which(dist.df.i_2$type == "storm"))/(0.5*t.sim)
    } 
    # Count the occurrences of fire the second half of simulations
    if("fire" %in% dist.df.i_2$type){
      out$freq.fire_half2[i] = length(which(dist.df.i_2$type == "fire"))/(0.5*t.sim)
    } 
    
    # Climate dataframe corresponding to simulation i
    clim.i = climate_dist_dflist[[out$plotcode[i]]][[out$ssp[i]]]$climate
    out$wai_change[i] = mean(clim.i$wai[c(1:10) + t.sim - 10]) - mean(clim.i$wai[c(1:10)])
    out$sgdd_change[i] = mean(clim.i$sgdd[c(1:10) + t.sim - 10]) - mean(clim.i$sgdd[c(1:10)])
  }
  
  # Data to plot the change in climatic conditions
  data.climchange = out %>%
    gather(key = "var", value = "diff", "sgdd_change", "wai_change") %>%
    group_by(climate, ssp, var) %>%
    summarize(pca1.mean = mean(pca1, na.rm = TRUE), 
              diff_mean = mean(diff, na.rm = TRUE), 
              diff_se = sd(diff, na.rm = TRUE)) %>%
    mutate(lwr = diff_mean - diff_se, 
           upr = diff_mean + diff_se, 
           var = gsub("\\_.+", "", var)) 
  
  # Finalize the formatting fo data to plot disturbances
  data.dist = out %>%
    group_by(climate, ssp) %>%
    summarize(pca1.mean = mean(pca1, na.rm = TRUE), 
              fire_half1 = mean(freq.fire_half1), 
              storm_half1 = mean(freq.storm_half1), 
              fire_half2 = mean(freq.fire_half2), 
              storm_half2 = mean(freq.storm_half2)) %>%
    gather(key = "key", value = "frequency", "storm_half1", "fire_half1", 
           "storm_half2", "fire_half2") %>%
    separate(col = "key", into = c("disturbance", "half"), sep = "\\_") %>%
    left_join(data.frame(half = c("half1", "half2"), 
                         period = c("2020-2060", "2061-2100")), by = "half")
  
  # Data with limits of each climate category
  data.clim = NFI_plots_selected %>%
    group_by(climate) %>%
    summarize(Hot = min(pca1, na.rm = TRUE), 
              Cold = max(pca1, na.rm = TRUE), 
              pca1.mean = mean(pca1, na.rm = TRUE)) %>%
    mutate(Hot = ifelse(Hot < min(data.climchange$pca1.mean), 
                        min(data.climchange$pca1.mean), Hot), 
           Cold = ifelse(Cold > max(data.climchange$pca1.mean), 
                         max(data.climchange$pca1.mean), Cold)) %>%
    mutate(upr.dist = max(data.dist$frequency), 
           lwr.dist = min(data.dist$frequency))
  
  # Plot the change in disturbance frequency 
  plot.dist = data.dist %>%
    ggplot(aes(x = pca1.mean, y = frequency, color = disturbance, linetype = ssp)) + 
    geom_rect(data = (data.clim %>% merge(data.frame(
      frequency = NA_real_, ssp = NA, period = unique(data.dist$period)))), 
      aes(fill = climate, xmin = Hot, xmax = Cold, ymin = lwr.dist, 
          ymax = upr.dist), color = "black", inherit.aes = TRUE, 
      alpha = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = color.vec) + 
    geom_line() + 
    facet_wrap(~ period) + 
    scale_color_manual(values = c('fire' = "red", 'storm' = "blue")) + 
    scale_linetype_manual(values = c('ssp126' = "dashed", 'ssp585' = "solid")) +
    xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
    ylab("Disturbance frequency") +
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank())
  
  # Plot the change in climatic conditions
  plot.clim = data.climchange %>%
    ggplot(aes(x = pca1.mean, y = diff_mean, linetype = ssp)) + 
    geom_rect(data = (data.clim %>% merge((
      data.climchange %>% group_by(var) %>% 
        summarize(min = min(lwr), max = max(upr)) %>% 
        mutate(ssp = NA, diff_mean = NA_real_)))), 
      aes(fill = climate, xmin = Hot, xmax = Cold, ymin = min,  ymax = max), 
      color = "black", inherit.aes = TRUE, alpha = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = color.vec) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), color = NA, fill = "grey", alpha = 0.5) +
    geom_line() + 
    facet_wrap(~ var, nrow = 1, scales = "free") + 
    xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
    ylab("Difference of average\n2090-2100 vs 2020-2030")+
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(face = "bold"), 
          legend.key = element_blank())
  
  
  
  
  ##%%%%%%%%%%%%%%%%%%%%%%%%
  ## -- Step 3 : assemble plots
  
  # Make the final plot
  plot.out = plot_grid(
    plot.map, plot_grid(plot.clim, plot.dist, ncol = 1, scale = 0.9,
                        labels = c("(b)", "(c)"), align = "v"), 
    nrow = 1, labels = c("(a)", ""))
  
  
  # Save the plot
  ggsave(file.out, plot.out, width = 25, height = 13 , units = "cm", 
         dpi = 600, bg = "white")
  
  
  # Return file saved
  return(file.out)
  
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Plot for statistical analyses ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








#' Function to plot the change in each sp.compo index between CC and ref
#' @param regional_pool regional pool per plotcode
#' @param sim_output_short dataframe containing output of the simulations
#' @param simul_list dataframe containing climate and plotcode per simulation
#' @param NFI_succession succession stage per plotcode
#' @param traits_compiled List with traits data, pca and functional axes
#' @param metric.ref Reference metric to use to quantify abundance ("BA" or "N")
#' @param file.out Name of the file to save, including path
plot_local_effect_permetric = function(
  regional_pool, sim_output_short, simul_list, NFI_succession, traits_compiled, 
  NFI_plots_selected, dist_occurence, metric.ref, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  
  # Calculate community shift per petric and per plotcode
  data = sim_output_short %>%
    gather(key = "metric", value = "weight", "N", "BA") %>%
    # Only focus on the second half of the simulations
    filter(time > floor(max(sim_output_short$time)*0.666)) %>%
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
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Calculate functional dispersion
    group_by(ID.simulation, metric, time) %>%
    mutate(p = weight/sum(weight), 
           plnp = p*log(p)) %>%
    summarize(FD = weighted.mean(z, w = weight), 
              H = -sum(plnp), 
              CWM1 = weighted.mean(GrSurv, w = weight), 
              CWM2 = weighted.mean(ShadeDrought, w = weight)) %>%
    # Average across timesteps
    ungroup() %>% group_by(ID.simulation, metric) %>%
    summarize(FD = mean(FD, na.rm = TRUE), 
              H = mean(H, na.rm = TRUE), 
              CWM1 = mean(CWM1, na.rm = TRUE), 
              CWM2 = mean(CWM2, na.rm = TRUE)) %>%
    # Calculate change per scenario
    gather(key = "sp.metric", value = "value", "FD", "H", "CWM1", "CWM2") %>%
    left_join(simul_list, by = "ID.simulation") %>%
    mutate(scenario = case_when(
      dist == "nodist" & ssp == "ssp126" ~ "ref", 
      dist == "dist" & ssp == "ssp585" ~ "CC", 
      TRUE ~ "other")) %>%
    filter(scenario != "other") %>% ungroup %>%
    dplyr::select(plotcode, metric, sp.metric, scenario, value) %>%
    spread(key = "scenario", value = "value") %>%
    # Remove plots where diversity is null
    mutate(Hnull = ifelse(sp.metric == "H" & CC == 0, 1, 0)) %>%
    group_by(plotcode) %>% mutate(Hnullmax = max(Hnull)) %>%
    ungroup() %>% filter(Hnull == 0) %>% dplyr::select(-Hnull, -Hnullmax) %>%
    # Calculate the difference between climate change and reference
    mutate(d = CC - ref) %>% dplyr::select(-CC, -ref)
  
  
  
  
  # Calculate functional diversity per plot
  data.FD = sim_output_short %>%
    # Only focus on the second half of the simulations
    filter(time == max(sim_output_short$time)) %>%
    # Only keep species and ID simulation
    dplyr::select(ID.simulation, species) %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Gather by functional axis
    gather(key = "axis", value = "trait_value", names(traits_compiled[[1]])) %>%
    # Calculate mean trait value for each axis at each simulation
    group_by(ID.simulation, axis) %>%
    mutate(c = mean(trait_value)) %>%
    ungroup() %>%
    # Calculate square distance of species to centroid along each axis
    mutate(zsq_per_axis = (trait_value - c)^2) %>%
    # Calculate the distance of each species to the centroid
    group_by(ID.simulation, species) %>%
    summarize(z = sqrt(sum(zsq_per_axis))) %>% ungroup() %>%
    # Calculate functional dispersion
    group_by(ID.simulation) %>%
    summarize(FD = mean(z), R = n()) %>%
    # Add plotcode
    left_join(simul_list[, c("ID.simulation", "plotcode")], by = "ID.simulation") %>%
    dplyr::select(plotcode, FD, R) %>% distinct()
  
  
  # Finish formatting by adding the different factors to data
  data.model = data %>%
    # Add climatic data
    left_join(NFI_plots_selected[, c("plotcode", "pca1", "climate")], 
              by = "plotcode") %>%
    # Add functional diversity
    left_join(data.FD, by = "plotcode") %>%
    # Add succession stage
    left_join((NFI_succession %>%
                 mutate(succession = case_when(
                   dqm_class == "succession1" ~ "early-succession", 
                   dqm_class == "succession2" ~ "mid-succession", 
                   dqm_class == "succession3" ~ "late-succession")) %>%
                 dplyr::select(plotcode, succession)), by = "plotcode") %>%
    # Add occurence of disturbances
    left_join(dist_occurence, by = "plotcode") %>%
    mutate(dist.bin = ifelse(fire.bin + storm.bin == 0, "undisturbed", 
                             "disturbed")) %>%
    # Remove plots where no composition change is possible
    filter(R > 1) %>%
    # Scale continuous variables
    mutate(succession = factor(succession, levels = paste0(
      c("early", "mid", "late"), "-succession"))) %>%
    as.data.frame()
  
  
  # Response variables 
  resp.var = unique(data.model$sp.metric)
  
  # Abundance variables 
  ab.var = unique(data.model$metric)
  
  # Diversity variables
  div.var = c("FD", "R")
  
  # Initialize the output list
  list.out = vector(mode = "list", length = length(resp.var))
  names(list.out) = resp.var
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(var = resp.var) %>%
    left_join(data.frame(
      var = c("H", "FD", "CWM1", "CWM2"), 
      title = c("Species_diversity", "Functional_diversity", 
                traits_compiled$title_axes$axis), 
      label = c("Species diversity\n(Shannon index)", "Functional diversity",
                traits_compiled$title_axes$title)
    ), by = "var")
  
  
  # Loop on all response variables
  for(r in 1:length(list.out)){
    
    
    # Initialize the list of outputs for list r
    list.out[[r]] = vector(mode = "list", length = length(ab.var))
    names(list.out[[r]]) = ab.var
    
    # Loop on all abundance metrics
    for(a in 1:length(ab.var)){
      
      # Subset data model for response variable r and abundance metric a
      data.model.r.a = subset(data.model, sp.metric == resp.var[r] & metric == ab.var[a])
      
      # Scale the continuous variables
      data.model.r.a_scaled = data.model.r.a
      data.model.r.a_scaled[, div.var] = scale(
        data.model.r.a[, div.var], center = TRUE, scale = TRUE)
      
      
      # Initialize the list of stats for each diversity metric
      list.out[[r]][[a]] = vector(mode = "list", length = length(div.var))
      names(list.out[[r]][[a]]) = div.var
      
      # Loop on all diversity metrics
      for(d in 1:length(div.var)){
        
        
        # Vector of simple effects
        vec.var = c(div.var[d], "succession", "dist.bin")
        # Vector of possible interactions
        vec.int = paste(combn(vec.var, 2)[1, ], combn(vec.var, 2)[2, ], sep = "*")
        # Initialize vector of formulas
        vec.form = paste(vec.var, collapse = " + ")
        # Loop on all possible number of interactions
        for(i in 1:length(vec.int)){
          # All combinations of i interactions
          comb.i = combn(vec.int, i)
          # Loop on all combinations
          for(j in 1:dim(comb.i)[2]){
            vec.form = c(vec.form, paste(c(vec.form[1], comb.i[, j]), collapse = " + "))
          }
        }
        
        # Initialize the list of models
        list.out[[r]][[a]][[d]]$mod.list = vector(
          mode = "list", length = length(vec.form))
        # Loop on all model formulas
        for(f in 1:length(vec.form)){
          eval(parse(text = paste0(
            "list.out[[r]][[a]][[d]]$mod.list[", f, "] = lmer(d ~ ", vec.form[f],
            " + (1|climate), data = data.model.r.a_scaled)")))
        }
        # Identify the best model
        mod.best = which.min(sapply(X = list.out[[r]][[a]][[d]]$mod.list, FUN = AIC))
        # Extract the best model
        list.out[[r]][[a]][[d]]$mod.best = list.out[[r]][[a]][[d]]$mod.list[[mod.best]]
        
        
      }
      
      # Identify the best model comparing the two diversity metrics
      # div.best = div.var[which.min(c(AIC(list.out[[r]][[a]][[1]]$mod.best), 
      #                                AIC(list.out[[r]][[a]][[2]]$mod.best)))]
      div.best = "R"
      list.out[[r]][[a]]$mod.best = list.out[[r]][[a]][[div.best]]$mod.best
      
      # Make predictions based on fixed effects
      # -- Extract fixed effects
      beta = fixef(list.out[[r]][[a]]$mod.best)
      # -- Extract variance covariance matrix
      vc = vcov(list.out[[r]][[a]]$mod.best)
      # -- Initialize data for predictions
      newdata <- expand.grid(
        R_unscaled = c(2, 3, 4), 
        succession = unique(data.model.r.a_scaled$succession)[order(unique(
          data.model.r.a_scaled$succession))], 
        dist.bin = unique(data.model.r.a_scaled$dist.bin)[order(unique(
          data.model.r.a_scaled$dist.bin))], 
        d = 0) %>%
        mutate(R = predict(lm(Y ~ R_unscaled, data.frame(
          Y = data.model.r.a_scaled$R, R_unscaled = data.model.r.a$R)),newdata = .))
      # -- Same formula without random effect
      form <- formula(paste0("d ~ ", paste(gsub("\\:", "\\*", rownames(
        Anova(list.out[[r]][[a]]$mod.best))), collapse = " + ")))
      # -- Generate matrix
      X <- model.matrix(form, newdata)
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
      list.out[[r]][[a]]$data.fit = newdata %>%
        dplyr::select(-d) %>%
        # Add predictions not in logit format
        mutate(fit = pred, lwr = lwr, upr = upr) %>%
        # Add the response variable 
        mutate(metric = ab.var[a], var = resp.var[r])
      
      # Add to a final dataframe
      if(a == 1 & r == 1) data.fit = list.out[[r]][[a]]$data.fit
      else data.fit = data.fit %>% rbind(list.out[[r]][[a]]$data.fit)
      
    }
    
    
  }
  
  # Make the plot
  plot.out = data.fit %>%
    filter(metric == metric.ref) %>%
    mutate(R_unscaled = case_when(
      succession == "early-succession" ~ R_unscaled - 0.1, 
      succession == "late-succession" ~ R_unscaled + 0.1, 
      TRUE ~ R_unscaled)) %>%
    mutate(R_unscaled = case_when(
      dist.bin == "disturbed" ~ R_unscaled - 0.025, 
      dist.bin == "undisturbed" ~ R_unscaled + 0.025)) %>%
    left_join(data.var, by = "var") %>%
    ggplot(aes(x = R_unscaled, y = fit, ymin = lwr, ymax = upr, 
               color = dist.bin, fill = succession)) + 
    geom_errorbar(width = 0.04) + 
    geom_point(shape = 21) + 
    scale_color_manual(values = c("black", "grey")) + 
    scale_fill_manual(values = c("#A7C957", "#6A994E", "#386641")) +
    scale_x_continuous(breaks = c(2:4)) +
    xlab("Richness of the local species pool") + 
    ylab("Community shift") + 
    facet_wrap(~ label, nrow = 2, scales = "free") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          legend.key = element_blank(), 
          legend.title = element_blank(), 
          strip.background = element_blank())
  

  
  # Save the plot
  ggsave(file.out, plot.out, width = 19, height = 12 , units = "cm", 
         dpi = 600, bg = "white")
  
  # Return file
  return(file.out)
  
  
}





#' Plot the effect of climate and disturbances on temporal change in different variables
#' @param sim_output_short formatted outputs of the simulations
#' @param NFI_succession Succession stage per climate and NFI plot
#' @param simul_list Informations on the simulations made
#' @param traits_compiled List of traits information
#' @param metric.ref "N" or "BA" : which metric should be used as abundance
#' @param file.out Name of the file to save, including path 
plot_biogeo_effect_per.metric = function(sim_output_short, NFI_succession, simul_list, 
                                      traits_compiled, metric.ref, file.out){
  
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Calculate functional diversity
  data.FD = sim_output_short %>%
    gather(key = "metric", value = "weight", "N", "BA") %>%
    # Only focus on the second half of the simulations
    filter(time > floor(max(sim_output_short$time)*0.666)) %>%
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
    filter(time > floor(max(sim_output_short$time)*0.666)) %>%
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
             dist == "dist" & ssp == "ssp126" ~ "Disturbance only", 
             dist == "dist" & ssp == "ssp585" ~ "Disturbance and climate change"
           )) %>%
    as.data.frame()
  
  
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
          left_join(data.frame(dqm_class = paste0("succession", c(1:3)), 
                               s = paste0(c("early", "mid", "late"), 
                                          "-succession")), by = "dqm_class") %>%
          mutate(s = factor(s, levels = paste0(
            c("early", "mid", "late"), "-succession")), 
            metric = metric.m) 
        
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
    # Baseline plot
    plotlist.out[[k]] = data.fit.i %>%
      filter(var == data.var$var[k]) %>%
      filter(dqm_class != "succession2") %>%
      filter(metric == metric.ref) %>%
      left_join(data.frame(dqm_class = paste0("succession", c(1:3)), 
                           s = paste0(c("early", "mid", "late"), 
                                      "-succession")), by = "dqm_class") %>%
      ggplot(aes(x = pca1, y = fit, group = scenario, 
                 color = scenario, fill = scenario, ymin = lwr, ymax = upr)) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_errorbar(data = (subset(data.points.i, var == data.var$var[k]) %>%
                              filter(dqm_class != "succession2" & metric == metric.ref)), 
                    aes(x = x.pos), inherit.aes = TRUE, alpha = 0.5) +
      geom_point(data = (subset(data.points.i, var == data.var$var[k]) %>%
                           filter(dqm_class != "succession2" & metric == metric.ref)), 
                 aes(x = x.pos), inherit.aes = TRUE, shape = 21) + 
      geom_line() + 
      geom_ribbon(alpha = 0.3, color = NA) + 
      xlab("Position along the climate axis (Hot-dry to cold-wet)") + 
      ylab(paste0("Change in \n", 
                  gsub("Ax.+\\)", "Axis", data.var$label[k]))) + 
      scale_color_manual(values = c('Disturbance only' = "#001219", 
                                    'Climate change only' = "#CA6702", 
                                    'Disturbance and climate change' = "#9B2226")) +
      scale_fill_manual(values = c('Disturbance only' = "#005F73", 
                                   'Climate change only' = "#EE9B00", 
                                   'Disturbance and climate change' = "#AE2012")) +
      facet_wrap( ~ s, nrow = 1) + 
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.key = element_blank(), 
            strip.background = element_blank(), 
            strip.text = element_text(size = 12))
    
    # Add diferent legend depending on the position
    if(k < dim(data.var)[1]){
      # Remove x axis title and ticks and remove legend
      plotlist.out[[k]] = plotlist.out[[k]] + 
        theme(legend.position = "none")
    } else{
      # Extract legend first
      plot.legend = get_legend(plotlist.out[[k]] + theme(
        legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 16)))
      # And then remove legend but leave x axis
      plotlist.out[[k]] = plotlist.out[[k]] + 
        theme(legend.position = "none")
    }
    
  }
  # Plot predictions of the model
  plot.out = plot_grid(
    plot_grid(plotlist = plotlist.out, ncol = 2, align = "hv", scale = 0.9,
              labels = paste0("(", letters[c(1:dim(data.var)[1])], ")")), 
    plot.legend, ncol = 1, rel_heights = c(1, 0.1)
  ) 
  
  # -- Save the plots
  ggsave(file.out, plot.out, width = 35, height = 16 , units = "cm", dpi = 600, bg = "white")
  
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






