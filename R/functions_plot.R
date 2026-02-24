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




#' Make a map of phi per species composition metric
#' @param NFI_data_sub subset of the NFI data with only plots selected
#' @param phi_per_scenario Metric phi for each plot and each scenario
#' @param metric.ref Do we use basal area ("BA') or density ("N") to wuantify abundance
#' @param phi.ref Do we show the absolute change in composition ("fixed") or the
#'                rate of change ("rate")
#' @param file.out Name of the file to save, including path
map_phi = function(NFI_data_sub, phi_per_scenario, metric.ref, phi.ref, file.out){
  
  # Bug fix since sf update
  sf::sf_use_s2(FALSE)
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("H", "FD", "cwm1", "cwm2"), 
    title = c("Species_diversity", "Functional_diversity", 
              "GrSurv", "ShadeDrought"), 
    label = c("Species\ndiversity", "Functional\ndiversity",
              "CWM on axis\nGrowth <-> Survival", 
              "CWM on axis\nShade tol. <-> Drought tol.")
  )
  
  # Choose the right metric based on reference defined
  # - If phi is epressed as a rate per year
  if(phi.ref == "rate"){
    # The final phi is the change in metric per decade 
    phi_per_scenario = phi_per_scenario %>%
      mutate(phi.final = phi.rate.percent*10)
    # Ajust the label
    phi.label = "\u03c6 (% of range\nobserved per decade)"
    phi.label.long = "\u03c6 : climate change effect on species composition\n(% of range observed per decade)"
  }
  # - If phi is epressed as an absolute change
  if(phi.ref == "fixed"){
    # Just change column name
    phi_per_scenario = phi_per_scenario %>%
      rename(phi.final = phi.percent)
    # Ajust the label
    phi.label = "\u03c6 (% of range\nobserved)"
    phi.label.long = "\u03c6 : climate change effect on species composition\n(% of range observed)"
  }
  
  # Prepare data for mapping
  data.map = phi_per_scenario %>%
    filter(metric == metric.ref & dist.scenario == "dist" & pool == "pool") %>%
    left_join((NFI_data_sub[, c("plotcode", "longitude", "latitude")] %>% distinct), 
              by = "plotcode") %>%
    ungroup() %>% 
    dplyr::select(plotcode, longitude, latitude, variable, phi.final) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  # Initialize the grid for plotting
  world_6933 <- st_transform(world, 4326) %>%
    st_make_grid(n = c(300, 300), what = 'polygons', square = FALSE,
                 flat_topped = TRUE) %>%
    st_as_sf() %>%
    mutate(hex = floor(as.numeric(rownames(.))))
  
  
  # Average variable k across each hexagon
  data.map_perhex = data.map %>%
    st_join(world_6933, join = st_within) %>%
    st_drop_geometry() %>%
    group_by(hex, variable) %>%
    summarize(phi = mean(phi.final, na.rm = TRUE), 
              n = n()) %>%
    filter(n > 7)
  
  # Add value of variable k to the grid. 
  data_plot = merge(world_6933, data.frame(variable = data.var$var)) %>%
    filter(hex %in% data.map_perhex$hex) %>%
    left_join(data.map_perhex, by = c("hex", "variable")) %>%
    left_join(data.var %>% rename(variable = var), by = "variable")
  
  # Initialize the list of plots
  plotlist.out = vector(mode = "list", length = dim(data.var)[1])
  
  # Loop on all variable to plot
  for(i in 1:dim(data.var)[1]){
    
    # Make the map
    map.i = ggplot() +
      geom_sf(data = ne_countries(scale = "medium", returnclass = "sf"), 
              aes(geometry = geometry),
              fill = "#343A40", color = "gray", show.legend = F, size = 0.2) +  
      geom_sf(data = subset(data_plot, variable == data.var$var[i]), aes(fill = phi)) +
      scale_fill_gradient2(
        low = '#1368AA', mid = 'white', high = '#CB1B16', midpoint = 0,
        name = phi.label, 
        guide = "colourbar") +
      coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) + 
      ggtitle(paste0("\u03c6(", data.var$var[i], "): Climate change effect on\n", 
                     data.var$label[i])) +
      theme(panel.background = element_rect(color = 'black', fill = 'white'), 
            panel.grid = element_blank(), 
            legend.key = element_blank(), 
            legend.position = "bottom",
            legend.title = element_text(hjust = 1, size = 12),
            axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks = element_blank(), 
            plot.title = element_text(hjust = 0.5, size = 10),
            strip.background = element_blank())
    
    # Make histogram for variable i
    hist.i = data.map_perhex %>%
      filter(variable == data.var$var[i]) %>%
      mutate(class = cut(phi, breaks = seq(from = min(.$phi, na.rm = TRUE), 
                                           to = max(.$phi, na.rm = TRUE), 
                                           length.out = 10))) %>%
      group_by(class) %>%
      summarize(n = n(), 
                phi.bin = median(phi, na.rm = TRUE)) %>%
      ungroup() %>%
      drop_na() %>%
      ggplot(aes(x = phi.bin, y = n)) +
      geom_bar(color = "gray", stat = "identity", aes(fill = phi.bin)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_vline(xintercept = mean(
        subset(data.map_perhex, variable == data.var$var[i])$phi, na.rm = TRUE), 
        color = "purple", linetype = "dashed") +
      scale_fill_gradient2(
        low = '#1368AA', mid = 'white', high = '#CB1B16', midpoint = 0) + 
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            axis.text.y = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_blank(), 
            legend.position = "none") 
    
    # Assemble to get plot i
    plotlist.out[[i]] = plot_grid((map.i + theme(legend.position = "none")), hist.i, 
                                  align = "v", rel_heights = c(1, 0.3), ncol = 1)
  }
  
  # Make a separate plot for xlab
  plot.xlab = ggdraw() + draw_label(phi.label.long, hjust = 0.5) 
  
  # Assemble all plots
  plot.out = plot_grid(plot_grid(plotlist = plotlist.out, nrow = 1, align = "hv", scale = 0.9), 
                       plot.xlab, ncol = 1, rel_heights = c(1, 0.1))
  
  # Save the plot
  ggsave(file.out, plot.out, width = 26, height = 15 , units = "cm", 
         bg = "white", dpi = 600)
  
  # Return file saved
  return(file.out)
  
} 

#' Plot the effect of climate and disturbances on temporal change in different variables
#' @param phi_per_scenario Metric phi for each plot and each scenario
#' @param metric.ref Do we use basal area ("BA') or density ("N") to wuantify abundance
#' @param phi.ref Do we show the absolute change in composition ("fixed") or the
#'                rate of change ("rate")
#' @param dir.out Directory where to save plot
plot_biogeo_effect = function(phi_per_scenario, metric.ref, phi.ref, dir.out){
  
  # Create output directory if needed
  create_dir_if_needed(paste0(dir.out, "/test"))
  
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("H", "FD", "cwm1", "cwm2"), 
    title = c("Species_diversity", "Functional_diversity", 
              "GrSurv", "ShadeDrought"), 
    label = c("Species\ndiversity", "Functional\ndiversity",
              "CWM on axis\nGrowth <-> Survival", 
              "CWM on axis\nShade tol. <-> Drought tol.")
  )
  
  # Choose the right metric based on reference defined
  # - If phi is epressed as a rate per year
  if(phi.ref == "rate"){
    # The final phi is the change in metric per decade 
    phi_per_scenario = phi_per_scenario %>%
      mutate(phi.final = phi.rate.percent*10)
    # Ajust the label
    phi.label = "\u03c6: climate change effect on species composition\n(% of range observed per decade)"
  }
  # - If phi is epressed as an absolute change
  if(phi.ref == "fixed"){
    # Just change column name
    phi_per_scenario = phi_per_scenario %>%
      rename(phi.final = phi.percent)
    # Ajust the label
    phi.label = "\u03c6: climate change effect on species composition\n(% of range observed)"
  }
  
  # Add the quadratic term of climate and remove NA's
  phi_per_scenario = phi_per_scenario %>% 
    # Keep metric of interest
    filter(metric == metric.ref) %>%
    rename(var = variable) %>%
    # Average phi per modality
    group_by(climate, dqm_class, var, dist.scenario, pool) %>%
    mutate(phi.final = as.numeric(phi.final)) %>%
    summarize(pca1 = mean(pca1, na.rm = TRUE),
              phi.mean = mean(phi.final, na.rm = TRUE),
              phi.se = sd(phi.final, na.rm = TRUE)/sqrt(n())) %>%
    mutate(pca1sq = pca1^2, 
           w = 1/phi.se) %>%
    drop_na()
  
  # Initialize list of residual plots
  resid.plotlist = vector(mode = "list", length = 4)
  
  
  # Loop on all metrics
  for(i in 1:dim(data.var)[1]){
    
    # Subset dataset with variable i
    data.i = subset(phi_per_scenario, var == data.var$var[i]) 
    
    # Full model with all interactions
    model.i.full = lmerTest::lmer(phi.mean ~ pca1 + pca1sq + dqm_class + 
                                    dqm_class*pca1 + dqm_class*pca1sq + dist.scenario + 
                                    dist.scenario*pca1 + dist.scenario*pca1sq + 
                                    pool + pool*pca1 + pool*pca1sq + 
                                    dqm_class*dist.scenario + dqm_class*pool + 
                                    pool*dist.scenario + 
                                    (1|climate), data = data.i, weights = w)
    
    # Temporarily attach data to global environment
    assign("data.i", data.i, envir = .GlobalEnv)
    
    # Reduce model with backward selection
    model.i.reduced = lmerTest::get_model(lmerTest::step(model.i.full, 
                                                         reduce.random = FALSE))
    
    # And then remove it from global environment
    rm("data.i", envir = .GlobalEnv)
    
    # Percentage of variance explained based on chi square
    var.explained.i = data.frame(var = rownames(Anova(model.i.reduced)), 
                                 chisq = Anova(model.i.reduced)$Chisq) %>%
      mutate(prop.var = chisq/sum(chisq)*100, 
             var.expl = data.var$var[i]) 
    
    # Make residuals plot
    plot.resid.i = data.frame(residuals = residuals(model.i.reduced), 
                              fitted = fitted(model.i.reduced), 
                              obs = data.i$phi.mean, 
                              pool = data.i$pool, 
                              dist.scenario = data.i$dist.scenario, 
                              dqm_class = data.i$dqm_class) %>%
      mutate(factor = paste(pool, dist.scenario, sep = " - ")) %>%
      ggplot(aes(x = fitted, y = obs)) + 
      geom_point(alpha = 0.2) + 
      geom_smooth(method = "loess", color = "red", se = FALSE) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") + 
      ggtitle(data.var$label[i]) + 
      facet_grid(dqm_class ~ factor) +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    # Add to the final residual list
    resid.plotlist[[i]] = plot.resid.i
    
    # Make predictions based on fixed effects
    # -- Extract fixed effects
    beta = fixef(model.i.reduced)
    # -- Extract variance vocariance matrix
    v = vcov(model.i.reduced)
    # -- Initialize data for predictions
    newdata.i <- expand.grid(
      pca1 = seq(from = min(data.i$pca1), 
                 to = max(data.i$pca1), length.out = 100), 
      pool = unique(data.i$pool)[order(unique(data.i$pool))], 
      dqm_class = unique(data.i$dqm_class)[order(unique(data.i$dqm_class))], 
      dist.scenario = unique(data.i$dist.scenario)[order(unique(data.i$dist.scenario))]) %>%
      mutate(pca1sq = pca1^2, phi.mean = 0)
    # -- Same formula without random plot
    form <- formula(paste0("phi.mean ~ ", paste(
      rownames(Anova(model.i.reduced)), collapse = " + ")))
    # -- Generate matrix
    X <- model.matrix(form, newdata.i)
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
    newdata.i = newdata.i %>% 
      mutate(fit = as.numeric(pred), lwr = as.numeric(lwr), upr = as.numeric(upr))
    # -- Add variable
    newdata.i$var = data.var$var[i]
    
    # Add to final datasets
    if(i == 1){
      data.predict = newdata.i
      var.explained = var.explained.i
    } else {
      data.predict = rbind(data.predict, newdata.i)
      var.explained = rbind(var.explained, var.explained.i)
    } 
    
  }
  
  # assemble the final plot of residuals
  plot.resid.out = plot_grid(plotlist = resid.plotlist, ncol = 2, scale = 0.95)
  
  # Add title and label to data.predict
  data.predict = data.predict %>% 
    left_join(data.var, by = "var") %>%
    mutate(pool.title = ifelse(pool == "pool", "Regional pool\nincluded", 
                               "Regional pool\nnot included"), 
           dist.title = ifelse(dist.scenario == "dist", 
                               "Disturbances\nincluded in\nsimulations", 
                               "Disturbances\nnot included in\nsimulations")) %>%
    mutate(label = factor(label, levels = data.var$label))
  
  # Make the plot of data and prediction
  plot.predict = phi_per_scenario %>%
    # Calculate lower and upper confidence interval
    mutate(lwr = phi.mean - phi.se, upr = phi.mean + phi.se) %>%
    # Ajust x position of points  
    mutate(pca1 = ifelse(dist.scenario == "dist", pca1 + 0.05, pca1 - 0.05)) %>%
    # Change title of labels
    mutate(pool.title = ifelse(pool == "pool", "Regional pool\nincluded", 
                               "Regional pool\nnot included"), 
           dist.title = ifelse(dist.scenario == "dist", 
                               "Disturbances\nincluded in\nsimulations", 
                               "Disturbances\nnot included in\nsimulations")) %>%
    left_join(data.var, by = "var") %>%
    mutate(label = factor(label, levels = data.var$label)) %>%
    # Make plot
    ggplot(aes(x = pca1, ymin = lwr, ymax = upr, color = pool.title, fill = pool.title), 
           group = interaction(pool.title, dist.title)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(data = data.predict, aes(alpha = dist.title), 
                inherit.aes = TRUE, color = NA) + 
    geom_line(data = data.predict, aes(linetype = dist.title, y = fit), 
              inherit.aes = TRUE) +
    geom_errorbar(width = 0) + 
    geom_point(aes(y = phi.mean, shape = dist.title, size = dist.title)) + 
    scale_color_manual(values = c("#386641", "#003049")) +
    scale_fill_manual(values = c("#6A994E", "#669BBC")) +
    scale_shape_manual(values = c(23, 21)) + 
    scale_size_manual(values = c(2, 2.5)) +
    scale_alpha_manual(values = c(0.7, 0.4)) +
    ggh4x::facet_grid2(dqm_class ~ label, independent = "y", scales = "free_y") + 
    xlab("Position along the climate axis\n(Hot-dry to cold-wet)") +
    ylab(phi.label) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          legend.key = element_blank(), 
          legend.title = element_blank(),
          strip.background = element_blank(), 
          legend.text = element_text(size = 12), 
          strip.text = element_text(size = 12), 
          axis.title = element_text(size = 12))
  
  # Plot percentage of variance explained
  # - Group quadratic terms
  data.variance = var.explained %>%
    # Add title per response variable
    left_join(data.var %>% rename(var.expl = var), by = "var.expl") %>%
    # Modify the name of variables
    mutate(var = gsub("sq", "", var), 
           var = gsub("pca1", "Mean climate", var), 
           var = gsub("pool", "Regional pool", var), 
           var = gsub("dqm\\_class", "Succession", var), 
           var = gsub("dist\\.scenario", "Disturbances", var), 
           var = gsub("\\:", "\\ x\\ ", var)) %>%
    # Sum the variance explained by pca1 and pca1sq
    group_by(var, label) %>%
    summarize(prop.cumul = sum(prop.var, na.rm = TRUE)) %>%
    ungroup() 
  # - Arrange by cumulated variance
  levels.variance = (data.variance %>% group_by(var) %>% 
                       summarize(m = mean(prop.cumul)) %>% arrange(m))$var
  # - Complete the plot
  plot.variance = data.variance %>%
    mutate(label = factor(label, levels = data.var$label), 
           var = factor(var, levels = levels.variance)) %>%
    ggplot(aes(x = var, y = prop.cumul, fill = prop.cumul)) + 
    geom_bar(stat = "identity", color = "black") +
    scale_fill_gradient(low = "white", high = "#C1121F") + 
    facet_wrap(~ label, nrow = 1) + 
    coord_flip() + 
    ylab("Relative proportion of variance explained (%)") + 
    xlab("") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none", 
          strip.text = element_text(size = 12), 
          axis.title = element_text(size = 12))
  
  # Assemble plots
  plot.out = plot_grid(
    plot_grid((ggplot() + theme_void()), plot.predict, nrow = 1, rel_widths = c(0.095, 1)), 
    (ggplot() + theme_void()), 
    plot_grid(plot.variance, (ggplot() + theme_void()), nrow = 1, rel_widths = c(1, 0.17)), 
    ncol = 1, rel_heights = c(1, 0.1, 0.4), labels = c("(a)", "", "(b)")) 
  
  # - Name plots
  file.predict = paste0(dir.out, "/predict_phi.jpg")
  file.resid = paste0(dir.out, "/fig_residuals_biogeo.pdf")
  
  # -- Save the plots
  ggsave(file.predict, plot.out, width = 35, height = 20 , units = "cm", bg = "white", dpi = 600)
  ggsave(file.resid, plot.resid.out, width = 29, height = 20 , units = "cm", bg = "white")
  
  # Return files generated
  return(c(file.predict, file.resid))
  
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
    geom_histogram(aes(fill = climate), color = NA, show.legend = FALSE) +
    scale_fill_manual(values = colorRampPalette(c("#FFBF69", "#C5D86D", "#AEB8FE"))(10)) +
    facet_grid(dqm_class ~ climate) + 
    theme_bw() + 
    xlab("Diameter at breast height (mm)")
  
  # Export the plot
  ggsave(file.out, plot.out, width = 22, height = 12, units = "cm", 
         dpi = 600, bg = "white")
  
  # Return the file saved
  return(file.out)
  
}

#' Function to plot predictions of disturbance coefficients from traits
#' @param intensity_file posterior of disturbance coefficients from Barrere et al. 2023
#' @param traits_compiled LIst of functional traits value
#' @param file.out Name of the file to save, including path
plot_disturb_coef = function(intensity_file, traits_compiled, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  
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
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    filter(disturbance == "fire") %>%
    ungroup() %>%
    group_by(species)  %>%
    summarise(a0 = mean(a0), b = mean(b), c = mean(c)) %>%
    left_join(traits_compiled$traits_imputed$ShadeDrought %>%
                dplyr::select("species", "BT" = "bark.thickness"), 
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
    left_join(traits_compiled$traits_imputed$GrSurv %>%
                dplyr::select("species", "WD" = "wood.density"), 
              by = "species") %>%
    drop_na()
  # - Fit model
  mod.storm = lm(cbind(a0, a1, b, c) ~ WD, data = data.wd.fit.storm)
  # - Scale parameters
  scale.storm = param %>% ungroup() %>%
    filter(disturbance == "storm") %>%
    dplyr::select(logratio.intercept, logratio.slope, dbh.intercept, dbh.slope) %>%
    distinct()
  
  
  # Prepare coefficients data to label plots with equation
  # - for fire
  data.coef.fire = coefficients(mod.fire) %>% as.data.frame() %>% 
    mutate(param = c("int", "slope")) %>% 
    gather(key = "parameter", value = "value", "a0", "b", "c") %>% 
    mutate(dist = "fire", value = round(value, digits = 3)) %>% 
    spread(key = "param", value = "value") %>% 
    mutate(label = paste0(parameter, "\ny = ", int, " + ", slope, "*x")) %>% 
    dplyr::select(-int, -slope)
  # - For storm
  data.coef.storm = coefficients(mod.storm) %>% as.data.frame() %>% 
    mutate(param = c("int", "slope")) %>% 
    gather(key = "parameter", value = "value", "a0", "a1", "b", "c") %>% 
    mutate(dist = "storm", value = round(value, digits = 3)) %>% 
    spread(key = "param", value = "value") %>% 
    mutate(label = paste0(parameter, "\ny = ", int, " + ", slope, "*x")) %>% 
    dplyr::select(-int, -slope)
  
  # Prepare data for plotting (with data and preditions)
  data = rbind((data.wd.fit.storm[, c("species", "WD")] %>%
                  mutate(dist = "storm") %>% 
                  cbind(predict(mod.storm, newdata = .)) %>%
                  rename(trait = WD) %>%
                  gather(key = "parameter", value = "fit", "a0", "a1", "b", "c")), 
               (data.bt.fit[, c("species", "BT")] %>%
                  mutate(dist = "fire") %>% 
                  cbind(predict(mod.fire, newdata = .)) %>%
                  rename(trait = BT) %>%
                  gather(key = "parameter", value = "fit", "a0", "b", "c"))) %>%
    left_join((bind_rows(list(storm = data.wd.fit.storm[, c("species", "a0", "a1", "b", "c")], 
                              fire = data.bt.fit[, c("species", "a0", "b", "c")]), 
                         .id = "dist") %>%
                 gather(key = "parameter", value = "value", "a0", "a1", "b", "c") %>%
                 drop_na()), 
              by = c("species", "dist", "parameter")) %>%
    left_join((rbind(data.coef.fire, data.coef.storm)), 
              by = c("dist", "parameter"))
  
  # Plot the results
  # - Plot the data and fit for storm
  plot.storm = data %>%
    filter(dist == "storm") %>% 
    ggplot(aes(x = trait)) + 
    geom_point(aes(y = value), shape = 21, color = "black", fill = "blue") + 
    geom_line(aes(y = fit), color = "blue") + 
    facet_wrap(~ label, nrow = 1, scales = "free") + 
    xlab("Wood density") + ylab("Parameter of storm\nsensitivity") + 
    theme_bw()
  # - Plot the data and fit for fire
  plot.fire = data %>%
    filter(dist == "fire") %>% 
    ggplot(aes(x = trait)) + 
    geom_point(aes(y = value), shape = 21, color = "black", fill = "red") + 
    geom_line(aes(y = fit), color = "red") + 
    facet_wrap(~ label, nrow = 1, scales = "free") + 
    xlab("Bark thickness") + ylab("Parameter of fire\nsensitivity") + 
    theme_bw()
  # Assemble the two plots
  plot.out = plot_grid(plot.storm, plot.fire, ncol = 1, align = "v", 
                       labels = c("(a)", "(b)"), scale = 0.9)
  
  # Save the plot
  ggsave(file.out, plot.out, width = 22, height = 13 , units = "cm", 
         dpi = 600, bg = "white")
  
  
  # Return file saved
  return(file.out)
  
}



#' Plot the distribution of species richness per climate
#' @param NFI_plots_selected plot level information on the NFI plots selected
#' @param file.out Name of the file to save, including path
plot_richness_distrib = function(NFI_plots_selected, file.out){
  
  # Create output directory if necessary
  create_dir_if_needed(file.out)
  
  # Initialize the data with information on richness distribution per climate
  data.clim = NFI_plots_selected %>% 
    dplyr::select(climate, lambda) %>% 
    distinct() %>%
    mutate(climate = factor(climate, levels = paste0("clim", c(1:dim(.)[1]))))
  
  # Loop on all climate
  for(i in 1:dim(data.clim)[1]){
    data.i = data.frame(climate = data.clim$climate[i], 
                        richness = rpois(1000, data.clim$lambda[i]))
    if(i == 1) data = data.i
    else data = rbind(data, data.i)
  }
  
  # Vector of color for plotting
  color.vec = colorRampPalette(c("#FFBF69", "#C5D86D", "#AEB8FE"))(
    length(unique(NFI_plots_selected$climate)))
  names(color.vec) = paste0("clim", c(1:length(color.vec)))
  
  # Make the plot
  plot.out = data %>%
    mutate(climate = factor(climate, levels = paste0("clim", c(1:dim(.)[1])))) %>%
    ggplot(aes(x = richness, fill = climate)) + 
    geom_histogram(color = "black", binwidth = 1, show.legend = FALSE) + 
    facet_wrap(~ climate, nrow = 2) + 
    scale_fill_manual(values = color.vec) + 
    theme_bw()
  
  # Save the plot
  ggsave(file.out, plot.out, width = 22, height = 12 , units = "cm", 
         dpi = 600, bg = "white")
  
  
  # Return file saved
  return(file.out)
  
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Export of tables ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to a table of the functional and climate position of each species
#' @param NFI_data_sub Subset of NFI data for the plots selected
#' @param NFI_plots_selected NFI plots that were selected for the analyses
#' @param traits_compiled list compiling traits data, pca and functional position
#' @param file.out Name of the file to save, including path
export_table_funclim_species = function(NFI_data_sub, NFI_plots_selected, 
                                        traits_compiled, file.out){
  
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Compile species level data
  data = NFI_data_sub %>%
    dplyr::select(plotcode, species) %>%
    distinct() %>%
    left_join(NFI_plots_selected[, c("plotcode", "pca1")], by = "plotcode") %>%
    group_by(species) %>%
    summarize(Cold = quantile(pca1, 0.975, na.rm = TRUE), 
              Optimum = mean(pca1, na.rm = TRUE), 
              Hot = quantile(pca1, 0.025, na.rm = TRUE)) %>%
    # Add coordinates on functional axes
    left_join((traits_compiled$species_coord %>%
                 mutate(species = gsub("\\_", "\\ ", species)) %>%
                 mutate(species = ifelse(species == "Betula", "Betula sp", species))), 
              by = "species") %>%
    # Add percentage of basal area
    left_join((NFI_data_sub %>% group_by(species) %>% summarize(ba = sum(ba_ha)) %>%
                 ungroup() %>% mutate(percent = round(ba/sum(ba)*100, digits = 1)) %>%
                 dplyr::select(-ba)), by = "species") %>%
    # Round each variable
    mutate(Cold = round(Cold, digits = 2), Optimum = round(Optimum, digits = 2), 
           Hot = round(Hot, digits = 2), GrSurv = round(GrSurv, digits = 2), 
           ShadeDrought = round(ShadeDrought, digits = 2)) %>%
    # Arrange by climatic optimum
    arrange(desc(Optimum))
  
  # Prepare the conversion in a table
  table = data.frame(
    col1 = c("", "", data$species), 
    col2 = c("Share in ", "basal area", paste0(data$percent, " %")), 
    col3 = c("Cold", "margin", data$Cold), 
    col4 = c("Climatic", "optimum", data$Optimum), 
    col5 = c("Hot", "margin", data$Hot), 
    col6 = c("GrSurv", "axis", data$GrSurv), 
    col7 = c("DroughtShade", "axis", data$ShadeDrought)
  )
  
  
  # Save tex file
  print(xtable(table, type = "latex", label = "table_species",
               caption = paste0("Share in basal area within the NFI dataset, ", 
                                "climatic niche (optimum, cold and hot margin ", 
                                "along the sgdd-wai pca axis) and trait value ", 
                                "along the growth-survival and shade - drought ", 
                                "tolerance trait axes for the 26 species included", 
                                " in the simulations"), 
               align = rep("c", dim(table)[2]+1)), 
        include.rownames=FALSE, hline.after = c(0, 2, dim(table)[1]), 
        include.colnames = FALSE, caption.placement = "top", size = "\\small",
        file = file.out)
  
  # Return the file saved
  return(file.out)
  
}




