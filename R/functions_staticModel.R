#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_chronosequence.R  
#' @description R script containing all functions relative to 
#               the static model analysis
#' @author Julien Barrere, Georges Kunstler
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to format climatic data for the SDM analysis
#' @param climate_dist_dflist climate and disturbance per plotcode
#' @param NFI_plots_selected plot level df of the plots selected for simulation
format_clim_SDM = function(climate_dist_dflist, NFI_plots_selected){
  
  # Prepare input dataset
  for(i in 1:length(names(climate_dist_dflist))){
    
    # Extract climatic data for each plot 
    df.i = rbind(data.frame(plotcode = names(climate_dist_dflist)[i], 
                            ssp = "ssp126", 
                            t = climate_dist_dflist[[i]]$rep1$ssp126$climate$t, 
                            sgdd = climate_dist_dflist[[i]]$rep1$ssp126$climate$sgdd, 
                            wai = climate_dist_dflist[[i]]$rep1$ssp126$climate$wai), 
                 data.frame(plotcode = names(climate_dist_dflist)[i], 
                            ssp = "ssp585", 
                            t = climate_dist_dflist[[i]]$rep1$ssp585$climate$t, 
                            sgdd = climate_dist_dflist[[i]]$rep1$ssp585$climate$sgdd, 
                            wai = climate_dist_dflist[[i]]$rep1$ssp585$climate$wai))
    
    # Add plot data to the output dataset
    if(i == 1) df = df.i
    else df = rbind(df, df.i)
    
  }
  
  # Average by climate
  out = df %>%
    left_join(NFI_plots_selected %>% dplyr::select(plotcode, climate), 
              by = "plotcode") %>%
    group_by(climate, ssp, t) %>%
    summarize(sgdd = mean(sgdd, na.rm = TRUE), 
              wai = mean(wai, na.rm = TRUE))
  
  # Return output dataset
  return(out)
  
}

#' Title
#'
#' @param df data base with plotcode species sgdd and wai with the climate of the plots and the species to predict in the plots
#' @param clim_pca pca coordinate 
#' @param coef_ba.reg.fec fitted models
#'
#' @return data.frame with same column as df but with the column ba.reg
#' @export
#'
#' @examples
predict.rel.ab <- function(df, clim_pca, coef_ba.reg.fec){
  # add pca1 pca2
  df$pca1 <- predict(clim_pca$pca, newdata = df)[, "PC1"]
  df$pca2 <- predict(clim_pca$pca, newdata = df)[, "PC2"]
  
  
  calc_fit <- function(x, y, a, b1, b2, c1, c2, rho) {
    term1 <- (x - b1)^2 / c1
    term2 <- (y - b2)^2 / c2
    term3 <- 2 * rho * (x - b1) * (y - b2) / sqrt(c1 * c2)
    a * exp(-0.5 * (1 / (1 - rho^2)) * (term1 + term2 - term3))
  }
  
  df_pred = df %>%
    left_join(coef_ba.reg.fec, by = "species") %>%
    mutate(fit = calc_fit(pca1, pca2,a, b1,b2, c1, c2, rho)) 
  return(df_pred)
}

#' Plot the change in sgdd and wai with time per scenario
#' @param data_clim_SDM climate data formated for SDm analysis
#' @param file.out Name of the file to save, including path
plot_climate_SDM = function(data_clim_SDM, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Plot the temporal change in sgdd and wai per climate 
  plot.out = data_clim_SDM %>% 
    pivot_longer(names_to = "var", values_to = "value", cols = c("sgdd", "wai")) %>% 
    mutate(climate = factor(climate, levels = paste0("clim", c(1:10)))) %>%
    ggplot(aes(x = t, y = value, color = ssp)) + 
    geom_line() + 
    facet_grid(var ~ climate, scales = "free") + 
    theme_bw()
  
  # Save the plot
  ggsave(file.out, plot.out, width = 25, height = 10 , units = "cm", 
         dpi = 600, bg = "white")
  
  # Return file saved
  return(file.out)
}

#' Function to calculate regional basal area from climate
#' @param clim_pca pca with sgdd and wai
#' @param data_clim_SDM climate df formatted for SDM analysis
#' @param coef_SDM coefficients of regional basal area model
#' @param species_list Dataframe of species present in each climate
get_reg_ba_SDM = function(clim_pca, data_clim_SDM, coef_SDM, species_list){
  
  # Merge climate data with species data
  data.in = data_clim_SDM %>%
    merge.data.frame(species_list[, c("climate", "species")], by = "climate") %>%
    mutate(species = gsub("\\_", "\\ ", species))
  
  # Apply GK function to calculate regional basal area per timestep
  out = predict.rel.ab(data.in, clim_pca, coef_SDM) %>%
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    dplyr::select(climate, ssp, t, pca1, sgdd, wai, species, ba.reg = fit)
  
  # Return output
  return(out)
}

#' Function to calculate species composition metrics for SDM analysis
#' @param traits_compiled list containing trait information per species
#' @param reg_ba_SDM df with regional basal area from SDM
get_sp.composition_SDM = function(traits_compiled, reg_ba_SDM){
  
  # Calculate functional diversity
  data.FD = reg_ba_SDM %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Gather by functional axis
    gather(key = "axis", value = "trait_value", names(traits_compiled[[1]])) %>%
    # Calculate cwm for each time step and each simulation
    group_by(climate, t, ssp, axis) %>%
    mutate(cwm = weighted.mean(trait_value, w = ba.reg)) %>%
    ungroup() %>%
    # Calculate square distance of species to centroid along each axis
    mutate(zsq_per_axis = (trait_value - cwm)^2) %>%
    # Calculate the distance of each species to the centroid
    group_by(climate, t, ssp, ba.reg, species) %>%
    summarize(z = sqrt(sum(zsq_per_axis))) %>% ungroup() %>%
    # Calculate functional dispersion
    group_by(climate, t, ssp) %>%
    summarize(FD = weighted.mean(z, w = ba.reg))
  
  
  # Extract simulation output on forest composition
  out = reg_ba_SDM %>%
    # Add traits data
    left_join(traits_compiled$species_coord, by = "species") %>%
    # Calculate composition metrics
    group_by(climate, t, ssp) %>%
    mutate(p = ba.reg/sum(ba.reg, na.rm = TRUE), 
           plnp = p*log(p)) %>%
    summarise(H = -sum(plnp, na.rm = TRUE), 
              cwm1 = weighted.mean(GrSurv, w = ba.reg), 
              cwm2 = weighted.mean(ShadeDrought, w = ba.reg)) %>%
    # Add functional diversity
    left_join(data.FD,  by = c("climate", "t", "ssp")) %>%
    # Final formatting
    ungroup() 
  
  
  
  # Return output
  return(out)
  
}

#' Plot change in species compositon metric across time per ssp and climate
#' @param sp.composition_SDM species composition at each timestep per scenario
#' @param file.out Name of the file to save, including path
plot_sp.composition_SDM = function(sp.composition_SDM, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Make the plot 
  plot.out = sp.composition_SDM %>%
    pivot_longer(names_to = "var", values_to = "value", cols = c("H", "FD", "cwm1", "cwm2")) %>% 
    mutate(climate = factor(climate, levels = paste0("clim", c(1:10))), 
           var = factor(var, levels = c("H", "FD", "cwm1", "cwm2"))) %>%
    ggplot(aes(x = t, y = value, color = ssp)) + 
    geom_line() + 
    facet_grid(var ~ climate, scales = "free") + 
    theme_bw() + 
    geom_vline(xintercept = 70, linetype = "dashed", color = "#344E41")
  
  # Save the plot
  ggsave(file.out, plot.out, width = 25, height = 10 , units = "cm", 
         dpi = 600, bg = "white")
  
  # Return file saved
  return(file.out)
}

#' Function to calculate phi for the SDM analysis
#' @param sp.composition_SDM species composition at each timestep per scenario
#' @param data_pool output of simulations with regional pool formatted 
get_delta_per_scenario_SDM = function(sp.composition_SDM, data_pool){
  
  # Calculate the range of each composition variable per climate 
  #     with same approach as in the simulations
  data_range = data_pool %>%
    filter(metric == "BA") %>%
    group_by(plotcode, variable, rep) %>%
    summarize(min = min(value, na.rm = TRUE), 
              max = max(value, na.rm = TRUE)) %>%
    mutate(diff = max - min) %>%
    ungroup() %>% group_by(variable) %>%
    summarize(range = mean(diff, na.rm = TRUE))
  
  # Calculate phi per scenario
  out = sp.composition_SDM %>%
    pivot_longer(names_to = "variable", values_to = "value", 
                 cols = c("H", "FD", "cwm1", "cwm2")) %>%
    filter(t %in% range(data_pool$time)) %>%
    mutate(time_cat = ifelse(t == min(data_pool$time), "tmin", "tmax")) %>%
    dplyr::select(- t) %>%
    pivot_wider(names_from = "time_cat", values_from = "value") %>%
    mutate(Delta = tmax - tmin) %>%
    dplyr::select(-tmin, -tmax) %>% 
    pivot_wider(names_from = ssp, values_from = Delta) %>%
    left_join(data_range, by = c("variable")) %>%
    mutate(delta = (ssp585 - ssp126)/range*100) %>%
    dplyr::select(-range, -ssp585, -ssp126)
  
  # Return output
  return(out)
  
}


#' Function to plot the comparison of delta in simulations and in SDM
#' @param delta_per_climate_SDM phi calculated from SDM 
#' @param delta phi calculated from simlulations
#' @param file.out Name of the file to save, including path
plot_delta_simulations_vs_SDM = function(
    delta_per_climate_SDM, delta, file.out){
  
  # Create output directory if needed
  create_dir_if_needed(file.out)
  
  # Plot of phi dereived from simulations
  plot.simulations = delta %>%
    filter(pool == "pool" & dist == "dist" & metric == "BA") %>%
    mutate(climate = factor(climate, levels = paste0("clim", c(1:10))), 
           variable = factor(variable, levels = c("H", "FD", "cwm1", "cwm2"))) %>%
    group_by(climate, variable) %>%
    summarize(mean = mean(delta, na.rm = TRUE), 
              se = sd(delta, na.rm = TRUE)/sqrt(n())) %>% 
    ggplot(aes(x = climate)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0) +
    geom_point(aes(y = mean), shape = 21, color = "black", fill = "grey") + 
    facet_wrap( ~ variable) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          axis.text.x = element_text(angle = 90)) + 
    ylab("\u03b4 (% of average\ncommunity change)") + xlab("")
  
  # Plot of phi derived from SDM
  plot.SDM = delta_per_climate_SDM %>% 
    mutate(climate = factor(climate, levels = paste0("clim", c(1:10))), 
           variable = factor(variable, levels = c("H", "FD", "cwm1", "cwm2"))) %>% 
    ggplot(aes(x = climate, y = delta)) + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    geom_point(shape = 21, color = "black", fill = "grey") + 
    facet_wrap( ~ variable) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          axis.text.x = element_text(angle = 90)) + 
    ylab("\u03b4 (% of average\ncommunity change)") + xlab("")
  
  # Assemble the two plots
  plot.out = plot_grid(plot.simulations, plot.SDM, align = "hv", nrow = 1,
                       labels = c("(a) Simulations", "(b) Pseudo-SDM"), 
                       scale = 0.9)
  
  # Save the plot
  ggsave(file.out, plot.out, width = 25, height = 15 , units = "cm", 
         dpi = 600, bg = "white")
  
  
  # Return file saved
  return(file.out)
}


