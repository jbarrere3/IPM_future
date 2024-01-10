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


#' Function to plot the relation between fire intensity and vpd
#' @param model_fire_vpd list containing the model and associated data
#' @param file.out name of the file to save, including path
plot_fire_vpd = function(model_fire_vpd, file.out){
  
  # Create directory if needed
  create_dir_if_needed(file.out)
  
  # Add predictions of the model to the data
  data = model_fire_vpd$data %>%
    ungroup() %>%
    mutate(fit.lwr = predict(model_fire_vpd$model, newdata = ., 
                             type = "quantile", at = 0.025), 
           fit.upr = predict(model_fire_vpd$model, newdata = ., 
                             type = "quantile", at = 0.975), 
           fit = predict(model_fire_vpd$model, newdata = ., type = "response"))
  
  # Plot data
  plot.out = data %>%
    ggplot(aes(x = vpd.mean, y = I.mean)) + 
    geom_ribbon(aes(ymin = fit.lwr, ymax = fit.upr), alpha = 0.5, color = NA, 
                fill = "#F8961E") +
    geom_point(aes(size = weight), shape = 21, color = "black", 
               fill = "#F9844A", alpha = 0.6) + 
    geom_line(aes(y = fit), color = "#F94144", size = 1) + 
    xlab("Mean summer VPD (kPa)") + ylab("Mean fire intensity") + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          legend.position = "none") + 
    ggtitle(paste0(
      "Chisq = ", as.numeric(round(Anova(model_fire_vpd$model)[2], digits = 2)), 
      " - ", pvalue(as.numeric(Anova(model_fire_vpd$model)[3]), add_p = TRUE)
    ))
  
  # Export the plot
  ggsave(file.out, plot.out, 
         width = 14, height = 10, units = "cm", dpi = 600, bg = "white")
}


#' Plot the effect of climate and disturbances on productivity
#' @param sim_output formatted outputs of the simulations
#' @param simul_list Informations on the simulations made
#' @param pca1_per_species coordinates on the trait pca of each species
#' @param dir.in Directory where to save figures
plot_climate_effect = function(sim_output, simul_list, pca1_per_species, dir.in){
  
  # Extract simulation output on forest structure and productivity
  data.str = sim_output %>%
    filter(species == "all") %>%
    group_by(ID.simulation) %>%
    summarize(prod = mean(prod, na.rm = TRUE), 
              dbh.mean = mean(dbh.mean, na.rm = TRUE), 
              dbh.var = mean(dbh.var, na.rm = TRUE)) %>%
    drop_na() 
  
  # Extract simulation output on forest composition
  data.sp = sim_output %>%
    filter(species != "all" & H > 0) %>%
    left_join(pca1_per_species, by = "species") %>%
    group_by(ID.simulation, time) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = weighted.mean(pca1, w = BA), 
              traits.pca1.var = weighted.var(pca1, w = BA)) %>%
    ungroup() %>% group_by(ID.simulation) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = mean(traits.pca1.mean, na.rm = TRUE), 
              traits.pca1.var = mean(traits.pca1.var, na.rm = TRUE)) %>%
    drop_na() 
  
  # Join all data together
  data = simul_list %>%
    left_join(data.str, by = "ID.simulation") %>%
    left_join(data.sp, by = "ID.simulation")
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("prod", "dbh.mean", "dbh.var", "H", "traits.pca1.mean", "traits.pca1.var"), 
    title = c("Productivity", "Mean_diameter", "Variance_diameter", 
              "Species_diversity", "Mean_traits", "Variance_traits"), 
    label = c("Productivity\n(m2/ha/year)", "Mean diameter\n(mm)", 
              "Variance of\ndiameter", "Species diversity\n(Shannon index)", 
              "Mean trait value\n(Fast growth -> high survival)", 
              "Functional diversity"), 
    transfo = c("log", "log", "log", "log", "none", "log")
  )
  
  # Initialize output
  out = c()
  
  # Loop on all variables for which to make analyses
  for(j in 1:dim(data.var)[1]){
    
    # Restrict the dataset to the variable j
    data.in = data %>%
      mutate(var = data[, data.var$var[j]]) %>%
      dplyr::select(var, ssp, dist, climate, plotcode, pca1) %>%
      drop_na() %>%
      # To double check
      filter(var > 0)
    
    # All climates for which to fit model
    climates.in = paste0("clim", c(1:length(unique(data.in$climate))))
    
    # Initialize the data containing the text to add to the plot
    data.text = data.frame(climate = climates.in, text = NA_character_)
    
    # Initialize plot list for the residuals
    plotlist.resid = vector(mode = "list", length = length(climates.in))
    
    # Loop on all climates
    for(i in 1:length(climates.in)){
      
      # Log transform if needed
      if(data.var$transfo[j] == "log"){
        mod.i = lmer(log(var) ~ ssp*dist + (1|plotcode), 
                     data = subset(data.in, climate == climates.in[i]))
      }else{
        mod.i = lmer(var ~ ssp*dist + (1|plotcode), 
                     data = subset(data.in, climate == climates.in[i]))
      }
      
      
      # Add the stat result to data text
      data.text$text[i] = paste(paste(rownames(Anova(mod.i)), 
                                      pvalue(Anova(mod.i)[, 3], add_p = TRUE, accuracy = 0.01), 
                                      sep = ": "), 
                                collapse = "\n")
      # Make a plot for the residuals
      plotlist.resid[[i]] = ggplot(augment(mod.i), aes(.fitted, .resid)) + 
        geom_point() +
        geom_hline(yintercept = 0, linetype = "dashed") + 
        geom_smooth(method = "loess") + 
        ggtitle(climates.in[i])
    }
    
    # Make the final residual plot
    plot.residuals = plot_grid(plotlist = plotlist.resid, scale = 0.9, 
                               nrow = 2)
    
    # Data for plotting
    data.plot = data.in %>%
      group_by(ssp, climate, dist) %>%
      summarize(pca.mean = mean(pca1, na.rm = TRUE), 
                mean = mean(var, na.rm = TRUE), 
                sd = sd(var, na.rm = TRUE), 
                lwr = quantile(var, 0.05, na.rm = TRUE), 
                upr = quantile(var, 0.95, na.rm = TRUE)) %>%
      mutate(x.pos = case_when(
        ssp == "ssp126" ~ pca.mean - 0.02*diff(range(.$pca.mean)), 
        ssp == "ssp585" ~ pca.mean + 0.02*diff(range(.$pca.mean)), 
        TRUE ~ pca.mean
      )) %>%
      mutate(x.pos = ifelse(dist == "dist", x.pos - 0.001*diff(range(.$pca.mean)), 
                            x.pos + 0.001*diff(range(.$pca.mean))))
    
    # Add to data text the x and y position
    data.text = data.text %>%
      left_join((data.plot %>% group_by(climate) %>% 
                   summarize(x.pos = min(x.pos)) %>% ungroup() %>% 
                   mutate(ssp = NA, dist = NA, mean = 1.5*max(data.plot$upr))), 
                by = "climate")
    
    # Plots with groups and facet
    plot.effect = data.plot %>%
      ggplot(aes(x = x.pos, y = mean, group = interaction(ssp, dist), 
                 color = ssp, fill = ssp)) + 
      geom_errorbar(aes(ymin = lwr, ymax = upr), 
                    width = 0) + 
      geom_line(aes(linetype = dist)) +
      geom_point(aes(shape = dist)) + 
      facet_wrap(~ factor(climate, levels = paste0(
        "clim", c(1:length(unique(simul_list$climate))))), 
        scales = "free_x", nrow = 2) + 
      scale_shape_manual(values = c('dist' = 22, 'nodist' = 23)) +
      scale_color_manual(values = c('ssp126' = "#005F73", 'ssp370' = "#CA6702", 
                                    'ssp585' = "#9B2226")) + 
      scale_fill_manual(values = c('ssp126' = "#0A9396", 'ssp370' = "#EE9B00", 
                                   'ssp585' = "#AE2012")) + 
      geom_text(data = data.text, aes(label = text), inherit.aes = TRUE, 
                hjust = "inward", vjust = "inward", size = 2.5, alpha = 0.8) + 
      xlab("") + ylab(data.var$label[j]) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            strip.background = element_blank(), 
            strip.text = element_text(face = "bold"), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.key = element_blank())
    
    # Name of the files to save
    file.residuals = paste0(dir.in, "/fig_residuals_", data.var$title[j], ".jpg")
    file.effect = paste0(dir.in, "/fig_effect_", data.var$title[j], ".jpg")
    
    # Create directory if needed
    create_dir_if_needed(file.effect)
    
    # Save plots
    ggsave(file.residuals, plot.residuals, width = 25, height = 10 , units = "cm", 
           dpi = 600, bg = "white")
    ggsave(file.effect, plot.effect, width = 25, height = 10 , units = "cm", 
           dpi = 600, bg = "white")
    
    # Add the saved plots to the output vector
    out = c(out, file.residuals, file.effect)
    
  }
  
  
  # Return the files saved
  return(out)
  
  
}




#' Plot the effect of climate and disturbances on temporal change in different variables
#' @param sim_output formatted outputs of the simulations
#' @param simul_list Informations on the simulations made
#' @param pca1_per_species coordinates on the trait pca of each species
#' @param dir.in Directory where to save figures
plot_temporalchange = function(sim_output, simul_list, pca1_per_species, dir.in){
  
  # Extract simulation output on forest structure and productivity
  data.str = sim_output %>%
    filter(species == "all") %>%
    mutate(timing = case_when(time <= 30 ~ "beg", 
                              time >= max(sim_output$time) - 30 ~ "end", 
                              TRUE ~ "mid")) %>%
    filter(timing != "mid") %>% group_by(ID.simulation, timing) %>%
    summarize(prod = mean(prod, na.rm = TRUE), 
              dbh.mean = mean(dbh.mean, na.rm = TRUE), 
              dbh.var = mean(dbh.var, na.rm = TRUE)) %>%
    drop_na() 
  
  # Extract simulation output on forest composition
  data.sp = sim_output %>%
    filter(species != "all" & H > 0) %>%
    left_join(pca1_per_species, by = "species") %>%
    group_by(ID.simulation, time) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = weighted.mean(pca1, w = BA), 
              traits.pca1.var = weighted.var(pca1, w = BA)) %>%
    mutate(timing = case_when(time <= 30 ~ "beg", 
                              time >= max(sim_output$time) - 30 ~ "end", 
                              TRUE ~ "mid")) %>%
    ungroup() %>% filter(timing != "mid") %>% group_by(ID.simulation, timing) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = mean(traits.pca1.mean, na.rm = TRUE), 
              traits.pca1.var = mean(traits.pca1.var, na.rm = TRUE)) %>%
    drop_na() 
  
  # Join all data together
  data = data.str %>%
    left_join(data.sp, by = c("ID.simulation", "timing")) %>%
    left_join(simul_list, by = "ID.simulation") %>%
    as.data.frame()
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("prod", "dbh.mean", "dbh.var", "H", "traits.pca1.mean", "traits.pca1.var"), 
    title = c("Productivity", "Mean_diameter", "Variance_diameter", 
              "Species_diversity", "Mean_traits", "Variance_traits"), 
    label = c("Productivity", "Mean diameter", 
              "Structural diversity", "Species diversity", 
              "Mean trait \n(Fast growth -> high survival)", 
              "Functional diversity"), 
    transfo = c("log", "log", "log", "log", "none", "log")
  )
  
  # Initialize output
  out = c()
  
  # FIX FUCKING BUG §§§§§
  # Loop on all variables for which to make analyses
  for(j in 1:dim(data.var)[1]){
    
    # Restrict the dataset to the variable j
    data.in = data %>%
      ungroup() %>%
      mutate(var = data[, data.var$var[j]], 
             pca1sq = pca1^2) %>%
      dplyr::select(var, climate, ssp, dist, timing, plotcode, pca1, pca1sq) %>%
      drop_na() %>%
      spread(key = "timing", value = "var")%>%
      mutate(var.change = log(end/beg), 
             ssp = factor(ssp, levels = c("ssp126", "ssp370", "ssp585"))) %>%
      filter(!(beg == 0 | end == 0))
    
    # Fit a first model
    mod = lmer(var.change ~ pca1*ssp*dist + pca1sq*ssp*dist + (1|plotcode), 
               data = data.in)
    # Initialize boolean to stop model selection and counter of interactions removed
    selection.over = FALSE; k=0
    # Start while loop
    while(!selection.over){
      print(k)
      # Simple and double interaction terms
      terms.all.in = rownames(Anova(mod))[grep(":", rownames(Anova(mod)))]
      terms.double.in = rownames(Anova(mod))[grep(":.+:", rownames(Anova(mod)))]
      terms.simple.in = terms.all.in[which(!(terms.all.in %in% terms.double.in))]
      # Initialize formula of the model
      form = paste(gsub("\\:", "\\*", terms.all.in), collapse = " + ")
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
          form = paste(gsub("\\:", "\\*", terms.all.in), collapse = " + ")
          # -- Fit model with the new formula 
          eval(parse(text = paste0(
            "mod = lmer(var.change ~ ", form, " + (1|plotcode), data = data.in)")))
        }
      }else{
        # Second configuration : no more double interactions but some simple are left
        if(length(terms.simple.in) > 1){
          # First sub-configuration : all simple interactions are significant
          if(!any(Anova(mod)[terms.simple.in, 3] > 0.05)){
            # We stop model selection, no interactions can be removed
            selection.over = TRUE
          }else{
            # Otherwise, we remove the least significant simple interaction 
            # -- increase counter
            k = k+1
            # -- remove from simple iteractions the least significant
            terms.simple.in = terms.simple.in[
              -which(Anova(mod)[terms.simple.in, 3] == max(Anova(mod)[terms.simple.in, 3]))]
            # -- update vector containing all terms
            terms.all.in = c(terms.simple.in, terms.double.in)
            # -- New formula
            form = paste(gsub("\\:", "\\*", terms.all.in), collapse = " + ")
            # -- Fit model with the new formula 
            eval(parse(text = paste0(
              "mod = lmer(var.change ~ ", form, " + (1|plotcode), data = data.in)")))
          }
          # Last configuration: only one simple interaction left: keep the model as it is
        }else{
          selection.over = TRUE
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
      pca1 = seq(from = quantile(data.in$pca1, 0.01), 
                 to = quantile(data.in$pca1, 0.99), length.out = 100), 
      ssp = unique(data.in$ssp), dist = unique(data.in$dist)) %>%
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
    
    # Plot residuals
    plot.residuals = ggplot(augment(mod), aes(.fitted, .resid)) + 
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_smooth(method = "loess")
    
    # Plot predictions of the model
    plot.effect = newdata %>%
      ggplot(aes(x = pca1, y = fit, group = interaction(ssp, dist), 
                 color = ssp, fill = ssp, linetype = dist)) + 
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
      geom_line() + 
      xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
      ylab(paste0(data.var$label[j], "\nchange (logratio)")) +
      scale_linetype_manual(values = c('dist' = "dashed", 'nodist' = "solid")) +
      scale_color_manual(values = c('ssp126' = "#005F73", 'ssp370' = "#CA6702", 
                                    'ssp585' = "#9B2226")) + 
      scale_fill_manual(values = c('ssp126' = "#0A9396", 'ssp370' = "#EE9B00", 
                                   'ssp585' = "#AE2012")) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.key = element_blank()) 
    
    # Add a horizontal line if zero is within the range of predicted values
    if(min(newdata$lwr) < 0 & max(newdata$upr) > 0) plot.effect = plot.effect + 
      geom_hline(yintercept = 0, linetype = "dashed")
    
    # Name of the files to save
    file.residuals = paste0(dir.in, "/fig_residuals_", data.var$title[j], ".jpg")
    file.effect = paste0(dir.in, "/fig_effect_", data.var$title[j], ".jpg")
    
    # Create directory if needed
    create_dir_if_needed(file.effect)
    
    # Save plots
    ggsave(file.residuals, plot.residuals, width = 17, height = 9 , units = "cm", 
           dpi = 600, bg = "white")
    ggsave(file.effect, plot.effect, width = 17, height = 9 , units = "cm", 
           dpi = 600, bg = "white")
    
    # Add the saved plots to the output vector
    out = c(out, file.residuals, file.effect)
    
  }
  
  
  # Return the files saved
  return(out)
  
}


#' Plot the effect of climate and disturbances on temporal change in different variables
#' @param sim_output formatted outputs of the simulations
#' @param simul_list Informations on the simulations made
#' @param pca1_per_species coordinates on the trait pca of each species
#' @param dir.in Directory where to save figures
plot_magnitudechange = function(sim_output, simul_list, pca1_per_species, dir.in){
  
  # Extract simulation output on forest structure and productivity
  data.str = sim_output %>%
    filter(species == "all" & time > floor(max(sim_output$time)/2)) %>%
    group_by(ID.simulation) %>%
    summarize(prod = mean(prod, na.rm = TRUE), 
              dbh.mean = mean(dbh.mean, na.rm = TRUE), 
              dbh.var = mean(dbh.var, na.rm = TRUE)) %>%
    drop_na() 
  
  # Extract simulation output on forest composition
  data.sp = sim_output %>%
    filter(species != "all" & H > 0) %>%
    filter(time > floor(max(sim_output$time)/2)) %>%
    left_join(pca1_per_species, by = "species") %>%
    group_by(ID.simulation, time) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = weighted.mean(pca1, w = BA), 
              traits.pca1.var = weighted.var(pca1, w = BA)) %>%
    ungroup() %>% group_by(ID.simulation) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = mean(traits.pca1.mean, na.rm = TRUE), 
              traits.pca1.var = mean(traits.pca1.var, na.rm = TRUE)) %>%
    drop_na() 
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("prod", "dbh.mean", "dbh.var", "H", "traits.pca1.mean", "traits.pca1.var"), 
    title = c("Productivity", "Mean_diameter", "Variance_diameter", 
              "Species_diversity", "Mean_traits", "Variance_traits"), 
    label = c("Productivity", "Mean diameter", 
              "Structural diversity", "Species diversity", 
              "Mean trait \n(Fast growth -> high survival)", 
              "Functional diversity"), 
    transfo = c("log", "log", "log", "log", "none", "log")
  )
  
  # Join all data together
  data = data.str %>%
    left_join(data.sp, by = "ID.simulation") %>%
    left_join(simul_list, by = "ID.simulation") %>%
    gather(key = "variable", value = "value", data.var$var) %>%
    drop_na()
  data = data %>%
    left_join((data %>% filter(dist == "nodist" & ssp == "ssp126") %>%
                 dplyr::select(plotcode, variable, value.ref = value)), 
              by = c("plotcode", "variable")) %>%
    mutate(var.change = 100*(value-value.ref)/value.ref) %>%
    as.data.frame()
  
  
  # Initialize output
  out = c()
  
  # Loop on all variables for which to make analyses
  for(j in 1:dim(data.var)[1]){
    
    # Restrict the dataset to the variable j
    data.in = data %>%
      ungroup() %>% filter(variable == data.var$var[j]) %>%
      mutate(pca1sq = pca1^2, 
             ssp = factor(ssp, levels = c("ssp126", "ssp370", "ssp585"))) %>%
      dplyr::select(plotcode, climate, ssp, dist, pca1, pca1sq, var.change) %>%
      filter(!is.infinite(var.change)) %>%
      filter(var.change < quantile(.$var.change, 0.99, na.rm = TRUE))
    
    # Fit a first model
    mod = lmer(var.change ~ pca1*ssp*dist + pca1sq*ssp*dist + (1|plotcode), 
               data = data.in)
    # Initialize boolean to stop model selection and counter of interactions removed
    selection.over = FALSE; k=0
    # Start while loop
    while(!selection.over){
      print(k)
      # Simple and double interaction terms
      terms.all.in = rownames(Anova(mod))[grep(":", rownames(Anova(mod)))]
      terms.double.in = rownames(Anova(mod))[grep(":.+:", rownames(Anova(mod)))]
      terms.simple.in = terms.all.in[which(!(terms.all.in %in% terms.double.in))]
      # Initialize formula of the model
      form = paste(gsub("\\:", "\\*", terms.all.in), collapse = " + ")
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
          form = paste(gsub("\\:", "\\*", terms.all.in), collapse = " + ")
          # -- Fit model with the new formula 
          eval(parse(text = paste0(
            "mod = lmer(var.change ~ ", form, " + (1|plotcode), data = data.in)")))
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
            form = paste(gsub("\\:", "\\*", terms.all.in), collapse = " + ")
            # -- Fit model with the new formula 
            eval(parse(text = paste0(
              "mod = lmer(var.change ~ ", form, " + (1|plotcode), data = data.in)")))
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
      pca1 = seq(from = quantile(data.in$pca1, 0.01), 
                 to = quantile(data.in$pca1, 0.99), length.out = 100), 
      ssp = unique(data.in$ssp), dist = unique(data.in$dist)) %>%
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
    
    # Prepare points data for plotting
    data.points = data.in %>%
      group_by(climate, ssp, dist) %>%
      summarize(pca1 = mean(pca1), 
                fit = mean(var.change, na.rm = TRUE), 
                lwr = quantile(var.change, 0.05, na.rm = TRUE), 
                upr = quantile(var.change, 0.95, na.rm = TRUE), 
                sd = sd(var.change, na.rm = TRUE)) %>%
      mutate(lwr = fit - sd, upr = fit + sd) %>%
      mutate(x.pos = case_when(
        ssp == "ssp126" ~ pca1 - 0.02*diff(range(.$pca1)), 
        ssp == "ssp585" ~ pca1 + 0.02*diff(range(.$pca1)), 
        TRUE ~ pca1
      )) %>%
      mutate(x.pos = ifelse(dist == "dist", x.pos - 0.0015*diff(range(.$pca1)), 
                            x.pos + 0.0015*diff(range(.$pca1))))
    
    # Plot residuals
    plot.residuals = ggplot(augment(mod), aes(.fitted, .resid)) + 
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_smooth(method = "loess")
    
    # Plot predictions of the model
    plot.effect = newdata %>%
      ggplot(aes(x = pca1, y = fit, group = interaction(ssp, dist), 
                 color = dist, fill = ssp, ymin = lwr, ymax = upr)) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_errorbar(data = data.points, aes(x = x.pos), 
                    inherit.aes = TRUE, alpha = 0.5) +
      geom_point(data = data.points, aes(x = x.pos), 
                 inherit.aes = TRUE, shape = 21) + 
      geom_line() + 
      geom_ribbon(alpha = 0.3, color = NA) + 
      xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
      ylab(paste0("Change in \n", data.var$label[j], " (%)")) + 
      scale_color_manual(values = c('dist' = "#343A40", 'nodist' = "#ADB5BD")) +
      scale_fill_manual(values = c('ssp126' = "#0A9396", 'ssp370' = "#EE9B00", 
                                   'ssp585' = "#AE2012")) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.key = element_blank())
    
    # Name of the files to save
    file.residuals = paste0(dir.in, "/fig_residuals_", data.var$title[j], ".jpg")
    file.effect = paste0(dir.in, "/fig_effect_", data.var$title[j], ".jpg")
    
    # Create directory if needed
    create_dir_if_needed(file.effect)
    
    # Save plots
    ggsave(file.residuals, plot.residuals, width = 17, height = 9 , units = "cm", 
           dpi = 600, bg = "white")
    ggsave(file.effect, plot.effect, width = 17, height = 9 , units = "cm", 
           dpi = 600, bg = "white")
    
    # Add the saved plots to the output vector
    out = c(out, file.residuals, file.effect)
    
  }
  
  
  # Return the files saved
  return(out)
  
  
}
