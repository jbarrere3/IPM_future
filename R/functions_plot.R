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
              dbh.var = mean(dbh.var, na.rm = TRUE), 
              stock = mean(BA, na.rm = TRUE)) %>%
    drop_na() 
  
  
  # Extract simulation output on forest composition
  data.sp = sim_output %>%
    filter(species != "all") %>%
    # Only focus on the second half of the simulations
    filter(time > floor(max(sim_output$time)/2)) %>%
    # Add traits data
    left_join(pca1_per_species, by = "species") %>%
    # Add climate margins data
    left_join((matreex::climate_species %>% filter(N == 2) %>%
                 dplyr::select(species = sp, optimum = PC1)), 
              by = "species") %>%
    group_by(ID.simulation, time) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = weighted.mean(pca1, w = BA), 
              traits.pca2.mean = weighted.mean(pca2, w = BA), 
              opticlim.mean = weighted.mean(optimum, w = BA)) %>%
    # Average metrics
    ungroup() %>% group_by(ID.simulation) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              cwm1 = mean(traits.pca1.mean, na.rm = TRUE), 
              cwm2 = mean(traits.pca2.mean, na.rm = TRUE), 
              opticlim = mean(opticlim.mean, na.rm = TRUE)) %>%
    drop_na() 
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("prod", "dbh.mean", "dbh.var", "H", "cwm1", 
            "cwm2", "stock", "opticlim"), 
    title = c("Productivity", "Mean_diameter", "Variance_diameter", 
              "Species_diversity", "CWM1", "CWM2", "Stocking", "Climate_optimum"), 
    label = c("Productivity", "Mean diameter", 
              "Structural diversity", "Species diversity", 
              "CWM1 \nHigh growth <--> High survival", 
              "CWM2 \nHigh <--> Low recruitment", 
              "Total Basal Area\n(m2/ha)", 
              "Mean species climatic optimum\n(hot-dry to cold-wet)"), 
    transfo = c("log", "log", "log", "log", "none", "none", "log", "none")
  )
  
  # Join all data together
  data = data.str %>%
    left_join(data.sp, by = "ID.simulation") %>%
    left_join(simul_list, by = "ID.simulation") %>%
    gather(key = "variable", value = "value", data.var$var) %>%
    drop_na()
  
  # Calculate the change relative to the reference scenario (ssp126 and no dist)
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
             ssp = factor(ssp, levels = c("ssp126", "ssp585"))) %>%
      dplyr::select(plotcode, climate, ssp, dist, pca1, pca1sq, var.change) %>%
      filter(!is.infinite(var.change)) %>%
      filter(var.change < quantile(.$var.change, 0.99, na.rm = TRUE) & 
               var.change > quantile(.$var.change, 0.01, na.rm = TRUE))
    
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
      scale_fill_manual(values = c('ssp126' = "#0A9396", 'ssp585' = "#AE2012")) +
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


#' Plot the effect of climate and disturbances on temporal change in different variables
#' @param sim_output formatted outputs of the simulations
#' @param simul_list Informations on the simulations made
#' @param pca1_per_species coordinates on the trait pca of each species
#' @param dir.in Directory where to save figures
plot_magnitudechange_2 = function(sim_output, simul_list, pca1_per_species, dir.in){
  
  # Extract simulation output on forest structure and productivity
  data.str = sim_output %>%
    filter(species == "all" & time > floor(max(sim_output$time)/2)) %>%
    group_by(ID.simulation) %>%
    summarize(prod = mean(prod, na.rm = TRUE), 
              dbh.mean = mean(dbh.mean, na.rm = TRUE), 
              dbh.var = mean(dbh.var, na.rm = TRUE), 
              stock = mean(BA, na.rm = TRUE)) %>%
    drop_na() 
  
  
  # Extract simulation output on forest composition
  data.sp = sim_output %>%
    filter(species != "all") %>%
    # Only focus on the second half of the simulations
    filter(time > floor(max(sim_output$time)/2)) %>%
    # Add traits data
    left_join(pca1_per_species, by = "species") %>%
    # Add climate margins data
    left_join((climate_species %>% filter(N == 2) %>%
                 dplyr::select(species = sp, optimum = PC1)), 
              by = "species") %>%
    group_by(ID.simulation, time) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              traits.pca1.mean = weighted.mean(pca1, w = BA), 
              traits.pca2.mean = weighted.mean(pca2, w = BA), 
              opticlim.mean = weighted.mean(optimum, w = BA)) %>%
    # Average metrics
    ungroup() %>% group_by(ID.simulation) %>%
    summarize(H = mean(H, na.rm = TRUE), 
              cwm1 = mean(traits.pca1.mean, na.rm = TRUE), 
              cwm2 = mean(traits.pca2.mean, na.rm = TRUE), 
              opticlim = mean(opticlim.mean, na.rm = TRUE)) %>%
    drop_na() 
  
  # Build data frame listing the variables to analyse
  data.var = data.frame(
    var = c("prod", "dbh.mean", "dbh.var", "H", "cwm1", 
            "cwm2", "stock", "opticlim"), 
    title = c("Productivity", "Mean_diameter", "Variance_diameter", 
              "Species_diversity", "CWM1", "CWM2", "Stocking", "Climate_optimum"), 
    label = c("Productivity", "Mean diameter", 
              "Structural diversity", "Species diversity", 
              "CWM1 \nHigh growth <--> High survival", 
              "CWM2 \nHigh <--> Low recruitment", 
              "Total Basal Area\n(m2/ha)", 
              "Mean species climatic optimum\n(hot-dry to cold-wet)"), 
    transfo = c("log", "log", "log", "log", "none", "none", "log", "none")
  )
  
  # Join all data together
  data = data.str %>%
    left_join(data.sp, by = "ID.simulation") %>%
    left_join(simul_list, by = "ID.simulation") %>%
    gather(key = "variable", value = "value", data.var$var) %>%
    drop_na()
  
  # Calculate the change relative to the reference scenario (ssp126 and no dist)
  data = data %>%
    left_join((data %>% filter(dist == "nodist" & ssp == "ssp126") %>%
                 dplyr::select(plotcode, variable, value.ref = value)), 
              by = c("plotcode", "variable")) %>%
    mutate(var.change = 100*(value-value.ref)/value.ref, 
           scenario = case_when(
             dist == "nodist" & ssp == "ssp126" ~ "reference", 
             dist == "nodist" & ssp == "ssp585" ~ "Climate change only", 
             dist == "dist" & ssp == "ssp126" ~ "Disturbance only", 
             dist == "dist" & ssp == "ssp585" ~ "Disturbance and climate change"
           )) %>%
    as.data.frame()
  
  
  # Initialize output
  out = c()
  
  # Loop on all variables for which to make analyses
  for(j in 1:dim(data.var)[1]){
    
    # Restrict the dataset to the variable j
    data.in = data %>%
      ungroup() %>% filter(variable == data.var$var[j]) %>%
      mutate(pca1sq = pca1^2, 
             ssp = factor(ssp, levels = c("ssp126", "ssp585"))) %>%
      dplyr::select(plotcode, climate, ssp, dist, pca1, pca1sq, var.change, scenario) %>%
      filter(!is.infinite(var.change)) %>%
      filter(var.change < quantile(.$var.change, 0.99, na.rm = TRUE) & 
               var.change > quantile(.$var.change, 0.01, na.rm = TRUE)) %>%
      filter(scenario != "reference")
    
    # Fit a first model
    mod = lmer(var.change ~ pca1*scenario + pca1sq*scenario + (1|plotcode), 
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
      scenario = unique(data.in$scenario)[order(unique(data.in$scenario))]) %>%
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
      group_by(climate, scenario) %>%
      summarize(pca1 = mean(pca1), 
                fit = mean(var.change, na.rm = TRUE), 
                lwr = quantile(var.change, 0.05, na.rm = TRUE), 
                upr = quantile(var.change, 0.95, na.rm = TRUE), 
                sd = sd(var.change, na.rm = TRUE)) %>%
      mutate(lwr = fit - sd, upr = fit + sd) %>%
      mutate(x.pos = case_when(
        scenario == "Disturbance only" ~ pca1 - 0.02*diff(range(.$pca1)), 
        scenario == "Disturbance and climate change" ~ pca1 + 0.02*diff(range(.$pca1)), 
        TRUE ~ pca1
      )) 
    
    # Plot residuals
    plot.residuals = ggplot(augment(mod), aes(.fitted, .resid)) + 
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_smooth(method = "loess")
    
    # Plot predictions of the model
    plot.effect = newdata %>%
      ggplot(aes(x = pca1, y = fit, group = scenario, 
                 color = scenario, fill = scenario, ymin = lwr, ymax = upr)) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_errorbar(data = data.points, aes(x = x.pos), 
                    inherit.aes = TRUE, alpha = 0.5) +
      geom_point(data = data.points, aes(x = x.pos), 
                 inherit.aes = TRUE, shape = 21) + 
      geom_line() + 
      geom_ribbon(alpha = 0.3, color = NA) + 
      xlab("Position along the climate axis\n(Hot-dry to cold-wet)") + 
      ylab(paste0("Change in \n", data.var$label[j], " (%)")) + 
      scale_color_manual(values = c('Disturbance only' = "#001219", 
                                    'Climate change only' = "#CA6702", 
                                    'Disturbance and climate change' = "#9B2226")) +
      scale_fill_manual(values = c('Disturbance only' = "#005F73", 
                                   'Climate change only' = "#EE9B00", 
                                   'Disturbance and climate change' = "#AE2012")) +
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






#' Plot traits PCA
#' @param traits dataframe containing trait value per species
#' @param sp.in.sim dataset generated by get_species_list() function
#' @param file.in name including path of the file to save
plot_traits_pca12 <- function(traits, traits_demo, sp.in.sim, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Compile the traits data
  data_traits = traits %>% 
    left_join(traits_demo, by = "species") %>%
    filter(species %in% sp.in.sim) %>%
    drop_na() 
  
  # Make the PCA
  pca <- prcomp((data_traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # Extract data for individuals
  data.ind = data.frame(species = data_traits$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2])  %>%
    filter(species %in% sp.in.sim) %>%
    mutate(species = gsub("\\_", "\\ ", species))
  
  # Range of pca1 and pca2
  range.pca1 = diff(range(data.ind$pca1))
  range.pca2 = diff(range(data.ind$pca2))
  
  # Extract data for variables
  data.var = data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2]) %>%
    mutate(var = gsub("\\.", "\\ ", var), 
           pca1 = pca1*(max(abs(data.ind$pca1))/max(abs(pca1)))*0.9, 
           pca2 = pca2*(max(abs(data.ind$pca2))/max(abs(pca2)))*0.9, 
           pca1.txt = pca1*1.1, pca2.txt = pca2*1.2) 
  
  # Space fraction between axis and text
  space.frac = 0.025
  
  # Range of x and y axis for plotting
  range.x = range(data.var$pca1.txt)*c(1.4, 1.5)
  range.y = range(data.var$pca2.txt)*c(1.3, 1.1)
  
  # Make the plot
  plot.out = data.ind %>%
    # Add coordinates for species along x axis
    left_join((data.ind %>%
                 arrange(pca1) %>%
                 mutate(pca1_x = seq(from = range.x[1]+0.1, to = range.x[2]-0.1, 
                                     length.out = dim(.)[1]), 
                        pca1_y = range.y[2] + space.frac*diff(range.y)) %>%
                 dplyr::select(species, pca1_x, pca1_y)), by = "species") %>%
    # Add coordinates for species along y axis
    left_join((data.ind %>%
                 arrange(pca2) %>%
                 mutate(pca2_y = seq(from = range.y[1]+0.1, to = range.y[2]-0.1, 
                                     length.out = dim(.)[1]), 
                        pca2_x = range.x[2] + space.frac*diff(range.x)) %>%
                 dplyr::select(species, pca2_x, pca2_y)), by = "species") %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_segment(data = data.var, aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.3, "cm"))) + 
    geom_point(size = 2, shape = 21, fill = "#22333B", color = "black", alpha = 0.5) +
    geom_text(data = data.var, aes(label = var, x = pca1.txt, y = pca2.txt)) +
    # Frame of the graph
    geom_rect(aes(xmin = range.x[1], xmax = range.x[2], ymin = range.y[1], ymax = range.y[2]), 
              color = "black", fill = NA) +
    # Zero horizontal and vertical lines
    geom_segment(x = range.x[1], xend = range.x[2], y=0, yend=0, linetype = "dashed") +
    geom_segment(x = 0, xend = 0, y = range.y[1], yend = range.y[2], linetype = "dashed") +
    # Segment to connect points to species name
    geom_segment(aes(xend = pca1_x), yend=range.y[2], linetype = "dotted", 
                 inherit.aes = TRUE, color = "#ADB5BD", size = 0.3) +
    geom_segment(aes(yend = pca2_y), xend=range.x[2], linetype = "dotted", 
                 inherit.aes = TRUE, color = "#ADB5BD", size = 0.3) +
    # Axis label
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)",
                "\nHigh growth <--> High survival")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)",
                "\nHigh recruitment <--> Low recruitment")) +
    # Text axis
    # -- x axis
    geom_text(data = (data.frame(pca1 = c(-10:10), pca2 = range.y[1] - space.frac*diff(range.y)) %>%
                        filter(pca1 > range.x[1] & pca1 < range.x[2])), 
              aes(label = pca1), color = "#6C757D", size = 3) +
    # -- y axis
    geom_text(data = (data.frame(pca1 = range.x[1] - space.frac*diff(range.x), pca2 = c(-10:10)) %>%
                        filter(pca2 > range.y[1] & pca2 < range.y[2])), 
              aes(label = pca2), color = "#6C757D", size = 3) +
    # text of the x top axis
    geom_text(angle = 90, size = 2, 
              aes(x = pca1_x, y = pca1_y, label = species), color = "#6C757D", 
              inherit.aes = TRUE, vjust = 0, hjust = 0, fontface = "italic") + 
    # text of the y top axis
    geom_text(size = 2, 
              aes(x = pca2_x, y = pca2_y, label = species), color = "#6C757D", 
              inherit.aes = TRUE, vjust = 0, hjust = 0, fontface = "italic") + 
    # Theme
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          axis.title = element_blank()) + 
    # Limits of the plot
    xlim(c(1.2, 1.5)*range.x) + ylim(c(1.2, 1.5)*range.y) + 
    # X title manual
    geom_text(x = mean(range.x), y = range.y[1] - 3*space.frac*diff(range.y), 
              label = paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)"), 
              size = 3) + 
    geom_text(x = mean(range.x), y = range.y[1] - 5*space.frac*diff(range.y), 
              label = "High growth <--> High survival", fontface = "bold", size = 3) + 
    # Y title manual
    geom_text(y = mean(range.y), x = range.x[1] - 5*space.frac*diff(range.x), 
              label = paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)"), 
              size = 3, angle = 90) + 
    geom_text(y = mean(range.y), x = range.x[1] - 3*space.frac*diff(range.x), angle = 90,
              label = "High recruitment <--> Low recruitment", fontface = "bold", size = 3)
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 12, height = 12, 
         units = "cm", dpi = 600, bg = "white")
  
  
  # return the name of all the plots made
  return(file.in)
}


#' Plot the forest structure and composition along the climatic gradient
#' @param NFI_climate Climate data for each NFI plot
#' @param NFI_data Individual tree data for each NFI plot
#' @param pca_demo_per_species coordinates of each species on trait pca axes
#' @param nclim Number of categorical climate
#' @param dir.out directory whee to save the generated plots
plot_str_compo_climate = function(NFI_climate, NFI_data, pca_demo_per_species, 
                                  nclim, dir.out){
  
  
  # Prepare a climate dataframe 
  # -- Get mean climate per plot
  meanclim = NFI_climate %>% dplyr::select(plotcode)
  meanclim$sgdd = rowMeans(NFI_climate %>% dplyr::select(grep("sgdd", colnames(.))))  
  meanclim$wai = rowMeans(NFI_climate %>% dplyr::select(grep("wai", colnames(.))))
  # -- Make a pca
  pca <- prcomp((meanclim %>% dplyr::select(-plotcode)), 
                center = T, scale = T)
  # -- Add pca coordinate to the dataframe
  meanclim$pca1 = get_pca_ind(pca)[[1]][, 1]
  # -- Group plots per climate
  meanclim = meanclim %>%
    mutate(climate = ntile(pca1, n = nclim)) %>%
    ungroup() %>% group_by(climate) %>%
    mutate(pca1_mean = mean(pca1)) %>% ungroup() %>%
    dplyr::select(plotcode, climate, pca1_mean)
  
  # Prepare dataset to plot species diversity and mean composition per climate
  data.sp = NFI_data %>%
    dplyr::select(plotcode, species, ba_ha) %>%
    group_by(plotcode, species) %>%
    summarize(ba_sp = sum(ba_ha, na.rm = TRUE)) %>%
    left_join((pca_demo_per_species %>% mutate(species = gsub("\\_", "\\ ", species))), 
              by = "species") %>%
    ungroup() %>% group_by(plotcode) %>%
    summarize(nsp = n(), 
              cwm1 = weighted.mean(pca1, w = ba_sp, na.rm = TRUE), 
              cwm2 = weighted.mean(pca2, w = ba_sp, na.rm = TRUE))
  
  # Prepare dataset with information on forest structure
  data.str = NFI_data %>%
    group_by(plotcode) %>%
    mutate(ba_ha = ba_ha/1000000) %>%
    summarize(dbh_mean = weighted.mean(dbh, w = Nha, na.rm = TRUE), 
              stock = sum(ba_ha, na.rm = TRUE), 
              density = sum(Nha, na.rm = TRUE))
  
  # Assemble data on composition and on structure
  data = data.sp %>%
    dplyr::select(-nsp) %>%
    left_join(data.str, by = "plotcode")
  
  # Dataframe containing variable names and associated caption
  data.caption = data.frame(
    var.name = c("cwm1", "cwm2", "dbh_mean", "stock", "density"), 
    caption = c("CWM1 \nHigh growth <--> High survival", 
                "CWM2 \nHigh <--> Low recruitment", 
                "Mean dbh\n(mm)", 
                "Total basal area\n(ba/ha)", 
                "Tree density\n(tree/ha)"))
  
  # Plot diversity per climate
  plot.diversity = meanclim %>%
    left_join((data.sp %>% dplyr::select(plotcode, nsp)), by = "plotcode") %>%
    mutate(N_species = paste0(nsp, " species")) %>%
    group_by(climate, N_species) %>%
    summarize(n = n()) %>%
    ungroup() %>% group_by(climate) %>%
    mutate(rate = n/sum(n)) %>%
    mutate(climate = factor(climate, levels = as.character(c(1:10)))) %>%
    ggplot(aes(x = climate, y = rate, fill = N_species)) +
    geom_bar(stat = "identity", color = "black")
  
  #¨File name for diversity plot
  file.div = paste0(dir.out, "/fig_diversity.jpg")
  # Create directory if needed
  create_dir_if_needed(file.div)
  # Save the plot
  ggsave(file.div, plot.diversity, width = 14, height = 10, units = "cm", 
         dpi = 600, bg = "white")
  # Initialize output with the file name
  out = c(file.div)
  
  # Loop on all plots to make
  for(i in 1:dim(data.caption)[1]){
    
    # Name a column containing variable i
    data$var.i = as.numeric(unlist(data[, data.caption$var.name[i]]))
    
    # Make the plot
    plot.i = data %>%
      left_join(meanclim, by = "plotcode") %>%
      mutate(climate = factor(climate, levels = as.character(c(1:10)))) %>%
      ggplot(aes(x = var.i, y = climate, fill = climate)) + 
      geom_density_ridges(alpha =  0.7) + 
      scale_fill_manual(values = colorRampPalette(c("orange", "blue"))(10)) +
      xlab(data.caption$caption[i]) + 
      ylab("Climate\n(Hot-dry to cold-wet)") +
      coord_flip() +
      theme(legend.position = "none", 
            panel.grid = element_blank(), 
            panel.background = element_rect(fill = "white", color = "black")) 
    
    # Name of file i
    file.i = paste0(dir.out, "/fig_", data.caption$var.name[i], ".jpg")
    
    # Save the plot
    ggsave(file.i, plot.i, width = 14, height = 10, units = "cm", 
           dpi = 600, bg = "white")
    
    # Add to the output list
    out = c(out, file.i)
  }
  
  # Return output
  return(out)
  
}

