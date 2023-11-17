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