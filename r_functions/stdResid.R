stdResid <- function(data, model, return.data = T, bound = 2.5, plot = F, show.bound = F, show.loess = F, ... ) {
  
  # Save model name for later
  main <- class(model)
  
  if (length( grep(paste(c('glmerMod', 'glm'), collapse="|"), main) ) > 0 ) {
    data <- data %>% 
      # Compute person residuals
      mutate(st_R = as.numeric(resid(model, type = 'pearson')),
             # Code cases as outliers if residual > bound
             Outlier = ifelse(abs(st_R) > bound, 1, 0))
  } else {
    data <- data %>% 
      # Compute scaled residuals
      mutate(st_R = as.numeric(scale(resid(model))),
             # Code cases as outliers if residual > bound
             Outlier = ifelse(abs(st_R) > bound, 1, 0))
  }
  
  
  # summary of outlying cases
  print(paste(round((sum(data$Outlier == 1) / nrow(data))*100, digits = 2), 
              '% are outliers',  sep = ' '))
  
  # Plot residuals against fitted values
  if (plot == T) {
    plot(fitted(model), data$st_R, 
         ylim = c(min(data$st_R)-.5, max(data$st_R)+.5), ...)
    
    # Draw red lines to mark boundaries
    if (show.bound == T) {
      abline(h=c(-bound, bound), col = 'red')
      abline(h=0, col='black', lty = 2)
    }
    
    # Draw loess fit (experimental, don't run)
    # if (show.loess ==  T) {
    #   lloess <- seq(min(fitted(model)), max(fitted(model)), length = nrow(model.frame(model))  )
    #   lines(lloess, 
    #         predict(loess(residuals(model) ~ fitted(model)), 
    #                 newdata = lloess),
    #         col = 'blue', lwd = 2)
    # }
    
  }
  
  # Return new data including residuals and outlier identifier 
  if (return.data == T) {
    return(data)
  }
  
}