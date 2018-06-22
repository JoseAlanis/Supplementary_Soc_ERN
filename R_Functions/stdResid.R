stdResid <- function(data, model, return.data = T, bound = 2.5, plot = F, show.bound = F, show.loess = F, ... ) {
  
  # Save model name for later
  main <- deparse(substitute(model))
  
  data <- data %>% 
    # Compute standaridised residuals
    mutate(st_R = as.numeric(scale(residuals( model, 'pearson' ))),
           # Code cases as outliers if residual > bound
           Outlier = ifelse(abs(st_R) > bound, 1, 0))
  
  # summary of outlying cases
  print(paste(round((sum(data$Outlier == 1) / sum(data$Outlier == 0))*100, digits = 2), 
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
