# Semi-partion R2
spR2 <- function(model) {
  anov_mod <- as.data.frame(model)
  
  rr <- ((anov_mod$NumDF / anov_mod$DenDF) * anov_mod[, grepl(names(anov_mod), pattern = 'F*value')] ) / 
    (1+((anov_mod$NumDF / anov_mod$DenDF) * anov_mod[, grepl(names(anov_mod), pattern = 'F*value')] ))
  
  print(rr)
  
}