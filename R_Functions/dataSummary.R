dataSummary <- function(x) {
  m <- median(x)
  ymin <- as.numeric(quantile(x, .25))
  ymax <- as.numeric(quantile(x, .75))
  return(c(y=m,ymin=ymin,ymax=ymax))
}
