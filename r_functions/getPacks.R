getPacks <- function(packages) {
  
  # list of packages missing
  not_installed <- which(!packages %in% installed.packages()[, 'Package'])
  
  # check wich packages are not intalled and install them
  if (length(not_installed)) {
    install.packages(packages[not_installed], dependencies = T)
  }
  
  # load all packages
  sapply(packages, require, character.only = T)
}
