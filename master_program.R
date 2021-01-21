library(pacman)
library(rmarkdown)

unload_pkgs <- function() p_unload(setdiff(p_loaded(),c("pacman","rmarkdown")), character.only = TRUE)

### process raw data

source('data_prep.R')

### knit all the notebooks specified below

rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("intro.Rmd")

rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("datadriven_unit.Rmd")
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("datadriven_summarization.Rmd") 
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("datadriven_normalization.Rmd") 
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("datadriven_DEA.Rmd") 

rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("modelbased_unit.Rmd")
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("modelbased_summarization.Rmd") 
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("modelbased_normalization.Rmd") 
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("modelbased_DEA.Rmd") 

rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("compare_defaults.Rmd") 
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("CONSTANd_vs_medianSweeping.Rmd") 
rm(list=setdiff(ls(),c('unload_pkgs'))); unload_pkgs(); render("datadriven_unit_rawratio.Rmd")
