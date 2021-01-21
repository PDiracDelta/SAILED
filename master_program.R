library(pacman)
library(rmarkdown)

unload_pkgs <- function() p_unload(setdiff(p_loaded(),c("pacman","rmarkdown")), character.only = TRUE)

knit_notebook <- function(name, notebook.params=NULL){
  render(input = paste0(name,".Rmd"), 
         output_file = file.path(dirname(paste0(name,".Rmd")), paste0(name, '_', data_name)),
         params=notebook.params)
}

# input data set name used as a suffix in notebook output files 
# (this variable is used in data_prep.R and all notebook .Rmd)
data_name='msstatstmt'

### process raw data
source('data_prep.R')

### knit all the notebooks specified below
rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook('intro')

#### data-driven notebooks ####
rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("datadriven_unit")

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("datadriven_summarization")

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("datadriven_normalization") 

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("datadriven_DEA") 

#### model-based notebooks ####
rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("modelbased_unit")

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("modelbased_summarization") 

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("modelbased_normalization") 

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("modelbased_DEA") 

#### other notebooks ####
rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("compare_defaults") 

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("CONSTANd_vs_medianSweeping") 

rm(list=setdiff(ls(),c('unload_pkgs', 'data_name', 'knit_notebook'))); unload_pkgs(); 
knit_notebook("datadriven_unit_rawratio")
