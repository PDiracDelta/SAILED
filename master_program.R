clean.session <- function(){
  # install remotes and download JGmisc from github
  # install.packages(remotes)
  # library(remotes)
  # remotes::install_github("jasongraf1/JGmisc")
  library(JGmisc)
  detachAllPackages() # detach packages
  library(rmarkdown)
}

### knit all the notebooks specified below
#st=Sys.time()

#clean.session(); rm(list = ls()[!ls()=='clean.session']); render("intro.Rmd") 
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("datadriven_unit.Rmd")
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("datadriven_summarization.Rmd")
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("datadriven_normalization.Rmd")
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("datadriven_DEA.Rmd")
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("CONSTANd_vs_medianSweeping.Rmd") 

clean.session(); rm(list = ls()[!ls()=='clean.session']); render("modelbased_unit.Rmd") 
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("modelbased_summarization.Rmd")
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("modelbased_normalization.Rmd")
clean.session(); rm(list = ls()[!ls()=='clean.session']); render("modelbased_DEA.Rmd")

clean.session(); rm(list = ls()[!ls()=='clean.session']); render("compare_defaults.Rmd")

#ed=Sys.time()
#ed-st # running time