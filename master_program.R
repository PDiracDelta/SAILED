clean.session <- function(){
  rm(list = ls()) # clean workspace
  # install remotes and download JGmisc from github
  # install.packages(remotes)
  # library(remotes)
  # remotes::install_github("jasongraf1/JGmisc")
  library(JGmisc)
  detachAllPackages() # detach packages
  library(rmarkdown)
}

### knit all notebooks specified below

clean.session(); render("intro.Rmd") 
# clean.session(); render("datadriven_unit.Rmd") 
# clean.session(); render("datadriven_summarization.Rmd") 
# clean.session(); render("datadriven_normalization.Rmd") 
# clean.session(); render("datadriven_DEA.Rmd")
# clean.session(); render("CONSTANd_vs_medianSweeping.Rmd") 

clean.session(); render("modelbased_unit.Rmd") 
# clean.session(); render("modelbased_summarization.Rmd")
# clean.session(); render("modelbased_normalization.Rmd")
# clean.session(); render("modelbased_DEA.Rmd")
