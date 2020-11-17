### knit all the notebooks specified below
library(rmarkdown)
st=Sys.time()
#render("intro.Rmd") 
render("datadriven_unit.Rmd", envir=new.env())
render("datadriven_summarization.Rmd", envir=new.env())
render("datadriven_normalization.Rmd", envir=new.env())
render("datadriven_DEA.Rmd", envir=new.env())
render("CONSTANd_vs_medianSweeping.Rmd", envir=new.env()) 

render("modelbased_unit.Rmd", envir=new.env()) 
render("modelbased_summarization.Rmd", envir=new.env())
render("modelbased_normalization.Rmd", envir=new.env())
render("modelbased_DEA.Rmd", envir=new.env())
ed=Sys.time()
ed-st # running time