library(ggplot2)
library(stringi)
library(grid)
library(gridExtra)
library(dendextend)
library(kableExtra)
library(psych)
library(tidyverse)
library(ggplotify)
library(ggpubr)

op <- par(no.readonly = TRUE)
source('manuscript_functions.R')
out_path <-'G:/My Drive/Isobaric labeling strategies/analyses output/' 
referenceCondition <- '0.5'
reportContrast <- '0.125'
# specify colours corresponding to biological conditions
condition.color <- tribble(
  ~Condition, ~Color,
  "0.125", 'black',
  "0.5", 'blue',
  "0.667", 'green',
  "1", 'red' )
data.list <- readRDS('input_data.rds')
dat.l <- data.list$dat.l 
spiked.proteins <- dat.l %>% distinct(Protein) %>% filter(stri_detect(Protein, fixed='ups')) %>% pull %>% as.character
sample.info <- get_sample_info(dat.l, condition.color)
channelNames <- remove_factors(unique(sample.info$Channel))
keep_objects=c('referenceCondition', 'reportContrast', 'condition.color', 'spiked.proteins','sample.info', 'channelNames', 'out_path','op')

### Conclusion 1 (main) ###
# M Class compare default LMM1 with default data-driven approach: scatter plots, violin plots, PCA plots: results very similar! 

load('compare_defaults_outdata.rda')
# recordPlot to save R-base graphs
dat <- list(dat.nonnorm.summ.l$data_driven
            ,dat.nonnorm.summ.l$model_based
            ,dat.norm.summ.l$model_based 
            ,dat.norm.summ.l$data_driven)
names(dat) <- c('Non-normalized,Data-driven','Non-normalized,Model-based',
                'Normalized,Data-driven', 'Normalized,Model-based')
# dat <- list(dat.norm.summ.l$model_based 
#             ,dat.norm.summ.l$data_driven)
# names(dat) <- c('Normalized,Data-driven', 'Normalized,Model-based')
p0 <- run_effect_plot_rep(dat)
par(mfrow=c(2,2))
p1 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$data_driven %>% select(-Protein), info=sample.info, 'Data-driven', show.legend = T,
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))
p2 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$model_based %>% select(-Protein), info=sample.info, 'Model-based',
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))
dea.cols <- colnames(dat.dea[[1]])
logFC.cols <- dea.cols[stri_detect_fixed(dea.cols, 'logFC')]
significance.cols <- dea.cols[stri_detect_fixed(dea.cols, 'q.mod')]
n.contrasts <- length(logFC.cols)
p3 <- scatterplot_rep(dat.dea, significance.cols, spiked.proteins, 
                      xlab.title='Data-driven (q-value)', ylab.title='Model-based (q-value)', reportContrast='0.125')
p4 <- scatterplot_rep(dat.dea, logFC.cols, spiked.proteins, 
                      xlab.title='Data-driven (logFC)', ylab.title='Model-based (logFC)', reportContrast='0.125')
p5 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, logFC.cols]), referenceCondition, reportContrast='0.125')

jpeg(paste0(out_path, "output_m1.jpeg"), width = 600, height = 500)
  ggarrange(p1,p2,p3,p4, labels="AUTO",ncol = 2, nrow = 2, heights=c(5,4))
dev.off()
jpeg(paste0(out_path, "output_s1.jpeg"), width = 600, height = 250)
  ggarrange(p0,p5, labels="AUTO",ncol = 2, nrow = 1)
dev.off()

### Conclusion 2 (main) ###
#N peptide-by-run is essential to remove run effect, though Hill and Oberg don't mention this

rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_normalization_outdata.rda')
variant.names=names(dat.norm.summ.w2)

jpeg(paste0(out_path, "output_m2.jpeg"), width = 500, height = 400)
par(mfrow=c(2,2))
pcaplot_rep(dat.nonnorm.summ.w2 %>% select(-'Protein'), info=sample.info, 'Raw', show.legend=T, par.list=list(oma=c(0,0,0,0), mar=c(2,2,1,1), mgp=c(1.2,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, paste('Variant',i,'normalization'))
}
dev.off()

dea.cols <- colnames(dat.dea[[1]])
logFC.cols <- dea.cols[stri_detect_fixed(dea.cols, 'logFC')]
significance.cols <- dea.cols[stri_detect_fixed(dea.cols, 'q.mod')]
n.contrasts <- length(logFC.cols)
names(dat.dea) <- c('Variant 1','Variant 2','Variant 3','Raw')
jpeg(paste0(out_path, "output_s2.jpeg"), width = 700, height = 350)
  scatterplot_rep2(dat.dea, significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()
jpeg(paste0(out_path, "output_s2bkp.jpeg"), width = 700, height = 350)
  scatterplot_rep2(dat.dea, logFC.cols, 'log2FC', spiked.proteins, reportContrast)
dev.off()


### Conclusion 3 (main) ###

### Conclusion 4 (main) ###