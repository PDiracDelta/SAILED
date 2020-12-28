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
library(cowplot)
library(viridis)

op <- par(no.readonly = TRUE)
source('manuscript_functions.R')
out_path <-'G:/My Drive/Isobaric labeling strategies/analyses output/manuscript_figures_tables/' 
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

# set standard plot width, height and resolution
# PCA plot
#1x3
width.plot1.a <- 1000
height.plot1.a <- 500
res.plot1.a <- 150

#2x2
width.plot1.b <- 1200
height.plot1.b <- 900
res.plot1.b <- 180

#2x3
width.plot1.c <- 1400
height.plot1.c <- 1000
res.plot1.c <- 200

# RE & violin plot
#1x2
width.plot2.a <- 1000
height.plot2.a <- 500
res.plot2.a <- 150

#1x1
width.plot2.b <- 800
height.plot2.b <- 500
res.plot2.b <- 120

# scatter plots
#3x3
width.plot3.a <- 1000
height.plot3.a <- 600
res.plot3.a <- 120

#4x4
width.plot3.b <- 1100
height.plot3.b <- 700
res.plot3.b <- 110

#6x6
width.plot3.c <- 1400
height.plot3.c <- 1000
res.plot3.c <- 150

keep_objects=c('referenceCondition', 'reportContrast', 'condition.color', 'spiked.proteins','sample.info', 'channelNames', 'out_path','op','keep_objects',
               'width.plot1.a','height.plot1.a','res.plot1.a', 
               'width.plot1.b','height.plot1.b','res.plot1.b',
               'width.plot1.c','height.plot1.c','res.plot1.c',
               'width.plot2.a','height.plot2.a','res.plot2.a', 
               'width.plot2.b','height.plot2.b','res.plot2.b',
               'width.plot3.a','height.plot3.a','res.plot3.a', 
               'width.plot3.b','height.plot3.b','res.plot3.b',
               'width.plot3.c','height.plot3.c','res.plot3.c')

### compare defaults ###
load('compare_defaults_outdata.rda')
# dat <- list(dat.nonnorm.summ.l$data_driven
#             ,dat.nonnorm.summ.l$model_based
#             ,dat.norm.summ.l$model_based 
#             ,dat.norm.summ.l$data_driven)
# names(dat) <- c('Non-normalized,Data-driven','Non-normalized,Model-based',
#                 'Normalized,Data-driven', 'Normalized,Model-based')
dat <- list(dat.norm.summ.l$model_based
            ,dat.norm.summ.l$data_driven)
names(dat) <- c('Data-driven', 'Model-based')
p0 <- run_effect_plot_rep(dat)
p1 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$data_driven %>% select(-Protein), info=sample.info, 'Data-driven', show.legend = T,legend.cex=0.9,legend.run.xy=c(-0.9,1.5),legend.cond.xy=c(-0.9,-1.8),legend.textwidth =1.5, 
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))
p2 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$model_based %>% select(-Protein), info=sample.info, 'Model-based', 
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))

dea.info <- get_dea_info(dat.dea)
p3 <- scatterplot_rep(dat.dea, dea.info$significance.cols, spiked.proteins, 
                      xlab.title='Data-driven (q-value)', ylab.title='Model-based (q-value)', reportContrast='0.125')
p4 <- scatterplot_rep(dat.dea, dea.info$logFC.cols, spiked.proteins, 
                      xlab.title='Data-driven (logFC)', ylab.title='Model-based (logFC)', reportContrast='0.125', legend_position = c(0.2, 0.75))
p5 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')

jpeg(paste0(out_path, "comparedefaults1.jpeg"), width = 1400, height = 1400, res=160)
ggarrange(p1,p2,p3,p4, labels="AUTO",ncol = 2, nrow = 2, heights=c(5,4))
dev.off()

jpeg(paste0(out_path, "comparedefaults2.jpeg"), width = width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p0,p5, labels="AUTO",ncol = 2, nrow = 1)
dev.off()

### raw ratios as good as normalized ratios ###
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_outdata.rda')
variant.names=names(dat.norm.summ.w2)

p1 <- as.ggplot(~pcaplot_rep(dat.nonnorm.summ.w2$Ratio %>% select(-Protein), info=sample.info, 'Raw ratios', show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.2,1.3),legend.cond.xy=c(-2.2,-0.5),legend.textwidth=3.5,
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))
p2 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$Ratio %>% select(-Protein), info=sample.info, 'Median-swept ratios',
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))
dat <- list(dat.nonnorm.summ.l$Ratio, dat.norm.summ.l$Ratio)
names(dat) <- c('Raw ratios', 'Median-swept ratios')
p3 <- run_effect_plot_rep(dat)
#width = 1400, height = 700, res=150
jpeg(paste0(out_path, "rawratio1.jpeg"), width = width.plot1.a , height = height.plot1.a+50 , res=res.plot1.a-20)
ggarrange(p1,p2,p3, labels="AUTO",ncol = 3, nrow = 1)
dev.off()

rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_rawratio_outdata.rda')
dat.dea[!(names(dat.dea) %in% c('log2Intensity', 'Ratio'))] <- NULL
dea.info <- get_dea_info(dat.dea)
dat.dea[!(names(dat.dea) %in% c('log2Intensity', 'Ratio'))] <- NULL
names(dat.dea)[2] <- 'Raw_ratios'

jpeg(paste0(out_path, "rawratio2.jpeg"), width = 1200, height = 600, res=110)
p1 <- as_grob(~scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast))
p2 <- as_grob(~scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast))
ggarrange(p1, p2)
dev.off()

### constand_vs_mediansweeping ###
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('constand_vs_mediansweeping_outdata.rda')
variant.names=names(dat.norm.summ.w2)
dea.info <- get_dea_info(dat.dea)

p1 <- as_grob(~scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast))
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
#width = 800, height = 500, res=100
jpeg(paste0(out_path, "constand_vs_mediansweeping1.jpeg"), width = 1000, height = 500, res=110)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1, widths=c(65,35))
dev.off()

#############################################
# standard plots                           #
# (same set of figures for each component)#
##########################################

### Unit ###

# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_unit_outdata.rda')
dat.norm.summ.w2[c("IntensityFix","RatioFix")]<- NULL
dat.dea[c("IntensityFix","RatioFix")]<- NULL
variant.names=names(dat.norm.summ.w2)
dea.info <- get_dea_info(dat.dea)

# EXTRA (PCA PLOT) of unnormalized data (doesn't matter if based on model-based or data-driven data)
jpeg(paste0(out_path, "unit1.jpeg"), width = width.plot1.a, height = height.plot1.a, res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(variant.names)){
  pca.scale=F; if (variant.names[i]=='Intensity') pca.scale=T
  if (i==1){
  pcaplot_rep(dat.nonnorm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i], show.legend = T,legend.cex=0.9,legend.run.xy=c(5,27),legend.cond.xy=c(5,-10),legend.textwidth=1.2)} else {
  pcaplot_rep(dat.nonnorm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i])}
}
dev.off()

# F1 (PCA PLOT)
jpeg(paste0(out_path, "mb_unit1.jpeg"), width = width.plot1.a, height = height.plot1.a, res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
  pca.scale=F; if (variant.names[i]=='Intensity') pca.scale=T
  if (i==1){
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i], show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.5,1.6),legend.cond.xy=c(-1.5,-1.8),legend.textwidth=1)} else {
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i])}
}
dev.off()

# F2 (REP + VIOLIN PLOT)
dat.norm.summ.l[c("IntensityFix","RatioFix")]<- NULL
p1 <- run_effect_plot_rep(dat.norm.summ.l)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
jpeg(paste0(out_path, "mb_unit2.jpeg"), width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
dat.dea[c('IntensityFix', 'RatioFix')] <- NULL
jpeg(paste0(out_path, "mb_unit3.jpeg"), width.plot3.a, height = height.plot3.a, res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "mb_unit4.jpeg"), width.plot3.a, height = height.plot3.a, res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_outdata.rda')
names(dat.norm.summ.w2)
variant.names=names(dat.norm.summ.w2)
dat.norm.summ.w2$Intensity_lateLog2 <- NULL
dea.info <- get_dea_info(dat.dea)

# F1 (PCA PLOT)
jpeg(paste0(out_path, "dd_unit1.jpeg"), width.plot1.a, height = height.plot1.a, res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
  #if (variant.names[i]=='Intensity') pca.scale=T
  if (i==1){
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i], show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.5,1.6),legend.cond.xy=c(-1.5,-1.65),legend.textwidth=1)} else {
      pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])
    }
}
dev.off()

# F2 (REP + VIOLIN PLOT)
p1 <- run_effect_plot_rep(dat.norm.summ.l)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
jpeg(paste0(out_path, "dd_unit2.jpeg"), width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
jpeg(paste0(out_path, "dd_unit3.jpeg"), width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "dd_unit4.jpeg"), width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

### Normalization ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_normalization_outdata.rda')
variant.names=names(dat.norm.summ.w2)
# F1 (PCA PLOT)
jpeg(paste0(out_path, "mb_norm1.jpeg"), width = width.plot1.b, height = height.plot1.b, res=res.plot1.b)
par(list(mfrow=c(2,2),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
pcaplot_rep(dat.nonnorm.summ.w2%>% select(-'Protein'), info=sample.info, title='Raw',show.legend = T,legend.cex=0.8,legend.run.xy=c(1,28),legend.cond.xy=c(1,-2))
for (i in 1:length(dat.norm.summ.w2)){
  pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])
}
dev.off()

# F2 (REP + VIOLIN PLOT)
dat <- list(dat.nonnorm.summ.l, dat.norm.summ.l$LMM1, dat.norm.summ.l$LMM2, dat.norm.summ.l$LMM3) 
names(dat) <- c('Raw', 'LMM1', 'LMM2', 'LMM3')
p1 <- run_effect_plot_rep(dat)
dea.info <- get_dea_info(dat.dea)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
jpeg(paste0(out_path, "mb_norm2.jpeg"), width = width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
# nams <- names(dat.dea)
# dat.dea <- dat.dea[c('raw',nams[!(nams=='raw')])]
# names(dat.dea) <- c('Raw', paste('Variant',1:length(variant.names)))

jpeg(paste0(out_path, "mb_norm3.jpeg"), width = width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "mb_norm4.jpeg"), width = width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_normalization_outdata.rda')
names(dat.norm.summ.w2)[names(dat.norm.summ.w2) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
variant.names=names(dat.norm.summ.w2)

# F1 (PCA PLOT)
jpeg(paste0(out_path, "dd_norm1.jpeg"), width = width.plot1.c , height = height.plot1.c, res=res.plot1.c)
par(list(mfrow=c(2,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
pcaplot_rep(dat.nonnorm.summ.w2%>% select(-'Protein'), info=sample.info, title='Raw',show.legend = T,legend.cex=0.8,legend.run.xy=c(1,28),legend.cond.xy=c(1,-2))
for (i in 1:length(dat.norm.summ.w2)){
  pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])
}
dev.off()

# F2 (REP + VIOLIN PLOT)
dat <- list(dat.nonnorm.summ.l, dat.norm.summ.l$medianSweeping, dat.norm.summ.l$CONSTANd, dat.norm.summ.l$NOMAD,
            dat.norm.summ.l$quantile, dat.norm.summ.l$quantileByCondition)
names(dat) <- c('Raw', names(dat.norm.summ.l))
names(dat)[names(dat) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
p1 <- run_effect_plot_rep(dat)
dea.info <- get_dea_info(dat.dea)
names(dat.dea)[names(dat.dea) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125', rot.xaxis.lab=30)
jpeg(paste0(out_path, "dd_norm2.jpeg"), width = width.plot2.a , height = height.plot2.a , res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
names(dat.dea)[names(dat.dea) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
names(dat.dea)[names(dat.dea)=='raw'] <- 'Raw'
jpeg(paste0(out_path, "dd_norm3.jpeg"), width = width.plot3.c , height = height.plot3.c , res=res.plot3.c)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "dd_norm4.jpeg"), width = width.plot3.c , height = height.plot3.c , res=res.plot3.c)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

### Summarization ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_summarization_outdata.rda')
variant.names=names(dat.norm.summ.w2)

# F1 (PCA PLOT)
jpeg(paste0(out_path, "mb_summ1.jpeg"), width = width.plot1.a , height = height.plot1.a , res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
  if (i==1){pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i], show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.5,1.6),legend.cond.xy=c(-1.5,-1.8),legend.textwidth=1)} else {
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])}
}
dev.off()

# F2 (REP + VIOLIN PLOT)
p1 <- run_effect_plot_rep(dat.norm.summ.l)
dea.info <- get_dea_info(dat.dea)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
jpeg(paste0(out_path, "mb_summ2.jpeg"), width = width.plot2.a , height = height.plot2.a , res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
jpeg(paste0(out_path, "mb_summ3.jpeg"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "mb_summ4.jpeg"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_summarization_outdata.rda')
variant.names=names(dat.norm.summ.w2)
# F1 (PCA PLOT)
jpeg(paste0(out_path, "dd_summ1.jpeg"), width = width.plot1.a , height = height.plot1.a , res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
  if (i==1){ 
   pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i], show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.5,1.6),legend.cond.xy=c(-1.5,-1.6),legend.textwidth=1)} else {
   pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])}
}
dev.off()

# F2 (REP + VIOLIN PLOT)
p1 <- run_effect_plot_rep(dat.norm.summ.l)
dea.info <- get_dea_info(dat.dea)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
jpeg(paste0(out_path, "dd_summ2.jpeg"), width = width.plot2.a , height = height.plot2.a , res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
jpeg(paste0(out_path, "dd_summ3.jpeg"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "dd_summ4.jpeg"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

### DEA ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_DEA_outdata.rda')
variant.names=names(dat.dea)
dea.info <- get_dea_info(dat.dea)

# F2 (REP + VIOLIN PLOT)
jpeg(paste0(out_path, "mb_DEA2.jpeg"), width = width.plot2.b , height = height.plot2.b , res=res.plot2.b)
violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
dev.off()

# F3 (2xSPLOM)
jpeg(paste0(out_path, "mb_DEA3.jpeg"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "mb_DEA4.jpeg"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_DEA_outdata.rda')
variant.names=names(dat.dea)
dea.info <- get_dea_info(dat.dea)
# F2 (REP + VIOLIN PLOT)
jpeg(paste0(out_path, "dd_DEA2.jpeg"), width = width.plot2.b , height = height.plot2.b , res=res.plot2.b)
violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
dev.off()

# F3 (2xSPLOM)
jpeg(paste0(out_path, "dd_DEA3.jpeg"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

jpeg(paste0(out_path, "dd_DEA4.jpeg"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()
