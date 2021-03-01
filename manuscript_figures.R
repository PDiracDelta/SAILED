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
data.list <- readRDS('input_data_msstatstmt.rds')
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
res.plot1.b <- 150

#2x3
width.plot1.c <- 1400
height.plot1.c <- 1000
res.plot1.c <- 200

# RE & violin plot
#1x2
width.plot2.a <- 900
height.plot2.a <- 350
res.plot2.a <- 120

#1x1
width.plot2.b <- 600
height.plot2.b <- 300
res.plot2.b <- 120

# scatter plots
#3x3
width.plot3.a <- 900
height.plot3.a <- 500
res.plot3.a <- 120

#4x4
width.plot3.b <- 1000
height.plot3.b <- 500
res.plot3.b <- 120

#6x6
width.plot3.c <- 1400
height.plot3.c <- 1000
res.plot3.c <- 120

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
load('compare_defaults_outdata_msstatstmt.rda')
# dat <- list(dat.nonnorm.summ.l$data_driven
#             ,dat.nonnorm.summ.l$model_based
#             ,dat.norm.summ.l$model_based 
#             ,dat.norm.summ.l$data_driven)
# names(dat) <- c('Non-normalized,Data-driven','Non-normalized,Model-based',
#                 'Normalized,Data-driven', 'Normalized,Model-based')
dat <- list(dat.norm.summ.l$model_based
            ,dat.norm.summ.l$data_driven)
names(dat) <- c('data_driven', 'model_based')
p0 <- run_effect_plot_rep(dat)
p1 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$data_driven %>% select(-Protein), info=sample.info, 'data_driven', show.legend = T,legend.cex=0.9,legend.run.xy=c(-0.9,1.5),legend.cond.xy=c(-0.9,-1.8),legend.textwidth =1.5, 
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))
p2 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$model_based %>% select(-Protein), info=sample.info, 'model_based', 
                             par.list=list(oma=c(2,1,1,1), mar=c(2,2,1,1), mgp=c(1.5,0.5,0))))

dea.info <- get_dea_info(dat.dea)
p3 <- scatterplot_rep(dat.dea, dea.info$significance.cols, spiked.proteins, 
                      xlab.title='data_driven (q-value)', ylab.title='model_based (q-value)', reportContrast='0.125')
p4 <- scatterplot_rep(dat.dea, dea.info$logFC.cols, spiked.proteins, 
                      xlab.title='data_driven (logFC)', ylab.title='model_based (logFC)', reportContrast='0.125', legend_position = c(0.2, 0.75))
p5 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')

png(paste0(out_path, "comparedefaults1.png"), width = 1400, height = 1400, res=160)
ggarrange(p1,p2,p3,p4, labels="AUTO",ncol = 2, nrow = 2, heights=c(5,4))
dev.off()

png(paste0(out_path, "comparedefaults2.png"), width = width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p0,p5, labels="AUTO",ncol = 2, nrow = 1)
dev.off()

### raw ratios as good as normalized ratios ###
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_outdata_msstatstmt.rda')
variant.names=names(dat.norm.summ.w2)

#c(2,1,1,1)
p1 <- as.ggplot(~pcaplot_rep(dat.nonnorm.summ.w2$ratio %>% select(-Protein), info=sample.info, 'raw ratio', show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.2,1.3),legend.cond.xy=c(-2.2,-0.5),legend.textwidth=3.5,
                             par.list=list(oma=c(1,1,1,1), mar=c(2,2,1,1), mgp=c(1,0.3,0)), cex.lab=0.8, cex.axis=0.8))
p2 <- as.ggplot(~pcaplot_rep(dat.norm.summ.w2$ratio %>% select(-Protein), info=sample.info, 'median-swept ratio',
                             par.list=list(oma=c(1,1,1,1), mar=c(2,2,1,1), mgp=c(1,0.3,0)), cex.lab=0.8, cex.axis=0.8))
dat <- list(dat.nonnorm.summ.l$ratio, dat.norm.summ.l$ratio)
names(dat) <- c('raw ratio', 'median-swept ratio')
p3 <- run_effect_plot_rep(dat)
#width = 1400, height = 700, res=150
png(paste0(out_path, "rawratio1.png"), width = width.plot1.a+10 , height = height.plot1.a+50 , res=res.plot1.a-30)
ggarrange(p1,p2,p3, labels="AUTO",ncol = 3, nrow = 1)
dev.off()

rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_rawratio_outdata_msstatstmt.rda')
dat.dea[!(names(dat.dea) %in% c('log2_intensity', 'ratio'))] <- NULL
dea.info <- get_dea_info(dat.dea)
dat.dea[!(names(dat.dea) %in% c('log2_intensity', 'ratio'))] <- NULL
names(dat.dea)[2] <- 'raw_ratio'

png(paste0(out_path, "rawratio2.png"), width = 1200, height = 600, res=110)
p1 <- as_grob(~scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast))
p2 <- as_grob(~scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast))
ggarrange(p1, p2)
dev.off()

### constand_vs_mediansweeping ###
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('constand_vs_mediansweeping_outdata_msstatstmt.rda')
variant.names=names(dat.norm.summ.w2)
dea.info <- get_dea_info(dat.dea)

p1 <- as_grob(~scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast))
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
#width = 800, height = 500, res=100
png(paste0(out_path, "constand_vs_mediansweeping1.png"), width = 1000, height = 500, res=110)
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
load('modelbased_unit_outdata_msstatstmt.rda')

dat.norm.summ.w2[c("intensity_fix","ratio_fix")]<- NULL
dat.dea[c("intensity_fix","ratio_fix")]<- NULL
variant.names=names(dat.norm.summ.w2)
dea.info <- get_dea_info(dat.dea)

# EXTRA (PCA PLOT) of unnormalized data (doesn't matter if based on model-based or data-driven data)
png(paste0(out_path, "unit1.png"), width = width.plot1.a, height = height.plot1.a, res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(variant.names)){
  pca.scale=F; if (variant.names[i]=='intensity') pca.scale=T
  if (i==1){
  pcaplot_rep(dat.nonnorm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i], show.legend = T,legend.cex=0.9,legend.run.xy=c(5,27),legend.cond.xy=c(5,-10),legend.textwidth=1.2)} else {
  pcaplot_rep(dat.nonnorm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i])}
}
dev.off()

# F1 (PCA PLOT)
png(paste0(out_path, "mb_unit1.png"), width = width.plot1.a, height = height.plot1.a, res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
  pca.scale=F; if (variant.names[i]=='intensity') pca.scale=T
  if (i==1){
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i], show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.5,1.6),legend.cond.xy=c(-1.5,-1.8),legend.textwidth=1)} else {
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, scale=pca.scale, title=variant.names[i])}
}
dev.off()

# F2 (REP + VIOLIN PLOT)
dat.norm.summ.l[c("intensity_fix","ratio_fix")]<- NULL
p1 <- run_effect_plot_rep(dat.norm.summ.l)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
png(paste0(out_path, "mb_unit2.png"), width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
dat.dea[c('intensity_fix', 'ratio_fix')] <- NULL
png(paste0(out_path, "mb_unit3.png"), width.plot3.a, height = height.plot3.a, res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "mb_unit4.png"), width.plot3.a, height = height.plot3.a, res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_outdata_msstatstmt.rda')
names(dat.norm.summ.w2)
variant.names=names(dat.norm.summ.w2)
dat.norm.summ.w2$intensity_lateLog2 <- NULL
dea.info <- get_dea_info(dat.dea)

# F1 (PCA PLOT)
png(paste0(out_path, "dd_unit1.png"), width.plot1.a, height = height.plot1.a, res=res.plot1.a)
par(list(mfrow=c(1,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
for (i in 1:length(dat.norm.summ.w2)){
  #if (variant.names[i]=='intensity') pca.scale=T
  if (i==1){
    pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i], show.legend = T,legend.cex=0.8,legend.run.xy=c(-1.5,1.6),legend.cond.xy=c(-1.5,-1.65),legend.textwidth=1)} else {
      pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])
    }
}
dev.off()

# F2 (REP + VIOLIN PLOT)
p1 <- run_effect_plot_rep(dat.norm.summ.l)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
png(paste0(out_path, "dd_unit2.png"), width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
png(paste0(out_path, "dd_unit3.png"), width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "dd_unit4.png"), width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

### Normalization ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_normalization_outdata_msstatstmt.rda')
variant.names=names(dat.norm.summ.w2)
# F1 (PCA PLOT)
png(paste0(out_path, "mb_norm1.png"), width = width.plot1.b, height = height.plot1.b, res=res.plot1.b)
par(list(mfrow=c(2,2),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
pcaplot_rep(dat.nonnorm.summ.w2%>% select(-'Protein'), info=sample.info, title='raw',show.legend = T,legend.cex=0.8,legend.run.xy=c(1,28),legend.cond.xy=c(1,-2))
for (i in 1:length(dat.norm.summ.w2)){
  pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])
}
dev.off()

# F2 (REP + VIOLIN PLOT)
dat <- list(dat.norm.summ.l$LMM1, dat.norm.summ.l$LMM2, dat.norm.summ.l$LMM3, dat.nonnorm.summ.l) 
names(dat) <- c('LMM1', 'LMM2', 'LMM3', 'raw')
p1 <- run_effect_plot_rep(dat)
dea.info <- get_dea_info(dat.dea)
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
png(paste0(out_path, "mb_norm2.png"), width = width.plot2.a, height = height.plot2.a, res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
# nams <- names(dat.dea)
# dat.dea <- dat.dea[c('raw',nams[!(nams=='raw')])]
# names(dat.dea) <- c('raw', paste('Variant',1:length(variant.names)))

png(paste0(out_path, "mb_norm3.png"), width = width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "mb_norm4.png"), width = width.plot3.b, height = height.plot3.b, res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_normalization_outdata_msstatstmt.rda')
#names(dat.norm.summ.w2)[names(dat.norm.summ.w2) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
variant.names=names(dat.norm.summ.w2)

# F1 (PCA PLOT)
png(paste0(out_path, "dd_norm1.png"), width = width.plot1.c , height = height.plot1.c, res=res.plot1.c)
par(list(mfrow=c(2,3),oma=rep(1,4), mar=c(2.5,2.5,1,1), mgp=c(1.5,0.5,0)))
pcaplot_rep(dat.nonnorm.summ.w2%>% select(-'Protein'), info=sample.info, title='raw',show.legend = T,legend.cex=0.8,legend.run.xy=c(1,28),legend.cond.xy=c(1,-2))
for (i in 1:length(dat.norm.summ.w2)){
  pcaplot_rep(dat.norm.summ.w2[[variant.names[i]]] %>% select(-'Protein'), info=sample.info, title=variant.names[i])
}
dev.off()

# F2 (REP + VIOLIN PLOT)
dat <- list(dat.norm.summ.l$median_sweeping, dat.norm.summ.l$CONSTANd, dat.norm.summ.l$NOMAD,
            dat.norm.summ.l$quantile1, dat.norm.summ.l$quantile2, dat.nonnorm.summ.l)
names(dat) <- c(names(dat.norm.summ.l), 'raw')
#names(dat)[names(dat) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
p1 <- run_effect_plot_rep(dat)
dea.info <- get_dea_info(dat.dea)
#names(dat.dea)[names(dat.dea) %in% c("quantile", "quantileByCondition")] <- c("quantile1", "quantile2")
p2 <- violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125', rot.xaxis.lab=30)
png(paste0(out_path, "dd_norm2.png"), width = width.plot2.a , height = height.plot2.a , res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
#names(dat.dea)[names(dat.dea) %in% c("quantile", "quantile2")] <- c("quantile1", "quantile2")
#names(dat.dea)[names(dat.dea)=='raw'] <- 'raw'
png(paste0(out_path, "dd_norm3.png"), width = width.plot3.c , height = height.plot3.c , res=res.plot3.c)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "dd_norm4.png"), width = width.plot3.c , height = height.plot3.c , res=res.plot3.c)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

### Summarization ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_summarization_outdata_msstatstmt.rda')
variant.names=names(dat.norm.summ.w2)

# F1 (PCA PLOT)
png(paste0(out_path, "mb_summ1.png"), width = width.plot1.a , height = height.plot1.a , res=res.plot1.a)
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
png(paste0(out_path, "mb_summ2.png"), width = width.plot2.a , height = height.plot2.a , res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
png(paste0(out_path, "mb_summ3.png"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "mb_summ4.png"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_summarization_outdata_msstatstmt.rda')
variant.names=names(dat.norm.summ.w2)
# F1 (PCA PLOT)
png(paste0(out_path, "dd_summ1.png"), width = width.plot1.a , height = height.plot1.a , res=res.plot1.a)
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
png(paste0(out_path, "dd_summ2.png"), width = width.plot2.a , height = height.plot2.a , res=res.plot2.a)
ggarrange(p1,p2,labels="AUTO",ncol = 2, nrow=1)
dev.off()

# F3 (2xSPLOM)
png(paste0(out_path, "dd_summ3.png"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "dd_summ4.png"), width = width.plot3.a , height = height.plot3.a , res=res.plot3.a)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

### DEA ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_DEA_outdata_msstatstmt.rda')
variant.names=names(dat.dea)
dea.info <- get_dea_info(dat.dea)

# F2 (REP + VIOLIN PLOT)
png(paste0(out_path, "mb_DEA2.png"), width = width.plot2.b , height = height.plot2.b , res=res.plot2.b)
violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
dev.off()

# F3 (2xSPLOM)
png(paste0(out_path, "mb_DEA3.png"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "mb_DEA4.png"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_DEA_outdata_msstatstmt.rda')
variant.names=names(dat.dea)
dea.info <- get_dea_info(dat.dea)
# F2 (REP + VIOLIN PLOT)
png(paste0(out_path, "dd_DEA2.png"), width = width.plot2.b , height = height.plot2.b , res=res.plot2.b)
violinplot_rep(lapply(dat.dea, function(x) x[spiked.proteins, dea.info$logFC.cols]), referenceCondition, reportContrast='0.125')
dev.off()

# F3 (2xSPLOM)
png(paste0(out_path, "dd_DEA3.png"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$significance.cols, 'q-values', spiked.proteins, reportContrast)
dev.off()

png(paste0(out_path, "dd_DEA4.png"), width = width.plot3.b , height = height.plot3.b , res=res.plot3.b)
scatterplot_rep2(dat.dea, dea.info$logFC.cols, 'logFC', spiked.proteins, reportContrast)
dev.off()
