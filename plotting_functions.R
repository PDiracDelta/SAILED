library(gridExtra)
library(dendextend)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# boxplot.ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# example function arguments' below for quicker development and debugging
# dat.beforenorm <- dat.summplot.l
# dat.afternorm <- dat.normplot.l
# comp <- 'unit'

# fix the issue with variable height of run labels

boxplot.ils <- function(dat, title, ...){  
    
  # first make sure that Run and Channel are factors 
  dat$Run <- as.factor(dat$Run)
  dat$Channel <- as.factor(dat$Channel)
  
  # prepare nice run and channel labels
  run.labels <- paste0('run', 1:length(levels(dat$Run)))
  n.run <- length(run.labels)
  channel.labels <- levels(dat$Channel)
  n.channel <- length(channel.labels)
  
  # prepare x coordinates of vertical lines separating MS runs
  x.lines <- rep(c(1:(n.run-1))*n.channel+0.5, each=3)
  x.lines[seq(3, length(x.lines), by=3)] <- NA
  
  # prepare  y coordinates of vertical lines separating MS runs
  y.range <- range(dat$response)
  y.range <- c(floor(y.range[1]-1), ceiling(y.range[2]+1))
  y.lines <- rep(c(y.range, NA), n.run-1)
  
  boxplot(response~factor(Run:Channel), data=dat, main=title, 
          notch=TRUE, xlab='sample', xaxt='n', ...)
  abline(h=mean(dat$response), col='red', lwd=1.5) # red vertical line
  lines(x.lines, y.lines, lty=1, lwd=2) # draw vertical lines after each run
  text(x = (1:(n.run))*(n.channel)-n.channel/2, # write run labels
       y = par('usr')[4]*1.02,
       labels = run.labels,
       xpd = NA,
       cex = 1.2)
  axis(side = 1, labels = FALSE, at=1:(n.run*n.channel)) # add x-axis tick marks
  text(x = 1:(n.run*n.channel), # add x-axis tick labels (channels)
       y = par("usr")[3] - 0.45,
       labels = channel.labels,
       xpd = NA,
       srt = 35, # rotate the labels by 35 degrees.
       cex = 1)
}

# use case (not run)
# boxplot.ils(dat, 'afternorm')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# boxplot.w
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

boxplot.w <- function(dat, study.design, title, ...){  
  # remove reference cols if still present
  study.design <- study.design[!(study.design$Channel %in% c('126', '131')),]
  # convert to long format again because ggplot for wide data is excruciating
  dat <- to_long_format(dat, study.design = study.design, merge_study_design = F)
  
  boxplot.ils(dat, title, ...)  
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# maplot.ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# investigate correctness of MA plots for model-based normalized data

# example function arguments' below for quicker development and debugging
# design.look <- dat.l %>% distinct(Run, Channel, Condition)
# dat <- dat.summplot.w[[1]] # before normalization
# samples.num <- 'Mixture1_1:130N' # condition1
# samples.denom <- 'Mixture2_2:128C' # condition 0.125
# scale <- 'log' #or raw
# title <- 'Before normalization'

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

maplot.ils <- function(dat, samples.num, samples.denom, scale, title){
  
  dat <- dat %>% drop_na(any_of(c(samples.num, samples.denom)))
  select.scale=match.arg(scale, c('log', 'raw'))
  num <- as.matrix(dat[, samples.num])
  denom <- as.matrix(dat[, samples.denom])

  if (select.scale=='log'){
    
    num <- apply(num, 1, function(x) mean(x, na.rm=TRUE))
    denom <- apply(denom, 1, function(x) mean(x, na.rm=TRUE))
    FC <- num-denom
    AVE <- (num+denom)/2
    
  } else if (select.scale=='raw') {
    
    num <- apply(num, 1, gm_mean)
    denom <- apply(denom, 1, gm_mean)
    FC <- log2(num/denom)
    AVE <- (log2(num)+log2(denom))/2
    
  }
  
  df <- data.frame(FC, AVE)
  
  ggplot(df, aes(x = AVE, y = FC)) +
    geom_point() +
    geom_smooth() +
    scale_y_continuous("log2FC") +
    scale_x_continuous("logAvg_intensity") +
    ggtitle(title) + 
    geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
    geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down
}

# use case (not run)
# samples.num <- dat.summplot.l[[1]] %>% filter(Condition=='1') %>% distinct(Run:Channel) %>% pull
# samples.denom <- dat.summplot.l[[1]] %>% filter(Condition=='0.125') %>% distinct(Run:Channel) %>% pull 
# p1 <- maplot.ils(dat.summplot.w[[1]], samples.num, samples.denom, scale='log', 'Before normalization')
# p2 <- maplot.ils(dat.normplot.w[[1]], samples.num, samples.denom, scale='log', 'After normalization')
# grid.arrange(p1, p2, ncol=2)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# pcaplot.ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# think of a better way to allocate legend
# conditions are not sorted properly

# example function arguments' below for quicker development and debugging
# dat <- dat.summplot.w[[1]]
# run.labels <- stri_replace(unlist(lapply(stri_split(colnames(dat), fixed=':'), function(x) x[1])), fixed='Mixture', 'Mix')
# condition.labels <- conditions.vec
# colour.labels <- cols.vec
# title <- 'Before normalization'

pcaplot.ils=function(dat, run.labels, condition.labels, colour.labels, title, scale=F){
  # drop NA values (they are due to proteins not detected in all runs.)
  pc.cr <- prcomp(t(dat %>% drop_na()), scale = scale, center = TRUE)
  sumpc.cr=summary(pc.cr)
  prop.var=paste('(',round(100*sumpc.cr$importance[2,1:2],2),' %)',sep='')
  axis.lab=paste(c('PC1','PC2'),prop.var)
  
  Ux1<-as.vector(pc.cr$x[,1])
  Ux2<-as.vector(pc.cr$x[,2])
  
  legend.map <- cbind(run=run.labels, run.numeric=as.numeric(as.factor(run.labels)), 
                     condition=condition.labels, colour=colour.labels) %>% as.data.frame 
  # points shape passed to 'pch' argument, starts from 15
  legend.map$run.numeric <- as.numeric(legend.map$run.numeric)+14 
  
  legend1 <- legend.map %>% distinct(run, run.numeric)
  legend2 <- legend.map %>% distinct(condition, colour)
  
  plot(Ux1,Ux2, col=colour.labels, pch=legend.map$run.numeric, 
       main=title, xlab=axis.lab[1], ylab=axis.lab[2], cex=1.5)
  legend('bottomleft', legend=legend1$run, pch=legend1$run.numeric, bty = "n", cex=1.1) # run
  legend('bottomright', legend=legend2$condition, text.col=as.character(legend2$colour), bty = "n", cex=1.1) # condition
}

# use case (not run)
# pcaplot.ils(dat, run.labels, condition.labels, colour.labels, title)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# dendrogram.ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# add legend

# example function arguments' below for quicker development and debugging
# dat <- dat.normplot.w[[1]][1:100, ]
# colour.labels <- cols.vec
# title <- 'Before normalization'

dendrogram.ils <- function(dat, colour.labels, title){
  par.mar.org <- par('mar')
  par(mar=c(3,4,1,7))
  dend_raw <- as.dendrogram(hclust(dist(t(dat %>% drop_na()))))
  labels_colors(dend_raw) <- colour.labels
  plot(dend_raw, horiz=TRUE, main=title)
  par(mar=par.mar.org) 
}

# use case (not run)
# dendrogram.ils(dat.normplot.w[[1]][1:100, ], cols.vec, 'Before normalization')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# volcanoplot.ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# example function arguments' below for quicker development and debugging
# dat <- dat.dea[[1]]
# logFC.cols
# q.cols
# variant.title <- variant.names[1]

volcanoplot.ils <- function(dat, variant.title){
  
  dat.cols <- colnames(dat)
  logFC.cols <- dat.cols[stri_detect(dat.cols, fixed='logFC')]
  q.cols <- dat.cols[stri_detect(dat.cols, fixed='q.mod')]
  
  # list for storing contrast specifc plots
  volcano.plots <- vector('list', length(logFC.cols))
  
  # contrast names
  contrast.names <- stri_replace(logFC.cols, fixed='logFC_', '')
  
  # iterate over contrasts  
  for (j in 1:length(logFC.cols)){
    df <- data.frame(logFC=dat[, logFC.cols[j]],
                     q.mod=dat[, q.cols[j] ])
    volcano.plots[[j]] <- ggplot(df, aes(x = logFC, y = -log10(q.mod))) +
      geom_point() +
      xlab("log2(FC)") +
      ylab("-log10(FDR)") +
      ggtitle(paste(contrast.names[j], variant.title, sep='_' )) +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") + 
      geom_vline(xintercept =  1, color = "black", linetype = "dashed") +
      geom_vline(xintercept = -1, color = "black", linetype = "dashed") 
  }
  
  # put all volcano plots in one grid
  do.call(grid.arrange, volcano.plots) # try to allocate them in columns or reshape the data and switch to facet_wrap/facet_grid  in ggplot
}

# use case (not run)
# volcanoplot.ils(dat.dea[[1]], variant.names[1])

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# violinplot.ils: violin plot for each condition and dashed line with expected fold change
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
violinplot.ils <- function(dat.spiked.logfc.l) {
  conditions.num <- sort(as.numeric(unique(dat.spiked.logfc.l[[1]]$condition)))
  segment_xy <- data.frame(xv=order(conditions.num), yv=log2(conditions.num/as.numeric(referenceCondition)))
  violin.plots <- lapply(dat.spiked.logfc.l, function(x) {
    p <- ggplot(x, aes(x=condition, y=logFC)) + geom_violin(draw_quantiles = TRUE) + 
      geom_segment(data=segment_xy, aes(x=xv-0.25, xend=xv+0.25, y=yv, yend=yv), 
                   color = 'red', linetype = 'dashed')
    return(p)})
  do.call(grid.arrange, violin.plots)
}
