#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# boxplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

boxplot_ils <- function(dat, title, ...){  
  # first make sure that Run and Channel are factors 
  dat$Run <- factor(dat$Run)
  dat$Channel <- factor(dat$Channel)
  dat$Sample <- dat$Run:dat$Channel
  
  # prepare nice run and channel labels
  run.labels <- stri_replace(levels(dat$Run), fixed='Mixture', replacement='Mix')
  n.run <- length(run.labels)
  channel.labels <- levels(dat$Channel)
  n.channel <- length(channel.labels)
  
  # prepare x coordinates of vertical lines separating MS runs
  x.lines <- rep(c(1:(n.run-1))*n.channel+0.5, each=3)
  x.lines[seq(3, length(x.lines), by=3)] <- NA
  
  # prepare  y coordinates of vertical lines separating MS runs
  y.range <- range(dat$response)
  y.range <- c(floor(y.range[1]), ceiling(y.range[2]))
  y.lines <- rep(c(y.range, NA), n.run-1)
  
  boxplot(response~Sample, data=dat, main=title, 
          notch=TRUE, xlab='sample', xaxt='n', ...)
  
  abline(h=mean(dat$response), col='red', lwd=1.5) # red vertical line
  lines(x.lines, y.lines, lty=1, lwd=2) # draw vertical lines after each run
  text(x = (1:(n.run))*(n.channel)-n.channel/2, # write run labels
       y = par('usr')[4]*1.04,
       labels = run.labels,
       xpd = NA,
       cex = 1.2)
  axis(side = 1, labels = FALSE, at=1:(n.run*n.channel)) # add x-axis tick marks
  text(x = 1:(n.run*n.channel), # add x-axis tick labels (channels)
       y = par("usr")[3],
       labels = channel.labels,
       xpd = NA,
       srt = 90, # rotate the labels by 35 degrees.
       cex = 1,
      adj=1.2)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# boxplot.w
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

boxplot_w <- function(dat, sample.info, title, ...){  
  # convert to long format again because ggplot for wide data is excruciating
  dat <- to_long_format(dat, sample.info = sample.info, merge_study_design = F)
  boxplot_ils(dat, title, ...)  
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# maplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# geometric mean function for computing logFC on raw scale

gm_mean = function(x, na.rm=TRUE){
  if (length(x)>1){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  } else x
}

# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
#           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

maplot_ils <- function(dat, samples.num, samples.denom, scale, title, spiked.proteins, xaxis.rank=TRUE){
  dat <- dat %>% drop_na(any_of(c(samples.num, samples.denom)))
  samples.num <- as.character(samples.num)
  samples.denom <- as.character(samples.denom)
  if (!('Protein' %in% colnames(dat))) dat <- dat %>% rownames_to_column('Protein')
  select.scale=match.arg(scale, c('log', 'raw'))
  num <- as.matrix(dat[, samples.num])
  denom <- as.matrix(dat[, samples.denom])
  if (select.scale=='log'){
    num <- apply(num, 1, mean)
    denom <- apply(denom, 1, mean)
    FC <- num-denom
    AVE <- (num+denom)/2
  } else if (select.scale=='raw') {
    num <- apply(num, 1, gm_mean)
    denom <- apply(denom, 1, gm_mean)
    FC <- log2(num/denom)
    AVE <- (log2(num)+log2(denom))/2
  }
  xaxis.title <- 'Avglog2Intensity'
  if (xaxis.rank){
    AVE=rank(AVE) # -length(AVE)+1  descending ranks
    xaxis.title <- 'rank(Avglog2Intensity)'}
  df <- data.frame(FC, AVE, protein.type=ifelse(dat$Protein %in% spiked.proteins, 'spiked-in', 'background'))
  ggplot(df, aes(x = AVE, y = FC)) +
    geom_point(aes(colour=protein.type)) +
    geom_smooth() +
    scale_y_continuous("log2FC") +
    scale_x_continuous(xaxis.title) +
    #scale_colour_manual(values=cbp2) +
    ggtitle(title) + 
    geom_hline(yintercept = 0.0, color = "black") + # one-to-one line
    geom_hline(yintercept = 1.0, color = "black", linetype = "dotted") + # 2-fold up
    geom_hline(yintercept = -1.0, color = "black", linetype = "dotted") # 2-fold down
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# pcaplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

pcaplot_ils=function(dat, info, title, scale=F, shape.vec=c(15,19,17,3)){
  # fix columns order (as in 'info' arg)
  dat <- dat[, match(remove_factors(sample.info$Sample), colnames(dat))]
  # drop NA values (they are due to proteins not detected in all runs.)
  pc.cr <- prcomp(t(dat %>% drop_na()), scale = scale, center = TRUE)
  sumpc.cr=summary(pc.cr)
  prop.var=paste('(',round(100*sumpc.cr$importance[2,1:2],2),' %)',sep='')
  axis.lab=paste(c('PC1','PC2'),prop.var)
  
  Ux1<-as.vector(pc.cr$x[,1])
  Ux2<-as.vector(pc.cr$x[,2])
  
  legend.run <- info %>% distinct(Run.short) %>% pull
  legend.cond <- info %>% distinct(Condition, Colour) %>% arrange(Condition)
  # if shape.vec is shorter than distinct Run values, then assign shapes starting from 0
  if (length(shape.vec)!=length(legend.run)) shape.vec <- 0:(length(legend.run)-1)

  plot(Ux1,Ux2, col=info$Colour, pch=shape.vec[match(info$Run.short, legend.run)], 
       main=title, xlab=axis.lab[1], ylab=axis.lab[2], cex=1.5)
  legend('bottomleft', legend=legend.run, pch=shape.vec[as.numeric(legend.run)], bty = "n", cex=1.1) # run legend
  legend('bottomright', legend=legend.cond$Condition, text.col=legend.cond$Colour, bty = "n", cex=1.1) # condition legend
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# dendrogram_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

dendrogram_ils <- function(dat, info, title){
  # fix columns order (as in 'info' arg)
  dat <- dat[, match(sample.info$Sample, colnames(dat))]
  # columns are already aligned with 'info', so we can use shorter labels from there 
  colnames(dat) <- sample.info$Sample.short 
  par.mar.org <- par('mar')
  par(mar=c(3,4,1,7))
  dend_raw <- as.dendrogram(hclust(dist(t(dat %>% drop_na()))))
  dendcolors = info$Colour[order.dendrogram(dend_raw)]
  labels_colors(dend_raw) <- dendcolors
  plot(dend_raw, horiz=TRUE, main=title)
  par(mar=par.mar.org) 
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# violinplot_ils: violin plot for each condition and dashed line with expected fold change
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

violinplot_ils_old <- function(dat.spiked.logfc.l) {
  conditions.num <- sort(as.numeric(unique(dat.spiked.logfc.l[[1]]$condition)))
  segment_xy <- data.frame(xv=order(conditions.num), yv=log2(conditions.num/as.numeric(referenceCondition)))
  violin.plots <- emptyList(names(dat.spiked.logfc.l))
  for (i in seq_along(violin.plots)) {
    violin.plots[[i]] <- ggplot(dat.spiked.logfc.l[[i]], aes(x=condition, y=logFC)) + 
      geom_violin(draw_quantiles = TRUE) + ggtitle(names(dat.spiked.logfc.l)[i]) + 
      geom_segment(data=segment_xy, aes(x=xv-0.25, xend=xv+0.25, y=yv, yend=yv), 
                   color = 'red', linetype = 'dashed') }
  do.call(grid.arrange, violin.plots)
}

violinplot_ils <- function(dat, referenceCondition){
  variant.names <- names(dat)
  dat.l <- lapply(dat, function(x) {
    x %>% rename_with(function(y) sapply(y, function(z) strsplit(z, '_')[[1]][2])) %>% pivot_longer(cols = everything(), names_to = 'condition', values_to = 'logFC') %>% add_column(Protein=rep(rownames(x), each=length(colnames(x)))) })
  dat.all <- bind_rows(dat.l, .id='variant')
  dat.all$variant <- factor(dat.all$variant, levels=variant.names)
  conditions.num <- sort(as.numeric(unique(dat.all$condition)))
  violin.plots <- vector('list', length(unique(dat.all$condition)))
  for (j in seq_along(violin.plots)) {
    plot.title <- paste0(conditions.num[j], ' vs ', referenceCondition, ' contrast')
    segment_xy <- data.frame(xv=1:length(variant.names), yv=log2(conditions.num[j]/as.numeric(referenceCondition)))    
    violin.plots[[j]] <- dat.all %>% filter(condition==conditions.num[j]) %>%
    ggplot(aes(x=variant, y=logFC)) + 
      geom_violin(draw_quantiles = TRUE) +
      ggtitle(plot.title) +
      geom_segment(data=segment_xy, aes(x=xv-0.25, xend=xv+0.25, y=yv, yend=yv), 
                   color = 'red', linetype = 'dashed') + xlab("")
      }
  do.call(grid.arrange, violin.plots)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# CVplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# feature.group arg can be set to 'Protein' or 'Peptide', depending on which level CV should be computed
# xaxis.group arg specifies the groups for which CV will be computed separately. It can be set to, e.g., 
# 'Condition' or 'Run' or 'Mixture' or ...

cvplot_ils <- function(dat, feature.group, xaxis.group, title, rmCVquan=0.95, ...){  
  dat <- dat %>% ungroup
  # compute CV per feature (Protein/Peptide) 
  CV.df <- dat %>% group_by(across(feature.group), across(xaxis.group)) %>% summarise(CV=sd(response, na.rm=TRUE)/mean(response, na.rm=TRUE))
  # compute CV quantile within the groups
  CV.quantiles.df <- CV.df %>% group_by(across(xaxis.group)) %>% summarise(quan=quantile(CV, rmCVquan, na.rm=TRUE))

  # avg of median CV within xaxis.group
  avgCV <- CV.df %>% group_by(across(xaxis.group)) %>% summarise(medCV=median(CV)) %>% pull %>% mean(na.rm=TRUE)
  avgCV <- round(avgCV, 5)
  
  # filter out features with CV larger greater than the quantile
  # this is done to get rid of outliers and improve visibility
  CV.df <- left_join(CV.df, CV.quantiles.df, by=xaxis.group) %>% filter(CV<quan)
  
  ff <- formula(paste0('CV~', xaxis.group))
  boxplot(ff, data=CV.df, main=paste0(title, ' (avg CV=', avgCV,')'), notch=T, xlab=xaxis.group, ...)  
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# scatterplot_ils: wrapper function on pairs.panels from 'psych' package
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# pairs.panels.my is a modified pairs.panels function such that the y=x identity line is plotted when lm=T
source('pairs_panels_idline.R')

scatterplot_ils <- function(dat, cols, stat, spiked.proteins){
  select.stat <- match.arg(stat, c('p-values', 'log2FC', 'q-values'))
  title <- paste("Spearman's correlation of", select.stat)
  contrast.names <- unlist(lapply(stri_split(cols, fixed='_'), function(x) x[2]))
  
  # align rows (proteins) between DEA variants
  if (length(unique(unlist(lapply(dat, nrow))))>1) stop('Unequal number of analysed proteins across different DEA methods')
  rw <- rownames(dat[[1]])
  dat=lapply(dat, function(x){
    ord <- match(rw, rownames(x))
    return(x[ord,])})
  rw <- rownames(dat[[1]])
  for (i in 1:length(cols)){
    df <- sapply(dat, function(x) x[, cols[i]]) %>% data.frame
    pairs.panels.idline(df, main=paste(title, paste0(contrast.names[i], ' vs ', referenceCondition, ' contrast'), sep='-')
                    , method='spearman', lm=T, ellipses = FALSE, 
                    ,pch=ifelse(rw %in% spiked.proteins, 'X', 'o') 
                    #,pch=ifelse(rw2 %in% spiked.proteins, 2, 1) 
                    ,col.points=ifelse(rw %in% spiked.proteins, '#E69F00', '#000000')
                    ,cex.points=ifelse(rw %in% spiked.proteins, 2, 1))
  }
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# volcanoplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

volcanoplot_ils <- function(dat, contrast.num, spiked.proteins){
  dat.cols <- colnames(dat[[1]])
  logFC.cols <- dat.cols[stri_detect(dat.cols, fixed='logFC')]
  significance.cols <- dat.cols[stri_detect(dat.cols, fixed='q.mod')]
  
  # list for storing variant specifc plots
  volcano.plots <- vector('list', length(dat))
  
  # contrast names
  contrast.names <- stri_replace(logFC.cols, fixed='logFC_', '')
  variant.title <- names(dat)
  
  # iterate over variants
  for (j in 1:length(dat)){
    df <- data.frame(logFC=dat[[j]][, logFC.cols[contrast.num]],
                     score.mod=dat[[j]][, significance.cols[contrast.num] ],
                     protein=ifelse(rownames(dat[[j]]) %in% spiked.proteins, 'spiked-in', 'background'))
    
    volcano.plots[[j]] <- ggplot(df, aes(x = logFC, y = -log10(score.mod), color=protein)) +
      geom_point() +
      xlab("log2(FC)") +
      ylab("-log10(q-value)") +
      ggtitle(paste(paste0(contrast.names[contrast.num], ' vs ', referenceCondition, ' contrast'), variant.title[j], sep='_' )) +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") + 
      geom_vline(xintercept =  1, color = "black", linetype = "dashed") +
      geom_vline(xintercept = -1, color = "black", linetype = "dashed") +
      theme(legend.position = "none")
  }
  #grid.arrange(grobs = volcano.plots, ncol=length(volcano.plots))
  grid.arrange(grobs = volcano.plots, ncol=2, nrow=ceiling(length(dat)/2))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# hist_ils - a very simple wrapper on hist
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

hist_ils <- function(x, ...) hist(x, xlab='Run effect p-value', breaks=15, xlim=c(0,1), ...)
