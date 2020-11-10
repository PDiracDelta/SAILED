#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# boxplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

boxplot_ils <- function(dat, title, ...){  
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
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

maplot_ils <- function(dat, samples.num, samples.denom, scale, title){
  dat <- dat %>% drop_na(any_of(c(samples.num, samples.denom)))
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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# pcaplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

pcaplot_ils=function(dat, info, title, scale=F){
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
  
  plot(Ux1,Ux2, col=info$Colour, pch=match(info$Run.short, legend.run)+15, 
       main=title, xlab=axis.lab[1], ylab=axis.lab[2], cex=1.5)
  
  legend('bottomleft', legend=legend.run, pch=as.numeric(legend.run)+15, bty = "n", cex=1.1) # run legend
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

violinplot_ils <- function(dat.spiked.logfc.l) {
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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# CVplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# feature.group arg can be set to 'Protein' or 'Peptide', depending on which level CV should be computed
# xaxis.group arg specifies the groups for which CV will be computed separately. It can be set to, e.g., 
# 'Condition' or 'Run' or 'Mixture' or ...

cvplot_ils <- function(dat, feature.group, xaxis.group, title, rmCVquan=0.95, abs.val=T, ...){  
  dat <- dat %>% ungroup
  # compute CV per feature (Protein/Peptide) 
  if (abs.val){
    CV.df <- dat %>% group_by(across(feature.group), across(xaxis.group)) %>% summarise(CV=sd(response, na.rm=TRUE)/abs(mean(response, na.rm=TRUE)))
  } else {
    CV.df <- dat %>% group_by(across(feature.group), across(xaxis.group)) %>% summarise(CV=sd(response, na.rm=TRUE)/mean(response, na.rm=TRUE))
  }
  # compute CV quantile within the groups
  CV.quantiles.df <- CV.df %>% group_by(across(xaxis.group)) %>% summarise(quan=quantile(CV, rmCVquan, na.rm=TRUE))
  
  # filter out features with CV larger greater than the quantile
  # this is done to get rid of outliers and improve visibility
  CV.df <- left_join(CV.df, CV.quantiles.df, by=xaxis.group) %>% filter(CV<quan)
  
  ff <- formula(paste0('CV~', xaxis.group))
  boxplot(ff, data=CV.df, main=title, notch=T, xlab=xaxis.group, ...)  
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# scatterplot_ils: wrapper function on pairs.panels from 'psych' package
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

scatterplot_ils <- function(dat, cols, stat){
  select.stat <- match.arg(stat, c('p-values', 'log2FC', 'q-values'))
  title <- paste("Spearman's correlation of", select.stat)
  contrast.names <- unlist(lapply(stri_split(cols, fixed='_'), function(x) x[2]))
  
  # align rows (proteins) between DEA variants
  if (length(unique(unlist(lapply(dat, nrow))))>1) stop('Unequal number of analysed proteins across different DEA methods')
  rw <- rownames(dat[[1]])
  dat=lapply(dat, function(x){
    ord <- match(rw, rownames(x))
    return(x[ord,])})
  
  for (i in 1:length(cols)){
    # names(dat) <- NULL # this line generates variant names on the plot
    df <- sapply(dat, function(x) x[, cols[i]]) %>% data.frame %>% drop_na()
    pairs.panels(df, main=paste(title, contrast.names[i], sep='_'), method='spearman', lm=T, pch=16, ellipses=F)
  }
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# volcanoplot_ils
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

volcanoplot_ils <- function(dat, contrast.num, spiked.proteins){
  dat.cols <- colnames(dat[[1]])
  logFC.cols <- dat.cols[stri_detect(dat.cols, fixed='logFC')]
  q.cols <- dat.cols[stri_detect(dat.cols, fixed='q.mod')]
  
  # list for storing variant specifc plots
  volcano.plots <- vector('list', length(dat))
  
  # contrast names
  contrast.names <- stri_replace(logFC.cols, fixed='logFC_', '')
  variant.title <- names(dat)
  
  # iterate over variants
  for (j in 1:length(dat)){
    df <- data.frame(logFC=dat[[j]][, logFC.cols[contrast.num]],
                     q.mod=dat[[j]][, q.cols[contrast.num] ],
                     protein=ifelse(rownames(dat[[j]]) %in% spiked.proteins, 'spiked-in', 'background'))
    
    volcano.plots[[j]] <- ggplot(df, aes(x = logFC, y = -log10(q.mod), color=protein)) +
      geom_point() +
      xlab("log2(FC)") +
      ylab("-log10(FDR)") +
      ggtitle(paste(contrast.names[contrast.num], variant.title[j], sep='_' )) +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") + 
      geom_vline(xintercept =  1, color = "black", linetype = "dashed") +
      geom_vline(xintercept = -1, color = "black", linetype = "dashed") +
      theme(legend.position = "none")
  }
  grid.arrange(grobs = volcano.plots, ncol=length(volcano.plots))
}

