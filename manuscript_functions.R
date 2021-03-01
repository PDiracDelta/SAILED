remove_factors <- function(x) {
  if (is.factor(x)) return(levels(x)[x]) else return(x)
}

get_sample_info <- function(dat, map){
  tmp <- dat %>% select(Run, Channel, Condition) %>% distinct %>% 
    left_join(map, by='Condition') %>% 
    mutate(Run.short=factor(stri_replace(Run, fixed='Mixture', 'Mix')), 
           Sample=Run:Channel,
           Sample.short=Run.short:Channel) %>% 
    relocate(Run.short, .after=Run) %>%
    relocate(c(Sample, Sample.short), .after=Channel) %>% arrange(Sample)
  return(tmp)
}

get_dea_info <- function(dat){
  dea.cols <- colnames(dat[[1]])
  logFC.cols <- dea.cols[stri_detect_fixed(dea.cols, 'logFC')]
  significance.cols <- dea.cols[stri_detect_fixed(dea.cols, 'q.mod')]
  n.contrasts <- length(logFC.cols)  
  return(list(logFC.cols=logFC.cols, significance.cols=significance.cols, n.contrasts=n.contrasts))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# pcaplot_rep
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

pcaplot_rep=function(dat, info, title, scale=F, shape.vec=c(15,19,17,3),show.legend=F, par.list=NULL, 
                     legend.cex=1,legend.run.xy=NULL, legend.cond.xy=NULL, legend.textwidth=1, cex.lab=0.9, cex.main=0.9, ...){
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
  legend.cond <- info %>% distinct(Condition, Color) %>% arrange(Condition)
  # if shape.vec is shorter than distinct Run values, then assign shapes starting from 0
  if (length(shape.vec)!=length(legend.run)) shape.vec <- 0:(length(legend.run)-1)
  if (!is.null(par.list)) par(par.list)
  plot(Ux1,Ux2, col=info$Color, pch=shape.vec[match(info$Run.short, legend.run)], xlab=axis.lab[1], ylab=axis.lab[2], cex.lab=cex.lab, ...)
  #plot(Ux1,Ux2, col=info$Color, pch=shape.vec[match(info$Run.short, legend.run)], xlab=axis.lab[1], ylab=axis.lab[2], cex.lab=0.9)
  grid()
  
  if (show.legend){
    text.run <- paste0('run',as.numeric(legend.run))
    text.cond <- legend.cond$Condition
    legend(legend.run.xy[1],legend.run.xy[2], text.width = strwidth(text.run)[1]*legend.textwidth, legend=text.run, pch=shape.vec[as.numeric(legend.run)], cex=legend.cex) # run legend
    legend(legend.cond.xy[1],legend.cond.xy[2], text.width = strwidth(text.cond)[1]*legend.textwidth, legend=text.cond, text.col=legend.cond$Color, cex=legend.cex) # condition legend
  }
  mtext(side=3, text=title, line=0, cex=cex.main, font=2)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# run_test - one-way ANOVA for testing the run effect
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# dat - data in long format

run_test <- function(dat){
  unique.proteins <- dat %>% group_by(Protein) %>% summarise(ndist=n_distinct(Run)) %>% filter(ndist>1) %>% pull(Protein) %>% as.character
  dat2 <- dat %>% filter(Protein %in% unique.proteins)
  pvalues=sapply(split(dat2, as.character(dat2$Protein)), function(x) summary(aov(x$response~x$Run, data=x))[[1]]$`Pr(>F)`[1])
  pvalues=as.data.frame(pvalues) %>% rownames_to_column('Protein')
  return(pvalues)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# run_effect_plot - Run effect p-value plot showing the distribution of p-values for
# all variants in one plot
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# dat is a named list

run_effect_plot_rep <- function(dat, main.title=''){
  if(is.null(names(dat))) stop('List names are missing')
  dat <- lapply(dat, run_test) %>% bind_rows(.id = "Variant")
  dat$Variant <- factor(dat$Variant, levels=(unique(dat$Variant)))
  ggplot(dat, aes(x=pvalues, group=Variant, colour=Variant)) +
    stat_density(aes(x=pvalues, y=..scaled..,color=Variant), position="dodge", geom="line", size=1.1)+
    ggtitle(main.title) +
    scale_color_viridis(discrete=TRUE) +
    xlab('p-value') +
    ylab('scaled density') + 
    theme_bw() +
    theme(legend.position="bottom", legend.title = element_blank())
  
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# scatterplot_rep
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

scatterplot_rep <- function(dat, cols, spiked.proteins, xlab.title, ylab.title, reportContrast, legend_position=NULL){
  variant.names <- names(dat)
  # align rows (proteins) between DEA variants
  if (length(unique(unlist(lapply(dat, nrow))))>1) stop('Unequal number of analysed proteins across different DEA methods')
  rw <- rownames(dat[[1]])
  dat=lapply(dat, function(x){
    ord <- match(rw, rownames(x))
    return(x[ord,])})
  rw <- rownames(dat[[1]])
  plots <- vector('list',length(cols))
  for (i in 1:length(cols)){
    # fix column order, variant1 is xx, variant2 is yy
    df <- sapply(dat, function(x) x[, cols[i]]) %>% data.frame 
    a <- match(variant.names[1], colnames(df))
    b <- match(variant.names[2], colnames(df))
    colnames(df)[a] <- 'xx'
    colnames(df)[b] <- 'yy'
    # show Pearson's correlation on the plot
    p.corr=round(cor(df$xx,df$yy),2)
    x.p.corr=range(df$xx)[2]
    y.p.corr=range(df$yy)[1]
    df$protein=ifelse(rw %in% spiked.proteins, 'spike-in', 'background')
    df$protein=relevel(as.factor(df$protein), ref='spike-in')
    if (stri_detect(cols[i], fixed=reportContrast)){
      plots[[i]] <- df %>% ggplot +
        geom_point(aes(x=xx, y=yy, colour=protein, shape=protein, size=protein))+
        geom_abline(intercept=0, slope=1, col='red') +
        scale_shape_manual(values=c('X', 'o'))+
        scale_color_manual(values=c('#E69F00','#000000'))+
        scale_size_manual(values=c(4,2))+
        xlab(xlab.title) + 
        ylab(ylab.title) + 
        #ggtitle(plot.title) + 
        theme_bw() +
        theme(legend.position = "none", legend.title = element_blank(), axis.title=element_text(size=8))+
        annotate("text", size=4, x = x.p.corr, y=y.p.corr, label = paste0('r=',p.corr), hjust=1)
      if (!is.null(legend_position)) plots[[i]] <- plots[[i]]+theme(legend.position=legend_position, legend.background=element_blank())
    } 
  }
  return(plots[[which(!is.null(plots))]])
}

scatterplot_rep2 <- function(dat, cols, stat, spiked.proteins, reportContrast){
  source('pairs_panels_idline.R')
  select.stat <- match.arg(stat, c('p-values', 'logFC', 'q-values'))
  title <- paste("Pearson's correlation of", select.stat)
  contrast.names <- unlist(lapply(stri_split(cols, fixed='_'), function(x) x[2]))
  
  # align rows (proteins) between DEA variants
  if (length(unique(unlist(lapply(dat, nrow))))>1) stop('Unequal number of analysed proteins across different DEA methods')
  rw <- rownames(dat[[1]])
  dat=lapply(dat, function(x){
    ord <- match(rw, rownames(x))
    return(x[ord,])})
  rw <- rownames(dat[[1]])
    i=which(stri_detect(cols, fixed=reportContrast))
    df <- sapply(dat, function(x) x[, cols[i]]) %>% data.frame
    pairs.panels.idline(df, main=title
                        , method='pearson', lm=T, ellipses = FALSE, 
                        ,pch=ifelse(rw %in% spiked.proteins, 'X', 'o') 
                        ,col.points=ifelse(rw %in% spiked.proteins, '#E69F00', '#000000')
                        ,cex.points=ifelse(rw %in% spiked.proteins, 2, 1)
                        ,digits=3)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# volcanoplot_rep
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#dat=lapply(dat.dea, function(x) x[spiked.proteins, logFC.cols])
volcanoplot_rep <- function(dat, contrast.num, spiked.proteins, refCond=referenceCondition){
  dat.cols <- colnames(dat[[1]])
  logFC.cols <- dat.cols[stri_detect(dat.cols, fixed='logFC')]
  significance.cols <- dat.cols[stri_detect(dat.cols, fixed='q.mod')]
  
  # list for storing variant specifc plots
  volcano.plots <- vector('list', length(dat))
  
  # contrast names
  contrast.names <- stri_replace(logFC.cols, fixed='logFC_', '')
  variant.title <- names(dat)
  
  # compute axis ranges
  # if it's possible to convert contrast names and refCond into numeric, plot true logFC 
  if (all(!is.na(as.numeric(contrast.names[contrast.num], refCond)))) {
    true.logFC <- log2(as.numeric(contrast.names[contrast.num])/as.numeric(refCond))} else true.logFC=0
  x.range <- range(c(true.logFC, unlist(lapply(dat, function(x) x[,logFC.cols[contrast.num]]))), na.rm=TRUE)
  y.max <- -log10(min(unlist(lapply(dat, function(x) x[,significance.cols[contrast.num]])), na.rm=TRUE))
  
  # iterate over variants
  for (j in 1:length(dat)){
    df <- data.frame(logFC=dat[[j]][, logFC.cols[contrast.num]],
                     score.mod=dat[[j]][, significance.cols[contrast.num] ],
                     protein=ifelse(rownames(dat[[j]]) %in% spiked.proteins, 'spike-in', 'background'))
    
    volcano.plots[[j]] <- ggplot(df, aes(x = logFC, y = -log10(score.mod), color=protein)) +
      geom_point() +
      xlab("log2(FC)") +
      ylab("-log10(q-value)") +
      ggtitle(paste(paste0(contrast.names[contrast.num], ' vs ', refCond, ' contrast'), variant.title[j], sep='_' )) +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") + 
      geom_vline(xintercept =  1, color = "black", linetype = "dashed") +
      geom_vline(xintercept = -1, color = "black", linetype = "dashed") +
      geom_vline(xintercept =  true.logFC, color = "violet", linetype = "dashed") +
      xlim(c(min(-1,x.range[1]), max(1,x.range[2])))+
      ylim(c(0,y.max))+
      theme(legend.position = "none")
  }
  #grid.arrange(grobs = volcano.plots, ncol=length(volcano.plots))
  grid.arrange(grobs = volcano.plots, ncol=2, nrow=ceiling(length(dat)/2))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# violinplot_rep: violin plot for each condition and dashed line with expected fold change
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

violinplot_rep <- function(dat, referenceCondition, reportContrast, rot.xaxis.lab=0){
  variant.names <- names(dat)
  dat.l <- lapply(dat, function(x) {
    x %>% rename_with(function(y) sapply(y, function(z) strsplit(z, '_')[[1]][2])) %>% pivot_longer(cols = everything(), names_to = 'condition', values_to = 'logFC') %>% add_column(Protein=rep(rownames(x), each=length(colnames(x)))) })
  dat.all <- bind_rows(dat.l, .id='variant') %>% filter(condition==reportContrast)
  dat.all$variant <- factor(dat.all$variant, levels=variant.names)
  conditions.num <- sort(as.numeric(unique(dat.all$condition)))
  plot.title <- paste0(conditions.num, ' vs ', referenceCondition, ' contrast')
  segment_xy <- data.frame(xv=1:length(variant.names), yv=log2(conditions.num/as.numeric(referenceCondition)))    
  violin.plots <- dat.all %>% filter(condition==conditions.num) %>%
    ggplot(aes(x=variant, y=logFC)) + 
    geom_violin(draw_quantiles = TRUE, fill='gray82') +
    #scale_color_viridis(discrete=TRUE) +
    geom_segment(data=segment_xy, aes(x=xv-0.25, xend=xv+0.25, y=yv, yend=yv), 
                 color = 'red', linetype = 'dashed') + xlab("")+
    theme_bw()+
    theme(legend.position='none',legend.title = element_blank(), axis.text.x = element_text(angle = rot.xaxis.lab))
  return(violin.plots)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# wrapper on confusionMatrix from caret package
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

conf_mat <- function(dat, score.var, significance.threshold, spiked.proteins){
  # extract contrast names
  dat.cols <- colnames(dat[[1]])
  logFC.cols <- dat.cols[stri_detect(dat.cols, fixed='logFC')]
  
  contrast.names <- stri_replace(logFC.cols, fixed='logFC_', '')
  n.contrasts <- length(contrast.names)
  
  nam.variants <- names(dat)
  n.variants <- length(dat)
  
  score.vars=paste(score.var, contrast.names, sep='_')  
  
  results <- vector('list', n.contrasts)
  names(results)=contrast.names
  
  tab <- vector('list', n.variants)
  names(tab)=names(dat)
  stats <- matrix(NA, ncol=n.variants, nrow=5)
  colnames(stats) <- nam.variants
  rownames(stats) <- c('Accuracy', 'Sensitivity', 'Specificity', 'PPV', 'NPV')
  
  for (i in 1:length(contrast.names)){
    tab[] <- NA
    stats[] <- NA
    for (j in 1:n.variants){
      pred.class <- factor(ifelse(dat[[j]][, score.vars[i]]<significance.threshold, 'DE','not_DE'), levels=c('not_DE', 'DE'))
      true.class<- factor(ifelse(rownames(dat[[j]]) %in% spiked.proteins, 'DE','not_DE'), levels=c('not_DE', 'DE'))
      tmp=caret::confusionMatrix(data=pred.class, reference=true.class, positive='DE')
      tab[[j]] <- tmp$table
      stats[, j] <- c(tmp$overall['Accuracy'], tmp$byClass[c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value')])
    }
    tmp <- matrix(unlist(lapply(tab, as.numeric)), nrow=2)
    colnames(tmp) <- rep(c('background', 'spiked'), n.variants)
    rownames(tmp) <- c('not_DE', 'DE')
    results[[i]] <- list(tab=tmp, stats=stats)
  }
  return(results)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# print conf_mat output
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#scale_font=999 will use 'scale_down' option to automatically decrease/increase font size to fit into the page
print_conf_mat <- function(dat, reportContrast, cap, scale_font=9){
  
  variant.names <- colnames(dat[[1]]$stats)
  myHeaderVariant <- c(1, rep(2, length(variant.names)))
  names(myHeaderVariant) <- c(" ",variant.names)
  # dat is a list of size equal to # of contrasts
  # output is presented by contrasts
  
  tabs <- vector('list', length(dat))
  
  for (i in 1:length(dat)){
    # print confusion table counts  
      if (scale_font==999){
        k1 <- kable(dat[[i]]$tab, format='latex', booktabs = T, caption=cap) %>%
          kable_styling(latex_options = "scale_down", "HOLD_position") %>%
          add_header_above(myHeaderVariant)
      } else {
        k1 <- kable(dat[[i]]$tab, format='latex', booktabs = T, caption=cap) %>%
          kable_styling(font_size = scale_font, "HOLD_position") %>%
          add_header_above(myHeaderVariant)}
    tabs[[i]] <- k1
  }
  tabs[[which(names(dat)==reportContrast)]]
}

# this version of 'print_conf_mat' saves each conf matrix as a separate graph in pdf 
print_conf_mat2 <- function(dat, reportContrast, fileName, scale_font=9){
  
  variant.names <- colnames(dat[[1]]$stats)
  myHeaderVariant <- c(1, rep(2, length(variant.names)))
  names(myHeaderVariant) <- c(" ",variant.names)
  # dat is a list of size equal to # of contrasts
  # output is presented by contrasts
  
  tabs <- vector('list', length(dat))
  
  for (i in 1:length(dat)){
    # print confusion table counts  
    if (scale_font==999){
      #, caption=cap
      k1 <- kable(dat[[i]]$tab, format='latex', booktabs = T) %>%
        kable_styling(latex_options = "scale_down", "HOLD_position") %>%
        add_header_above(myHeaderVariant)  #%>% save_kable(file=fileName)
    } else {
      k1 <- kable(dat[[i]]$tab, format='latex', booktabs = T) %>%
        kable_styling(font_size = scale_font, "HOLD_position") %>%
        add_header_above(myHeaderVariant)} #%>% save_kable(file=fileName)
    tabs[[i]] <- k1
  }
  tabs[[which(names(dat)==reportContrast)]] %>% save_kable(file=fileName)
}


