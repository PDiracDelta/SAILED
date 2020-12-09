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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# pcaplot_rep
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# dat=dat.nonnorm.summ.w2 %>% select(-'Protein')
# info=sample.info
# title='aa'
# scale=F;
# shape.vec=c(15,19,17,3)
# par.list=list(oma=c(3,2,2,2), mar=c(3,3,2,1), mgp=c(1.5,0.5,0))

pcaplot_rep=function(dat, info, title, scale=F, shape.vec=c(15,19,17,3),show.legend=F, par.list=NULL,...){
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
  #plot(Ux1,Ux2, col=info$Color, pch=shape.vec[match(info$Run.short, legend.run)], xlab=axis.lab[1], ylab=axis.lab[2],...)
  plot(Ux1,Ux2, col=info$Color, pch=shape.vec[match(info$Run.short, legend.run)], xlab=axis.lab[1], ylab=axis.lab[2])
  grid()
  if (show.legend){
    legend('top', legend=paste0('run',as.numeric(legend.run)), pch=shape.vec[as.numeric(legend.run)]) # run legend
    legend('bottom', legend=legend.cond$Condition, text.col=legend.cond$Color) # condition legend
  }
  #bty = "n"
  # mtext(side = 1, text = axis.lab[1], line = 2, cex=0.8)
  # mtext(side = 2, text = axis.lab[2], line = 2, cex=0.8)
  mtext(side=3, text=title, line=0, font=2)
  #par(op)
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
  dat<- lapply(dat, run_test) %>% bind_rows(.id = "Variant")
  ggplot(dat, aes(x=pvalues, group=Variant, colour=Variant)) +
    stat_density(aes(x=pvalues, y=..scaled..,color=Variant), position="dodge", geom="line", size=1.5)+
    ggtitle(main.title) +
    xlab('p-value') +
    ylab('scaled density') + 
    theme(legend.title = element_blank())+
    theme(legend.position="bottom")
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# scatterplot_rep
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# library(GGally)
# dat=dat.dea$data_driven %>% rownames_to_column('Protein') %>% select(Protein, logFC_0.125) %>%
#   left_join(dat.dea$model_based %>% rownames_to_column('Protein') %>% select(Protein, logFC_0.125), by='Protein') %>% select(-Protein)
# m <- ggpairs(dat)
# for (i in 2:m$nrow) {
#   for (j in 1:(i-1)) {
#     m[i,j] = m[i,j] + geom_abline(intercept=0,slope=1)
#   }
# }

# dat=dat.dea
# cols=significance.cols
# stats='q-values'
# spiked.proteins=spiked.proteins
# plot.title='Aaa'
# reportContrast='0.125'
# xlab.title='a'
# ylab.title='b'

scatterplot_rep <- function(dat, cols, spiked.proteins, xlab.title, ylab.title, reportContrast){
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
    df$protein=ifelse(rw %in% spiked.proteins, 'spiked-in', 'background')
    if (stri_detect(cols[i], fixed=reportContrast)){
      plots[[i]] <- df %>% ggplot +
        geom_point(aes(x=xx, y=yy, colour=protein))+
        geom_abline(intercept=0, slope=1, col='black') +
        xlab(xlab.title) + 
        ylab(ylab.title) + 
        #ggtitle(plot.title) + 
        theme_bw() +
        theme(legend.position="bottom", legend.title = element_blank(), axis.title=element_text(size=8))+
        annotate("text", size=4, x = x.p.corr, y=y.p.corr, label = paste0('r=',p.corr), hjust=1)
    }     
  }
  return(plots[[which(!is.null(plots))]])
}

# dat=dat.dea;
# cols=logFC.cols
# stat='log2FC'
# spiked.proteins

scatterplot_rep2 <- function(dat, cols, stat, spiked.proteins, reportContrast){
  source('pairs_panels_idline.R')
  select.stat <- match.arg(stat, c('p-values', 'log2FC', 'q-values'))
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
    pairs.panels.idline(df, main=paste(title, paste0(contrast.names[i], ' vs ', referenceCondition, ' contrast'), sep='-')
                        , method='pearson', lm=T, ellipses = FALSE, 
                        ,pch=ifelse(rw %in% spiked.proteins, 'X', 'o') 
                        #,pch=ifelse(rw2 %in% spiked.proteins, 2, 1) 
                        ,col.points=ifelse(rw %in% spiked.proteins, '#E69F00', '#000000')
                        ,cex.points=ifelse(rw %in% spiked.proteins, 2, 1))
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
                     protein=ifelse(rownames(dat[[j]]) %in% spiked.proteins, 'spiked-in', 'background'))
    
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

violinplot_rep <- function(dat, referenceCondition, reportContrast){
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
    geom_violin(draw_quantiles = TRUE) +
    #ggtitle(plot.title) +
    geom_segment(data=segment_xy, aes(x=xv-0.25, xend=xv+0.25, y=yv, yend=yv), 
                 color = 'red', linetype = 'dashed') + xlab("")
  return(violin.plots)
}
