library(tidyverse)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# make empty list with names pre-defined.
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
emptyList <- function(names) {
  return(sapply(names,function(x) NULL))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# convert data to long format (assume Run, Mixture columns already present)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

to_long_format<-function(x, study.design, merge_study_design=T) {
  # remove reference cols if still present
  study.design <- study.design[!(study.design$Channel %in% c('126', '131')),]
  quan.cols = unique(study.design$Channel)
  x <- x %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='response', values_drop_na=FALSE)
  
  # merge Condition, TechRepMixture, BioReplicate variables from study.design
  # remove cols to add if already present to avoid errors.
  if (merge_study_design) {
    x <- x[,!(colnames(x) %in% c('Condition', 'TechRepMixture', 'BioReplicate'))]
    x <- left_join(x, study.design, by=c('Mixture', 'Run', 'Channel')) %>%
    relocate(TechRepMixture, .after=Mixture) %>%
    relocate(Condition, .after=TechRepMixture) %>%
    relocate(BioReplicate, .after=Condition) %>% select(-X)
  }
  return(x)
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function performing mean or median aggregation of variables specified in var.names argument
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

aggFunc=function(dat, var.names, agg.method='mean'){
  
  library(dtplyr)
  select.met=match.arg(agg.method, c('mean', 'median', 'sum'))
  
  dat2 <- lazy_dt(dat)
  out.dat=dat2 %>%
      group_by(Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate, Protein, Peptide) %>%
      summarize_at(var.names, eval(parse(text=select.met))) %>%
    as_tibble()
    
  return(out.dat)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function for mixed models DEA (without empirical bayes moderation)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

mixed.model.dea <- function(dat, mod.formula, conditions){
  
  response <- stri_trim_both( stri_split_fixed(mod.formula, pattern='~')[[1]][1]) # extract response var name
  mod.formula <- as.formula(mod.formula)
  dat <- dat[dat$Condition %in% conditions, ] # filter out other conditions
  
  # proteins with at least one not missing quantification value
  # in both conditions
  proteins <- dat %>% group_by(Protein, Condition) %>% 
    summarize(n=sum(!is.na(response))) %>% filter(n>=1) %>%
    group_by(Protein) %>% summarize(nn=n()) %>%
    filter(nn==2) %>% pull(Protein)
  
  nproteins <- length(proteins)
  
  results <- matrix(NA, nrow=length(proteins), ncol=7)
  colnames(results) <- c('logFC', 't.ord', 't.mod', 'p.ord', 'p.mod', 'q.ord', 'q.mod')
  rownames(results) <- proteins
  
  for (i in seq_along(proteins)) {
    dat.sub <- dat[dat$Protein==proteins[i], ]
    mod=lmer(mod.formula, data=dat.sub, control=lmerControl(check.nobs.vs.nlev="warning", check.nobs.vs.nRE="warning"))
    sum.mod=summary(mod)
    logFC <- sum.mod$coefficients[2,1] # estimate of the log2-fold-change corresponding to the effect size
    df.r <- NA # residual degrees of freedom assiciated with ordinary t-statistic and p-value
    df.0 <- NA # degrees of freedom associated with s2.0
    s2.0 <- NA # estimated prior value for the variance
    s2 <- vcov(mod)[2, 2] # sample variance
    s2.post <- s2 # posterior value for the variance
    t.ord <- sum.mod$coefficients[2,4] # vector of ordinary t-statistic: using sigma vs s2.post
    t.mod <- sum.mod$coefficients[2,4] # moderated t-statistic
    p.ord <- sum.mod$coefficients[2,5] # ordinary p-value corresonding to the ordinary t-statistic
    p.mod <- sum.mod$coefficients[2,5] # moderated p-value corresonding to the moderated t-statistic
    
    results[i, ] <- c(logFC, t.ord, t.mod, p.ord, p.mod, NA, NA)
    #print(i)
  } 
  
  results <- as.data.frame(results)
  if (nproteins>1) results$q.ord <- p.adjust(results$p.ord , method='BH') # ordinary q-value corresponding to the ordinary t-statistic
  else results$q.ord <- results$p.ord
  if(nproteins>1) results$q.mod <- p.adjust(results$p.mod , method='BH') # moderated q-value corresponding to the moderated t-statistic
  else results$q.mod <- results$p.mod
  
  # results <- results %>% mutate(candidate=factor(case_when(q.mod <.01 ~ 'high',
  #                                            q.mod >.01 & q.mod<.05 ~ 'med',
  #                                            q.mod >.05 & q.mod<.1 ~ 'low',
  #                                            q.mod >.1 ~ 'no'), levels=c('high', 'med', 'low', 'no')))
  
  return(results)
} 

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function for printing # of up/down/not regulated proteins
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

regulated.proteins <- function(dea.mat, score.var, conditions, cut.off){
  cat('numerator condition:', conditions[1],', ' , 'denominator condition:', conditions[2], '\n')
  c(`Up`=sum(dea.mat[, score.var]<cut.off & dea.mat[, 'logFC']>0),
    `Down`=sum(dea.mat[, score.var]<cut.off & dea.mat[, 'logFC']<0),
    `Not Sig`=sum(dea.mat[, score.var]>cut.off))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# wrapper on confusionMatrix from caret package
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# example function arguments below for quicker development and debugging
# dat=dat.dea[[1]]
# score.var='q.mod'
# cut.off <- 0.05
# i <- 1

conf.mat <- function(dat, score.var,cut.off){
  
  dat.cols <- colnames(dat)
  logFC.cols <- dat.cols[stri_detect(dat.cols, fixed='logFC')]
  contrast.names <- stri_replace(logFC.cols, fixed='logFC_', '')
  
  score.vars=paste(score.var, contrast.names, sep='_')  
  conf.list <- vector('list', length(contrast.names))
  conf.metrics <- matrix(NA, ncol=length(score.vars), nrow=5)
  colnames(conf.metrics) <- contrast.names
  rownames(conf.metrics) <- c('Accuracy', 'Sensitivity', 'Specificity', 'PPV', 'NPV')
  
  for (i in 1:length(contrast.names)){
    pred.class <- factor(ifelse(dat[, score.vars[i]]<cut.off, 'DE','not DE'), levels=c('not DE', 'DE'))
    true.class<- factor(ifelse(rownames(dat) %in% spiked.proteins, 'DE','not DE'), levels=c('not DE', 'DE'))
    tmp=caret::confusionMatrix(data=pred.class, reference=true.class, positive='DE')
    tmp.tab <- tmp$table
    conf.list[[i]] <- tmp$table
    conf.metrics[, i] <- c(tmp$overall['Accuracy'], tmp$byClass[c('Sensitivity', 'Specificity', 'Pos Pred Value', 'Neg Pred Value')])
  }
  
  conf.tab <- matrix(unlist(conf.list), ncol=2, byrow=TRUE)
  conf.tab <- cbind(rep(contrast.names, each=2), conf.tab)
  colnames(conf.tab) <- c('contrast', 'background','spiked')
  rownames(conf.tab) <- rep(c('not DEA', 'DEA'), length(contrast.names))
  return(list(tab=conf.tab, metrics=conf.metrics))
}

# use case (not run)
# cm <- conf.mat(dat.dea[[1]], 'q.mod', contrast.names, 0.05)
# xtable(cm$tab, type='html')
# xtable(cm$metrics, type='html')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# get design matrix for use in moderated t-test
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

get_design_matrix <- function(referenceCondition, study.design) {
  # ANOVA-like design matrix for use in moderated_ttest, indicating group 
  # (condition) membership of each entry in all_channels.
  
  # remove internal reference cols if still present
  study.design <- study.design[!(study.design$Channel %in% c('126', '131')),]
  
  otherConditions = setdiff(unique(study.design$Condition), referenceCondition)
  study.design.unitedchannel <- study.design %>% unite(Channel, c(Channel,Run))
  all_channels = study.design.unitedchannel %>% select(Channel) %>% pull
  N_channels = length(all_channels) 
  N_conditions = 1+length(otherConditions)
  design = matrix(rep(0,N_channels*N_conditions), c(N_channels, N_conditions))
  design[, 1] <- 1  # reference gets 1 everywhere
  rownames(design) <- all_channels
  colnames(design) <- c(referenceCondition, otherConditions)
  for (i in 1:N_channels) {  # for each channel in each condition, put a "1" in the design matrix.
    design[study.design.unitedchannel$Channel[i], study.design.unitedchannel$Condition[i]] = 1
  }
  return(design)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# moderated t-test for data-driven approach
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

moderated_ttest <- function(dat, design, scale) {
  library(limma)
  # This version makes MAXIMAL use of limma. It is an extension of eb.fit, but 
  # with terminology from moderated_ttest_extracted, which returns 
  # output for more than 1 sample comparison.
  design <- design[match(colnames(dat), rownames(design)),]  # fix column order
  ngenes <- dim(dat)[1]
  fit <- eBayes(lmFit(dat, design))
  logFC <- fit$coefficients # estimate of the log2-fold-change corresponding to the effect size
  reference_condition <- colnames(design)[colSums(design) == nrow(design)]
  reference_averages <- fit$coefficients[,reference_condition]
  if (scale!='log') {  # the inputs are RAW and NOT LOG2 transformed.
    # therefore, coefficients are NOT log fold changes! --> Adjust
    # fitted_response_matrix <- fit$coefficients %*% t(fit$design)
    logFC <- log2((logFC+reference_averages)/reference_averages)
  }
  # independent of `scale`, the reference condition averages should be subtracted to obtain 0 fold change
  logFC[,reference_condition] <- logFC[,reference_condition] - reference_averages
  df.r <- fit$df.residual # residual degrees of freedom assiciated with ordinary t-statistic and p-value
  df.0 <- rep(fit$df.prior, ngenes) # degrees of freedom associated with s2.0
  s2.0 <- rep(fit$s2.prior, ngenes) # estimated prior value for the variance
  s2 <- (fit$sigma)^2 # sample variance
  s2.post <- fit$s2.post # posterior value for the variance
  t.ord <- fit$coefficients / fit$stdev.unscaled / fit$sigma # vector of ordinary t-statistic: using sigma vs s2.post
  t.mod <- fit$t # moderated t-statistic
  p.ord <- 2*pt(-abs(t.ord), fit$df.residual) # ordinary p-value corresonding to the ordinary t-statistic
  p.mod <- fit$p.value # moderated p-value corresonding to the moderated t-statistic
  if(ngenes>1) q.ord <- apply(X = p.ord, MARGIN = 2, FUN = p.adjust, method='BH') # ordinary q-value corresponding to the ordinary t-statistic
  # if(ngenes>1) q.ord <- qvalue(p.ord)$q#, pi0=1)$q # avoid qvalue library when using BH correction
  else q.ord <- p.ord
  if(ngenes>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') # moderated q-value corresponding to the moderated t-statistic
  # if(ngenes>1) q.mod <- qvalue(p.mod)$q#, pi0=1)$q # avoid qvalue library when using BH correction
  else q.mod <- p.mod
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', colnames(logFC))
  colnames(t.ord) <- paste0('t.ord', '_', colnames(t.ord))
  colnames(t.mod) <- paste0('t.mod', '_', colnames(t.mod))
  colnames(p.ord) <- paste0('p.ord', '_', colnames(p.ord))
  colnames(p.mod) <- paste0('p.mod', '_', colnames(p.mod))
  colnames(q.ord) <- paste0('q.ord', '_', colnames(q.ord))
  colnames(q.mod) <- paste0('q.mod', '_', colnames(q.mod))
  results <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  rownames(results) <- rownames(dat)
  return(results)
}
