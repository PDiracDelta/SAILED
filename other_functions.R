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