library(tidyverse)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# make empty list with names pre-defined.
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
emptyList <- function(names) {
  return(sapply(names,function(x) NULL))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# interpret factors as a vector of character arrays. Factors are an invention of Satan.
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
remove_factors <- function(x) {
  return(levels(x)[x])
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
    relocate(BioReplicate, .after=Condition)  
    # %>%  select(-X) # is this needed?
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

# dat <- dat.norm.l[[2]]
# mod.formula='response ~ Condition + (1|Run:Channel)'
# referenceCondition

mixed.model.dea <- function(dat, mod.formula, referenceCondition){
  
  mod.formula <- as.formula(mod.formula)
  
  proteins <- dat %>% distinct(Protein) %>% pull(Protein) %>% as.character
  nproteins <- length(proteins)
  
  dat$Condition <- relevel(dat$Condition, referenceCondition)
  condition.levels <- levels(droplevels(dat$Condition))
  n.conditions <- length(condition.levels)
  
  contrast.names <- condition.levels[-match(referenceCondition,condition.levels)]
  
  possib_mod <- possibly(function(x) lmer(mod.formula, data=dat[dat$Protein==x, ], control=lmerControl(check.nobs.vs.nlev="warning", check.nobs.vs.nRE="warning")), otherwise=NULL)
  
  logFC <- matrix(NA, nrow=nproteins, ncol=n.conditions-1)
  t.mod <- p.mod <- logFC
  i=3
  for (i in seq_along(proteins)){
    mod <- possib_mod(proteins[i])
    if (!is.null(mod)){
      sum.mod <- summary(mod)
      logFC[i,]=sum.mod$coefficients[-1,1] # estimate of the log2-fold-change corresponding to the effect size
      t.mod[i,]=sum.mod$coefficients[-1,4] # moderated t-statistic
      p.mod[i,]=sum.mod$coefficients[-1,5] # moderated p-value corresonding to the moderated t-statistic
    } else {
      logFC[i,] <- rep(NA, n.conditions-1)
      t.mod[i,] <- rep(NA, n.conditions-1)
      p.mod[i,] <- rep(NA, n.conditions-1)
    }
  }
  
  if(nproteins>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') else {
    q.mod <- p.mod
  } # moderated q-value corresponding to the moderated t-statistic
  
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', contrast.names)
  colnames(t.mod) <- paste0('t.mod', '_', contrast.names)
  colnames(p.mod) <- paste0('p.mod', '_', contrast.names)
  colnames(q.mod) <- paste0('q.mod', '_', contrast.names)
  results <- data.frame(logFC, t.mod, p.mod, q.mod)
  rownames(results) <- proteins
  
  # results <- results %>% drop_na() # removing proteins with missing values may affect
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
# dat=dat.dea
# score.var='q.mod'
# cut.off <- 0.05

conf.mat <- function(dat, score.var, cut.off, spiked.proteins){
  
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
      pred.class <- factor(ifelse(dat[[j]][, score.vars[i]]<cut.off, 'DE','not_DE'), levels=c('not_DE', 'DE'))
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

# use case (not run)
# cm <- conf.mat(dat.dea[[1]], 'q.mod', contrast.names, 0.05)
# xtable(cm$tab, type='html')
# xtable(cm$metrics, type='html')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# print conf.mat output
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

print.conf.mat <- function(dat){
  
  myHeader1 <- c(1, rep(2, length(variant.names)))
  names(myHeader1) <- c(" ",variant.names)
  
  # dat is a list of size equal to # of contrasts
  # output is presented by contrasts
  for (i in 1:length(dat)){
    
    myHeader2 <- c(1, 2*length(variant.names))
    names(myHeader2) <- c(" ", names(dat)[i])
    
    myHeader3 <- c(1, length(variant.names))
    names(myHeader3) <- c(" ", names(dat)[i])
    
    # print confusion table counts  
    print(
      kable(dat[[i]]$tab) %>%
        kable_styling(bootstrap_options = c("striped", "hover"), full_width=F) %>%
        add_header_above(myHeader1) %>%
        add_header_above(myHeader2)
    )
    
    # print confusion table statistics
    print(
      kable(dat[[i]]$stats, digits=3) %>%
        kable_styling(bootstrap_options = c("striped", "hover"), full_width=F) %>%
        add_header_above(myHeader3)
    )
  }
}

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
    logFC[,reference_condition] <- logFC[,reference_condition] - 1  # correct reference
  } else {   # still have to correct reference
    logFC[,reference_condition] <- logFC[,reference_condition] - reference_averages }
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
  # remove referenceCondition values; they are irrelevant
  return(results %>% select(-contains(reference_condition)))
}
