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
  if (is.factor(x)) return(levels(x)[x]) else return(x)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# # create data frame with sample info (distinct Run,Channel, Sample, Condition, Colour)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

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
# return a character vector with peptide sequences identified in every MS run
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

get_inner_peptides <- function(dat){
  unique.pep=dat %>% 
    group_by(Run) %>%
    distinct(Peptide) %>% 
    mutate(val=1)
  unique.pep <- xtabs(val~Peptide+Run, data=unique.pep)
  tmp <- apply(unique.pep, 1, function(x) all(x==1))
  return(rownames(unique.pep)[tmp])
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# convert data to long format (assume Run, Mixture columns already present)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

to_long_format<-function(x, sample.info, merge_study_design=T) {
  quan.cols = unique(sample.info$Channel)
  x <- x %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='response', values_drop_na=FALSE)
  
  # merge Condition from sample.info
  # remove cols to add if already present to avoid errors.
  if (merge_study_design) {
    x <- x[,!(colnames(x) %in% c('Condition'))]
    x <- left_join(x, sample.info, by=c('Run', 'Channel')) }
  return(x)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function to switch name order of quantification columns,e.g.
# 127C:Mixture1_1 changed to Mixture1_1:127C, or vice-versa
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

flip_colnames=function(dat, fixed.cols, sep=':'){
  col.dat <- colnames(dat)
  ind <- !(col.dat %in% fixed.cols)
  x <- stri_split(col.dat[ind], fixed=sep)
  unlist(lapply(x, function(y) paste(y[2], y[1], sep=sep)))  
  colnames(dat)[ind]=unlist(lapply(x, function(y) paste(y[2], y[1], sep=sep)))  
  return(dat)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function performing mean or median aggregation on variables specified in var.names argument
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

aggFunc=function(dat, var.names, group.vars, agg.method='mean'){
  select.method=match.arg(agg.method, c('mean', 'median', 'sum'))
  
  dat2 <- lazy_dt(dat)
  out.dat=dat %>%
    group_by(across(all_of(group.vars))) %>%
    summarize_at(var.names, eval(parse(text=select.method))) %>%
    as_tibble()
    
  return(out.dat)
}

median_sweep <- function(dat, margin, fun){
  return(sweep(dat, margin, apply(dat[,cols], margin, median, na.rm=T), FUN=fun))}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function for mixed models DEA (without empirical bayes moderation)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

lmm_dea <- function(dat, mod.formula, referenceCondition, scale){
  
  mod.formula <- as.formula(mod.formula)
  
  proteins <- dat %>% distinct(Protein) %>% pull(Protein) %>% as.character
  nproteins <- length(proteins)
  
  dat$Condition <- relevel(dat$Condition, referenceCondition)
  condition.levels <- levels(droplevels(dat$Condition))
  n.conditions <- length(condition.levels)
  
  contrast.names <- condition.levels[-match(referenceCondition,condition.levels)]
  
  # wrap up lmer call with 'possibly' function to continue the for loop despite errors
  lmer_safe <- possibly(function(x) lmer(mod.formula, data=dat[dat$Protein==x, ],
                                         control=lmerControl(check.nobs.vs.nlev="warning", check.nobs.vs.nRE="warning", check.nlev.gtr.1 = "warning")), 
                        otherwise=NULL)
  logFC <- matrix(NA, nrow=nproteins, ncol=n.conditions-1)
  t.mod <- p.mod <- logFC
  n.obs <- s2.RE <- s2.resid <- matrix(NA, nrow=nproteins, ncol=1)
  
  for (i in seq_along(proteins)){
    mod <- lmer_safe(proteins[i])
    if (!is.null(mod)){
      sum.mod <- summary(mod)
      # estimate of the log2-fold-change corresponding to the effect size
      if (scale=='log') logFC[i,]=sum.mod$coefficients[-1,1] else if (scale=='raw'){
        rat <- (sum.mod$coefficients[-1,1]+sum.mod$coefficients[1,1])/sum.mod$coefficients[1,1]
        logFC[i, rat>0]=log2(rat[rat>0]); logFC[i, rat<=0] <- NA # some of the calculated ratios may be negative!
        }
      t.mod[i,]=sum.mod$coefficients[-1,4] # t-statistic
      p.mod[i,]=sum.mod$coefficients[-1,5] # p-value
      var.mod <- as.data.frame(VarCorr(mod))
      s2.resid[i,] <- var.mod[var.mod$grp=='Residual','vcov']
      s2.RE[i,] <- var.mod[var.mod$grp!='Residual','vcov']
      n.obs[i,] <- nrow(mod@frame)
      
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
  results <- data.frame(logFC, t.mod, p.mod, q.mod, n.obs=n.obs, s2.RE, s2.resid)
  rownames(results) <- proteins
  
  # results <- results %>% drop_na() # removing proteins with missing values may affect
  return(results)
} 

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# function for printing # of up/down/not regulated proteins
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

regulated_proteins <- function(dea.mat, score.var, conditions, cut.off){
  cat('numerator condition:', conditions[1],', ' , 'denominator condition:', conditions[2], '\n')
  c(`Up`=sum(dea.mat[, score.var]<cut.off & dea.mat[, 'logFC']>0),
    `Down`=sum(dea.mat[, score.var]<cut.off & dea.mat[, 'logFC']<0),
    `Not Sig`=sum(dea.mat[, score.var]>cut.off))
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

print_conf_mat <- function(dat, referenceCondition){
  variant.names <- colnames(dat[[1]]$stats)
  myHeaderVariant <- c(1, rep(2, length(variant.names)))
  names(myHeaderVariant) <- c(" ",variant.names)
  # dat is a list of size equal to # of contrasts
  # output is presented by contrasts
  
  tabs <- vector('list', length(dat))
  
  for (i in 1:length(dat)){
    myHeader2 <- c(1, 2*length(variant.names))
    names(myHeader2) <- c(" ", paste0(names(dat)[i], ' vs ', referenceCondition, ' contrast'))
    myHeader3 <- c(1, length(variant.names))
    names(myHeader3) <- c(" ", paste0(names(dat)[i], ' vs ', referenceCondition, ' contrast'))
    # print confusion table counts  
    k1 <- kable(dat[[i]]$tab) %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width=F) %>%
      add_header_above(myHeaderVariant) %>%
      add_header_above(myHeader2)
    # print confusion table statistics
    k2 <- kable(dat[[i]]$stats, digits=3) %>%
      kable_styling(bootstrap_options = c("striped", "hover"), full_width=F) %>%
      add_header_above(myHeader3) 
    tabs[[i]] <- list(htmltools::HTML(k1), htmltools::HTML(k2))
  }
  return(htmltools::tagList(tabs))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# get design matrix for use in moderated t-test
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

get_design_matrix <- function(referenceCondition, sample.info) {
  # ANOVA-like design matrix for use in moderated_ttest, indicating group 
  # (condition) membership of each entry in all_channels.
  otherConditions = setdiff(unique(sample.info$Condition), referenceCondition)
  all_channels = sample.info$Sample
  N_channels = length(all_channels) 
  N_conditions = 1+length(otherConditions)
  design = matrix(rep(0,N_channels*N_conditions), c(N_channels, N_conditions))
  design[, 1] <- 1  # reference gets 1 everywhere
  rownames(design) <- all_channels
  colnames(design) <- c(referenceCondition, otherConditions)
  for (i in 1:N_channels) {  # for each channel in each condition, put a "1" in the design matrix.
    design[sample.info$Sample[i], sample.info$Condition[i]] <- 1 }
  return(design)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# moderated t-test for data-driven approach
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

moderated_ttest <- function(dat, design, scale) {
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

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Wilcoxon rank-sum test
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Wilcoxon rank-sum test (equivalent to Mann–Whitney U test) comparing outcomes between 
# two independent groups of samples. For each protein separately, the test is applied to each otherCondition w.r.t. referenceCondition.
# Wilcoxon rank-sum test employs ranks of quantification values, therefore logFC must be computed manually:
#   if logFC.method='diff': difference of arithmetic means computed on log2 scale
#   if logFC.method='ratio': log2 ratio of arithmetic means computed on raw scale

wilcoxon_test <- function(dat, sample.info, referenceCondition, otherConditions, logFC.method='ratio'){
  proteins <- dat %>% distinct(Protein) %>% pull(Protein) %>% as.character
  nproteins <- length(proteins)
  n.conditions <- length(otherConditions)
  refCondCols <- sample.info %>% filter(Condition==referenceCondition) %>% 
    distinct(Run, Channel) %>% mutate(sample=paste(Run, Channel, sep=':')) %>% pull(sample)
  
  logFC <- matrix(NA, nrow=nproteins, ncol=n.conditions)
  t.mod <- p.mod <- logFC
  
  for (i in 1:n.conditions){
    CondCols <- sample.info %>% filter(Condition==otherConditions[i]) %>% 
      distinct(Run, Channel) %>% mutate(sample=paste(Run, Channel, sep=':')) %>% pull(sample)
    
    wilcox.results=row_wilcoxon_twosample(x=dat[,CondCols], y=dat[,refCondCols])
    
    if (logFC.method=='diff') {logFC[,i] <- rowMeans(dat[,CondCols], na.rm = T)-rowMeans(dat[,refCondCols], na.rm = T)} 
    else if (logFC.method=='ratio') {
      num <- rowMeans(2^dat[,CondCols], na.rm = T)
      denom <- rowMeans(2^dat[,refCondCols], na.rm = T)
      logFC[,i] <- ifelse(denom!=0, log2(num/denom), NA) 
    }
    t.mod[,i] <- wilcox.results$statistic
    p.mod[,i] <- wilcox.results$pvalue
  }
  
  if(nproteins>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') else {
    q.mod <- p.mod
  } # moderated q-value corresponding to the moderated t-statistic
  
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', otherConditions)
  colnames(t.mod) <- paste0('t.mod', '_', otherConditions)
  colnames(p.mod) <- paste0('p.mod', '_', otherConditions)
  colnames(q.mod) <- paste0('q.mod', '_', otherConditions)
  results <- data.frame(logFC, t.mod, p.mod, q.mod)
  rownames(results) <- proteins
  return(results)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# permutation test
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Two-group protein-by-protein permutation test based on the difference in arithmetic means. 
# Columns of the m x n quantification matrix (m proteins and n biological samples) 
# are B times randomly shuffled (drawing n column indices without replacement) in order to generate the null distribution 
# of the test statistic. P-values are computed by comparing the observed test statistic value with the null distribution.

permutation_test_manual <- function(dat, sample.info, referenceCondition, otherConditions, B=1000, seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  proteins <- dat %>% distinct(Protein) %>% pull(Protein) %>% as.character
  nproteins <- length(proteins)
  n.conditions <- length(otherConditions)
  #allCols <- colnames(dat)[-1]
  refCondCols <- sample.info %>% filter(Condition==referenceCondition) %>% 
    distinct(Run, Channel) %>% mutate(sample=paste(Run, Channel, sep=':')) %>% pull(sample)
  refCondColsLen <- length(refCondCols)
  
  logFC <- matrix(NA, nrow=nproteins, ncol=n.conditions)
  t.mod <- p.mod <- logFC

  for (i in 1:n.conditions){
    
    CondCols <- sample.info %>% filter(Condition==otherConditions[i]) %>% 
      distinct(Run, Channel) %>% mutate(sample=paste(Run, Channel, sep=':')) %>% pull(sample)
    # the observed test statistic value (difference in arithmetic means)
    obsStat <- rowMeans(dat[,CondCols], na.rm = T)-rowMeans(dat[,refCondCols], na.rm = T)
    
    permStat <- matrix(NA, nrow=nproteins, ncol=B)
    for (b in 1:B){
      shuffledChannels <- sample(c(refCondCols, CondCols))
      datRef <- dat[, shuffledChannels[1:refCondColsLen]]
      datCond <- dat[, shuffledChannels[(refCondColsLen+1):length(shuffledChannels)]]
      # test statistic value (difference in arithmetic means) computed on permuted data
      permStat[,b] <- rowMeans(datCond, na.rm=TRUE)-rowMeans(datRef, na.rm=TRUE)
    }
    
    logFC[,i] <- obsStat 
    t.mod[,i] <- NA
    sum.tmp <- apply(cbind(abs(obsStat), abs(permStat)), 1, function(x) sum(x[-1]>x[1], na.rm=T))
    p.mod[,i] <- (1+sum.tmp)/(1+B)
  }
  
  if(nproteins>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') else {
    q.mod <- p.mod
  } # moderated q-value corresponding to the moderated t-statistic
  
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', otherConditions)
  colnames(t.mod) <- paste0('t.mod', '_', otherConditions)
  colnames(p.mod) <- paste0('p.mod', '_', otherConditions)
  colnames(q.mod) <- paste0('q.mod', '_', otherConditions)
  results <- data.frame(logFC, t.mod, p.mod, q.mod)
  rownames(results) <- proteins
  return(results)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Fisher-Pitman permutation test from 'coin' package.
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# dat should be data in long format
# ... pass arguments to oneway_test function from 'coin' package. this will be mainly used
# to specify how to compute p-values:
# distribution='exact'
# distribution=approximate(nresample=10000)
# distribution='asymptotic'

permutation_test <- function(dat, referenceCondition, otherConditions, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  proteins <- dat %>% distinct(Protein) %>% pull(Protein) %>% as.character
  nproteins <- length(proteins)
  n.conditions <- length(otherConditions)
  logFC <- matrix(NA, nrow=nproteins, ncol=n.conditions); rownames(logFC)=proteins
  p.mod <- logFC
  dat$Condition <- relevel(as.factor(dat$Condition), referenceCondition)
  dat$Protein <- remove_factors(dat$Protein)
  
  for (i in 1:n.conditions){
    sdat <- dat %>% filter(Condition %in% c(referenceCondition, otherConditions[i]))  
    sdat.split <- split(sdat, sdat$Protein)  
    
    p.tmp <- sapply(sdat.split, function(x){
      return(pvalue(oneway_test(response~Condition, data=x, ...)))
    })
    p.mod[names(p.tmp),i] <- p.tmp
    logFC.tmp <- sapply(sdat.split, function(x){
      a <- mean(x$response[x$Condition==referenceCondition], na.rm=TRUE)
      b <- mean(x$response[x$Condition==otherConditions[i]], na.rm=TRUE)
      return(b-a)})
    logFC[names(logFC.tmp),i] <- logFC.tmp
  }
  if(nproteins>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') else {
    q.mod <- p.mod
  } # moderated q-value corresponding to the moderated t-statistic
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', otherConditions)
  colnames(p.mod) <- paste0('p.mod', '_', otherConditions)
  colnames(q.mod) <- paste0('q.mod', '_', otherConditions)
  results <- data.frame(logFC, p.mod, q.mod)
  rownames(results) <- proteins
  return(results)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# ROTS
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

rots_test <- function(dat, referenceCondition, otherConditions, seed=NULL, ...){
  if (!is.null(seed)) set.seed(seed)
  proteins <- dat %>% distinct(Protein) %>% pull(Protein) %>% as.character
  nproteins <- length(proteins)
  n.conditions <- length(otherConditions)
  refCondCols <- sample.info %>% filter(Condition==referenceCondition) %>% 
    distinct(Run, Channel) %>% mutate(sample=paste(Run, Channel, sep=':')) %>% pull(sample)
  refCondColsLen <- length(refCondCols)
  
  logFC <- matrix(NA, nrow=nproteins, ncol=n.conditions)
  t.mod <- p.mod <- logFC
  
  for (i in 1:n.conditions){
    CondCols <- sample.info %>% filter(Condition==otherConditions[i]) %>% 
      distinct(Run, Channel) %>% mutate(sample=paste(Run, Channel, sep=':')) %>% pull(sample)
    CondColsLen <- length(CondCols)
    
    results = ROTS(dat[,c(refCondCols,CondCols)] 
                   # the reference group should get smaller group indicator value
                  ,groups = c(rep(1, refCondColsLen), rep(0, CondColsLen)), seed = seed, ...)
    
    logFC[,i] <- results$logfc
    t.mod[,i] <- results$d
    p.mod[,i] <- results$pvalue
  }
  
  if(nproteins>1) q.mod <- apply(X = p.mod, MARGIN = 2, FUN = p.adjust, method='BH') else {
    q.mod <- p.mod
  } # moderated q-value corresponding to the moderated t-statistic
  
  # incorporate entity type into colnames to overwrite identical factor names
  colnames(logFC) <- paste0('logFC_', otherConditions)
  colnames(t.mod) <- paste0('t.mod', '_', otherConditions)
  colnames(p.mod) <- paste0('p.mod', '_', otherConditions)
  colnames(q.mod) <- paste0('q.mod', '_', otherConditions)
  results <- data.frame(logFC, t.mod, p.mod, q.mod)
  rownames(results) <- proteins
  return(results)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# DEqMS
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# wrapper function on the DEqMS method from the 'DEqMS' package.
# psm.counts should be a dataframe with two columns: the first one containing protein names and 
# the second one containing corresponding PSM counts

deqms_test <- function(dat, design, scale, psm.counts) {
  design <- design[match(colnames(dat), rownames(design)),]  # fix column order
  ngenes <- dim(dat)[1]
  fit <- eBayes(lmFit(dat, design))
  ord=match(rownames(fit$coefficients), psm.counts[,1])
  fit$count = psm.counts[ord,2]
  fit = spectraCounteBayes(fit)
  
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
  df.0 <- rep(fit$sca.dfprior, ngenes) # degrees of freedom associated with s2.0
  s2.0 <- fit$sca.priorvar # estimated prior value for the variance
  s2 <- (fit$sigma)^2 # sample variance
  s2.post <- fit$sca.postvar	# posterior value for the variance
  PSMcount <- fit$count
  t.ord <- fit$coefficients / fit$stdev.unscaled / fit$sigma # vector of ordinary t-statistic: using sigma vs s2.post
  t.mod <- fit$sca.t	 # moderated t-statistic
  p.ord <- 2*pt(-abs(t.ord), fit$df.residual) # ordinary p-value corresonding to the ordinary t-statistic
  p.mod <- fit$sca.p # moderated p-value corresonding to the moderated t-statistic
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
  results <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post, PSMcount)
  rownames(results) <- rownames(dat)
  # remove referenceCondition values; they are irrelevant
  return(results %>% select(-contains(reference_condition)))
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# get_anova
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# extracts ordinary linear model output from the moderated_ttest call and renames all columns 
# containing .ord to .mod (for consistency - .mod columns are used in confusion matrix, scatter plots, 
# volcano plots, violin plots)

get_anova <- function(dat, design, scale){
  mod <- moderated_ttest(dat, design, scale)
  allcols <- colnames(mod)
  ind <- stri_detect(allcols, fixed='.mod')
  # select columns corresponding to ordinary linear model only
  mod <- mod[,!ind & !(colnames(mod) %in% c("df.r", "df.0", "s2.0", "s2.post"))]
  # rename '.ord' to '.mod'
  colnames(mod) <- stri_replace(colnames(mod), fixed='.ord', replacement='.mod')
  return(mod)
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# run_test - one-way ANOVA for testing the run effect
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# dat - data in long format

run_test <- function(dat){
  unique.proteins <- dat %>% group_by(Protein) %>% summarise(ndist=n_distinct(Run)) %>% filter(ndist>1) %>% pull(Protein) %>% as.character
  dat2 <- dat %>% filter(Protein %in% unique.proteins)
  pvalues=sapply(split(dat2, as.character(dat2$Protein)), function(x) summary(aov(x$response~x$Run, data=x))[[1]]$`Pr(>F)`[1])
  pvalues=as.data.frame(pvalues)
  return(pvalues)
}
