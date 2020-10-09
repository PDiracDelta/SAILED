source("http://bioconductor.org/biocLite.R")
biocLite()
library(limma)
library(qvalue)
# read artificial iTRAQ data set
dat <- read.csv("http://www.biostat.jhsph.edu/~kkammers/software/eupa/example_iTRAQ.csv") 
dim(dat)
str(dat)
cha <- c("X113", "X114", "X115", "X116", "X117", "X118", "X119", "X121")
# data preprocessing, load all functions
source("http://www.biostat.jhsph.edu/~kkammers/software/eupa/source.functions.r")
dat <- read.peptides(dat, cha) 
dim(dat) # 10396 out of 10396 peptides left (no missing values in this artificial data set)
# identify proteins from peptide spectra
dat <- quantify.proteins(dat, cha) # 371 proteins identified
# find "one-hit wonders"
dat.onehit <- subset(dat, dat$n.peptides == 1) 
dim(dat.onehit) # 64 proteins are identified by one peptide only
# eliminate "one-hit wonders"
dat <- subset(dat, dat$n.peptides != 1)
dim(dat) # 307 proteins are identified by at least two peptides
# boxplot: intensities of all eight channels after data preprocessing and normalization
par(mfrow=c(1,1), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
boxplot(dat[, 1:length(cha)],  ylim = c(-3, 3), main="Boxplot normalized Intensities")
# define treatment and control groups for two group comparison, assuming 4 cases and 4 controls
tr <- c("X113", "X114", "X115", "X116")
ct <- c("X117", "X118", "X119", "X121")

# define design according to syntax of limma package
design <- model.matrix(~factor(c(2,2,2,2,1,1,1,1)))
design
colnames(design) <- c("Intercept", "Diff")
res.eb <- eb.fit(dat[, c(tr,ct)], design)
head(res.eb)

# volcano plots for ordinary and moderated p-values
rx <- c(-1, 1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, ceiling(max(-log10(res.eb$p.ord), -log10(res.eb$p.mod))))

par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
par(las=1, xaxs="i", yaxs="i")

plot(res.eb$logFC, -log10(res.eb$p.ord), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of ordinary p-values")

plot(res.eb$logFC, -log10(res.eb$p.mod), pch=21, bg="lightgrey", cex=0.9,
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,6,1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of moderated p-values")

# read three simulated iTRAQ experiments
dat1 <- read.csv("http://www.biostat.jhsph.edu/~kkammers/software/eupa/example1_iTRAQ.csv") 
dat2 <- read.csv("http://www.biostat.jhsph.edu/~kkammers/software/eupa/example2_iTRAQ.csv") 
dat3 <- read.csv("http://www.biostat.jhsph.edu/~kkammers/software/eupa/example3_iTRAQ.csv") 

# data preoprocessing as described in the previous section for each experiment separately
cha <- c("X113", "X114", "X115", "X116", "X117", "X118", "X119", "X121")

dat1 <- read.peptides(dat1, cha) 
dat2 <- read.peptides(dat2, cha) 
dat3 <- read.peptides(dat3, cha) 

dat1 <- quantify.proteins(dat1, cha)
dat2 <- quantify.proteins(dat2, cha)
dat3 <- quantify.proteins(dat3, cha)

# eliminate "one-hit wonders"
dat1 <- subset(dat1, dat1$n.peptides != 1)
dat2 <- subset(dat2, dat2$n.peptides != 1)
dat3 <- subset(dat3, dat3$n.peptides != 1)

# limma, factorial design; the first 8 columns of ech data set contain the channel measurements
dat <- data.frame(dat1[, 1:8], dat2[, 1:8], dat3[, 1:8])
tr <- as.factor(rep(c(2,2,2,2,1,1,1,1), 3))
ex <- as.factor(c(rep(1,8), rep(2,8), rep(3,8)))
design <- model.matrix(~ ex + tr)
res.eb.mult <- eb.fit.mult(dat, design)
head(res.eb.mult)

# volcano plots for ordinary and moderated p-values
rx <- c(-1, 1)*max(abs(res.eb$logFC))*1.1
ry <- c(0, ceiling(max(-log10(res.eb$p.ord), -log10(res.eb$p.mod))))

par(mfrow=c(1,2), font.lab=2, cex.lab=1.2, font.axis=2, cex.axis=1.2)
par(las=1, xaxs="i", yaxs="i")

plot(res.eb.mult$logFC, -log10(res.eb.mult$p.ord), pch=21, bg="lightgrey", cex=0.9, 
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,ry[2],1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of ordinary p-values")

plot(res.eb.mult$logFC, -log10(res.eb.mult$p.mod), pch=21, bg="lightgrey", cex=0.9,
     xlim=rx, ylim=ry, xaxt="n",
     xlab="fold change", ylab="-log10  p-value")
abline(v=seq(-2,2,1), col="lightgray", lty="dotted")
abline(h=seq(0,6,1), col="lightgray", lty="dotted")
axis(1, seq(-2,2,1), paste(c("1/4","1/2","1/1","2/1","4/1")))
title("volcano plot of moderated p-values")

# Appendix: R code for source functions
read.peptides <- function(dat, cha){
  output <- NULL
  
  dat$Sequence <- as.character(dat$Sequence)
  dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
  dat$Quan.Usage <- as.character(dat$Quan.Usage)
  dat$Quan.Info <- as.character(dat$Quan.Info)
  dat$Isolation.Interference <- as.numeric(as.character(dat$Isolation.Interference))
  
  dat <- subset(dat, Isolation.Interference<=30)  
  dat <- subset(dat, Quan.Usage=="Used")
  dat <- subset(dat, Protein.Group.Accessions!="")
  dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))  
}

quantify.proteins <- function(dat, cha){
  e.function <- function(x, seq) tapply(x, seq, median)
  output <- NULL
  
  dat$Sequence <- toupper(dat$Sequence) # Capital letters
  accessions <- as.character(unique(dat$Protein.Group.Accessions))
  n.proteins <- length(accessions)
  n.cha <- length(cha)
  
  for(k in 1:n.proteins){
    id <- accessions[k]
    sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
    sdat[cha] <- log2(sdat[cha])
    sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
    pdat <- sdat[, -1]
    n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
    temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])          
    n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)    
    if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
    pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
    output <- rbind(output, pdat)
  }
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
  output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
  output[,1:n.cha] <- round(output[,1:n.cha],3)
  row.names(output) <- accessions
  output <- as.data.frame(output)
  return(output)
}

eb.fit <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coefficients[, 2]
  df.r <- fit.eb$df.residual
  df.0 <- rep(fit.eb$df.prior, n)
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
  t.mod <- fit.eb$t[, 2]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, 2]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
  results.eb <- results.eb[order(results.eb$p.mod), ]
  return(results.eb)
}

eb.fit.mult <- function(dat, design){
  n <- dim(dat)[1]
  fit <- lmFit(dat, design)
  fit.eb <- eBayes(fit)
  logFC <- fit.eb$coef[, "tr2"]
  df.0 <- rep(fit.eb$df.prior, n)
  df.r <- fit.eb$df.residual
  s2.0 <- rep(fit.eb$s2.prior, n)
  s2 <- (fit.eb$sigma)^2
  s2.post <- fit.eb$s2.post
  t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
  t.mod <- fit.eb$t[, "tr2"]
  p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
  p.mod <- fit.eb$p.value[, "tr2"]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  results.eb.mult <- data.frame(logFC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.0, df.r, s2.0, s2, s2.post)
  results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
  return(results.eb.mult)
}