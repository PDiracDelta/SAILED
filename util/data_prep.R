library(tidyverse)
library(stringi)
source('util/other_functions.R')

# suffix used when saving the processed data set ('input_data_<data_name>.rds')
dat.raw <- read.delim('data/PSMs.csv', sep = '\t')  # create symlink
dat.raw.org <- dat.raw

# read in the study design data frame
study.design <- read.delim('data/msstatstmt_studydesign.csv', sep=',')  # create symlink

# rename quantification columns
tmp.fun <- function(x){
  stri_replace(x, fixed='Abundance..', replacement='')
}
dat.raw <- dat.raw %>% rename_with(.fn=tmp.fun, .cols=starts_with('Abundance..'))

# drop channels with reference samples and rename some more variables
dat.raw <- dat.raw %>% select(-c('126','131')) %>% 
  rename(Protein='Protein.Accessions', Peptide='Annotated.Sequence', RT='RT..min.', PTM='Modifications')

# save names of quantification columns for later use
quan.cols <- paste0(rep(127:130, each=2), c('C','N'))

# focus only on single proteins, not on protein groups (as in the MSstatTMT publication)
dat.raw <- dat.raw %>% filter(X..Proteins==1)

# check for empty Protein variable - there should be none of them
dat.raw %>% filter(Protein=='') %>% nrow

# create Run variable
mix.loc=stri_locate(str=dat.raw$Spectrum.File, regex='Mixture')[1,]
st=mix.loc[1]
ed=stri_locate(str=dat.raw$Spectrum.File, regex='.raw')[1,1]-1

dat.raw <- dat.raw %>%
  mutate(Run=stri_replace_first(stri_sub(Spectrum.File, from=st, to=ed), '', fixed='0'),
         Mixture=stri_sub(Spectrum.File, from=mix.loc[1], to=mix.loc[2]+1))

# generate noise from Y_mtcb~N(0,sigma.err), where m-mixture; t-techrep; c-condition; b-biological replicate
sigma.err <- 0.2
noise.df <- study.design %>% filter(!(Channel %in% c('126', '131'))) %>% select('Run', 'Channel')
set.seed(1991)
noise.df$err <- 2^rnorm(n=nrow(noise.df), mean=0,sd=sigma.err)
noise.df <- pivot_wider(noise.df, id_cols=c('Run'), names_from=Channel,values_from=err, names_prefix='err' )
dat.raw <- left_join(dat.raw, noise.df, by='Run')
dat.raw[,quan.cols] <- dat.raw[,quan.cols]*dat.raw[,paste0('err',quan.cols)]

# construct approximation of Intensity column, which for some reason is missing... but necessary for iPQF
dat.raw <- dat.raw %>% rename(DeltaMZ='Deltam.z..Da.')  %>% mutate(TotalIntensity=rowSums(.[,quan.cols], na.rm = T),)

# note that for this data only Mascot was used
table(dat.raw$Identifying.Node)

# remove PSM redundancy due to multiple PSM engine score e.g. one peptide with 2 rows with identical RT,
# charge, PTM, and abundance values, but different 'Identifying.Node' and 'DeltaScore' values: Mascot (A2) and Mascot (A4).
# in such case, keep only one record (doesn't matter which one)
dat.raw <- dat.raw %>% distinct_at(vars(Mixture, Run, Protein, Peptide, RT, Charge, PTM, Isolation.Interference....,quan.cols, TotalIntensity, Ions.Score, DeltaMZ))

# create flags for PSMs: 
#   having isolation interference <=30% or NA [isoInterOk]; 
#   having 0 missing quantification values [noNAs]
#   representing 'one-hit' wonders, i.e., proteins identified by only one peptide
#   representing shared peptides (note that shared peptides will not be removed until summarization and DEA)

dat.raw$isoInterOk=ifelse(dat.raw$Isolation.Interference....<=30 | is.na(dat.raw$Isolation.Interference....), T, F)

dat.raw$noNAs=ifelse(complete.cases(dat.raw[, quan.cols]), T, F)

# this flag should be calculated regardless of MS run
# otherwise one protein could be marked as 'one-hit wonder' in one run and as not 'one-hit wonder' in another run
# then if you would remove the one-hit wonder case, you would create more missing values affecting further analyses
onehit.df <- dat.raw %>% group_by(Protein) %>% 
  summarize(ndist=n_distinct(Peptide), onehit.protein=ifelse(ndist==1, T, F)) %>% select(-c(ndist))

# in theory, there should be no shared peptides anymore (after excluding protein groups)
# but let's calculate a flag for this anyway
# this flag should be calculated regardless of MS run
shared.peptide.df <- dat.raw %>% group_by(Peptide) %>% 
  summarize(ndist=n_distinct(Protein), 
            shared.peptide=ifelse(ndist>1, T, F)) %>% select(-c(ndist)) 

dat.raw <- inner_join(dat.raw, onehit.df, by=c('Protein'))
dat.raw <- inner_join(dat.raw, shared.peptide.df, by=c('Peptide'))

# select only useful columns
dat.raw <- dat.raw %>%
  select(Mixture, Run, Protein, Peptide, RT, Charge, PTM, quan.cols, isoInterOk, noNAs, onehit.protein, shared.peptide, TotalIntensity, Ions.Score, DeltaMZ) 

# turn the data into long format (useful for modeling)
dat.l <- dat.raw %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='intensity', values_drop_na=FALSE) %>%
  select(Mixture, Run, Channel, Protein, Peptide, RT, Charge, PTM, intensity, TotalIntensity, Ions.Score, DeltaMZ, isoInterOk, noNAs, onehit.protein, shared.peptide)

# merge Condition, TechRepMixture, BioReplicate variables from study.design
dat.l <- left_join(dat.l, study.design, by=c('Mixture', 'Run', 'Channel')) %>%
  relocate(TechRepMixture, .after=Mixture) %>%
  relocate(Condition, .after=TechRepMixture) %>%
  relocate(BioReplicate, .after=Condition)

# convert character variables to factors and drop unused factor levels
dat.l <- dat.l %>% mutate(across(c(Mixture:Peptide, Charge, PTM ), .fns=as.factor)) %>% droplevels

### finally, remove the proteins overlapped between spiked-in proteins and background proteins
# extract the list of spiked-in proteins
ups.prot <- unique(dat.l[grepl("ups", dat.l$Protein), "Protein"]) %>% pull %>% as.character
# extract the list of background proteins
bg.prot <- unique(dat.l[!grepl("ups", dat.l$Protein), "Protein"]) %>% pull %>% as.character
# overlapped proteins between spiked-in proteins and background proteins
inter <- ups.prot[which(gsub("ups", "", ups.prot) %in% bg.prot)]
# generate the list of proteins to remove
protein.remove <- c(inter, gsub("ups", "", inter))
# remove the overlapped proteins 
dat.l <- dat.l %>% filter(!(Protein %in% protein.remove))

# remove shared peptides
dat.l <- dat.l %>% filter(!shared.peptide)

# keep spectra with (isolation interference <=30 or NA) and no missing quantification channels
dat.l <- dat.l %>% filter(isoInterOk & noNAs)

# remember that you can also exclude proteins with just one peptide (onehit.protein flag), though we decided not to do that

# and now return to semi-wide format (wide only within runs)
dat.w <- dat.l %>% pivot_wider(id_cols=-one_of(c('Condition', 'BioReplicate')), names_from=Channel, values_from=intensity)

# specify parameters used in each notebook:
  referenceCondition <- '0.5'
  # specify colours corresponding to biological conditions
  condition.color <- tribble(
    ~Condition, ~Color,
    "0.125", 'black',
    "0.5", 'blue',
    "0.667", 'green',
    "1", 'red' )
  # quantification channels ordered according to the study design
  channelsOrdered <- c("127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C")
  # specify samples and conditions for ma plots
  ma.onesample.num <- 'Mixture2_1:127C'
  ma.onesample.denom <- 'Mixture1_2:129N'
  ma.allsamples.num <- '0.5'
  ma.allsamples.denom <- '0.125'
params <- list(referenceCondition=referenceCondition,
               condition.color=condition.color, 
               channelsOrdered=channelsOrdered, 
               ma.onesample.num=ma.onesample.num, 
               ma.onesample.denom=ma.onesample.denom, 
               ma.allsamples.num=ma.allsamples.num, 
               ma.allsamples.denom=ma.allsamples.denom)

# save data in wide and long format
if ('X' %in% colnames(dat.l)) { dat.l$X <- NULL }
saveRDS(list(dat.l=dat.l, dat.w=dat.w, data.params=params), paste0('data/input_data', '.rds'))  # make symlink
