library(tidyverse)
library(ggplot2)
library(stringi)

path.data='G:/My Drive/Isobaric labeling strategies/data'
dat.raw.fullvar <- read.delim(paste0(path.data, '/PSMs.csv'))
dat.raw <- read.delim(paste0(path.data, '/PSMs_cleaned.csv'))
dat.raw <- dat.raw.fullvar[, c(colnames(dat.raw), 'Identifying.Node', 'Spectrum.File')]
dat.raw.org <- dat.raw

# read in the study design data frame
study.design=read.delim(paste0(path.data, './msstatstmt_studydesign.csv'))

# rename quantification columns
tmp.fun <- function(x){
  stri_replace(x, fixed='Abundance..', replacement='')
}
dat.raw <- dat.raw %>% rename_with(.fn=tmp.fun, .cols=starts_with('Abundance..'))

# drop channels with reference samples and rename some more variables
dat.raw <- dat.raw %>% select(-c('126','131')) %>% 
  rename(Protein=Master.Protein.Accessions, Peptide=Annotated.Sequence, RT=RT..min., PTM=Modifications)

# save names of quantification columns for later use
quan.cols <- paste0(rep(127:130, each=2), c('C','N'))

# remove records with empty Protein column
dat.raw <- dat.raw %>% filter(Protein!='')

# remove PSM redundancy due to multiple PSM engine score
table(dat.raw$Identifying.Node)
check1=dat.raw %>% group_by(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference...., across(quan.cols)) %>% 
  summarize(ndist=n_distinct(Identifying.Node)) %>% filter(ndist>1)

dat.raw <- dat.raw %>% distinct_at(vars(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference....,quan.cols))

# create flags for PSMs: 
#   having isolation interference <30% or NA [isoInterOk]; 
#   having 0 missing quantification values [noNAs]
#   representing 'one-hit' wonders, i.e., proteins identified by only one peptide
#   representing shared peptides (note that shared peptides will not be removed until summarization and DEA)

dat.raw$isoInterOk=ifelse(dat.raw$Isolation.Interference....<=30 | is.na(dat.raw$Isolation.Interference....), 'Y', 'N')

dat.raw$noNAs=ifelse(complete.cases(dat.raw[, quan.cols]), 'Y', 'N')

onehit.df <- dat.raw %>% group_by(Spectrum.File, Protein) %>% 
  summarize(ndist=n_distinct(Peptide), onehit.protein=ifelse(ndist==1, 'Y', 'N')) %>%
  select(-c(ndist))

shared.peptide.df <- dat.raw %>% group_by(Spectrum.File, Peptide) %>% 
  summarize(ndist=n_distinct(Protein), 
            shared.peptide=ifelse(ndist>1, 'Y', 'N')) %>%
  select(-c(ndist))

dat.raw <- inner_join(dat.raw, onehit.df, by=c('Spectrum.File', 'Protein'))
dat.raw <- inner_join(dat.raw, shared.peptide.df, by=c('Spectrum.File', 'Peptide'))

# select only useful columns
dat.raw <- dat.raw %>%
  select(Spectrum.File, Protein, Peptide, RT, Charge, PTM, quan.cols, isoInterOk, noNAs, onehit.protein, shared.peptide) 

# turn the dat into long format (useful for modeling) and create Run variable
mix.loc=stri_locate(str=dat.raw$Spectrum.File, regex='Mixture')[1,]
st=mix.loc[1]
ed=stri_locate(str=dat.raw$Spectrum.File, regex='.raw')[1,1]-1

dat.l <- dat.raw %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='Intensity', values_drop_na=FALSE) %>%
  mutate(Run=stri_replace_first(stri_sub(Spectrum.File, from=st, to=ed), '', fixed='0'),
         Mixture=stri_sub(Spectrum.File, from=mix.loc[1], to=mix.loc[2]+1)) %>%
  select(Mixture, Run, Channel, Protein, Peptide, RT, Charge, PTM, Intensity, isoInterOk, noNAs, onehit.protein, shared.peptide)

# merge Condition, TechRepMixture, BioReplicate variables from study.design
dat.l <- left_join(dat.l, study.design, by=c('Mixture', 'Run', 'Channel')) %>%
  relocate(TechRepMixture, .after=Mixture) %>%
  relocate(Condition, .after=TechRepMixture) %>%
  relocate(BioReplicate, .after=Condition)

# create factors 
dat.l <- dat.l %>% mutate(across(c(Mixture:Peptide, Charge, PTM ), .fns=as.factor))

# and now return to wide format (useful for some visualizations)
dat.w <- dat.l %>% pivot_wider(id_cols=-one_of(c('Condition', 'BioReplicate')), names_from=Channel, values_from=Intensity)

# save data in wide and long format
saveRDS(list(dat.l=dat.l, dat.w=dat.w), paste0(path.data, '/input_data.rds'))
