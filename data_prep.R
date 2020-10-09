library(tidyverse)
library(ggplot2)
library(stringi)

path.data='G:/My Drive/Isobaric labeling strategies/data'
data.raw.fullvar <- read.delim('./PSMs.csv')
data.raw <- read.delim('./PSMs_cleaned.csv')
data.raw <- data.raw.fullvar[, c(colnames(data.raw), 'Identifying.Node', 'Spectrum.File')]
data.raw.org <- data.raw

# read in the study design data frame
study.design=read.delim('./msstatstmt_studydesign.csv')

# rename quantification columns
tmp.fun <- function(x){
  stri_replace(x, fixed='Abundance..', replacement='')
}
data.raw <- data.raw %>% rename_with(.fn=tmp.fun, .cols=starts_with('Abundance..'))

# drop channels with reference samples and rename some more variables
data.raw <- data.raw %>% select(-c('126','131')) %>% 
  rename(Protein=Master.Protein.Accessions, Peptide=Annotated.Sequence, RT=RT..min., PTM=Modifications)

# save names of quantification columns for later use
quan.cols <- paste0(rep(127:130, each=2), c('C','N'))

# remove PSM redundancy due to multiple PSM engine score
table(data.raw$Identifying.Node)
check1=data.raw %>% group_by(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference...., across(quan.cols)) %>% 
  summarize(ndist=n_distinct(Identifying.Node)) %>% filter(ndist>1)

data.raw <- data.raw %>% distinct_at(vars(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference....,quan.cols))

# create flags for PSMs: 
#   having isolation interference <30% or NA [isoInterOk]; 
#   having 0 missing quantification values [noNAs]
#   representing 'one-hit' wonders, i.e., proteins identified by only one peptide
#   representing shared peptides (note that shared peptides will not be removed until summarization and DEA)

data.raw$isoInterOk=ifelse(data.raw$Isolation.Interference....<=30 | is.na(data.raw$Isolation.Interference....), 'Y', 'N')

data.raw$noNAs=ifelse(complete.cases(data.raw[, quan.cols]), 'Y', 'N')

onehit.df <- data.raw %>% group_by(Spectrum.File, Protein) %>% 
  summarize(ndist=n_distinct(Peptide), onehit.protein=ifelse(ndist==1, 'Y', 'N')) %>%
  select(-c(ndist))

shared.peptide.df <- data.raw %>% group_by(Spectrum.File, Peptide) %>% 
  summarize(ndist=n_distinct(Protein), 
            shared.peptide=ifelse(ndist>1, 'Y', 'N')) %>%
  select(-c(ndist))

data.raw <- inner_join(data.raw, onehit.df, by=c('Spectrum.File', 'Protein'))
data.raw <- inner_join(data.raw, shared.peptide.df, by=c('Spectrum.File', 'Peptide'))

# summary of flags
table(data.raw$isoInterOk)
table(data.raw$noNAs)
table(data.raw$onehit.protein)
table(data.raw$shared.peptide)

# select only useful columns
data.raw <- data.raw %>%
  select(Spectrum.File, Protein, Peptide, RT, Charge, PTM, quan.cols, isoInterOk, noNAs, onehit.protein, shared.peptide) 

# turn the data into long format (useful for modeling) and create Run variable
mix.loc=stri_locate(str=data.raw$Spectrum.File, regex='Mixture')[1,]
st=mix.loc[1]
ed=stri_locate(str=data.raw$Spectrum.File, regex='.raw')[1,1]-1

data.long <- data.raw %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='Intensity', values_drop_na=FALSE) %>%
  mutate(Run=stri_replace_first(stri_sub(Spectrum.File, from=st, to=ed), '', fixed='0'),
         Mixture=stri_sub(Spectrum.File, from=mix.loc[1], to=mix.loc[2]+1)) %>%
  select(Mixture, Run, Channel, Protein, Peptide, RT, Charge, PTM, Intensity, isoInterOk, noNAs, onehit.protein, shared.peptide)

# merge Condition and TechRep variables from study.design
data.long <- left_join(data.long, study.design[, c('Run', 'Channel', 'Condition', 'TechRep')], by=c('Run', 'Channel')) %>%
  relocate(TechRep, .after=Mixture) %>%
  relocate(Condition, .after=Channel)

# and then return to wide format (useful for visualization)
data.wide <- data.long %>% pivot_wider(names_from=Channel, values_from=Intensity)

# save untouched data and data in wide and long format
saveRDS(list(data.raw=data.raw.org, data.long=data.long, data.wide=data.wide), 'input_data.rds')
