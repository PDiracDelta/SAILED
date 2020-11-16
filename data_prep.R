library(tidyverse)
library(stringi)

dat.raw <- read.delim('PSMs.csv', sep = '\t')  # create symlink
dat.raw.org <- dat.raw

# read in the study design data frame
study.design=read.delim('msstatstmt_studydesign.csv')  # create symlink

# rename quantification columns
tmp.fun <- function(x){
  stri_replace(x, fixed='Abundance..', replacement='')
}
dat.raw <- dat.raw %>% rename_with(.fn=tmp.fun, .cols=starts_with('Abundance..'))

# rename some more variables
dat.raw <- dat.raw %>% rename(Protein='Protein.Accessions', Peptide='Annotated.Sequence', RT='RT..min.', PTM='Modifications')

# save names of quantification columns for later use
quan.cols <- c(126, paste0(rep(127:130, each=2), c('C','N')), 131)

# focus only on single proteins, not on protein groups (as in the MSstatTMT publication)
dat.raw <- dat.raw %>% filter(X..Proteins==1)

# check for empty Protein variable - there should be none of them
dat.raw %>% filter(Protein=='') %>% nrow

# remove PSM redundancy due to multiple PSM engine score
table(dat.raw$Identifying.Node)
check1=dat.raw %>% group_by(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference...., across(quan.cols)) %>% 
  summarize(ndist=n_distinct(Identifying.Node)) %>% filter(ndist>1)

# construct approximation of Intensity column, which for some reason is missing... but necessary for iPQF
dat.raw <- dat.raw %>% rename(DeltaMZ='Deltam.z..Da.')  %>% mutate(TotalIntensity=rowSums(.[,quan.cols], na.rm = T),)

dat.raw <- dat.raw %>% distinct_at(vars(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference....,quan.cols, TotalIntensity, Ions.Score, DeltaMZ))

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
  summarize(ndist=n_distinct(Peptide), onehit.protein=ifelse(ndist==1, T, F)) %>%
  select(-c(ndist))

# this flag should be calculated regardless of MS run
shared.peptide.df <- dat.raw %>% group_by(Peptide) %>% 
  summarize(ndist=n_distinct(Protein), 
            shared.peptide=ifelse(ndist>1, T, F)) %>%
  select(-c(ndist))

dat.raw <- inner_join(dat.raw, onehit.df, by=c('Protein'))
dat.raw <- inner_join(dat.raw, shared.peptide.df, by=c('Peptide'))

# select only useful columns
dat.raw <- dat.raw %>%
  select(Spectrum.File, Protein, Peptide, RT, Charge, PTM, quan.cols, isoInterOk, noNAs, onehit.protein, shared.peptide, TotalIntensity, Ions.Score, DeltaMZ) 

# turn the dat into long format (useful for modeling) and create Run variable
mix.loc=stri_locate(str=dat.raw$Spectrum.File, regex='Mixture')[1,]
st=mix.loc[1]
ed=stri_locate(str=dat.raw$Spectrum.File, regex='.raw')[1,1]-1

dat.l <- dat.raw %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='Intensity', values_drop_na=FALSE) %>%
  mutate(Run=stri_replace_first(stri_sub(Spectrum.File, from=st, to=ed), '', fixed='0'),
         Mixture=stri_sub(Spectrum.File, from=mix.loc[1], to=mix.loc[2]+1)) %>%
  select(Mixture, Run, Channel, Protein, Peptide, RT, Charge, PTM, Intensity, TotalIntensity, Ions.Score, DeltaMZ, isoInterOk, noNAs, onehit.protein, shared.peptide)

# merge Condition, TechRepMixture, BioReplicate variables from study.design
dat.l <- left_join(dat.l, study.design, by=c('Mixture', 'Run', 'Channel')) %>%
  relocate(TechRepMixture, .after=Mixture) %>%
  relocate(Condition, .after=TechRepMixture) %>%
  relocate(BioReplicate, .after=Condition)

# convert character variables to factors and drop unused factor levels
dat.l <- dat.l %>% mutate(across(c(Mixture:Peptide, Charge, PTM ), .fns=as.factor)) %>% droplevels

# create 'Sample' variable, which is Run by Channel interaction
# dat.l <- dat.l %>% mutate(Sample=Run:Channel) %>% relocate(Sample, .after=Channel)
 
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

# and now return to semi-wide format (wide only within runs)
dat.w <- dat.l %>% pivot_wider(id_cols=-one_of(c('Condition', 'BioReplicate')), names_from=Channel, values_from=Intensity)

# save data in wide and long format
if ('X' %in% colnames(dat.l)) { dat.l$X <- NULL }
saveRDS(list(dat.l=dat.l, dat.w=dat.w), 'input_data.rds')  # make symlink
