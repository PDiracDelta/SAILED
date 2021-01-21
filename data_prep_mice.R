library(tidyverse)
library(stringi)

tmt1=read.csv('B_a.txt',sep='\t', stringsAsFactors=FALSE); tmt1$Spectrum.File='Mixture1'
tmt2=read.csv('G_a.txt',sep='\t', stringsAsFactors=FALSE); tmt2$Spectrum.File='Mixture2'
tmt3=read.csv('R_a.txt',sep='\t', stringsAsFactors=FALSE); tmt3$Spectrum.File='Mixture3'
dat.raw <- rbind(tmt1,tmt2,tmt3)

#BGR -run123

# study design
des.mat <- matrix(c('BM3','PM3','TAM4','BM4','TAM3','PM4',
                    'PM5','TAM6','BM5','TAM5','PM6','BM6',
                    'TAM1','BM1','PM1','PM2','BM2','TAM2'),
                  byrow=TRUE,ncol=6)
sample=as.vector(t(des.mat))
study.design <- data.frame( 
  TechRepMixture='1',
  Channel=as.character(rep(126:(126+ncol(des.mat)-1))),
  Condition=case_when(
    str_detect(sample,'BM') ~ 'BM',
    str_detect(sample,'PM') ~ 'PM',
    str_detect(sample,'TAM') ~ 'TAM'))
study.design <- study.design %>% mutate(
  Mixture=paste0('Mixture',rep(1:3, each=ncol(des.mat))),
  BioReplicate=paste(Mixture,Condition,sep='_'),
  Run=paste(Mixture,rep(1:3, each=ncol(des.mat)), sep='_'))

###############################
# dat.raw <- read.delim('PSMs.csv', sep = '\t')  # create symlink
# dat.raw.org <- dat.raw

# read in the study design data frame
#study.design=read.delim('msstatstmt_studydesign.csv')  # create symlink

# rename quantification columns
colnames(dat.raw)[colnames(dat.raw) %in% paste0('X',126:131)] <- 126:131

# rename some more variables
dat.raw <- dat.raw %>% 
  rename(Protein='Protein.Accessions', Peptide='Annotated.Sequence', RT='RT..min.', PTM='Modifications')

# save names of quantification columns for later use
quan.cols <- as.character(126:131)

# focus only on single proteins, not on protein groups (as in the MSstatTMT publication)
dat.raw <- dat.raw %>% filter(X..Proteins==1)

# check for empty Protein variable - there should be none of them
dat.raw %>% filter(Protein=='') %>% nrow

# remove PSM redundancy due to multiple PSM engine score
table(dat.raw$Identifying.Node)
# First Mascot, then Sequest
tmp1=dat.raw %>%
  group_by(Spectrum.File, Protein, Peptide, RT, Charge, PTM) %>%
  summarize(ndist=n()) %>%
  filter(ndist>1) %>%
  select (-c(ndist))
tmp2=dat.raw %>%
  inner_join(tmp1,by=c('Spectrum.File', 'Protein', 'Peptide', 'RT', 'Charge', 'PTM')) %>%
  filter(Identifying.Node.Type=='Mascot')
tmp3=dat.raw %>%
  anti_join(tmp1,by=c('Spectrum.File', 'Protein', 'Peptide', 'RT', 'Charge', 'PTM'))
dat.raw=union(tmp2,tmp3)

check1=dat.raw %>% group_by(Spectrum.File, Protein, Peptide, RT, Charge, PTM, Isolation.Interference...., across(quan.cols)) %>% 
  summarize(ndist=n_distinct(Identifying.Node)) %>% filter(ndist>1)

# construct approximation of Intensity column, which for some reason is missing... but necessary for iPQF
#dat.raw <- dat.raw %>% rename(DeltaMZ='Deltam.z..Da.')  %>% mutate(TotalIntensity=rowSums(.[,quan.cols], na.rm = T),)
dat.raw <- dat.raw %>% rename(DeltaMZ='Deltam.z..Da.', TotalIntensity='Intensity')

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

# turn the data into long format (useful for modeling) and create Run variable
mix.loc=stri_locate(str=dat.raw$Spectrum.File, regex='Mixture')[1,]
st=mix.loc[1]
ed=stri_locate(str=dat.raw$Spectrum.File, regex='.raw')[1,1]-1

dat.l <- dat.raw %>% pivot_longer(cols=quan.cols, names_to='Channel', values_to='intensity', values_drop_na=FALSE) %>%
  mutate(Run=paste(Spectrum.File,'1',sep='_'),
         Mixture=Spectrum.File) %>%
  select(Mixture, Run, Channel, Protein, Peptide, RT, Charge, PTM, intensity, TotalIntensity, Ions.Score, DeltaMZ, isoInterOk, noNAs, onehit.protein, shared.peptide)

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
#ups.prot <- unique(dat.l[grepl("ups", dat.l$Protein), "Protein"]) %>% pull %>% as.character
# extract the list of background proteins
#bg.prot <- unique(dat.l[!grepl("ups", dat.l$Protein), "Protein"]) %>% pull %>% as.character
# overlapped proteins between spiked-in proteins and background proteins
#inter <- ups.prot[which(gsub("ups", "", ups.prot) %in% bg.prot)]
# generate the list of proteins to remove
#protein.remove <- c(inter, gsub("ups", "", inter))
# remove the overlapped proteins 
#dat.l <- dat.l %>% filter(!(Protein %in% protein.remove))

# remove shared peptides
dat.l <- dat.l %>% filter(!shared.peptide)

# keep spectra with (isolation interference <=30 or NA) and no missing quantification channels
dat.l <- dat.l %>% filter(isoInterOk & noNAs)

# and now return to semi-wide format (wide only within runs)
dat.w <- dat.l %>% pivot_wider(id_cols=-one_of(c('Condition', 'BioReplicate')), names_from=Channel, values_from=intensity)

# save data in wide and long format
if ('X' %in% colnames(dat.l)) { dat.l$X <- NULL }
saveRDS(list(dat.l=dat.l, dat.w=dat.w), 'input_data_mice.rds')  # make symlink
