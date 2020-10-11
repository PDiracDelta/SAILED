# create the study design data frame (only Mixture 1,2 and Technical Replicates 1,2)
library(MSstatsTMT)
path.data <- 'G:/My Drive/Isobaric labeling strategies/data'
study.design <- annotation.pd %>%
  filter(Mixture %in% c('Mixture1', 'Mixture2') & TechRepMixture %in% c(1,2)) %>%
  mutate(Run=paste(Mixture, TechRepMixture, sep='_')) %>%
  select(-Fraction)
write.table(study.design, paste0(path.data, '/msstatstmt_studydesign.csv'), sep='\t')
