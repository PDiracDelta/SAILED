library(tidyverse)
library(kableExtra)
library(stringi)
source('manuscript_functions.R')
out_path <-'G:/My Drive/Isobaric labeling strategies/analyses output/' 
referenceCondition <- '0.5'
reportContrast <- '0.125'
# specify colours corresponding to biological conditions
condition.color <- tribble(
  ~Condition, ~Color,
  "0.125", 'black',
  "0.5", 'blue',
  "0.667", 'green',
  "1", 'red' )
data.list <- readRDS('input_data.rds')
dat.l <- data.list$dat.l 
spiked.proteins <- dat.l %>% distinct(Protein) %>% filter(stri_detect(Protein, fixed='ups')) %>% pull %>% as.character
sample.info <- get_sample_info(dat.l, condition.color)
channelNames <- remove_factors(unique(sample.info$Channel))
keep_objects=c('referenceCondition', 'reportContrast', 'condition.color', 'spiked.proteins','sample.info', 'channelNames', 'out_path','op','keep_objects')

### Unit scale ###

# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_unit_outdata.rda')
dea.info <- get_dea_info(dat.dea)
dat.dea$intensity_fix <- NULL; dat.dea$ratio_fix <- NULL;
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='mb_unit_conf.pdf')

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_unit_outdata.rda')
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='dd_unit_conf.pdf')

### Summarization ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_summarization_outdata.rda')
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='mb_summ_conf.pdf')

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_summarization_outdata.rda')
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='dd_summ_conf.pdf')

### Normalization ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_normalization_outdata.rda')
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='mb_norm_conf.pdf')

# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_normalization_outdata.rda')
# names(dat.dea)[names(dat.dea)=='quantile'] <- 'quantile1'
# names(dat.dea)[names(dat.dea)=='quantileByCondition'] <- 'quantile2'
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='dd_norm_conf.pdf')

### DEA ###
# MODEL-BASED
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('modelbased_DEA_outdata.rda')
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='mb_DEA_conf.pdf')
# DATA-DRIVEN
rm(list = ls()[!(ls() %in% keep_objects)])
source('manuscript_functions.R')
load('datadriven_DEA_outdata.rda')
dea.info <- get_dea_info(dat.dea)
cm <- conf_mat(dat.dea, 'q.mod', 0.05, spiked.proteins)
print_conf_mat2(cm, '0.125',  scale_font=999, fileName='dd_DEA_conf.pdf')