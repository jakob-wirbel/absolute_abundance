# ##############################################################################
#
## Absolute abundance paper code
#
# ##############################################################################

library("tidyverse")
library("mlr3")
library("mlr3learners")
library("mlr3tuning")

set.seed(2024)

# ##############################################################################
# data exploration

# read in data
df.data <- read_tsv('./data/alloHCT_data.tsv', 
                    col_types = cols())

# look at correlation
df.data %>% 
  select_if(is.numeric) %>% 
  pivot_longer(-log_copies_per_extraction) %>% 
  ggplot(aes(y=log_copies_per_extraction, x=value)) + 
  geom_point() + 
  facet_wrap(~name, scales = 'free_x')

g <- df.data %>% 
  select_if(is.numeric) %>% 
  pivot_longer(-log_copies_per_extraction) %>% 
  filter(str_detect(name, 'k__')) %>% 
  ggplot(aes(y=log_copies_per_extraction, x=value)) + 
    geom_point() + 
    facet_wrap(~name, scales = 'free_x') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab('Relative abundance') + ylab('16S copies') + 
    geom_smooth(formula='y~x', method='lm')
ggsave(g, filename='./figures/correlation_taxa.pdf',
       width = 8, height = 4, useDingbats=FALSE)  

df.data %>% 
  select_if(is.numeric) %>% 
  pivot_longer(-log_copies_per_extraction) %>% 
  group_by(name) %>% 
  summarise(c=cor(log_copies_per_extraction, value, method='spearman'),
            p=cor.test(log_copies_per_extraction, value, 
                       method='spearman', exact = FALSE)$p.value)
# # A tibble: 7 × 3
# name                    c         p
# <chr>               <dbl>     <dbl>
# 1 DNA_concentration  0.916  1.41e-206
# 2 alpha              0.339  2.10e- 15
# 3 bac_frac           0.267  6.56e- 10
# 4 k__Archaea        -0.0406 3.56e-  1
# 5 k__Bacteria        0.0387 3.79e-  1
# 6 k__Eukaryota      -0.171  8.82e-  5
# 7 unclassified       0.0314 4.76e-  1

g <- df.data %>% 
  ggplot(aes(x=DNA_concentration, y=log_copies_per_extraction)) + 
  geom_point() + 
  geom_smooth(method='lm', formula='y~log(x)') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  xlab('DNA concentration [ng/µL]') + 
  ylab('log10(copies per extraction)')
ggsave(g, filename='./figures/correlation_DNA.pdf',
       width=4, height = 4, useDingbats=FALSE)
g <- df.data %>% 
  ggplot(aes(x=alpha, y=log_copies_per_extraction)) + 
  geom_point() + 
  geom_smooth(method='lm', formula='y~x') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  xlab('Alpha diversity [Shannon]') + 
  ylab('log10(copies per extraction)')
ggsave(g, filename='./figures/correlation_alpha.pdf',
       width=4, height = 4, useDingbats=FALSE)

# ##############################################################################
# train model

# trained the original models with `DNA_conc` instead of `DNA_concentration` 
# as the column name, so for comparability, we will have to rename
colnames(df.data)[2] <- 'DNA_conc'

# Full model
learner <- lrn("regr.ranger", mtry= to_tune(1,7),
               num.trees=to_tune(50, 500), 
               importance='impurity')

task <- as_task_regr(df.data %>% 
                       select(-c(Sample_ID, Plate)), 
                     target = "log_copies_per_extraction")
tnr_grid_search = tnr("grid_search", resolution = 5, batch_size = 5)
at = auto_tuner(
  tuner = tnr_grid_search,
  learner = learner,
  resampling = rsmp("cv", folds = 3),
  measure = msr('regr.rsq'),
)
rsmp_cv3 <- rsmp("repeated_cv", repeats = 10, folds = 10)

a <- Sys.time()
rr_full <- invisible(resample(task, at, rsmp_cv3, store_models = TRUE))
Sys.time() - a
# Time difference of 12.90 mins

# DNA-only model
learner <- lrn("regr.ranger", mtry= 1,
               num.trees=to_tune(50, 500), 
               importance='impurity')
task <- as_task_regr(df.data %>% 
                       select(log_copies_per_extraction, DNA_conc),
                     target = "log_copies_per_extraction")
tnr_grid_search = tnr("grid_search", resolution = 5, batch_size = 5)
at = auto_tuner(
  tuner = tnr_grid_search,
  learner = learner,
  resampling = rsmp("cv", folds = 3),
  measure = msr('regr.rsq'),
)
rsmp_cv3 <- rsmp("repeated_cv", repeats = 10, folds = 10)

a <- Sys.time()
rr_dna <- invisible(resample(task, at, rsmp_cv3, store_models = TRUE))
Sys.time() - a
# Time difference of 1.970487 mins

## save
save(rr_dna, file='./models/dna_model.Rdata')
save(rr_full, file='./models/full_model.Rdata')

# get predictions for evaluation
rr_full$aggregate(measures = c(msr('regr.mse'), msr('regr.mae'), 
                               msr('regr.rsq'), msr('regr.sae')))
# regr.mse    regr.mae    regr.rsq    regr.sae 
# 0.08270097  0.20571318  0.86123495 10.65579589 
rr_dna$aggregate(measures = c(msr('regr.mse'), msr('regr.mae'), 
                              msr('regr.rsq'), msr('regr.sae')))
# regr.mse   regr.mae   regr.rsq   regr.sae 
# 0.1039912  0.2307719  0.8229848 11.9535049 

df.cor.full <- tibble(id=rr_full$prediction()$row_ids,
                      truth=rr_full$prediction()$truth,
                      response=rr_full$prediction()$response) %>% 
  group_by(id) %>% 
  summarise(truth=unique(truth), response=mean(response))
cor(df.cor.full$truth, df.cor.full$response)                        # 0.9355706
cor(df.cor.full$truth, df.cor.full$response, method = 'spearman')   # 0.9136876
x <- epiR::epi.ccc(df.cor.full$truth, df.cor.full$response)
x$rho.c                                                             # 0.9330186 0.921178 0.9431332
g.full <- df.cor.full %>% 
  ggplot(aes(x=truth, y=response)) + geom_point() + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_bw()
ggsave(g.full, filename='./figures/full_model_predictions.pdf',
       width = 4, height = 4, useDingbats=FALSE)

# dna only
df.cor.dna <- tibble(id=rr_dna$prediction()$row_ids,
                     truth=rr_dna$prediction()$truth,
                     response=rr_dna$prediction()$response) %>% 
  group_by(id) %>% 
  summarise(truth=unique(truth), response=mean(response))
cor(df.cor.dna$truth, df.cor.dna$response)                          # 0.9174385
cor(df.cor.dna$truth, df.cor.dna$response, method = 'spearman')     # 0.8898321
x <- epiR::epi.ccc(df.cor.dna$truth, df.cor.dna$response) 
x$rho.c                                                             # 0.9164678 0.9016164 0.9291608
g.dna <- df.cor.dna %>% 
  ggplot(aes(x=truth, y=response)) + geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw()
ggsave(g.dna, filename='./figures/dna_model_predictions.pdf',
       width = 4, height = 4, useDingbats=FALSE)

# test which model is better?
df.cor.all <- df.cor.full %>% 
  rename(response_full=response) %>% 
  full_join(df.cor.dna %>% rename(response_dna=response),
            by=c('id', 'truth')) %>% 
  mutate(diff_full=abs(truth-response_full)) %>% 
  mutate(diff_dna=abs(truth-response_dna))
t.test(df.cor.all$diff_dna, df.cor.all$diff_full, paired = TRUE)

# Paired t-test
# 
# data:  df.cor.all$diff_dna and df.cor.all$diff_full
# t = 3.654, df = 517, p-value = 0.0002846
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   0.01107163 0.03682157
# sample estimates:
#   mean difference 
# 0.0239466 



# get model importance
df.importance <- map(seq_len(100), .f=function(x){
  enframe(rr_full$learners[[1]]$learner$importance()) %>% mutate(iter=x)}) %>% 
  bind_rows()

g <- df.importance %>% 
  group_by(iter) %>% mutate(full.size=sum(value)) %>% 
  mutate(rel_weight=value/full.size) %>% 
  group_by(name) %>% 
  summarize(r=mean(rel_weight), s=sd(rel_weight), .groups='drop') %>% 
  arrange(r) %>% 
  mutate(name=factor(name, levels=name)) %>% 
  ggplot(aes(x=name, y=r)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=r-s, ymax=r+s), width=0.1) + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.y = element_blank()) + 
  coord_flip() + 
  ylab('Mean relative model weight') + xlab('')
ggsave(g, filename='./figures/model_weight.pdf', 
       width=5, height = 4, useDingbats=FALSE)

# leave-one-plate-out
df.lopo <- list()
for (p in unique(df.data$Plate)){
  test <- df.data %>% 
    filter(Plate==p)
  train <- df.data %>% 
    filter(Plate!=p)
  # train
  learner <- lrn("regr.ranger", mtry= to_tune(1,7),
                 num.trees=to_tune(50, 500), 
                 importance='impurity')
  task <- as_task_regr(train %>% 
                         select(-c(Sample_ID, Plate)), 
                       target = "log_copies_per_extraction")
  tnr_grid_search = tnr("grid_search", resolution = 5, batch_size = 5)
  at = auto_tuner(
    tuner = tnr_grid_search,
    learner = learner,
    resampling = rsmp("cv", folds = 3),
    measure = msr('regr.rsq'),
  )
  rsmp_cv3 = rsmp("repeated_cv", repeats = 10, folds = 10)
  rr_plate = resample(task, at, rsmp_cv3, store_models = TRUE)
  
  # predict
  task_pred <- as_task_regr(test, target='log_copies_per_extraction')
  pred.matrix <- matrix(NA, nrow=nrow(test), ncol=100)
  rownames(pred.matrix) <- test$Sample_ID
  colnames(pred.matrix) <- paste0('model_', seq_len(100))  
  pred.matrix.full <- pred.matrix
  pb <- progress::progress_bar$new(total=100)
  for (i in seq_len(100)){
    pred <- rr_plate$learners[[i]]$predict(task=task_pred)
    pred.matrix.full[,i] <- pred$response
    pb$tick()
  }
  df.lopo[[p]] <- enframe(rowMeans(pred.matrix.full), name='Sample_ID', 
                      value='pred_full')
}

df.lopo %>% 
  bind_rows() %>% 
  full_join(df.data, by='Sample_ID') %>% 
  group_by(Plate)  %>% 
  mutate(se.full=(log_copies_per_extraction-pred_full)^2) %>% 
  summarise(mse.full=mean(se.full),
            med.full=median(se.full),
            pr.full=cor(log_copies_per_extraction, pred_full),
            r2.full=mlr3measures::rsq(truth=log_copies_per_extraction, response = pred_full),
            rho.full=cor(log_copies_per_extraction, pred_full, method='spearman'),
            ccc.full=epiR::epi.ccc(log_copies_per_extraction, pred_full)$rho.c$est[1])

## A tibble: 6 × 7
# Plate mse.full med.full pr.full r2.full rho.full ccc.full$est
# <chr>    <dbl>    <dbl>   <dbl>   <dbl>    <dbl>        <dbl>
# 1 P01     0.0703   0.0313   0.967   0.911    0.956        0.952
# 2 P03     0.0485   0.0161   0.959   0.913    0.915        0.954
# 3 P05     0.0557   0.0113   0.931   0.865    0.926        0.927
# 4 P07     0.151    0.0541   0.916   0.763    0.873        0.885
# 5 P10     0.0976   0.0286   0.938   0.860    0.914        0.932
# 6 P13     0.0935   0.0233   0.933   0.842    0.895        0.900

g <- df.lopo %>% 
  bind_rows() %>% 
  full_join(df.data, by='Sample_ID') %>% 
  ggplot(aes(x=log_copies_per_extraction, y=pred_full, col=Plate)) + 
  geom_abline(slope = 1, intercept = 0, lty=2) + 
  geom_point() + 
  facet_wrap(~Plate) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) + 
  xlab('Real log10(16S copies per extraction)') + 
  ylab('Predicted log10(16S copies per extraction)') + 
  ggthemes::scale_colour_tableau(guide='none') + 
  theme(strip.background = element_blank())
ggsave(g, filename='./figures/lopo.pdf',
       width = 6, height = 4, useDingbats=FALSE)
