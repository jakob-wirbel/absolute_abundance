# ##############################################################################
#
## Check if we see the same with the PD data
#
# ##############################################################################

library("tidyverse")
library("mlr3")

ref.data <- read_tsv('./data/alloHCT_data.tsv')
feat.external <- read_tsv('./data/PD_data.tsv')

# ##############################################################################
# check that the DNA concentration relationship holds
g <- feat.external %>% 
  bind_rows(ref.data %>% mutate(type='reference')) %>% 
  arrange(desc(type)) %>% 
  ggplot(aes(x=DNA_concentration, y=log_copies_per_extraction, col=type)) + 
    geom_point() + 
    geom_smooth(method='lm', formula='y~log(x)') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    scale_colour_manual(values=c('#E98300', '#7F7776'))
ggsave(g, filename='./figures/PD_external_correlation.pdf',
       width = 5, height = 4, useDingbats=FALSE)

cor.test(feat.external$log_copies_per_extraction, 
         feat.external$DNA_concentration, method='spearman')

#' Spearman's rank correlation rho
#' 
#' data:  feat.external$log_copies_per_extraction and feat.external$DNA_concentration
#' S = 24556, p-value < 2.2e-16
#' alternative hypothesis: true rho is not equal to 0
#' sample estimates:
#'       rho 
#' 0.9778263 

# ##############################################################################
# external predictions

# load the model
load('./models/full_model.Rdata')

# prepare the prediction task
pred.external <- as_task_regr(
  feat.external %>% rename(DNA_conc=DNA_concentration), 
  target = "log_copies_per_extraction")

# make predictions
pred.mat <- matrix(data=NA, nrow=nrow(feat.external), ncol=100,
                   dimnames = list(feat.external$Sample_ID, 
                                   paste0('model_', seq(100))))
for (i in seq(100)){
  p <- rr_full$learners[[i]]$learner$predict(pred.external)
  pred.mat[,i] <- p$response
}

# compare
pred.all <- enframe(rowMeans(pred.mat), 
                    name='Sample_ID', 
                    value='pred_copies')
df.eval <- feat.external %>% 
  full_join(pred.all, by='Sample_ID') %>% 
  filter(!is.na(log_copies_per_extraction))

# pearson
cor(df.eval$log_copies_per_extraction, 
    df.eval$pred_copies) # 0.9608967
# spearman
cor(df.eval$log_copies_per_extraction, 
    df.eval$pred_copies, method='spearman') # 0.9562386
# CCC
x <- epiR::epi.ccc(df.eval$log_copies_per_extraction, 
                   df.eval$pred_copies) # rho: 0.9606072
# mse
msr('regr.mse')$fun(response=df.eval$pred_copies, 
                    truth=df.eval$log_copies_per_extraction) # 0.02175194
# mae
msr('regr.mae')$fun(response=df.eval$pred_copies, 
                    truth=df.eval$log_copies_per_extraction) # 0.1025855
# rsq
msr('regr.rsq')$fun(response=df.eval$pred_copies, 
                    truth=df.eval$log_copies_per_extraction) # 0.9209634
# sae
msr('regr.sae')$fun(response=df.eval$pred_copies, 
                    truth=df.eval$log_copies_per_extraction) # 19.28608
# plot
g <- df.eval %>% 
  ggplot(aes(x=log_copies_per_extraction, y=pred_copies)) + 
    geom_abline(slope = 1, intercept = 0, lty=2, colour='darkgrey') + 
    geom_point() + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab('Measured log10(16S copies per extraction)') +
    ylab('Predicted log10(16S copies per extraction)') 
ggsave(g, filename='./figures/PD_external_pred.pdf',
       width = 4, height = 4, useDingbats=FALSE)
