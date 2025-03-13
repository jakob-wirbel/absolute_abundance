# ##############################################################################
#
## Apply the model to your own data
#
# ##############################################################################


.f_apply_model <- function(df, model){
  require("tidyverse")
  require("mlr3")  
  require("progress")
  
  # checks and balances
  # Todo: could do more checks that will be done later within mlr3 to give some
  # better error messages... Maybe let's see how many people will actually use 
  # this script
  
  # needed columns
  needed.colnames <- c('DNA_concentration', 
                       'unclassified', 
                       'k__Bacteria', 'k__Archaea', 'k__Eukaryota',
                       'bac_frac', 'alpha', 'Sample_type')
  ov.colnames <- setdiff(needed.colnames, colnames(df))
  if (length(ov.colnames) > 0){
    stop("Some required columns are missing!\n\tMissing are: ", 
         paste(ov.colnames, collapse = ', '))
    
  }
  # NA values
  no.nas <- map(colnames(df), .f=function(x){sum(is.na(df[[x]]))}) %>% 
    unlist()
  if (any(no.nas > 0)){
    stop("Some columns contain missing data! Please remove NAs or ",
         "add in some dummy variables\n\tColumns with NAs are: ", 
         paste(colnames(df)[which(no.nas > 0)], collapse = ', '))
  }
  
  # create prediction task
  pred.external <- as_task_regr(
    df %>% mutate('dummy_target'=0) %>% 
      rename(DNA_conc=DNA_concentration),
    target = "dummy_target")
  
  # create data frame for predictions
  pred.mat <- matrix(data=NA, nrow=nrow(df), ncol=100,
                     dimnames = list(NULL,
                                     paste0('model_', seq(100))))
  pb <- progress::progress_bar$new(total=100)
  # predict
  for (i in seq(100)){
    p <- model$learners[[i]]$learner$predict(pred.external)
    pred.mat[,i] <- p$response
    pb$tick()
  }
  
  # collate
  df.out <- df %>% 
    mutate(pred_copies=rowMeans(pred.mat))
  
  # return
  return(df.out)
}
