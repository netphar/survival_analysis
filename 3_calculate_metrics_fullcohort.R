rm(list=ls(all=TRUE))
options(scipen=999)
setwd('/Users/zagidull/Documents/fimm_files/survival_all/survival_loukavaara/data')

library(pacman)
pacman::p_load(
    survival,
    Hmisc,
    rms,
    pec,
    riskRegression,
    timeROC,
    tidyverse,
    SurvMetrics
)
set.seed(2111410)

D <- read_csv(file = '~/Desktop/D_june14.csv', 
              col_types = cols(
                  .default = col_factor(),
                  id = col_double(),
                  death = col_double(),
                  time = col_double(),
                  age = col_double(),
                  bmi = col_double(),
                  hemoglobin = col_double(),
                  hematocrit = col_double(),
                  leukocytes = col_double(),
                  thrombocytes = col_double(),
                  creatinine = col_double(),
                  Ca125 = col_double())) %>%
    as.data.frame

# Function definitions
get_metrics_flexi <- function(df, formula, mat_tree=NULL){
    # setup 
    surv_obj <- Surv(df$time, df$death)
    times <- sort(unique(df$time))
    times_12 <- times[times < 13]
    times_36 <- times[times < 37]
    times_60 <- times[times < 61]
    med_ind <- floor(median(1:length(times)))
    
    if(is.null(mat_tree)){#for cox models
        dat <- df[,all.vars(formula)]
        dat <- droplevels(dat) # in case some factor vars after resample end up missing a level
        cox <- coxph(formula, data = dat, x = T, y = T)
        
        # get cox survival probabilities
        mat_all <- predictSurvProb(cox, dat, times)
        mat_12 <- predictSurvProb(cox, dat, times_12)
        mat_36 <- predictSurvProb(cox, dat, times_36)
        mat_60 <- predictSurvProb(cox, dat, times_60)
    }else{# when feeding matrix of proba from tree models. matrix of proba is n_patients x number of unique timepoints for a given group
        #times_ind <- match(times, )
        mat_all <- mat_tree
        mat_12 <- mat_all[,which(times < 13)] # 1 year
        mat_36 <- mat_all[,which(times < 37)] # 3 year
        mat_60 <- mat_all[,which(times < 61)] # 5 years
    }
    # get IBS 
    ibs_all <- IBS(surv_obj, sp_matrix = mat_all, times)
    ibs_12 <- IBS(surv_obj, sp_matrix = mat_12, times_12)
    ibs_36 <- IBS(surv_obj, sp_matrix = mat_36, times_36)
    ibs_60 <- IBS(surv_obj, sp_matrix = mat_60, times_60)
    names(ibs_all) <- 'overall IBS'
    names(ibs_12) <- 'IBS at 1y'
    names(ibs_36) <- 'IBS at 3y'
    names(ibs_60) <- 'IBS at 5y'
    
    # get c-index
    cind <- Cindex(surv_obj, mat_all[,med_ind])
    names(cind) <- 'Harrell C'
    return(c(cind,ibs_12,ibs_36,ibs_60,ibs_all))
}

bootstrap_metrics <- function(df_vars, formula, mat_tree=NULL, bs_iter=120){
    holder <- list()
    for(i in c(1:bs_iter)){
        index_boot <- sample(1:nrow(df_vars), size = nrow(df_vars), replace = TRUE) # get resampled index
        df_boot <- df_vars[index_boot,]
        
        if(is.null(mat_tree)){
            res_bs <- tryCatch(expr = {get_metrics_flexi(df_boot, formula=formula, mat_tree=mat_tree)}, error = function(e) {return(NA)} )
        } else {
            index_col_boot <- match( sort(unique(df_boot$time)), sort(unique(df_vars$time)) ) # when bootstrapping we end up not having some time points, this accounts for it by subsetting the matrix of probabilities
            mat_boot <- mat_tree[index_boot, index_col_boot]
            #browser()
            res_bs <- get_metrics_flexi(df_boot, formula=formula, mat_tree=mat_boot)
        }
        holder[[i]] <- res_bs
    }
    bs_out <- do.call(rbind, holder)
    return(bs_out)
}

get_ci <- function(df_bs){
    cis <- base::apply(df_bs, MARGIN = 2, function (x) quantile(x, na.rm = TRUE, probs=c(0.025, 0.975))) # empirical bootstrap
    means <- base::apply(df_bs, MARGIN = 2, function (x) mean(x, na.rm = TRUE))
    out <- rbind(upper=cis[2,], means, lower=cis[1,])
    return(round(out, digits = 4))
}

#####
# trad for all 
# results calculated on cluster
load(file = 'boot_cox_trad_all.Rdata')
load(file = 'boot_tree_trad_all.Rdata')
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
trad <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))

# ext for all 
load(file = 'boot_cox_ext_all.Rdata')
load(file = 'boot_tree_ext_all.Rdata')
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
ext <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))

load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/boots_results_all_30june.Rdata")
holder_all[['trad']] <- trad
holder_all[['ext']] <- ext
save(holder_all, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/boots_results_all_1july.Rdata")

