rm(list=ls(all=TRUE))

set.seed(2111410)
options(scipen=999)
library(pacman)
pacman::p_load(survival, Hmisc, rms, pec, riskRegression, timeROC, tidyverse, SurvMetrics)
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

#####
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
# Prepare objects to calculate and store bootstrapped metrics (IBS at 1y-2y-5y-all and Harrell's cindex)
holder = list()
num_iter <- 120 # num bootstrap iterations
fx_trad <-Surv(time, death) ~ age+stage+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm # trad variable set
fx_ext <-Surv(time, death) ~ age+stage+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm+diameter_more_5cm+cytology+ER+CD171 # extended variable set

#####
# TRAD by ProMisE
# mmrd - ext
df <- D %>%
    filter(ProMisE=='MMRd') %>%
    as.data.frame
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_mmrd_29june.Rdata")
df_cox <- bootstrap_metrics(df, formula=fx_trad, bs_iter = num_iter)
df_tree <- bootstrap_metrics(df, formula=fx_trad, mat_tree = mat_trad_mmrd, bs_iter = num_iter)
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
trad_mmrd <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))
holder[['trad_mmrd']] <-  trad_mmrd

# nsmp - tradd
df <- D %>%
    filter(ProMisE=='NSMP') %>%
    as.data.frame
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_nsmp_29june.Rdata")
df_cox <- bootstrap_metrics(df, formula=fx_trad, bs_iter = num_iter)
df_tree <- bootstrap_metrics(df, formula=fx_trad, mat_tree = mat_trad_nsmp, bs_iter = num_iter)
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
trad_nsmp <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))
holder[['trad_nsmp']] <-  trad_nsmp

# p53ab - trad
df <- D %>%
    filter(ProMisE=='p53ab') %>%
    as.data.frame
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_p53ab_29june.Rdata")
df_cox <- bootstrap_metrics(df, formula=fx_trad, bs_iter = num_iter)
df_tree <- bootstrap_metrics(df, formula=fx_trad, mat_tree = mat_trad_p53ab, bs_iter = num_iter)
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
trad_p53 <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))
holder[['trad_p53']] <-  trad_p53

holder_trad <- holder
save(holder_trad, file = '~/Documents/fimm_files/survival_all/survival_loukavaara/data/boots_results_trad_29june.Rdata')

#####
# EXT 
# mmrd ext
df <- D %>%
    filter(ProMisE=='MMRd') %>%
    as.data.frame
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_mmrd_29june.Rdata")
df_cox <- bootstrap_metrics(df, formula=fx_ext, bs_iter = num_iter)
df_tree <- bootstrap_metrics(df, formula=fx_ext, mat_tree = mat_ext_mmrd, bs_iter = num_iter)
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
ext_mmrd <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))
holder[['ext_mmrd']] <-  ext_mmrd

# nsmp - ext
df <- D %>%
    filter(ProMisE=='NSMP') %>%
    as.data.frame
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_nsmp_29june.Rdata")
df_cox <- bootstrap_metrics(df, formula=fx_ext, bs_iter = num_iter)
df_tree <- bootstrap_metrics(df, formula=fx_ext, mat_tree = mat_ext_nsmp, bs_iter = num_iter)
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
ext_nsmp <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))
holder[['ext_nsmp']] <-  ext_nsmp

# p53ab - ext
df <- D %>%
    filter(ProMisE=='p53ab') %>%
    as.data.frame
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_p53ab_29june.Rdata")
df_cox <- bootstrap_metrics(df, formula=fx_ext, bs_iter = num_iter)
df_tree <- bootstrap_metrics(df, formula=fx_ext, mat_tree = mat_ext_p53ab, bs_iter = num_iter)
ci_cox <- get_ci(df_cox)
ci_tree <- get_ci(df_tree)
ext_p53 <-list('cox'=ci_cox,'tree'=ci_tree, 'bootstrap'=list('bs_cox'=df_cox,'bs_tree'=df_tree))
holder[['ext_p53']] <-  ext_p53

holder_all <- holder
save(holder_all, file = '~/Documents/fimm_files/survival_all/survival_loukavaara/data/boots_results_all_30june.Rdata')
