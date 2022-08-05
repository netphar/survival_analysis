rm(list=ls(all=TRUE))
setwd('/Users/zagidull/Documents/fimm_files/survival_all/survival_loukavaara/data')

library(pacman)
pacman::p_load(tidyverse, iai, survival, pec, SurvMetrics, riskRegression)

# setup global vars for iai
num_threads = 5
Sys.setenv(JULIA_NUM_THREADS=num_threads)

# global vars for modelling
bs_iter = 1000
random_seed = 1
fx_trad <-Surv(time, death) ~ age+stage+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm
fx_ext <-Surv(time, death) ~ age+stage+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm+diameter_more_5cm+cytology+ER+CD171 # extended variable set

# func defs
make_factor_simple <- function(x) {
    attr(x, "spec") <- NULL
    for (i in names(x)) {
        if ((i != 'death') &&
            (i != 'id') &&
            (i != 'time') &&
            (i != 'age') &&
            (i != 'bmi') &&
            (i != 'postop_helsinki_score') &&
            (i != 'risk_factors') &&
            (i != 'hemoglobin') &&
            (i != 'hematocrit') &&
            (i != 'leukocytes') &&
            (i != 'thrombocytes') &&
            (i != 'creatinine') &&
            (i != 'hcg') &&
            (i != 'hcgbv') &&
            (i != 'Ca125')) {
            x[, i] <- factor(x[, i], ordered = FALSE)
        }
    }
    x
}
predict_proba_iai_flexi <- function(df_input, learner, vars){ # this is for cohorts *statified* by ProMisE
    times <- sort(unique(df_input$time)) 
    train_X <- df_input[,vars] 
    
    mat_tree <- 1-as.matrix(do.call(cbind, lapply(times, function(x) iai::predict(learner, df_input, t=x))))
    return(mat_tree)
}
get_ci <- function(df_bs){
    cis <- base::apply(df_bs, MARGIN = 2, function (x) quantile(x, na.rm = TRUE, probs=c(0.025, 0.975))) # empirical bootstrap
    means <- base::apply(df_bs, MARGIN = 2, function (x) mean(x, na.rm = TRUE))
    out <- rbind(upper=cis[2,], means, lower=cis[1,])
    return(round(out, digits = 4))
} 
get_metrics_cind <- function(df_input, mat_tree, fx, bs_iter){
    times <- sort(unique(df_input$time))
    vars <- all.vars(fx)
    surv_obj <- Surv(df_input$time, df_input$death)
    med_ind <- floor(median(1:length(times))) # harrell's Cind at median ~= Uno's time-dependent Cind
    cox <- coxph(fx, data = df_input[, vars], x = T, y = T)
    mat_cox <- predictSurvProb(cox, df_input[, vars], times)
    
    holder_c <- vector("list", length = bs_iter) # loops are my friends
    #browser()
    for(i in c(1:bs_iter)){
        ## cindex pec COX and TREE with bs
        index_boot <- sample(1:nrow(df_input), size = nrow(df_input), replace = TRUE)  
        cind <- pec::cindex(
            list("tree"=matrix(mat_tree[index_boot, med_ind]), 
                 "cox"=matrix(mat_cox[index_boot, med_ind])),
            formula=fx,
            data=df_input[index_boot, vars],
            eval.times=times[med_ind])$AppCindex
        holder_c[[i]] <- unlist(cind)
    }
    cis <- get_ci(do.call(rbind, holder_c))
    return(cis)
}
get_metrics_ibs <- function(df_input, mat_tree, fx, bs_iter, cutoff_time=NULL){
    if(is.null(cutoff_time)){
        times <- sort(unique(df_input$time))
    } else {
        times <- sort(unique(df_input$time))
        cutoff_horizon_idx <- which(times <= cutoff_time) # here we use cutoff_time to get IBS at some timepoint
        times <- times[cutoff_horizon_idx]
        mat_tree <- mat_tree[,cutoff_horizon_idx]
    }
    vars <- all.vars(fx)
    surv_obj <- Surv(df_input$time, df_input$death)
    med_ind <- floor(median(1:length(times))) # harrell's Cind at median ~= Uno's time-dependent Cind
    cox <- coxph(fx, data = df_input[, vars], x = T, y = T)
    mat_cox <- predictSurvProb(cox, df_input[, vars], times)
    
    holder_i <- vector("list", length = bs_iter)
    for(i in c(1:bs_iter)){
        ## cindex pec COX and TREE with bs
        index_boot <- sample(1:nrow(df_input), size = nrow(df_input), replace = TRUE)  
        
        ## ibs pec COX and TREE with bs
        times_bs <- unique(sort(df_input[index_boot, 'time']))
        index_time <- match(times_bs, times) # need to get timepoint idx again, since we bsing with replacement
        index_time <- index_time[ !(is.na(index_time))]
        mat_tree_bs <- cbind(rep(1, nrow(df_input)), as.matrix(mat_tree[index_boot, index_time]))
        mat_cox_bs <- cbind(rep(1, nrow(df_input)), as.matrix(mat_cox[index_boot, index_time]))
        #browser()
        ibs <- pec::crps(
            pec::pec(
                list("tree"=mat_tree_bs, "cox"=mat_cox_bs),
                formula=Surv(time, death)~1,
                data=df_input[index_boot, vars],
                times= times[index_time],
                exact=FALSE,
                cens.model="marginal",
                splitMethod="none",
                B=0)
        )
        holder_i[[i]] <- ibs
        
    }
    ibs <- get_ci(t(do.call(cbind, holder_i)))
    return(ibs)
}

# load data
## with missinig
dataset <- read_delim(
    "endometrial_carcinoma_clean.csv",
    na = "",
    delim=",",
    quote = '"',
    col_types = cols(
        .default = col_character(),
        id = col_double(),
        death = col_double(),
        time = col_double(),
        age = col_double(),
        bmi = col_double(),
        postop_helsinki_score = col_double(),
        risk_factors = col_double(),
        hemoglobin = col_double(),
        hematocrit = col_double(),
        leukocytes = col_double(),
        thrombocytes = col_double(),
        creatinine = col_double(),
        hcg = col_double(),
        hcgbv = col_double(),
        Ca125 = col_double())) %>%
    select(-starts_with("transform")) %>%
    #filter(death!=2) %>%
    as.data.frame

dataset <- dataset %>%
    make_factor_simple() %>%
    as.data.frame

## imputed == no missing
D <- read_csv(file = 'D_june14.csv', 
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

# predict_proba_iai <- function(df, learner, trad=T){ # this is for a full cohort model, ie *not* stratified by ProMisE
#     times <- sort(unique(df$time)) 
#     if(trad==T){
#         train_X <- df[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')] 
#     }else{
#         train_X <- df[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
#                          'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')] 
#     }
#     
#     mat_tree <- 1-as.matrix(do.call(cbind, lapply(times, function(x) iai::predict(learner, df, t=x))))
#     return(mat_tree)
# }

########################################
# ALL
########################################
df <- D 
df$ER <- relevel(df$ER, ref='Negative')
df$stage <- factor(df$stage, ordered = T, levels=c('I','II','III','IV'))

########
# TRAD #
########
## lets use real data for training and imputed data for valid
ids <- dataset[ !(is.na(dataset$ProMisE)),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm','ProMisE')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

## train OST
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=1000,
        death_minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 6:10
)
iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_trad_all <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_trad_all_new <- predict_proba_iai_flexi(df, best_learner_trad_all, vars=vars) # now predicted

## get metrics with 1000 iterations of bootstrap. Harrell's C-index, IBS at 1y/2y/5y/all
fx_trad_all <- update(fx_trad, ~.+ProMisE)
cind <- get_metrics_cind(df, mat_trad_all_new, fx_trad_all, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_trad_all_new, fx_trad_all, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_trad_all_new, fx_trad_all, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_trad_all_new, fx_trad_all, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_trad_all_new, fx_trad_all, bs_iter = bs_iter)
all_trad <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)

############
# EXTENDED #
############
## lets use real data for training and imputed data for valid
ids <- dataset[ !(is.na(dataset$ProMisE)),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm',
          'diameter_more_5cm','cytology','ER','CD171')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

## train OST
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=1000
    ),
    criterion = 'loglikelihood',
    max_depth = 6:10,
    death_minbucket = 3
)
iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_ext_all <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_ext_all_new <- predict_proba_iai_flexi(df, best_learner_ext_all, vars=vars) # now predicted

## get metrics with 1000 iterations of bootstrap. Harrell's C-index, IBS at 1y/2y/5y/all
fx_ext_all <- update(fx_ext, ~.+ProMisE)
cind <- get_metrics_cind(df, mat_ext_all_new, fx_ext_all, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_ext_all_new, fx_ext_all, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_ext_all_new, fx_ext_all, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_ext_all_new, fx_ext_all, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_ext_all_new, fx_ext_all, bs_iter = bs_iter)
all_ext <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)

########################################
# p53ab
########################################
promise_label <- 'p53ab'
df <- D %>%
    filter(ProMisE==promise_label) %>%
    as.data.frame
df$ER <- relevel(df$ER, ref='Negative')
df$stage <- factor(df$stage, ordered = T, levels=c('I','II','III','IV'))

########
# TRAD #
########
## lets use real data for training and imputed data for valid
ids <- dataset[which(dataset$ProMisE==promise_label),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

## train OST
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=2500,
        death_minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 4:9
)
iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_trad_p53ab <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_trad_p53ab_new <- predict_proba_iai_flexi(df, best_learner_trad_p53ab, vars=vars) # now predicted

## get metrics with 1000 iterations of bootstrap. Harrell's C-index, IBS at 1y/2y/5y/all
cind <- get_metrics_cind(df, mat_trad_p53ab_new, fx_trad, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_trad_p53ab_new, fx_trad, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_trad_p53ab_new, fx_trad, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_trad_p53ab_new, fx_trad, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_trad_p53ab_new, fx_trad, bs_iter = bs_iter)
p53ab_trad <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)
# load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/individual_survival_curves_OSTmodel/old/mat_ext_p53ab_29june.Rdata") # used then

# lets look at the p53ab with CPH using FSII
summary(coxph(fx_ext, data = df[, all.vars(fx_ext)], x = T, y = T))

############
# EXTENDED #
############
## lets use real data for training and imputed data for valid
ids <- dataset[which(dataset$ProMisE==promise_label),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm',
          'diameter_more_5cm','cytology','ER','CD171')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

## train OST
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=2500,
        death_minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 4:9,

)
iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_ext_p53ab <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_ext_p53ab_new <- predict_proba_iai_flexi(df, best_learner_ext_p53ab, vars=vars) # now predicted

## get metrics with 1000 iterations of bootstrap. Harrell's C-index, IBS at 1y/2y/5y/all
cind <- get_metrics_cind(df, mat_ext_p53ab_new, fx_ext, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_ext_p53ab_new, fx_ext, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_ext_p53ab_new, fx_ext, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_ext_p53ab_new, fx_ext, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_ext_p53ab_new, fx_ext, bs_iter = bs_iter)
p53ab_ext <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)
# load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/individual_survival_curves_OSTmodel/old/mat_ext_p53ab_29june.Rdata") # used then

########################################
# NSMP
########################################
promise_label <- 'NSMP'
df <- D %>%
    filter(ProMisE==promise_label) %>%
    as.data.frame
df$ER <- relevel(df$ER, ref='Negative')
df$stage <- factor(df$stage, ordered = T, levels=c('I','II','III','IV'))

########
# TRAD #
########
## lets use real data for training and imputed data for valid
ids <- dataset[which(dataset$ProMisE==promise_label),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=2500,
        death_minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 4:12
)
iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_trad_nsmp <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_trad_nsmp_new <- predict_proba_iai_flexi(df, best_learner_trad_nsmp, vars=vars) # now predicted

## get metrics with 1000 iterations of bootstrap. Harrell's C-index, IBS at 1y/2y/5y/all
cind <- get_metrics_cind(df, mat_trad_nsmp_new, fx_trad, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_trad_nsmp_new, fx_trad, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_trad_nsmp_new, fx_trad, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_trad_nsmp_new, fx_trad, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_trad_nsmp_new, fx_trad, bs_iter = bs_iter)
nsmp_trad <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)
# load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/individual_survival_curves_OSTmodel/old/mat_ext_p53ab_29june.Rdata") # used then

############
# EXTENDED #
############
ids <- dataset[which(dataset$ProMisE==promise_label),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm',
          'diameter_more_5cm','cytology','ER','CD171')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=2500,
        death_minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 4:12,
    
)
iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_ext_nsmp <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_ext_nsmp_new <- predict_proba_iai_flexi(df, best_learner_ext_nsmp, vars=vars) # now predicted

## get metrics with 1000 iterations of bootstrap. Harrell's C-index, IBS at 1y/2y/5y/all
cind <- get_metrics_cind(df, mat_ext_nsmp_new, fx_ext, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_ext_nsmp_new, fx_ext, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_ext_nsmp_new, fx_ext, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_ext_nsmp_new, fx_ext, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_ext_nsmp_new, fx_ext, bs_iter = bs_iter)
nsmp_ext <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)
# load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/individual_survival_curves_OSTmodel/old/mat_ext_p53ab_29june.Rdata") # used then

########################################
# MMRd
########################################
promise_label <- 'MMRd'
df <- D %>%
    filter(ProMisE==promise_label) %>%
    as.data.frame
df$ER <- relevel(df$ER, ref='Negative')
df$stage <- factor(df$stage, ordered = T, levels=c('I','II','III','IV'))

########
# TRAD #
########
## lets use real data for training and imputed data for valid
ids <- dataset[which(dataset$ProMisE==promise_label),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=500,
        minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 4:10
)

iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_trad_mmrd <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_trad_mmrd_new <- predict_proba_iai_flexi(df, best_learner_trad_mmrd, vars=vars) # now predicted

cind <- get_metrics_cind(df, mat_trad_mmrd_new, fx_trad, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_trad_mmrd_new, fx_trad, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_trad_mmrd_new, fx_trad, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_trad_mmrd_new, fx_trad, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_trad_mmrd_new, fx_trad, bs_iter = bs_iter)
mmrd_trad <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)
# load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/individual_survival_curves_OSTmodel/old/mat_ext_p53ab_29june.Rdata") # used then

############
# EXTENDED #
############
ids <- dataset[which(dataset$ProMisE==promise_label),'id']
train_X <- df[df$id %in% ids,]
test_X <- df[ !(df$id %in% ids),]

## prep temp data
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm',
          'diameter_more_5cm','cytology','ER','CD171')

train_died <- as.logical(train_X[,'death'])
train_times <- train_X[,'time']
train_X <- train_X[,vars] 

test_died <- as.logical(test_X[,'death'])
test_times <- test_X[,'time']
test_X <- test_X[,vars] 

all_died  <- as.logical(df[,'death'])
all_times <- df[,'time']
all_X <- df[,vars]

grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed=random_seed,
        show_progress=F,
        smooth_survival_curves=TRUE,
        ls_num_tree_restarts=500,
        minbucket=2
    ),
    criterion = 'loglikelihood',
    max_depth = 4:10
)

iai::fit(grid, train_X, train_died, train_times, test_X, test_died, test_times, validation_criterion='harrell_c_statistic')
best_learner_ext_mmrd <- iai::get_learner(grid)
rm(grid)

## predict indi survival curves
mat_ext_mmrd_new <- predict_proba_iai_flexi(df, best_learner_ext_mmrd, vars=vars) # now predicted

cind <- get_metrics_cind(df, mat_ext_mmrd_new, fx_ext, bs_iter = bs_iter)
ibs_1y <- get_metrics_ibs(df, mat_ext_mmrd_new, fx_ext, bs_iter = bs_iter, cutoff_time = 12)
ibs_2y <- get_metrics_ibs(df, mat_ext_mmrd_new, fx_ext, bs_iter = bs_iter, cutoff_time = 24)
ibs_5y <- get_metrics_ibs(df, mat_ext_mmrd_new, fx_ext, bs_iter = bs_iter, cutoff_time = 60)
ibs_all <- get_metrics_ibs(df, mat_ext_mmrd_new, fx_ext, bs_iter = bs_iter)
mmrd_ext <- list('cind'=cind, 'ibs_1y'=ibs_1y, 'ibs_2y'=ibs_2y, 'ibs_5y'=ibs_5y, 'ibs_all'=ibs_all)
# load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/individual_survival_curves_OSTmodel/old/mat_ext_p53ab_29june.Rdata") # used then

# save bs results 
holder <- list('all_trad'=all_trad,
               'all_ext'=all_ext,
               'mmrd_trad'=mmrd_trad,
               'mmrd_ext'=mmrd_ext,
               'nsmp_trad'=nsmp_trad,
               'nsmp_ext'=nsmp_ext,
               'p53ab_trad'=p53ab_trad,
               'p53ab_ext'=p53ab_ext)
save(holder, file = "output/boots_results_all_4aug.Rdata")

# save indi surv curves
MASS::write.matrix(mat_trad_p53ab_new, 
                   file = "individual_survival_curves_OSTmodel/p53ab_trad.csv",
                   sep = ',')
MASS::write.matrix(mat_ext_p53ab_new, 
                   file = "individual_survival_curves_OSTmodel/p53ab_ext.csv",
                   sep = ',')
MASS::write.matrix(mat_trad_nsmp_new, 
                   file = "individual_survival_curves_OSTmodel/nsmp_trad.csv",
                   sep = ',')
MASS::write.matrix(mat_ext_nsmp_new, 
                   file = "individual_survival_curves_OSTmodel/nsmp_ext.csv",
                   sep = ',')
MASS::write.matrix(mat_trad_mmrd_new, 
                   file = "individual_survival_curves_OSTmodel/mmrd_trad.csv",
                   sep = ',')
MASS::write.matrix(mat_ext_mmrd_new, 
                   file = "individual_survival_curves_OSTmodel/mmrd_ext.csv",
                   sep = ',')
MASS::write.matrix(mat_trad_all_new, 
                   file = "individual_survival_curves_OSTmodel/all_trad.csv",
                   sep = ',')
MASS::write.matrix(mat_ext_all_new, 
                   file = "individual_survival_curves_OSTmodel/all_ext.csv",
                   sep = ',')

# get tree
## make fig 3
iai::write_dot('output/fig3_p53_ext.dot', best_learner_ext_p53ab)
# afterwards change dot file manually by adding fontsize=20 into node/edge properties. Also add two spaces after every "expected survival: xx.xx"
# then render as svg with graphviz dot -Tfile.dot -o file.svg
# then use inkscape to modify text if necessary. 
# then export as png with dpi=300 setting

######
# CC
# D_cc <- read_csv(file = '~/Desktop/D_cc_june20.csv', 
#                  col_types = cols(
#                      .default = col_factor(),
#                      id = col_double(),
#                      death = col_double(),
#                      time = col_double(),
#                      age = col_double(),
#                      bmi = col_double(),
#                      hemoglobin = col_double(),
#                      hematocrit = col_double(),
#                      leukocytes = col_double(),
#                      thrombocytes = col_double(),
#                      creatinine = col_double(),
#                      Ca125 = col_double())) %>%
#     as.data.frame
# 
# # traditional optimal survival forest CC
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed+1,
#         missingdatamode = "separate_class",
#         num_threads = 1,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:11,
# )
# train_died <- as.logical(D_cc[,'death'])
# train_times <- D_cc[,'time']
# train_X <- D_cc[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_trad_cc <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_trad_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # 0.8635481 with depth 8
# 
# mat_trad_cc <- predict_proba_iai(D_cc, best_learner_trad_cc)
# save(mat_trad_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_cc.Rdata")
# 
# # extended optimal survival forest CC
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed-1,
#         missingdatamode = "separate_class",
#         num_threads = 1,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:11,
# )
# train_died <- as.logical(D_cc[,'death'])
# train_times <- D_cc[,'time']
# train_X <- D_cc[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
#                    'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')] 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_ext_cc <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_ext_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # 0.8365357 with depth 10
# 
# mat_ext_cc <- predict_proba_iai(D_cc, best_learner_ext_cc)
# save(mat_ext_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_cc.Rdata")
# 
# ## CC
# ### TRAIN WITHOUT PROMISE
# predict_proba_iai_flexi <- function(df, learner, vars){
#     times <- sort(unique(df$time)) # get times where events happened
#     train_X <- df[,vars] 
#     
#     mat_tree <- 1-as.matrix(do.call(cbind, lapply(times, function(x) iai::predict(learner, df, t=x))))
#     return(mat_tree)
# }
# 
# # TRAD - cc
# vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')
# # MMRD - cc
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed,
#         missingdatamode = "separate_class",
#         num_threads = 3,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:12,
# )
# df <- D_cc %>%
#     filter(ProMisE=='MMRd') %>%
#     as.data.frame
# train_died <- as.logical(df[,'death'])
# train_times <- df[,'time']
# train_X <- df[,vars] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_trad_mmrd_cc <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_trad_mmrd_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # cc 0.8693982 with depth 11
# mat_trad_mmrd_cc <- predict_proba_iai_flexi(df, best_learner_trad_mmrd_cc, vars=vars)
# save(mat_trad_mmrd, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_mmrd_cc.Rdata")
# 
# # NSMP cc
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed,
#         missingdatamode = "separate_class",
#         num_threads = 3,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:12,
# )
# df <- D_cc %>%
#     filter(ProMisE=='NSMP') %>%
#     as.data.frame
# train_died <- as.logical(df[,'death'])
# train_times <- df[,'time']
# train_X <- df[,vars] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_trad_nsmp_cc <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_trad_nsmp_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8879802 with depth 9
# mat_trad_nsmp_cc <- predict_proba_iai_flexi(df, best_learner_trad_nsmp_cc, vars=vars)
# save(mat_trad_nsmp_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_nsmp_cc.Rdata")
# 
# # p53ab cc
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed,
#         missingdatamode = "separate_class",
#         num_threads = 3,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:12,
# )
# df <- D_cc %>%
#     filter(ProMisE=='p53ab') %>%
#     as.data.frame
# train_died <- as.logical(df[,'death'])
# train_times <- df[,'time']
# train_X <- df[,vars] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_trad_p53ab_cc <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_trad_p53ab_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.758492 with depth 4
# mat_trad_p53ab_cc <- predict_proba_iai_flexi(df, best_learner_trad_p53ab_cc, vars=vars)
# save(mat_trad_p53ab_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_p53ab_cc.Rdata")
# 
# # EXT - cc
# vars <- c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
#           'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')
# # MMRD - cc
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed+1,
#         missingdatamode = "separate_class",
#         num_threads = 3,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:11,
# )
# df <- D_cc %>%
#     filter(ProMisE=='MMRd') %>%
#     as.data.frame
# train_died <- as.logical(df[,'death'])
# train_times <- df[,'time']
# train_X <- df[,vars] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_ext_mmrd_cc <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_ext_mmrd_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.7713416 with depth 3
# mat_ext_mmrd <- predict_proba_iai_flexi(df, best_learner_var_mmrd, vars=vars)
# save(mat_ext_mmrd, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_mmrd.Rdata")
# 
# # NSMP impt
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed,
#         missingdatamode = "separate_class",
#         num_threads = 3,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:12,
# )
# df <- D %>%
#     filter(ProMisE=='NSMP') %>%
#     as.data.frame
# train_died <- as.logical(df[,'death'])
# train_times <- df[,'time']
# train_X <- df[,vars] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_ext_nsmp <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_ext_nsmp, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8455827 with depth 4
# mat_ext_nsmp <- predict_proba_iai_flexi(df, best_learner_ext_nsmp, vars=vars)
# save(mat_ext_nsmp, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_nsmp.Rdata")
# 
# # p53ab impt
# grid <- iai::grid_search(
#     iai::optimal_tree_survival_learner(
#         random_seed = seed,
#         missingdatamode = "separate_class",
#         num_threads = 3,
#         show_progress=TRUE,
#         smooth_survival_curves=TRUE
#     ),
#     criterion = 'loglikelihood',
#     max_depth = 1:12,
# )
# df <- D %>%
#     filter(ProMisE=='p53ab') %>%
#     as.data.frame
# train_died <- as.logical(df[,'death'])
# train_times <- df[,'time']
# train_X <- df[,vars] 
# 
# iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
# best_learner_ext_p53ab <- iai::get_learner(grid)
# iai::get_best_params(grid)
# rm(grid)
# iai::score(best_learner_ext_p53ab, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8434117 with depth 4
# mat_ext_p53ab <- predict_proba_iai_flexi(df, best_learner_ext_p53ab, vars=vars)
# save(mat_ext_p53ab, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_p53ab.Rdata")

### ###
# end #
### ###