# limit use Julia and R in the same notebook
rm(list=ls(all=TRUE))
setwd('/Users/zagidull/Documents/fimm_files/survival_all/survival_loukavaara/data')

library(pacman)
pacman::p_load(tidyverse, survival, iai)
num_threads = 3
Sys.setenv(JULIA_NUM_THREADS=num_threads)
seed = 897978

#####
# function defs
predict_proba_iai <- function(df, learner, trad=T){ # this is for a full cohort model, ie *not* stratified by ProMisE
    times <- sort(unique(df$time)) 
    if(trad==T){
        train_X <- df[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')] 
    }else{
        train_X <- df[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
                        'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')] 
    }
    
    mat_tree <- 1-as.matrix(do.call(cbind, lapply(times, function(x) iai::predict(learner, df, t=x))))
    return(mat_tree)
}

predict_proba_iai_flexi <- function(df, learner, vars){ # this is for cohorts *statified* by ProMisE
    times <- sort(unique(df$time)) 
    train_X <- df[,vars] 
    
    mat_tree <- 1-as.matrix(do.call(cbind, lapply(times, function(x) iai::predict(learner, df, t=x))))
    return(mat_tree)
}

# load dataset
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

####
# optimal survival forest training
# changing seed makes a *big difference* in terms of validation performance (cindex)

#####
# traditional full cohort
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 897978,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 6:14,
)
train_died <- as.logical(D[,'death'])
train_times <- D[,'time']
train_X <- D[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad <- iai::get_learner(grid)
iai::get_params(best_learner_trad)$max_depth
rm(grid)
iai::score(best_learner_trad, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8620492 depth=10 seed 897978 D_june14

mat_trad <- predict_proba_iai(D, best_learner_trad)
save(mat_trad, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_28june.Rdata")

#####
# extended full cohort
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 44898,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 5:12,
)
train_died <- as.logical(D[,'death'])
train_times <- D[,'time']
train_X <- D[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
                'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')] 
# seed 56123 cindex 0.8817923, seed 44898 cindex 0.9111172
# for (i in 16:20){
#     s = floor(seed / i)
#     grid <- iai::grid_search(
#         iai::optimal_tree_survival_learner(
#             random_seed = s,
#             missingdatamode = "separate_class",
#             num_threads = 3,
#             show_progress=F,
#             smooth_survival_curves=TRUE
#         ),
#         criterion = 'loglikelihood',
#         max_depth = 3:12,
#     )
#     iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
#     best_learner_ext <- iai::get_learner(grid)
#     score <- iai::score(best_learner_ext, train_X, train_died, train_times, criterion='harrell_c_statistic')
#     print(list(seed=s, cind=score)) # trad 0.8484507 /// depth 6 and cindex 0.8601328)
# }

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext <- iai::get_learner(grid)
iai::get_params(best_learner_ext)$max_depth
rm(grid)
iai::score(best_learner_ext, train_X, train_died, train_times, criterion='harrell_c_statistic') # seed 44898 cindex 0.9111172 depth 11
mat_ext <- predict_proba_iai(D, best_learner_ext)
save(mat_ext, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_28june.Rdata")

#####
# stratify by PROMISE before training trees

# TRAD - definitions
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')
######
# MMRD - impt
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 897978,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='MMRd') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_mmrd <- iai::get_learner(grid)
iai::get_params(best_learner_trad_mmrd)$max_depth
iai::score(best_learner_trad_mmrd, train_X, train_died, train_times, criterion='harrell_c_statistic') # cindex 0.8698862 depth 11 seed 897978
rm(grid)
mat_trad_mmrd <- predict_proba_iai_flexi(df, best_learner_trad_mmrd, vars=vars)
save(mat_trad_mmrd, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_mmrd_29june.Rdata")

#####
# NSMP impt
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 897976,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='NSMP') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_nsmp <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_trad_nsmp, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.9164425 with depth 11 seed 897976
mat_trad_nsmp <- predict_proba_iai_flexi(df, best_learner_trad_nsmp, vars=vars)
save(mat_trad_nsmp, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_nsmp_29june.Rdata")

#####
# p53ab impt
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed-8,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='p53ab') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_p53ab <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_trad_p53ab, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.793953 with depth 7 seed 897978 or 897970?
mat_trad_p53ab <- predict_proba_iai_flexi(df, best_learner_trad_p53ab, vars=vars)
save(mat_trad_p53ab, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_p53ab_29june.Rdata")

# lets look at the promise cat with CPH
fx_trad_promise <-Surv(time, death) ~ age+stage+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm
dat <- df[,all.vars(fx_trad_promise)]
dat <- droplevels(dat)
summary(coxph(fx_trad_promise, data = dat, x = T, y = T))


#####
# EXT - definitions
vars <- c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
          'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')
#####
# MMRD ext
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 897958,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='MMRd') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 
iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_mmrd <- iai::get_learner(grid)
iai::get_best_params(grid)
iai::score(best_learner_ext_mmrd, train_X, train_died, train_times, criterion='harrell_c_statistic') # depth 7 seed 897958 cindex 0.9033342// ext 0.8637999 with depth 7 seed 897976
rm(grid)
mat_ext_mmrd <- predict_proba_iai_flexi(df, best_learner_ext_mmrd, vars=vars)
save(mat_ext_mmrd, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_mmrd_29june.Rdata")

#####
# NSMP ext
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 897964,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='NSMP') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_nsmp <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_ext_nsmp, train_X, train_died, train_times, criterion='harrell_c_statistic') # ext 0.9213386 with depth 7 seed 897964
mat_ext_nsmp <- predict_proba_iai_flexi(df, best_learner_ext_nsmp, vars=vars)
save(mat_ext_nsmp, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_nsmp_29june.Rdata")

#####
# p53ab ext
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = 897958,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=F,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='p53ab') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_p53ab <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_ext_p53ab, train_X, train_died, train_times, criterion='harrell_c_statistic') # ext 0.8990295 with depth 8 andd seed 897958.  trad 0.8434117 with depth 4 // 897970 or 897972, 897967
mat_ext_p53ab <- predict_proba_iai_flexi(df, best_learner_ext_p53ab, vars=vars)
save(mat_ext_p53ab, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_p53ab_29june.Rdata")

# lets look at the promise cat with CPH
fx_ext_promise <-Surv(time, death) ~ age+stage+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm+diameter_more_5cm+cytology+ER+CD171
dat <- df[,all.vars(fx_ext_promise)]
dat$ER <- relevel(dat$ER, ref='Negative')
dat <- dat %>%
    mutate(stage = fct_relevel(stage, c('I','II','III','IV') )) %>%
    as.data.frame()

summary(coxph(fx_ext_promise, data = dat, x = T, y = T))

######
# CC
D_cc <- read_csv(file = '~/Desktop/D_cc_june20.csv', 
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

# traditional optimal survival forest CC
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed+1,
        missingdatamode = "separate_class",
        num_threads = 1,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:11,
)
train_died <- as.logical(D_cc[,'death'])
train_times <- D_cc[,'time']
train_X <- D_cc[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_cc <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_trad_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # 0.8635481 with depth 8

mat_trad_cc <- predict_proba_iai(D_cc, best_learner_trad_cc)
save(mat_trad_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_cc.Rdata")

# extended optimal survival forest CC
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed-1,
        missingdatamode = "separate_class",
        num_threads = 1,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:11,
)
train_died <- as.logical(D_cc[,'death'])
train_times <- D_cc[,'time']
train_X <- D_cc[,c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
                   'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')] 
iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_cc <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_ext_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # 0.8365357 with depth 10

mat_ext_cc <- predict_proba_iai(D_cc, best_learner_ext_cc)
save(mat_ext_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_cc.Rdata")

## CC
### TRAIN WITHOUT PROMISE
predict_proba_iai_flexi <- function(df, learner, vars){
    times <- sort(unique(df$time)) # get times where events happened
    train_X <- df[,vars] 
    
    mat_tree <- 1-as.matrix(do.call(cbind, lapply(times, function(x) iai::predict(learner, df, t=x))))
    return(mat_tree)
}

# TRAD - cc
vars <- c('age', 'stage', 'histological_subgroup','deep_myometrial_invasion','lymphovascular_invasion','diameter_more_3cm')
# MMRD - cc
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D_cc %>%
    filter(ProMisE=='MMRd') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_mmrd_cc <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_trad_mmrd_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # cc 0.8693982 with depth 11
mat_trad_mmrd_cc <- predict_proba_iai_flexi(df, best_learner_trad_mmrd_cc, vars=vars)
save(mat_trad_mmrd, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_mmrd_cc.Rdata")

# NSMP cc
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D_cc %>%
    filter(ProMisE=='NSMP') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_nsmp_cc <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_trad_nsmp_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8879802 with depth 9
mat_trad_nsmp_cc <- predict_proba_iai_flexi(df, best_learner_trad_nsmp_cc, vars=vars)
save(mat_trad_nsmp_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_nsmp_cc.Rdata")

# p53ab cc
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D_cc %>%
    filter(ProMisE=='p53ab') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_trad_p53ab_cc <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_trad_p53ab_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.758492 with depth 4
mat_trad_p53ab_cc <- predict_proba_iai_flexi(df, best_learner_trad_p53ab_cc, vars=vars)
save(mat_trad_p53ab_cc, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_trad_p53ab_cc.Rdata")

# EXT - cc
vars <- c('age', 'stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
          'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')
# MMRD - cc
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed+1,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:11,
)
df <- D_cc %>%
    filter(ProMisE=='MMRd') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_mmrd_cc <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_ext_mmrd_cc, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.7713416 with depth 3
mat_ext_mmrd <- predict_proba_iai_flexi(df, best_learner_var_mmrd, vars=vars)
save(mat_ext_mmrd, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_mmrd.Rdata")

# NSMP impt
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='NSMP') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_nsmp <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_ext_nsmp, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8455827 with depth 4
mat_ext_nsmp <- predict_proba_iai_flexi(df, best_learner_ext_nsmp, vars=vars)
save(mat_ext_nsmp, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_nsmp.Rdata")

# p53ab impt
grid <- iai::grid_search(
    iai::optimal_tree_survival_learner(
        random_seed = seed,
        missingdatamode = "separate_class",
        num_threads = 3,
        show_progress=TRUE,
        smooth_survival_curves=TRUE
    ),
    criterion = 'loglikelihood',
    max_depth = 1:12,
)
df <- D %>%
    filter(ProMisE=='p53ab') %>%
    as.data.frame
train_died <- as.logical(df[,'death'])
train_times <- df[,'time']
train_X <- df[,vars] 

iai::fit(grid, train_X, train_died, train_times, validation_criterion='harrell_c_statistic')
best_learner_ext_p53ab <- iai::get_learner(grid)
iai::get_best_params(grid)
rm(grid)
iai::score(best_learner_ext_p53ab, train_X, train_died, train_times, criterion='harrell_c_statistic') # trad 0.8434117 with depth 4
mat_ext_p53ab <- predict_proba_iai_flexi(df, best_learner_ext_p53ab, vars=vars)
save(mat_ext_p53ab, file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/mat_ext_p53ab.Rdata")




### ###
# end #
### ###
