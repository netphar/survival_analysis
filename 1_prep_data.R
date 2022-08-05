setwd('/Users/zagidull/Documents/fimm_files/survival_all/survival_loukavaara/data')
library(tidyverse)
library(mice)
library(ranger)
library(sjmisc)

set.seed(21183410)

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
std <- function(x) {
    # "division by 2 is because that puts them on the same scale and binary response vars.
    # If p is the probability of observing one, the maximum variance for a Bernoulli random variable 
    # is achieved when p = 1 / 2 and it is equal to p (1 - p) = 1/ 4
    # dividing by 2 standard deviation it will make the feature have a variance of 1/4.
    # It is a trick I learned from this book https://books.google.fi/books?id=lV3DIdV0F9AC"
    (x - mean(x)) / (2 * sd(x))
}

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

# recoding
## recode stages.
D <- dataset %>%
    mutate(stage = recode_factor(stage, 
                                 'IB'='I',
                                 'IA'='I',
                                 'II'='II',
                                 'IIIa'='III',
                                 'IIIb'='III',
                                 'IIIC1'='III',
                                 'IIIC2'='III',
                                 'IVb'='IV',
                                 .ordered = TRUE)) %>%
    as.data.frame()

## recode adjuvant_therapy the only thing we do is combine Brachy&Chemo into Chemo
D <- D %>%
    mutate(adjuvant_therapy = recode_factor(adjuvant_therapy, 
                                            'Brachy and Chemo'='Chemo',
                                            .ordered = FALSE
    )) %>%
    as.data.frame()

## recode postop_mayo_criteria
D <- D %>%
    mutate(postop_mayo_criteria = recode_factor(postop_mayo_criteria, 
                                                'High-Intermediate'='High',
                                                'Low-Intermediate'='High',
                                                .ordered = FALSE
    )) %>%
    as.data.frame()

## order myometrial_invasion to low-mid-high. make ordered
D <- D %>%
    mutate(myometrial_invasion = recode_factor(myometrial_invasion, 
                                               'MI <= 33%' = 'Low',
                                               '33% < MI <= 66%'='Mid',
                                               'MI > 66%'='High',
                                               .ordered = TRUE
    )) %>%
    as.data.frame()

D <- D %>%
    select(-c('postop_milwaukee_score','postop_helsinki_score', 'TCGAoma', 'risk_group', 'risk_factors')) %>%
    as.data.frame()

## save for complete case
# windsorize outliers at 0.99 as per https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal_RMD.md
cols <- colnames(D)[sapply(D, is.numeric)][6:13] # # hemoglobin hematocrit leukocytes thrombocytes creatinine hcg hcgbv Ca125
for (col in cols){
    var <- D[[col]]
    q <- quantile(var[!is.na(var)], .99, type = 1)
    D[[col]][!is.na(var)] <-  pmin(var[!is.na(var)], q)
}

# prepare CC file
D_cc <- D[!is.na(D$ProMisE), ]
tmp_cc <- D_cc[, 2:ncol(D_cc)] # removing ID column
# impute CC file
tmp_imp_cc <- mice::parlmice(tmp_cc, m = 120, n.core = 6, maxit = 10, n.imp.core = 20,
                          defaultMethod = c("rf", "rf", "polyreg", "polr"), visitSequence='monotone', rfPackage='ranger', 
                          print = FALSE, cl.type = "FORK", cluster.seed=343242)
tmp_imp_combined_cc <- sjmisc::merge_imputations(dat=tmp_cc, imp=tmp_imp_cc)
D_imp_cc <- cbind(D_cc[!(names(D_cc) %in% names(tmp_imp_combined_cc))], tmp_imp_combined_cc)
D_cc <- D_imp_cc %>%
    filter(death!=2) %>%
    filter(ProMisE!='POLE') %>%
    as.data.frame
write_csv(D_cc, file = 'output/D_cc_june20.csv')

# prepare imputed data (not CC)
# When the incomplete data are covariates in the analysis model, the analysis model outcome must be used to predict the missing covariate values. 
# Although this practice may seem counter-intuitive, it is in fact essential. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2998703/ 
tmp <- D[, 2:ncol(D)] # removing ID column
tmp_imp <- mice::parlmice(tmp, m = 120, n.core = 6, maxit = 10, n.imp.core = 20, ntree=100,
                          defaultMethod = c("rf", "rf", "polyreg", "polr"), visitSequence='monotone', rfPackage='ranger', 
                          print = FALSE, cl.type = "FORK", cluster.seed=343242)

### RF imputation does not predict anything for POLE class, while polyreg predicts some. so we go with polyreg, since it feels more realistic
# RF predicts less of NSMP and more of MMRd -> closer to reality. 
# NSMP class grows the most after imputation


# check how cont imputed values look like, e.g.
stripplot(tmp_imp, Ca125, xlab = "Imputation number")
# check how 5 random imputed rows look like for ProMisE
how_many_times = 20
holder <- vector('list', length = how_many_times)
for (times in 1:how_many_times){
    row = sample(1:nrow(tmp_imp$imp$ProMisE), size = 1)
    holder[[times]] <- table(unlist(tmp_imp$imp$ProMisE[row,], use.names = F))
}
do.call(rbind, holder)

# combine imputation results using mode or mean for cat or numeric. 
tmp_imp_combined <- sjmisc::merge_imputations(dat=tmp, imp=tmp_imp)
D_imp <- cbind(D[!(names(D) %in% names(tmp_imp_combined))], tmp_imp_combined)

# to check count tables
flat_table(D, ProMisE, death)
flat_table(D_imp, ProMisE, death)

# we gonna remove death=2 and POLE since noone died there
D <- D_imp %>%
    filter(death!=2) %>%
    filter(ProMisE!='POLE') %>%
    as.data.frame

write_csv(D, file = 'output/D_june20.csv')

# CONT VARS START
tmp <- D %>%
    select(-c('time','id')) %>%
    as.data.frame()

# compare mean and st.dev between death=0/death=1/death=2 groups
tmp %>%
    group_by(death) %>%
    summarise(across(where(is.numeric), list(mean=mean,sd=sd), na.rm = TRUE, .names="{.col}_{.fn}")) %>%
    round(., 2) %>%
    as.data.frame

# ttest on numeric columns start
f <- function(df, colname){
    p<-t.test( df[ (df$death==1), colname], df[ (df$death==0), colname])$p.value
    print(c(colname, p))
}
cols <- colnames(tmp[,sapply(tmp, is.numeric)])[-1] # get colnames for all numeric features and remove death, since it is our grouping var
for (c in cols){
    f(tmp,c)
}
# age, hemoglobin, hematocrit, leukocytes, thrombocytes, hcgbv and Ca125 are distributed differently between death=0/1 groups

# correlations between numeric columns, ie check collinearity. report >0.2
# age&creatinine=0.28, age&hcg=0.26, age&hcgbv=0.23, 
# hemoglobin&hematocrit=0.96, hemoglobin&thrombocytes=-0.28,hemoglobin&hcgbv=-0.2 
# hematocrit&thrombocytes=-0.23, hematocrit&hcgbv=-0.21,
# leukocytes&thrombocytes=0.36, leukocytes&Ca125=0.25,
# thrombocytes&Ca125=0.25, hcg&hcgbv=0.22, 
round(cor(tmp[!(tmp$death==2),][, sapply(tmp[!(tmp$death==2),], is.numeric) ], use='pairwise.complete.obs'), 2)
    ## care when adding both of the vars to a linear model

# i wanna figure out which of colinear cols are better for use in linear coxph,
# e.g. hcgbv *gonadotropin b* is better than hcg. hcgbv is preferred
mylogit <- glm(death ~ leukocytes, data=tmp[!(tmp$death==2),], family = "binomial")
mylogit1 <- glm(death ~ thrombocytes, data=tmp[!(tmp$death==2),], family = "binomial")
mylogit2 <- glm(death ~ leukocytes+thrombocytes, data=tmp[!(tmp$death==2),], family = "binomial") 
mylogit3 <- glm(death ~ leukocytes+thrombocytes+leukocytes*thrombocytes, data=tmp[!(tmp$death==2),], family = "binomial") 
# after checking logreg hcgbv is preferred to hcg *gonadotropin b*, hemoglobin is preferred to hematocrit, leukocytes is preferred to thrombocytes

##
# we ll probably proceed with age, hemoglobin, leukocytes, hcgbv, Ca125
## 

tmp <- df

# lets look at binary var associations. So factor variables with 2 levels
f<-function(x) {
    c(mean(x == "Yes" | x == "Positive" | x == "High risk", na.rm = TRUE),
      mean(x == "No" | x == "Negative" | x == "Low risk", na.rm = TRUE))
}

res0<- tmp %>%
    filter(death==0) %>%
    select(where(is.factor)) %>%
    select(where(~ nlevels(.x) <3)) %>%
    apply(., 2, f) %>%
    round(., 2) %>%
    data.frame(., row.names = c('0_yes','0_no'))

res1<- tmp %>%
    filter(death==1) %>%
    select(where(is.factor)) %>%
    select(where(~ nlevels(.x) <3)) %>%
    apply(., 2, f) %>%
    round(., 2) %>%
    data.frame(., row.names = c('1_yes','1_no'))
res <- t(rbind(res0, res1))

res[which(abs(res[,1]-res[,3]) > 0.1),] # difference between means
# large differences (>0.1, 0-1 max) between groups for
# CD171
# p53_ab
# lymphovascular_invasion
# cytology
# ER 
# PR
# relapse
# hyperplasia 
# diameter more than 2/3/5 cm
# preop_histology
# deep_myometrial_invasion
# TIL mod_or_abundant


# Associations between multilevel nominal/categorical data (nominal data is categorical, where cats are ordered), ordinal data is where cats are unordered
# following Measures of Association How to Choose? 10.1177/8756479308317006
res0<- tmp %>%
    filter(death==0) %>%
    select(where(is.factor)) %>%
    select(where(~ nlevels(.x) >2))
rho0<-matrix(nrow=ncol(res0), ncol=ncol(res0))

for (i in 1:(ncol(rho0))) {
    for (j in 1:(nrow(rho0))) {
        rho0[i,j]=DescTools::Lambda(res0[,i], res0[,j]) # kruskall lambda 0 to 1
    }
}

idx0 <- which(rho0>0.35, arr.ind = TRUE)

for (i in 1:length(colnames(res0)[idx0[,'row']])){
    if (colnames(res0)[idx0[i,'row']] != colnames(res0)[idx0[i,'col']])
    {print( c(colnames(res0)[idx0[i,'row']],colnames(res0)[idx0[i,'col']]))}
}

res1<- tmp %>%
    filter(death==1) %>%
    select(where(is.factor)) %>%
    select(where(~ nlevels(.x) >2))
rho1<-matrix(nrow=ncol(res1), ncol=ncol(res1))

for (i in 1:(ncol(rho1))) {
    for (j in 1:(nrow(rho1))) {
        rho1[i,j]=DescTools::Lambda(res1[,i], res1[,j]) # kruskall lambda 0 to 1
    }
}

idx1 <- which(rho1>0.35, arr.ind = TRUE)

for (i in 1:length(colnames(res1)[idx1[,'row']])){
    if (colnames(res1)[idx1[i,'row']] != colnames(res1)[idx1[i,'col']])
    {print( c(colnames(res1)[idx1[i,'row']],colnames(res1)[idx1[i,'col']]))}
}
# in death=0 and death=1
# strong (defined as above 0.35) association between
# "histological_subgroup" "final_histology"  
# "adjuvant_therapy" "stage"  

# remove POLE, since noone died from cancer. only other reason (death==2)
addmargins(table(D[which(D$death==0), 'ProMisE'], useNA = 'always'))
addmargins(table(D[which(D$death==1), 'ProMisE'], useNA = 'always'))

###
# imputation OLD - scrambled order of removing death=2 and POLE, no proper imp tuning, using std-scaled dataset
###
# we have not removed any data yet. 
# candidates for removal are death==2 (competing risks) and ProMisE==POLE (noone died)
# Q: in which order to impute 
#    on all data
#    or after removing death==2
#    or after removing POLE
#    or after removing both death==2 and POLE

# for now i will impute on all the data and then remove death==2 and POLE group
# tmp_DnP <- mice::mice(D, m = 2, maxit = 25, visitSequence='monotone', defaultMethod=c('rf', 'rf','rf','polr'))
# 
# tmp_DnP <- mice::parlmice(D, m = 5, n.core=1, maxit = 25, defaultMethod=c('rf', 'rf','logreg','polr'), print = FALSE, cl.type = "FORK", cluster.seed=343242)
# D_DnP <- sjmisc::merge_imputations(D, tmp_DnP)
# D_imp <- cbind(D[!(names(D) %in% names(D_DnP))], D_DnP)
# D_imputedDefault <- D_imp %>%
#     filter(death!=2) %>%
#     filter(ProMisE!='POLE') %>%
#     as.data.frame
# 
# D_std <- D_imputedDefault %>%
#     mutate(
#         age = std(age),
#         bmi = std(bmi),
#         hemoglobin = std(hemoglobin),
#         hematocrit = std(hematocrit),
#         leukocytes = std(leukocytes),
#         thrombocytes = std(thrombocytes),
#         creatinine = std(creatinine),
#         Ca125 = std(Ca125)
#     ) %>%
#     mutate(ProMisE = droplevels(ProMisE)) %>%
#     mutate(ProMisE = relevel(ProMisE, ref='NSMP')) %>%
#     as.data.frame
# 
# #write_csv(D_std, file = 'D_std_March11.csv')

