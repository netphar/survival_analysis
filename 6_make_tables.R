setwd('/Users/zagidull/Documents/fimm_files/survival_all/survival_loukavaara/data')

library(pacman)
pacman::p_load(tidyverse, broom, tableone, survival)

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

df <- D[,c('age', 'death','stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
           'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')]
vars <- c('age','stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
          'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','CD171','cytology','ER')
df$ER <- relevel(df$ER, ref='Negative')
out <- print(CreateTableOne(vars=vars, strata=c('death'), data=df, test=F), 
             quote=F, 
             explain=T, 
             nonnormal=c('age'), 
             minmax=T,
             showAllLevels=T,
             noSpaces=T,
             dropEqual=T)
write.csv(out, file = "myTable1.csv")

df <- D[,c('age','time', 'death','stage', 'histological_subgroup','ProMisE','deep_myometrial_invasion',
           'lymphovascular_invasion','diameter_more_3cm','diameter_more_5cm','cytology','ER','CD171')]
df$ER <- relevel(df$ER, ref='Negative')
df <- df %>%
    mutate(stage = fct_relevel(stage, c('I','II','III','IV') )) %>%
    as.data.frame()
df <- df %>%
    mutate(ProMisE = fct_relevel(ProMisE, c('NSMP','MMRd','p53ab') )) %>%
    as.data.frame()
fx_trad <-Surv(time, death) ~ age+stage+ProMisE+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm # trad variable set
fx_ext <-Surv(time, death) ~ age+stage+ProMisE+histological_subgroup+deep_myometrial_invasion+lymphovascular_invasion+diameter_more_3cm+diameter_more_5cm+cytology+ER+CD171 # extended variable set


t <- tidy(coxph(fx_ext, data = df))[,c(1,2,5)]
s <- summary(coxph(fx_ext, data = df))
s <- s$conf.int[,c(1,3,4)]
rownames(s) <- NULL
a <- cbind(t,s)
colnames(a)[2] <- c('coef')
rownames(a) <- NULL
a <- a[,c(1,2,4,5,6,3)]
a[,c(2,3,4,5)] <- round(a[,c(2,3,4,5)], digits = 2)
write.csv(a, file = "cox_regr.csv")

t <- tidy(coxph(fx_trad, data = df))[,c(1,2,5)]
s <- summary(coxph(fx_trad, data = df))
s <- s$conf.int[,c(1,3,4)]
rownames(s) <- NULL
a <- cbind(t,s)
colnames(a)[2] <- c('coef')
rownames(a) <- NULL
a <- a[,c(1,2,4,5,6,3)]
a[,c(2,3,4,5)] <- round(a[,c(2,3,4,5)], digits = 2)
write.csv(a, file = "cox_regr_trad.csv")
