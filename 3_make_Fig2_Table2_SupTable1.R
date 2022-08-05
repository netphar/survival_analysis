rm(list=ls(all=TRUE))
setwd('/Users/zagidull/Documents/fimm_files/survival_all/survival_loukavaara/data')
library(pacman)

pacman::p_load(survival, tidyverse)
load(file = "output/boots_results_all_4aug.Rdata")

# func defs
get_cind_new <- function(l_input){
    ns <- names(l_input)
    h <- list()
    for(n in ns){
        df <- as.data.frame(t(l_input[[n]]$cind))
        r1 <- unlist(strsplit(n,"_"))[1]
        df$row.names <- rep(r1, 2)
        type_pre <- rownames(df)
        type_post <- unlist(strsplit(n,"_"))[2]
        df$type <- paste(type_pre, type_post, sep = '_')
        h[[n]] <- df
    }
    out <- do.call(rbind,h)
    rownames(out) <- NULL
    out
}
get_IBS_new <- function(l_input, name){
    type_post <- unlist(strsplit(name, '_'))[1]
    ns <- names(l_input[[name]])[2:length(names(l_input[[name]]))]
    h <- list()
    for(n in ns){
        df <- as.data.frame(t(l_input[[name]][[n]]))
        #r1 <- unlist(strsplit(n,"_"))[1]
        df$row.names <- rep(n, 3)
        type_pre <- rownames(df)
        
        df$type <- paste(type_pre, type_post, sep = '_')
        h[[n]] <- df
    }
    out <- do.call(rbind,h)
    rownames(out) <- NULL
    out
}
make_barplot <- function(df_t, df_e, name, legend_position='none'){
    test1 <- df_t[str_detect(df_t$type, name), ]
    test1$type <- base::sapply(as.character(test1$type), 
                               function(x){ 
                                   paste(str_to_title(unlist(strsplit(x, '_'))[1]), 'FSI', sep=' ') 
                               }, 
                               USE.NAMES = F)
    test2 <- df_e[str_detect(df_e$type, name), ]
    test2$type <- base::sapply(as.character(test2$type), 
                               function(x){ 
                                   paste(str_to_title(unlist(strsplit(x, '_'))[1]), 'FSII', sep=' ') 
                               }, 
                               USE.NAMES = F)
    test <- rbind(test1,  test2)
    test <- test[!(test$type=='Reference FSII'),] # remove duplicated ref scores. They are identical for both FSI and FSII since we use K-M estimator to get them
    test[test$type=='Reference FSI','type'] <- rep('KM', 4)
    test$type <- factor(test$type, levels = c('KM','Cox FSI','Tree FSI','Cox FSII','Tree FSII'))
    test$row.names <- str_replace_all(test$row.names, "ibs_", "")
    rownames(test) <- NULL
    if(name=='all'){
        title <- 'Full Cohort'
    } else if(name=='mmrd') {
        title <- 'MMRd'
    } else if(name=='nsmp'){
        title <- 'NSMP'
    } else {
        title <- 'p53ab'
    }
    p <- ggplot(test, aes(x=row.names, y=means, fill=type)) +
        geom_bar(stat="identity", color="black", position=position_dodge()) +
        geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
        scale_fill_manual(values=c(
            'red',
            "#FCEEC5",
            "#409893",
            "#F1A319",
            "#577A81"
        )) +
        ylim(0, 0.25) +
        labs(title=title, x='', y='IBS') +
        theme_minimal() +
        theme(legend.position = legend_position, legend.title= element_blank(), text = element_text(size=20))
    p
}

# get cindex
df_cind <- get_cind_new(holder)

# get IBS
h_trad <- list()
h_ext <- list()
for(n in names(holder)){
    if(str_detect(n,'_trad')){
        h_trad[[n]] <- get_IBS_new(holder, name=n)
    } else {
        h_ext[[n]] <- get_IBS_new(holder, name=n)
    }
}
df_trad <- do.call(rbind, h_trad)
df_ext <- do.call(rbind, h_ext)

rownames(df_trad) <- NULL
rownames(df_ext) <- NULL

df_trad$type<- factor(df_trad$type, 
                      levels = c('Reference_all','cox_all','tree_all',
                                 'Reference_mmrd','cox_mmrd','tree_mmrd',
                                 'Reference_nsmp','cox_nsmp','tree_nsmp',
                                 'Reference_p53ab','cox_p53ab','tree_p53ab'))
df_ext$type<- factor(df_ext$type, 
                     levels = c('Reference_all','cox_all','tree_all',
                                'Reference_mmrd','cox_mmrd','tree_mmrd',
                                'Reference_nsmp','cox_nsmp','tree_nsmp',
                                'Reference_p53ab','cox_p53ab','tree_p53ab'))
rm(h_ext)
rm(h_trad)

# lets make a panel of four barplots 
# all/p53/mmrd/nsmp with FSI/FSII var sets for cox/ost and a single ref (negative control) in each
p_full <- make_barplot(df_trad, df_ext, name='all', legend_position = 'bottom')
p_mmrd <- make_barplot(df_trad, df_ext, name='mmrd')    
p_nsmp <- make_barplot(df_trad, df_ext, name='nsmp')        
p_p53ab <- make_barplot(df_trad, df_ext, name='p53ab')        
legend <- cowplot::get_legend(p_full)
p_full <- p_full + theme(legend.position="none")    
    
p_all <- gridExtra::grid.arrange(p_full, p_mmrd, p_nsmp, p_p53ab, legend, ncol=2, nrow = 3, 
             layout_matrix = rbind(c(1,2), c(3,4), c(5,5)),
             widths = c(2.7, 2.7), heights = c(2.5, 2.5, 0.2))

ggsave("figure 2.png", device = 'png', p_all)

####
# make tables
## IBS
temp <- cbind(df_trad[,c('means','type')], df_ext[,c('means', 'row.names')])
temp <- temp[!str_detect(temp$type, 'Reference'), ]
rownames(temp) <- NULL
temp <- temp[,c(4,2,1,3)]
colnames(temp)[c(3,4)] <- c('IBS FSI', 'IBS FSII')
temp$diff_perc <- round((temp$`IBS FSI`-temp$`IBS FSII`)/temp$`IBS FSI` *100, 1) # get % improvement of the IBS. We take the initial (FSI) IBS as the basis

temp[str_detect(temp$type, 'tree'), ] %>%
    group_by(row.names) %>%
    summarize(mean_diff = mean(diff_perc)) # trees improve in % by ibs_1y 15.4, ibs_2y 21.6, ibs_5y 21.0, ibs_all 16.2
temp[str_detect(temp$type, 'cox'), ] %>%
    group_by(row.names) %>%
    summarize(mean_diff = mean(diff_perc)) # cox improves ini % by ibs_1y 5.65, ibs_2y 7.5 , ibs_5y 3.85, ibs_all 4.4 

write.csv(temp, file = "tables/ibs_change_FSIvsFSII.csv")

## cindex
temp <- rbind(df_cind[df_cind$type=='cox_trad', c('type','row.names','means', 'upper','lower')], 
              df_cind[df_cind$type=='cox_ext',c('type','row.names','means','upper','lower')])
temp['ci95'] <- temp$upper-temp$lower
rownames(temp) <- seq(1,8)
temp <- temp[c(4,1,2,3,8,5,6,7), -c(4,5)]
rownames(temp) <- seq(1,8)
temp$diff_EXTbetterTRAD <- c(rep(NA, 4), 
                    round((temp[temp$type=='cox_ext','means']-temp[temp$type=='cox_trad','means'])/temp[temp$type=='cox_ext','means']*100, 1))
cind_change_cox <- temp

temp <- rbind(df_cind[df_cind$type=='tree_trad', c('type','row.names','means', 'upper','lower')], 
              df_cind[df_cind$type=='tree_ext',c('type','row.names','means','upper','lower')])
temp['ci95'] <- temp$upper-temp$lower
rownames(temp) <- seq(1,8)
temp <- temp[c(4,1,2,3,8,5,6,7), -c(4,5)]
rownames(temp) <- seq(1,8)
temp$diff_EXTbetterTRAD <- c(rep(NA, 4), 
                    round((temp[temp$type=='tree_ext','means']-temp[temp$type=='tree_trad','means'])/temp[temp$type=='tree_ext','means']*100, 1))
cind_change_all <- rbind(cind_change_cox, temp)
cind_change_all$diff_TREEbetterCOX <- c(
    rep(NA, 4),
    rep(NA, 4),
    round((cind_change_all[cind_change_all$type =='tree_trad', 'means']-cind_change_all[cind_change_all$type =='cox_trad', 'means'])/cind_change_all[cind_change_all$type =='tree_trad', 'means']*100, 1),
    round((cind_change_all[cind_change_all$type =='tree_ext', 'means']-cind_change_all[cind_change_all$type =='cox_ext', 'means'])/cind_change_all[cind_change_all$type =='tree_ext', 'means']*100, 1)
)


write.csv(cind_change_all, file = "tables/cind_change.csv")



##########################################
# OLD for first draft (23 july)
#########################################

#load(file = "boots_results_all_1july.Rdata")
# get_cind <- function(temp1, name){
#     df1 <- as.data.frame(t(temp1$cox[,1]))
#     df1$row.names <- colnames(temp1$cox)[1]
#     df1$type <- paste('cox',name,sep='_')
#     df2 <- as.data.frame(t(temp1$tree[,1]))
#     df2$row.names <- colnames(temp1$tree)[1]
#     df2$type <- paste('tree',name,sep='_')
#     return(rbind(df1, df2))
# }
# get_IBS <- function(temp1, name){
#     df1 <- as.data.frame(t(temp1$cox[,2:5]))
#     df1$row.names <- rownames(df1)
#     df1$type <- paste('cox',name,sep='_')
#     df2 <- as.data.frame(t(temp1$tree[,2:5]))
#     df2$row.names <- rownames(df2)
#     df2$type <- paste('tree',name,sep='_')
#     return(rbind(df1, df2))
# }

# PLOT IBS
# p_trad <- ggplot(df_trad, aes(x=row.names, y=means, fill=type)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     scale_fill_manual(values=c(
#         'red',
#         "#FCEEC5",
#         "#137A81",
#         'red',
#         "#FCEEC5",
#         "#137A81",
#         'red',
#         "#FCEEC5",
#         "#137A81",
#         'red',
#         "#FCEEC5",
#         "#137A81"
#     )) +
#     ylim(0, 0.25) +
#     labs(title='Feature set I', x='', y='IBS') +
#     theme_minimal() +
#     theme(text = element_text(size=20))
# 
# p_ext <- ggplot(df_ext, aes(x=row.names, y=means, fill=type)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     scale_fill_manual(values=c(
#         'red',
#         "#FCEEC5",
#         "#137A81",
#         'red',
#         "#FCEEC5",
#         "#137A81",
#         'red',
#         "#FCEEC5",
#         "#137A81",
#         'red',
#         "#FCEEC5",
#         "#137A81"
#     )) +
#     ylim(0, 0.25) +
#     labs(title='Feature set II', x='', y='IBS') +
#     theme_minimal() +
#     theme(legend.position = 'none', text = element_text(size=20))
# #p_ext
# prow <- cowplot::plot_grid(p_trad+theme(legend.position = 'none'), p_ext)
# #prow
# legend <- cowplot::get_legend(
#     # create some space to the left of the legend
#     p_trad + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# cowplot::plot_grid(prow, legend, rel_widths = c(2, .4))
# 
# # PLOT
# # cindex trad
# df1 <- get_cind(holder_all$trad_mmrd, name='trad')
# df2 <- get_cind(holder_all$trad_nsmp, name='trad')
# df3 <- get_cind(holder_all$trad_p53, name='trad')
# df4 <- get_cind(holder_all$trad, name='trad')
# df_trad <- rbind(df1, df2, df3, df4)
# #df_trad$type <- paste('trad', df_trad$type, sep="_")
# 
# # cindex ext
# df1 <- get_cind(holder_all$ext_mmrd, name='ext')
# df2 <- get_cind(holder_all$ext_nsmp, name='ext')
# df3 <- get_cind(holder_all$ext_p53, name='ext')
# df4 <- get_cind(holder_all$ext, name='ext')
# df_ext <- rbind(df1, df2, df3, df4)
# #df_ext$type <- paste('ext', df_ext$type, sep="_")
# 
# # trad + ext
# df_cind <- rbind(df_trad, df_ext)
# df_cind$row.names <- c('mmrd', 'mmrd', 'nsmp', 'nsmp', 'p53ab', 'p53ab', 'all', 'all',
#                        'mmrd', 'mmrd', 'nsmp', 'nsmp', 'p53ab', 'p53ab', 'all', 'all'
#                        )
# # mmrd
# temp <- df_cind[df_cind$row.names=='mmrd',]
# temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
# p_mmrd <- ggplot(temp, aes(x=type, y=means)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     ylim(0, 1) +
#     labs(title='MMRD', x='', y='C-ind') +
#     theme_minimal()
# # nsmp
# temp <- df_cind[df_cind$row.names=='nsmp',]
# temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
# p_nsmp <- ggplot(temp, aes(x=type, y=means)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     ylim(0, 1) +
#     labs(title='NSMP', x='', y='C-ind') +
#     theme_minimal()
# # p53ab
# temp <- df_cind[df_cind$row.names=='p53ab',]
# temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
# p_p53ab <- ggplot(temp, aes(x=type, y=means)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     ylim(0, 1) +
#     labs(title='p53ab', x='', y='C-ind') +
#     theme_minimal()
# # all
# temp <- df_cind[df_cind$row.names=='all',]
# temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
# p_all <- ggplot(temp, aes(x=type, y=means)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     ylim(0, 1) +
#     labs(title='all', x='', y='C-ind') +
#     theme_minimal()
# cowplot::plot_grid(p_mmrd, p_nsmp, p_p53ab, p_all)
# 
# 
# 
# # IBS trad by ProMisE and for full cohort
# df1 <- get_IBS(holder_all$trad_mmrd, name='mmrd')
# df2 <- get_IBS(holder_all$trad_nsmp, name='nsmp')
# df3 <- get_IBS(holder_all$trad_p53, name='p53ab')
# df4 <- get_IBS(holder_all$trad, name='all')
# df_trad <- rbind(df1, df2, df3, df4)
# 
# df_trad$type <- factor(df_trad$type, levels = c('cox_all','tree_all','cox_mmrd','tree_mmrd','cox_nsmp','tree_nsmp','cox_p53ab','tree_p53ab'))
# 
# p_trad <- ggplot(df_trad, aes(x=row.names, y=means, fill=type)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     scale_fill_manual(values=c(
#         'red',
#         "#C58E03",
#         "#137A81",
#         'red',
#         "#F1A319",
#         "#409893",
#         'red',
#         "#F7D200",
#         "#73C6B6",
#         'red',
#         "#FCEEC5",
#         "#D2F6E6"
#     )) +
#     ylim(0, 0.35) +
#     labs(title='Feature set I', x='', y='IBS') +
#     theme_minimal() +
#     theme(text = element_text(size=20))
# #p_trad
# 
# # IBS ext by ProMisE and for full cohort
# df1 <- get_IBS(holder_all$ext_mmrd, name='mmrd')
# df2 <- get_IBS(holder_all$ext_nsmp, name='nsmp')
# df3 <- get_IBS(holder_all$ext_p53, name='p53ab')
# df4 <- get_IBS(holder_all$ext, name='all')
# 
# df_ext <- rbind(df1, df2, df3, df4)
# df_ext$type <- factor(df_ext$type, levels = c('cox_all', 'tree_all','cox_mmrd','tree_mmrd','cox_nsmp','tree_nsmp','cox_p53ab','tree_p53ab'))
# p_ext <- ggplot(df_ext, aes(x=row.names, y=means, fill=type)) +
#     geom_bar(stat="identity", color="black", position=position_dodge()) +
#     geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
#     scale_fill_manual(values=c(
#         'red',
#         "#C58E03",
#         "#137A81",
#         'red',
#         "#F1A319",
#         "#409893",
#         'red',
#         "#F7D200",
#         "#73C6B6",
#         'red',
#         "#FCEEC5",
#         "#D2F6E6"
#     )) +
#     ylim(0, 0.35) +
#     labs(title='Feature set II', x='', y='IBS') +
#     theme_minimal() +
#     theme(legend.position = 'none', text = element_text(size=20))
# #p_ext
# prow <- cowplot::plot_grid(p_trad+theme(legend.position = 'none'), p_ext)
# #prow
# legend <- cowplot::get_legend(
#     # create some space to the left of the legend
#     p_trad + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# cowplot::plot_grid(prow, legend, rel_widths = c(2, .4))
# 
# 
# ####
# 
# p + labs(title='MMRd traditional vars', x='IBS', y='timepoint') +
#     theme_minimal()
# 
# cind <-c((1-perf_all[2,1]/perf_all[1,1])*100,
#          (1-perf_all[4,1]/perf_all[3,1])*100,
#          (1-perf_all[6,1]/perf_all[5,1])*100,
#          (1-perf_all[8,1]/perf_all[7,1])*100)
# 
# data <- rbind(
#     (1-perf_all[2,2:5]/perf_all[1,2:5])*100,
#     (1-perf_all[4,2:5]/perf_all[3,2:5])*100,
#     (1-perf_all[6,2:5]/perf_all[5,2:5])*100,
#     (1-perf_all[8,2:5]/perf_all[7,2:5])*100)
# data <- cbind(data, cind)
# data <- as.data.frame(data)
# rownames(data) <- c('cc trad Cox vs IAI',
#                     'cc ext Cox vs IAI', 
#                     'imp trad Cox vs IAI', 
#                     'imp ext Cox vs IAI')
# data <- reshape2::melt(as.matrix(data))
# 
# ggplot(data, aes(fill=Var1, y=value, x=Var2)) + 
#     geom_bar(position="dodge", stat="identity") + 
#     scale_fill_brewer(palette="Paired") + 
#     labs(x ="", y = "change in %", fill = '', title = 'overall') + 
#     theme_bw() +
#     theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5))
# ggsave('~/Desktop/compare_all.png',height = 10, units = 'in')
