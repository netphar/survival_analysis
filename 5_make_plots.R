library(pacman)

pacman::p_load(survival, tidyverse)
load(file = "~/Documents/fimm_files/survival_all/survival_loukavaara/data/boots_results_all_1july.Rdata")

# f defs
get_cind <- function(temp1, name){
    df1 <- as.data.frame(t(temp1$cox[,1]))
    df1$row.names <- colnames(temp1$cox)[1]
    df1$type <- paste('cox',name,sep='_')
    df2 <- as.data.frame(t(temp1$tree[,1]))
    df2$row.names <- colnames(temp1$tree)[1]
    df2$type <- paste('tree',name,sep='_')
    return(rbind(df1, df2))
}

get_IBS <- function(temp1, name){
    df1 <- as.data.frame(t(temp1$cox[,2:5]))
    df1$row.names <- rownames(df1)
    df1$type <- paste('cox',name,sep='_')
    df2 <- as.data.frame(t(temp1$tree[,2:5]))
    df2$row.names <- rownames(df2)
    df2$type <- paste('tree',name,sep='_')
    return(rbind(df1, df2))
}

# PLOT
# cindex trad
df1 <- get_cind(holder_all$trad_mmrd, name='trad')
df2 <- get_cind(holder_all$trad_nsmp, name='trad')
df3 <- get_cind(holder_all$trad_p53, name='trad')
df4 <- get_cind(holder_all$trad, name='trad')
df_trad <- rbind(df1, df2, df3, df4)
#df_trad$type <- paste('trad', df_trad$type, sep="_")

# cindex ext
df1 <- get_cind(holder_all$ext_mmrd, name='ext')
df2 <- get_cind(holder_all$ext_nsmp, name='ext')
df3 <- get_cind(holder_all$ext_p53, name='ext')
df4 <- get_cind(holder_all$ext, name='ext')
df_ext <- rbind(df1, df2, df3, df4)
#df_ext$type <- paste('ext', df_ext$type, sep="_")

temp <- rbind(df_cind[df_cind$type=='cox_trad', c('type','row.names','means', 'upper','lower')], 
              df_cind[df_cind$type=='cox_ext',c('type','row.names','means','upper','lower')])
temp['diff'] <- temp$upper-temp$lower
rownames(temp) <- seq(1,8)
temp <- temp[c(4,1,2,3,8,5,6,7), -c(4,5)]
rownames(temp) <- seq(1,8)
write.csv(temp, file = "perf_cox.csv")
temp <- rbind(df_cind[df_cind$type=='tree_trad', c('type','row.names','means', 'upper','lower')], 
              df_cind[df_cind$type=='tree_ext',c('type','row.names','means','upper','lower')])
temp['diff'] <- temp$upper-temp$lower
rownames(temp) <- seq(1,8)
temp <- temp[c(4,1,2,3,8,5,6,7), -c(4,5)]
rownames(temp) <- seq(1,8)
write.csv(temp, file = "perf_tree.csv")

# trad + ext
df_cind <- rbind(df_trad, df_ext)
df_cind$row.names <- c('mmrd', 'mmrd', 'nsmp', 'nsmp', 'p53ab', 'p53ab', 'all', 'all',
                       'mmrd', 'mmrd', 'nsmp', 'nsmp', 'p53ab', 'p53ab', 'all', 'all'
                       )
# mmrd
temp <- df_cind[df_cind$row.names=='mmrd',]
temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
p_mmrd <- ggplot(temp, aes(x=type, y=means)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
    ylim(0, 1) +
    labs(title='MMRD', x='', y='C-ind') +
    theme_minimal()
# nsmp
temp <- df_cind[df_cind$row.names=='nsmp',]
temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
p_nsmp <- ggplot(temp, aes(x=type, y=means)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
    ylim(0, 1) +
    labs(title='NSMP', x='', y='C-ind') +
    theme_minimal()
# p53ab
temp <- df_cind[df_cind$row.names=='p53ab',]
temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
p_p53ab <- ggplot(temp, aes(x=type, y=means)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
    ylim(0, 1) +
    labs(title='p53ab', x='', y='C-ind') +
    theme_minimal()
# all
temp <- df_cind[df_cind$row.names=='all',]
temp$type <- factor(temp$type, levels = c('cox_trad','tree_trad','cox_ext','tree_ext'))
p_all <- ggplot(temp, aes(x=type, y=means)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(x=type, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
    ylim(0, 1) +
    labs(title='all', x='', y='C-ind') +
    theme_minimal()
cowplot::plot_grid(p_mmrd, p_nsmp, p_p53ab, p_all)



# IBS trad by ProMisE and for full cohort

df1 <- get_IBS(holder_all$trad_mmrd, name='mmrd')
df2 <- get_IBS(holder_all$trad_nsmp, name='nsmp')
df3 <- get_IBS(holder_all$trad_p53, name='p53ab')
df4 <- get_IBS(holder_all$trad, name='all')
df_trad <- rbind(df1, df2, df3, df4)

df_trad$type <- factor(df_trad$type, levels = c('cox_all','tree_all','cox_mmrd','tree_mmrd','cox_nsmp','tree_nsmp','cox_p53ab','tree_p53ab'))

p_trad <- ggplot(df_trad, aes(x=row.names, y=means, fill=type)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
    scale_fill_manual(values=c(
        "#C58E03",
        "#137A81",
        "#F1A319",
        "#409893",
        "#F7D200",
        "#73C6B6",
        "#FCEEC5",
        "#D2F6E6"
    )) +
    ylim(0, 0.2) +
    labs(title='Feature set I', x='', y='IBS') +
    theme_minimal() +
    theme(text = element_text(size=20))
#p_trad

# IBS ext by ProMisE and for full cohort
df1 <- get_IBS(holder_all$ext_mmrd, name='mmrd')
df2 <- get_IBS(holder_all$ext_nsmp, name='nsmp')
df3 <- get_IBS(holder_all$ext_p53, name='p53ab')
df4 <- get_IBS(holder_all$ext, name='all')

df_ext <- rbind(df1, df2, df3, df4)
df_ext$type <- factor(df_ext$type, levels = c('cox_all', 'tree_all','cox_mmrd','tree_mmrd','cox_nsmp','tree_nsmp','cox_p53ab','tree_p53ab'))
p_ext <- ggplot(df_ext, aes(x=row.names, y=means, fill=type)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    geom_errorbar(aes(x=row.names, ymin=lower, ymax=upper), width=0.4, position=position_dodge(.9)) + 
    scale_fill_manual(values=c(
        "#C58E03",
        "#137A81",
        "#F1A319",
        "#409893",
        "#F7D200",
        "#73C6B6",
        "#FCEEC5",
        "#D2F6E6"
    )) +
    ylim(0, 0.2) +
    labs(title='Feature set II', x='', y='IBS') +
    theme_minimal() +
    theme(legend.position = 'none', text = element_text(size=20))
#p_ext
prow <- cowplot::plot_grid(p_trad+theme(legend.position = 'none'), p_ext)
#prow
legend <- cowplot::get_legend(
    # create some space to the left of the legend
    p_trad + theme(legend.box.margin = margin(0, 0, 0, 12))
)
cowplot::plot_grid(prow, legend, rel_widths = c(2, .4))

temp <- cbind(df_trad[,c('means','type')], df_ext[,c('means', 'row.names')])
rownames(temp) <- NULL
temp <- temp[,c(4,2,1,3)]
colnames(temp)[c(3,4)] <- c('IBS FSI', 'IBS FSII')
temp$diff_perc <- round((temp$`IBS FSI`-temp$`IBS FSII`)/temp$`IBS FSI` *100, 3)
write.csv(temp, file = "ibs.csv")
####

p + labs(title='MMRd traditional vars', x='IBS', y='timepoint') +
    theme_minimal()

cind <-c((1-perf_all[2,1]/perf_all[1,1])*100,
         (1-perf_all[4,1]/perf_all[3,1])*100,
         (1-perf_all[6,1]/perf_all[5,1])*100,
         (1-perf_all[8,1]/perf_all[7,1])*100)

data <- rbind(
    (1-perf_all[2,2:5]/perf_all[1,2:5])*100,
    (1-perf_all[4,2:5]/perf_all[3,2:5])*100,
    (1-perf_all[6,2:5]/perf_all[5,2:5])*100,
    (1-perf_all[8,2:5]/perf_all[7,2:5])*100)
data <- cbind(data, cind)
data <- as.data.frame(data)
rownames(data) <- c('cc trad Cox vs IAI',
                    'cc ext Cox vs IAI', 
                    'imp trad Cox vs IAI', 
                    'imp ext Cox vs IAI')
data <- reshape2::melt(as.matrix(data))

ggplot(data, aes(fill=Var1, y=value, x=Var2)) + 
    geom_bar(position="dodge", stat="identity") + 
    scale_fill_brewer(palette="Paired") + 
    labs(x ="", y = "change in %", fill = '', title = 'overall') + 
    theme_bw() +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5))
ggsave('~/Desktop/compare_all.png',height = 10, units = 'in')
