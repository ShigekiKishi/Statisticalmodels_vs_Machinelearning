library(tidyverse)
library(ggsci)


kishi <- read.table('kishi/rmse_data201001.csv', sep = ',', header = TRUE)
colnames(kishi) <- c('id', 'crop', 'disease', 'disease_name', 'norm_mu', 'null_mu', 'norm_sd', 'null_sd', 'z_score', 'samplesize')
kishi$disease <- paste(kishi$crop, kishi$disease, sep = '__')



jsun.1 <- read.table('model/results/cucumber/valid_summary.tsv', sep = '\t', header = FALSE)
jsun.2 <- read.table('model/results/eggplant/valid_summary.tsv', sep = '\t', header = FALSE)
jsun.3 <- read.table('model/results/tomato/valid_summary.tsv', sep = '\t', header = FALSE)
jsun.4 <- read.table('model/results/strawberry/valid_summary.tsv', sep = '\t', header = FALSE)
jsun <- rbind(jsun.1, jsun.2, jsun.3, jsun.4)
colnames(jsun) <- c('disease', 'norm_mu', 'norm_sd', 'null_mu',  'null_sd',  'z_score')
jsun$disease <- tolower(jsun$disease)


d <- dplyr::full_join(jsun, kishi, by = 'disease', suffix = c('.sun', '.kishi'))
d$crop <- sapply(strsplit(d$disease, '__'), '[', 1)
d$disease


# Z-score vs sample size
df <- rbind(data.frame(disease = d$disease, zscore = d$z_score.sun, rmse = d$norm_mu.sun,
                       n = d$samplesize, crop = d$crop, method = 'NN'),
            data.frame(disease = d$disease, zscore = d$z_score.kishi, rmse = d$norm_mu.kishi,
                       n = d$samplesize, crop = d$crop, method = 'Bayes'))

pl <- ggplot(df, aes(x = n, y = zscore, color = crop)) + 
            geom_point() + scale_color_nejm() +
            facet_grid(~ method)
png('kishi/compresult/zscores_scatter.png', 1700, 800, res = 250)
print(pl)
dev.off()
 
pl <- ggplot(df, aes(x = rmse, y = zscore, color = crop)) + 
            geom_point() + scale_color_nejm() +
            facet_grid(~ method)
png('kishi/compresult/zscores_scatter2.png', 1700, 800, res = 250)
print(pl)
dev.off()
 


df <- data.frame(disease = d$disease, crop = d$crop,
                 difference = d$norm_mu.kishi - d$norm_mu.sun,
                 n = d$samplesize)
df <- df[!is.na(df$difference), ]
df <- df[order(df$difference), ]
df$disease <- factor(df$disease, levels = df$disease)
px <- ggplot(df, aes(x = n, y = difference, color = crop)) +
      geom_point() + scale_color_nejm() +
      ylab('RMSE (Bayes) − RMSE (NN)') 
png('kishi/compresult/rmsediff_scatter.png', 1200, 1000, res = 250)
print(px)
dev.off()


df <- data.frame(disease = d$disease, crop = d$crop,
                 difference = d$null_mu.kishi - d$null_mu.sun,
                 n = d$samplesize)
df <- df[!is.na(df$difference), ]
df <- df[order(df$difference), ]
df$disease <- factor(df$disease, levels = df$disease)
px <- ggplot(df, aes(x = n, y = difference, color = crop)) +
      geom_point() + scale_color_nejm() +
      ylab('RMSE (Bayes) − RMSE (NN)') 
png('kishi/compresult/rmsediff_scatter-null.png', 1200, 1000, res = 250)
print(px)
dev.off()


df <- data.frame(disease = d$disease, crop = d$crop,
                 null_diff = d$null_mu.kishi - d$null_mu.sun,
                 norm_diff = d$norm_mu.kishi - d$norm_mu.sun,
                 n = d$samplesize)
df$disease <- factor(df$disease, levels = df$disease)
px <- ggplot(df, aes(x = null_diff, y = norm_diff, color = crop)) +
      geom_point() + scale_color_nejm() +
      ylab('RMSE (Bayes) − RMSE (NN) [normal]') + xlab('RMSE (Bayes) - RMSE (NN) [null]')
png('kishi/compresult/rmsediff_scatter-normnull.png', 1200, 1000, res = 250)
print(px)
dev.off()








px <- ggplot(d, aes(x = norm_mu.sun, y = norm_mu.kishi, color = samplesize)) +
        geom_point() + coord_fixed() + 
        xlab('RMSE (NN)') + ylab('RMSE (Bayes)')
png('kishi/compresult/rmse_scatter.png', 1200, 1000, res = 250)
print(px)
dev.off()
px <- ggplot(d, aes(x = norm_mu.sun, y = norm_mu.kishi, color = samplesize)) +
        geom_point() + coord_fixed() + scale_x_log10() + scale_y_log10()+
        xlab('RMSE (NN)') + ylab('RMSE (Bayes)')
png('kishi/compresult/rmse_scatter_log10.png', 1200, 1000, res = 250)
print(px)
dev.off()









for (crop in c('cucumber', 'tomato', 'strawberry', 'eggplant')) {
    figprefix = paste0('kishi/compresult/', crop, '__')

    d_crop <- d %>% filter(crop == !!crop)
    
   
    df <- data.frame(disease = d_crop$disease,
                     difference = d_crop$norm_mu.kishi - d_crop$norm_mu.sun,
                     n = d_crop$samplesize)
    df <- df[!is.na(df$difference), ]
    df <- df[order(df$difference), ]
    df$disease <- factor(df$disease, levels = df$disease)
    pl <- ggplot(df, aes(x = disease, y = difference)) +
            geom_bar(stat = 'identity') +
            ylab('RMSE (Bayes) − RMSE (NN)') +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
            coord_flip()
    
    png(paste0(figprefix, 'rmsediff.png'), 1400, 1200, res = 200)
    print(pl)
    dev.off()

    px <- ggplot(df, aes(x = n, y = difference)) +
            geom_point() +
            ylab('RMSE (Bayes) − RMSE (NN)') 
    png(paste0(figprefix, 'rmsediff_samplesize.png'), 1400, 1200, res = 200)
    print(px)
    dev.off()

}







