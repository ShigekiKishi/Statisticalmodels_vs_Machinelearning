library(tidyverse)
library(ggsci)

## set up path to a summary file
cv_sum_fpath <- 'cv_results/summary.tsv'


## set up feature type: `category` or `decimal`
##   category: month and prefecture name were converted to one-hot vectors during modeling.
##   decimal:  use temperature, precipitation, longitude, latitude instead of month and prefecture.
feature_type <- 'category'
#feature_type <- 'decimal'



## main process
d <- read_tsv(cv_sum_fpath)
d$id <- apply(d[, c(1, 2, 4)], 1, paste, sep = '', collapse = '____')

d_null <- d %>% filter(feature_type == UQ(feature_type) & data_type == 'random')
d_altr <- d %>% filter(feature_type == UQ(feature_type) & data_type == 'shuffle')
#### number of rows of d_null and d_altr are not equal,
#### since some processs were aborted during modeling of randomized datasets or shuffled datasets.

## merge d_null and d_altr for calculating z-score simply
df <- full_join(d_null, d_altr, by = 'id', suffix = c('_null', '_altr'))
df <- df[!is.na(df$model_altr), ]       # removed data that is only appeared in shuffled dataset
df <- df[!is.na(df$model_null), ]  # removed data that is only appeared in randomzed dataset

## calculate the difference of mean RMSE between shuffled and randomized datasets and z-score
df$mean_diff <- df$cv_mean_altr - df$cv_mean_null
df$z         <- (df$cv_mean_altr - df$cv_mean_null) / df$cv_var_null


#df <- df[df$disease_null != 'cucumber__tsuruwarebyo__hompohasseimenseki', ]


fig1 <- df %>%
    ggplot(aes(y = disease_null, x = z, color = model_null)) +
        geom_point() +
        scale_color_locuszoom() +
        ylab('disease') + xlab('Z-score') + xlim(-5, 5) +
        facet_wrap(~ crop_null, ncol = 2)

fig2 <- df %>%
    ggplot(aes(y = disease_null, x = z, fill = model_null)) +
        geom_bar(stat = 'identity') +
        facet_grid(crop_null ~ model_null, scales = 'free_x') +
        scale_fill_locuszoom() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ylab('disease') + xlab('Z-score') + xlim(-5, 5)



png(paste0('cv_results/', feature_type, '_fig1.png'), 2000, 2000, res = 220)
print(fig1)
dev.off()
png(paste0('cv_results/', feature_type, '_fig2.png'), 2000, 3000, res = 220)
print(fig2)
dev.off()



