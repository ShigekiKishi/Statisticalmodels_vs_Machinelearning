library(tidyverse)
library(ggExtra)

summarise_category <- function(dpath) {
    cate <- list.dirs(dpath, full.names = FALSE, recursive = FALSE)
    cate_design <- data.frame(
        dpath = paste0(dpath, '/', cate),
        crop = str_split(cate, pattern = '__', simplify = TRUE)[, 1],
        disease = str_split(cate, pattern = '__', simplify = TRUE)[, 2],
        datatype = str_split(cate, pattern = '__', simplify = TRUE)[, 3]
    )
    
    cate_design
}


summarise_permutation_test <- function(dpath) {
    
    .get_best <- function(fpath) {
        x <- read.table(fpath, header = TRUE)
        mean(x$valid)
    }
    
    .sum_p <- function(p) {
        model_archs <- c('L1', 'L2')
        mse <- matrix(NA, nrow = 100, ncol = length(model_archs))
        colnames(mse) <- model_archs
    
        for (fpath in sort(list.files(p, pattern = 'tsv', full.names = TRUE))) {
            model_arch <- str_split(basename(fpath), pattern = '_', simplify = TRUE)[,1]
            i <- as.integer(gsub('.pth', '', str_split(basename(fpath), pattern = '_', simplify = TRUE)[,2]))
            mse[i, model_arch] <- .get_best(fpath)
        }
        
        mse
    }
    
    mse_p1 <- .sum_p(paste0(dpath, '/weights_permutest_1'))
    mse_p2 <- .sum_p(paste0(dpath, '/weights_permutest_2'))
    
    mse_L1 <- data.frame(row_shuffle = mse_p1[, 'L1'], random = mse_p2[, 'L1'])
    mse_L2 <- data.frame(row_shuffle = mse_p1[, 'L2'], random = mse_p2[, 'L2'])
    
    list(L1 = mse_L1, L2 = mse_L2)
}

summarise_valid <- function(dpath) {
    .get_best <- function(fpath) {
        x <- read.table(fpath, header = TRUE)
        mean(x$valid)
    }
    
    mse_L1 <- .get_best(paste0(dpath, '/weights/L1_best.pth_stats.tsv'))
    mse_L2 <- .get_best(paste0(dpath, '/weights/L2_best.pth_stats.tsv'))
    
    list(L1 = mse_L1, L2 = mse_L2)
}




if (TRUE) {
    
cate_design <- summarise_category('formatted_data_dw')
cate_design$model <- NA
cate_design$rmse <- NA
cate_design$pvalue <- NA

for (cate in 1:nrow(cate_design)) {
    
    # whole dataset
    x <- read.table(paste0(cate_design$dpath[cate], '/data.tsv'), header = TRUE)
    x <- x[!is.na(x$incidence), ]
    x$date <- as.Date(paste(x$year, x$month, '01', sep = '-'))
    hist_0 <- ggplot(x, aes(x = incidence)) +
                geom_histogram(bins = 20) +
                ggtitle(paste0('n = ', nrow(x)))
    png(paste0(cate_design$dpath[cate], '/data.hist.png'), 800, 500, res = 250)
    print(hist_0)
    dev.off()
    scat_0 <- ggplot(x, aes(x = date, y = incidence)) +
                geom_point() +
                ggtitle(paste0(paste(cate_design[cate, 2:4], collapse = ' '), '\n(n = ', nrow(x), ')')) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

    png(paste0(cate_design$dpath[cate], '/data.scatter.png'), round(100 * length(unique(x$year))), 700, res = 250)
    print(scat_0)
    dev.off()


    # permutation results
    mse <- try(summarise_permutation_test(cate_design$dpath[cate]), silent = TRUE)
    if (class(mse) != 'try-error') {
        mse_df <- rbind(
            pivot_longer(mse$L1, everything(), names_to = 'test', values_to = 'mse') %>% mutate(model = 'L1'),
            pivot_longer(mse$L2, everything(), names_to = 'test', values_to = 'mse') %>% mutate(model = 'L2')
        )
        mse_df$rmse <- sqrt(mse_df$mse)
    } else {
        mse <- NULL
        mse_df <- NULL
    }
    
    # best result
    mse_best <- try(summarise_valid(cate_design$dpath[cate]), silent = TRUE)
    if (class(mse_best) != 'try-error') {
        rmse_best <- lapply(mse_best, sqrt)
    } else {
        rmse_best <- NULL
    }
    
    if (!is.null(mse) && !is.null(rmse_best)) {
        hist_L1 <- mse_df %>% filter(model == 'L1') %>%
                ggplot(aes(x = rmse)) +
                geom_histogram(bins = 20) +
                geom_vline(xintercept = rmse_best$L1, linetype = 'dashed', color = '#C71000') +
                facet_grid(test ~ .) +
                ggtitle(paste(cate_design[cate, 2:5], collapse = ' '))
    
        hist_L2 <- mse_df %>% filter(model == 'L2') %>%
                ggplot(aes(x = rmse)) +
                geom_histogram(bins = 20) +
                geom_vline(xintercept = rmse_best$L2, linetype = 'dashed', color = '#C71000') +
                facet_grid(test ~ .) +
                ggtitle(paste(cate_design[cate, 2:5], collapse = ' '))
    
        png(paste0(cate_design$dpath[cate], '/permutation_test.L1.png'), 800, 1000, res = 250)
        print(hist_L1)
        dev.off()
    
        png(paste0(cate_design$dpath[cate], '/permutation_test.L2.png'), 800, 1000, res = 250)
        print(hist_L2)
        dev.off()
    
        if (rmse_best$L1 > rmse_best$L2) {
            cate_design$model[cate] <- 'L2'
            cate_design$rmse[cate] <- rmse_best$L2
            cate_design$pvalue[cate] <- t.test(sqrt(mse$L2$row_shuffle), sqrt(mse$L2$random))$p.value
        } else {
            cate_design$model[cate] <- 'L1'
            cate_design$rmse[cate] <- rmse_best$L1
            cate_design$pvalue[cate] <- t.test(sqrt(mse$L1$row_shuffle), sqrt(mse$L1$random))$p.value
        }
    }
    
    
}


write.table(cate_design, file = 'formatted_data_dw/summary.xls',
            sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


}




