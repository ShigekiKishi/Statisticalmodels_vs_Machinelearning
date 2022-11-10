library(tidyverse)
library(ggsci)


GENERATE_DATASET   <- FALSE
NORM_SPLIT_DATASET <- TRUE






# hash for changing japanese prefecture name into english
pref.name.ja <- c('北海道', '青森県', '岩手県', '宮城県', '秋田県', '山形県',
                  '福島県', '茨城県', '栃木県', '群馬県', '埼玉県', '千葉県',
                  '東京都', '神奈川県', '新潟県', '富山県', '石川県', '福井県',
                  '山梨県', '長野県', '岐阜県', '静岡県', '愛知県', '三重県',
                  '滋賀県', '京都府', '大阪府', '兵庫県', '奈良県', '和歌山県',
                  '鳥取県', '島根県', '岡山県', '広島県', '山口県', '徳島県',
                  '香川県', '愛媛県', '高知県', '福岡県', '佐賀県', '長崎県',
                  '熊本県', '大分県', '宮崎県', '鹿児島県', '沖縄県')
pref.name.en <- c('hokkaido', 'aomori', 'iwate', 'miyagi', 'akita', 'yamagata',
                  'fukushima', 'ibaraki', 'tochigi', 'gunma', 'saitama', 'chiba',
                  'tokyo', 'kanagawa', 'niigata', 'toyama', 'ishikawa', 'fukui',
                  'yamanashi', 'nagano', 'gifu', 'shizuoka', 'aichi', 'mie',
                  'shiga', 'kyoto', 'osaka', 'hyogo', 'nara', 'wakayama',
                  'tottori', 'shimane', 'okayama', 'hiroshima', 'yamaguchi', 'tokushima',
                  'kagawa', 'ehime', 'kochi', 'fukuoka', 'saga', 'nagasaki',
                  'kumamoto', 'oita', 'miyazaki', 'kagoshima', 'okinawa')
ja2en <- pref.name.en
names(ja2en) <- pref.name.ja








read_jppnet_csv <- function(dpath = './dat/現況報告') {
    dat <- NULL
    for (fpath in list.files(dpath, full.names = TRUE)) {
        if (is.null(dat)) {
            dat <- read_csv(fpath, local = locale(encoding = 'CP932'))
        } else {
            dat <- bind_rows(dat, read_csv(fpath, local = locale(encoding = 'CP932')))
        }
    }
    
    dat <- dat %>% dplyr::rename(date = 調査年月日, prefecture = 都道府県名称,
                                 crop = 作物名称, disease = 病害虫名称, method = 調査名,
                                 level = 程度, value = データ, unit = 単位) %>%
                   mutate(date = as.Date(date),
                          prefecture = as.character(ja2en[prefecture]),
                          level = ifelse(is.na(level), NA, recode(level,
                                        `高`=5, `多`=5, `やや高`=4, `やや多`=4, `並`=3,
                                        `やや少` = 2, `やや低`=2, `少`=1, `低`=1, .default =0)),
                          value = as.numeric(value),
                          year = str_split(date, pattern = '-', simplify = T)[, 1],
                          month = str_split(date, pattern = '-', simplify = T)[, 2],
                          day = str_split(date, pattern = '-', simplify = T)[, 3]) %>%
                   select(date, year, month, day, prefecture, crop, disease, method, level, value, unit)
    
    dat
}


read_weather_data <- function(wcode = 'temp') {
    
    if (wcode == 'temp') {
        fpath <- './dat/airtemp200127.csv'
    } else {
        fpath <- './dat/precip200127.csv'
    }
    
    dat <- read_csv(fpath) %>%
                pivot_longer(-c(year, month), names_to = 'prefecture', values_to = 'value')
    years <- sort(unique(dat$year))
        
    dlist <- vector('list', length(years))
    names(dlist) <- as.character(years)
        
    for (y in names(dlist)) {
        dfy <- as.data.frame(dat %>% filter(year == y) %>%
                    select(month, prefecture, value) %>%
                    pivot_wider(names_from = month, values_from = value))
        df <- dfy[, -1]
        rownames(df) <- dfy[, 1]
        dlist[[y]] <- as.matrix(df[pref.name.en, ])
    }
      
    dlist
}



make_incidence_matrix <- function(dat, replace.na = NA) {
    year <- sort(unique(dat$year))
    
    incidence_matrices <- vector('list', length = length(year))
    names(incidence_matrices) <- year
    
    for (y in year) {
        mat <- matrix(replace.na, ncol = 12, nrow = length(ja2en))
        rownames(mat) <- ja2en
        colnames(mat) <- 1:12
        
        dat_y <- dat %>% filter(year == y) %>%
                    select(month, prefecture, value) %>%
                    group_by(month, prefecture) %>%
                    summarize(value = mean(value)) %>%
                    ungroup() %>%
                    spread(key = month, value = value)
        for (m in colnames(dat_y)) {
            if (m != 'prefecture') {
                mat[dat_y$prefecture, as.integer(m)] <- dat_y[[m]]
            }
        }
        mat[is.na(mat)] <- replace.na
        incidence_matrices[[y]] <- mat[pref.name.en, ]
    }
    
    incidence_matrices
}



generate_dataset_from_mesh <- function(incidence_matrices = NULL) {
    
    # make emtpy matrix for generating test dataset
    if (is.null(incidence_matrices)) {
        mat <- matrix(-1, ncol = 12, nrow = length(ja2en))
        rownames(mat) <- ja2en
        colnames(mat) <- 1:12
        incidence_matrices <- list('2022' = mat)
    }
    
    # prefecture to gis
    pre2gis <- read_tsv('./dat/pre2gis.tsv')
    pre2lng <- pre2gis$lng
    pre2lat <- pre2gis$lat
    names(pre2lng) <- names(pre2lat) <- pre2gis$prefecture
    
    # prefecture to mesh code
    mesh2env <- read_csv('./dat/mesh2Env.csv')
    pre2mesh <- c(NA)
    for (i in 1:nrow(pre2gis)) {
        dd <- sqrt((mesh2env$lat - pre2gis$lat[i])^2 + (mesh2env$lon - pre2gis$lng[i])^2)
        pre2mesh <- c(pre2mesh, mesh2env$mesh2[dd == min(dd)][1])
    }
    pre2mesh <- pre2mesh[-1]
    names(pre2mesh) <- pre2gis$prefecture
    
    # mesh code to temperature and rain matrix
    rain <- temp <- matrix(NA, nrow = length(pre2mesh), ncol = 12)
    rownames(rain) <- rownames(temp) <- names(pre2mesh)
    for (i in 1:12) {
        for (m in 1:length(pre2mesh)) {
            temp[m, i] <- as.numeric(mesh2env[mesh2env$mesh2 == pre2mesh[m], i + 3 + 12])
            rain[m, i] <- as.numeric(mesh2env[mesh2env$mesh2 == pre2mesh[m], i + 3])
        }
    }
    
    
    # dataset generation
    mat2env <- NULL
    mat <- as.data.frame(incidence_matrices[[1]])
    mat$prefecture <- rownames(incidence_matrices[[1]])
    mat <- mat %>% gather(key = month, value = value, - prefecture) %>%
                    mutate(month = as.integer(month))
    for (i in 1:nrow(mat)) {
        df <- data.frame(year = '2020', month = mat$month[i],
                   prefecture = mat$prefecture[i],
                   temperature = temp[mat$prefecture[i], mat$month[i]],
                   rain = rain[mat$prefecture[i], mat$month[i]],
                   lat = pre2lat[mat$prefecture[i]],
                   lng = pre2lng[mat$prefecture[i]])
        if (is.null(mat2env)) {
            mat2env <- df
        } else {
            mat2env <- rbind(mat2env, df)
        }
    }
    rownames(mat2env) <- NULL
    mat2env <- mat2env[, -c(1, 2, 3)]


    dat <- NULL
    for (y in names(incidence_matrices)) {
        mat <- as.data.frame(incidence_matrices[[y]])
        mat$prefecture <- rownames(incidence_matrices[[y]])
        mat <- mat %>% gather(key = month, value = value, - prefecture) %>%
                    mutate(year = y, month = as.integer(month))
        mat <- cbind(mat, mat2env) %>%
                select(year, month, prefecture, temperature, rain, lat, lng, value)
        if (is.null(dat)) {
            dat <- mat
        } else {
            dat <- rbind(dat, mat)
        }
    }
    rownames(dat) <- NULL
    
    dat
}



generate_dataset <- function(incidence_matrices = NULL) {
    
    # make emtpy matrix for generating test dataset
    if (is.null(incidence_matrices)) {
        mat <- matrix(NA, ncol = 12, nrow = length(ja2en))
        rownames(mat) <- ja2en
        colnames(mat) <- 1:12
        incidence_matrices <- list('2022' = mat)
    }
    
    
    # weather data
    wdat_temp <- read_weather_data('temp')
    wdat_rain <- read_weather_data('rain')
    
    
    # GIS data
    pre2gis <- read_tsv('./dat/pre2gis.tsv')
    pre2lng <- pre2gis$lng
    pre2lat <- pre2gis$lat
    names(pre2lng) <- names(pre2lat) <- pre2gis$prefecture
    
    
    # data generation
    train_dat <- NULL
    .calc_matlist_mean <- function(dlist) {
        mat <- NULL
        for (i in seq(dlist)) {
            if (is.null(mat)) {
                mat <- dlist[[i]]
            } else {
                mat <- mat + dlist[[i]]
            }
        }
        mat <- mat / length(dlist)
        mat
    }
    
    for (y in names(incidence_matrices)) {
        if (y != '2022') {
            mat_  <- incidence_matrices[[y]]
            temp_ <- wdat_temp[[y]]
            rain_ <- wdat_rain[[y]]
            lat_  <- matrix(pre2lat[rownames(mat_)], ncol = ncol(mat_), nrow = nrow(mat_))
            lng_  <- matrix(pre2lng[rownames(mat_)], ncol = ncol(mat_), nrow = nrow(mat_))
        } else {
            mat_ <- incidence_matrices[[y]]
            temp_ <- .calc_matlist_mean(wdat_temp)
            rain_ <- .calc_matlist_mean(wdat_rain)
            lat_  <- matrix(pre2lat[rownames(mat_)], ncol = ncol(mat_), nrow = nrow(mat_))
            lng_  <- matrix(pre2lng[rownames(mat_)], ncol = ncol(mat_), nrow = nrow(mat_))
        }
        df_ <- data.frame(incidence = as.numeric(mat_), year = as.integer(y),
                          month = rep(as.integer(colnames(mat_)), each = nrow(mat_)),
                          prefecture = as.character(rownames(mat_)),
                          temp = as.numeric(temp_), rain = as.numeric(rain_),
                          lat = as.numeric(lat_), lng = as.numeric(lng_))
        train_dat <- rbind(train_dat, df_)
    }
    
    train_dat
}






calc_nf <- function(dat_records) {
    nf.temp <- list(mean = mean(dat_records$temp), sd = sd(dat_records$temp))
    nf.rain <- list(mean = mean(dat_records$rain), sd = sd(dat_records$rain))
    nf.lat  <- list(mean = mean(dat_records$lat), sd = sd(dat_records$lat))
    nf.lng  <- list(mean = mean(dat_records$lng), sd = sd(dat_records$lng))
    nf <- list(temp = nf.temp, rain = nf.rain, lat = nf.lat, lng = nf.lng)
    nf
}

normalize_records <- function(dat_records, nf) {
    dat_records$temp <- (dat_records$temp - nf$temp$mean) / nf$temp$sd
    dat_records$rain <- (dat_records$rain - nf$rain$mean) / nf$rain$sd
    dat_records$lat <- (dat_records$lat - nf$lat$mean) / nf$lat$sd
    dat_records$lng <- (dat_records$lng - nf$lng$mean) / nf$lng$sd
    if ('value' %in% colnames(dat_records)) {
        dat_records$value <- dat_records$value
    }
    dat_records
}




arrange_dataset <- function(data_dpath, cutoff = 2013) {
    
    # load formatted data
    dat_records <- read_tsv(paste0(data_dpath, '/data.tsv'))
    
    # split the original data to trainning set and validation set, and normalize both sets
    dat_records_train <- dat_records %>% filter(year <= cutoff)
    dat_records_valid <- dat_records %>% filter(year > cutoff)
    nf_train <- calc_nf(dat_records_train)
    dat_records_train_std <- normalize_records(dat_records_train, nf_train)
    dat_records_valid_std <- normalize_records(dat_records_valid, nf_train)
    
    # generate whole dataset for fixing model
    dat_records_test <- generate_dataset()
    nf_all <- calc_nf(dat_records)
    dat_records_std <- normalize_records(dat_records, nf_all)
    dat_records_test_std <- normalize_records(dat_records_test, nf_all)
    
    # remove NA
    dat_records_std <- dat_records_std[!is.na(dat_records_std$incidence), ]
    dat_records_train_std <- dat_records_train_std[!is.na(dat_records_train_std$incidence), ]
    dat_records_valid_std <- dat_records_valid_std[!is.na(dat_records_valid_std$incidence), ]
    
    # write records
    # file name must be data_std.tsv, train_std.tsv, valid_std.tsv, and test_std.tsv
    write_delim(dat_records_std,       path = paste0(data_dpath, '/data_std.tsv'), delim = '\t', na = 'NA')
    write_delim(dat_records_test_std,  path = paste0(data_dpath, '/test_std.tsv'), delim = '\t', na = 'NA')
    if (nrow(dat_records_train_std) > 0 && nrow(dat_records_valid_std) > 0) {
        write_delim(dat_records_train_std, path = paste0(data_dpath, '/train_std.tsv'), delim = '\t', na = 'NA')
        write_delim(dat_records_valid_std, path = paste0(data_dpath, '/valid_std.tsv'), delim = '\t', na = 'NA')
    }
    
    
    # generate train (train) data and train (valid) for cross-validation
    dir.create(paste0(data_dpath, '/cv'), showWarnings = FALSE)
    set.seed(2020)
    n_cv <- 10
    if (nrow(dat_records_train_std) > n_cv) {
        n_ <- floor(nrow(dat_records_train_std) / n_cv)
        k_ <- nrow(dat_records_train_std) %% n_cv
        if (k_ == 0) {
            idx_ <- sample(c(rep(1:n_cv, each = n_)))
        } else {
            idx_ <- sample(c(rep(1:n_cv, each = n_), 1:k_))
        }
        for (i in 1:n_cv) {
            dat_cv_train <- dat_records_train_std[idx_ != i, ]
            dat_cv_valid <- dat_records_train_std[idx_ == i, ]
            write_delim(dat_cv_train, path = paste0(data_dpath, '/cv/train_std_', i, '.tsv'), delim = '\t', na = 'NA')
            write_delim(dat_cv_valid, path = paste0(data_dpath, '/cv/valid_std_', i, '.tsv'), delim = '\t', na = 'NA')
        }
    }
    
}






randomize_dataset <- function(dat_records, dpath, method = 'norm') {
    # shuffle datasets (rows) and generate 100 datasets.
    if (nrow(dat_records) <= 10) {
        return(FALSE)
    }
    
    dir.create(dpath, showWarnings = FALSE)
    
    for (i in 1:100) {
        dpath_i <- paste0(dpath, '/', sprintf('%03d', i))
        dir.create(dpath_i, showWarnings = FALSE)
        
        set.seed(2020 + i)
        dat_records_rand <- dat_records[sample(1:nrow(dat_records)), ]
        
        if (method == 'random') {
            for (di in 1:ncol(dat_records_rand)) {
                dat_records_rand[[di]] <- sample(dat_records_rand[[di]])
            }
        }
        
        train_raw <- dat_records_rand[1:round(nrow(dat_records_rand) * 0.8), ]
        valid_raw <- dat_records_rand[(round(nrow(dat_records_rand) * 0.8) + 1):nrow(dat_records_rand), ]
        
        nf_train <- calc_nf(train_raw)
        train_std <- normalize_records(train_raw, nf_train)
        valid_std <- normalize_records(valid_raw, nf_train)
        
        write_delim(train_std, path = paste0(dpath_i, '/train_std.tsv'),
                    delim = '\t', na = 'NA')
        write_delim(valid_std, path = paste0(dpath_i, '/valid_std.tsv'),
                    delim = '\t', na = 'NA')
        
        # generate cv
        n_cv <- 10
        
        dpath_cv <- paste0(dpath_i, '/cv')
        dir.create(dpath_cv, showWarnings = FALSE)
        
        if (nrow(train_std) > n_cv) {
            n_ <- floor(nrow(train_std) / n_cv)
            k_ <- nrow(train_std) %% n_cv
            if (k_ == 0) {
                idx_ <- sample(c(rep(1:n_cv, each = n_)))
            } else {
                idx_ <- sample(c(rep(1:n_cv, each = n_), 1:k_))
            }
            for (j in 1:n_cv) {
                train_std_cv <- train_std[idx_ != j, ]
                valid_std_cv <- train_std[idx_ == j, ]
                write_delim(train_std_cv, path = paste0(dpath_cv, '/train_std_', j, '.tsv'), delim = '\t', na = 'NA')
                write_delim(valid_std_cv, path = paste0(dpath_cv, '/valid_std_', j, '.tsv'), delim = '\t', na = 'NA')
            }
        }
        
    }

}







if (GENERATE_DATASET) {
    # exmaple code for generating 'Cucumber__udonkobyo__hompohatsubyoyoritsu/data.tsv'
    jppnet_dat <- read_jppnet_csv()
    dat <- jppnet_dat %>%
                filter(str_detect(crop, 'キュウリ')) %>%
                filter(disease == 'うどんこ病') %>%
                filter(method == '本圃発病葉率') %>%
                filter(!is.na(value))
    incidence_matrices <- make_incidence_matrix(dat, replace.na = NA)
    dat_records <- generate_dataset(incidence_matrices)
    
    # folder name should be `crop__disease__datatype`, file name should be `data.tsv`
    write_delim(dat_records, path = 'formatted_data/Cucumber__udonkobyo__hompohatsubyoyoritsu/data.tsv')
}
 



if (NORM_SPLIT_DATASET) {
    project_data_dpath <- './formatted_data'
    
    for (data_dpath in list.dirs(project_data_dpath, recursive = FALSE)) {
        print(data_dpath)
        dat_records <- read_tsv(paste0(data_dpath, '/data.tsv')) %>% filter(!is.na(incidence))
        
        # data to generate observed data distribution
        norm_dpath <- paste0(data_dpath, '/norm_datasets')
     #   if (!file.exists(norm_dpath)) {
            randomize_dataset(dat_records, norm_dpath, 'norm')
     #   }
        
        # data to generate null distribution
        null_dpath <- paste0(data_dpath, '/null_datasets')
     #   if (!file.exists(null_dpath)) {
            randomize_dataset(dat_records, null_dpath, 'random')
     #   }
    }
}



