### シンプル重回帰するRスクリプト

library(tidyverse)
library(magrittr)
library(MASS)
library(Metrics)
#library(furrr)
#plan(multiprocess) # for スパコン用
library(rsample)
# install.packages("glmnet)
library(glmnet) #repos = "https://cran.ism.ac.jp/")
## 普通はplan(multisession) or plan(multicore)

# データの読み込み
# Data
# process_files.Rから渡されるのはfile_name ("~~.csv")

## glm

factor_name <- c("Month", "longi", "lati", "airtemp", "precip")


##/////////////////////////////////////////////////////////
# 2-fold cross validation to choice explanatory variables -----------------
##/////////////////////////////////////////////////////////

calc_2cv <- function(train_d){ ## 関数化

  # 2分割
  spl_2dats <- train_d %>% rsample::vfold_cv(v=2)
  train_d1 <- analysis(spl_2dats$splits[[1]])
  train_d2 <- analysis(spl_2dats$splits[[2]])

  #train_dat %>% head()

  # train_d1とtrain_d2の行数をそろえる
  if(nrow(train_d1) != nrow(train_d2)){
    if(nrow(train_d1) > nrow(train_d2)){
      train_d1 <- train_d1[1:nrow(train_d2),]
    }else{
      train_d2 <- train_d2[1:nrow(train_d1),]
    }
  }

  rmse_table <- NULL
  ## no factors
  res <- glm(newY~1, data=train_d1)
  pred_res <- predict(res)
  rmse1 <- rmse(train_d1$newY, pred_res)
  res <- glm(newY~1, data=train_d2)
  pred_res <- predict(res)  
  rmse2 <- rmse(train_d2$newY, pred_res)
  rmse_m <- mean(c(rmse1, rmse2))
  # factor 0 or 1
  term_factor <- c(0,0,0,0,0)
  # 積み重ねる
  rmse_table <- rbind(rmse_table, c(term_factor, rmse_m)) 
  
  ## more than 0 factors
  for(i in 1:5){ ## glmnetは説明変数が2以上
    factor_mat <- combn(factor_name, i)
    
    for(j in 1:ncol(factor_mat)){
      glm_factors <- factor_mat[,j]
      # train_d1
      #factor_df <- train_d1 %>% dplyr::select(all_of(glm_factors))
      
      if(sum(train_d1$newY)==0){  ## when all records = 0 !
        rmse2<- NA
      }
      else{
        factor_form <- paste(glm_factors, collapse="+")
        eval(parse(text = paste0("res <- glm(newY~", factor_form, ", family='gaussian', data=train_d1)")))
        pred_res <- predict(res, newdata = train_d2) #%>% dplyr::select(all_of(glm_factors)) %>% as.matrix())
        rmse2 <- rmse(train_d2$newY, pred_res)
      }
      
      # train_d2
      #factor_df <- train_d2 %>% dplyr::select(all_of(glm_factors))
      
      if(sum(train_d2$newY)==0){  ## when all records = 0 !
        rmse1 <- NA
      }
      else{
        factor_form <- paste(glm_factors, collapse="+")
        eval(parse(text = paste0("res <- glm(newY~", factor_form, ", family='gaussian', data=train_d2)")))
        pred_res <- predict(res, newdata = train_d1)# newx = train_d1 %>% dplyr::select(all_of(glm_factors)) %>% as.matrix())
        rmse1 <- rmse(train_d1$newY, pred_res)
      }
      
      # rmse mean
      rmse_m <- mean(c(rmse1, rmse2))
      # factor 0 or 1
      term_factor <- ifelse(match(factor_name, glm_factors) %>% is.na(), 0, 1)
      # 積み重ねる
      rmse_table <- rbind(rmse_table, c(term_factor, rmse_m))
      
    }
  }

  min_row <- rmse_table %>% as.data.frame() %>% set_colnames(c(factor_name, "RMSE")) %>% slice(which.min(RMSE))

  return(min_row)

} # end of calc_2cv() function

#


##OK
#ここで選ばれた変数をつかって検証データでRMSEを求める。
RMSE_calc <- function(train_d, valid_d){
  # choice parameters
  params <- calc_2cv(train_d)
  use_factors <- params[1:5]

  valid_dat <- valid_d
  glm_factors <- factor_name[ifelse(use_factors==1, TRUE, FALSE)] # 説明変数のベクター

  if(length(glm_factors) == 0){ 
    res <- glm(newY~1, data=train_d)
    pred_res <- predict(res, newdata = valid_dat)
    rmse3 <- rmse(valid_dat$newY, pred_res)
    # factor 0 or 1
    term_factor <- c(0,0,0,0,0)
  }else{ # factors > 0
    #glm_factors
    factor_form <- paste(glm_factors, collapse="+")
    eval(parse(text = paste0("res <- glm(newY~", factor_form, ", family='gaussian', data=train_d)")))
    
    pred_res <- predict(res, newdata = valid_dat)
    rmse3 <- rmse(valid_dat$newY, pred_res)
    term_factor <- ifelse(match(factor_name, glm_factors) %>% is.na(), 0, 1)

  }


  return(c(term_factor,rmse3))

} ## end of 2-fold cross validation function

##/////////////////////////////////////////////////////////
# z_score calculation by 10-fold cross validation -----------------
##/////////////////////////////////////////////////////////

calc_10cv <- function(x, y){ # x:train_list, y: valid_list
  ## calculating 10 normal RMSEs
  rmse_normal <- mapply(RMSE_calc, x, y) #return matrix
  
  ## calculating 10 randomized RMSEs
  ## making null model (train_list and valid_list)
  train_nlist <- lapply(x, function(a){
    newx <- a
    newx[,c(factor_name, "newY")] %<>% apply(2, sample)
    newx %<>% as.data.frame()
    return(newx)
  })
  valid_nlist <- lapply(y, function(a){
    newx <- a
    newx[,c(factor_name, "newY")] %<>% apply(2, sample)
    newx %<>% as.data.frame()
    return(newx)
  })
  ## strawberry(-Season)にも対応できるよう
  rmse_null <-  mapply(RMSE_calc, train_nlist, valid_nlist)
  
  # calculating z score
  
  z_score <- (mean(rmse_normal[6,])-mean(rmse_null[6,]))/sd(rmse_null[6,])
  
  return(c(z_score, rmse_normal[6,], rmse_null[6,]))
} ## end of 10-fold cross validation


##/////////////////////////////////////////////////////////
# applying 10-fold cv function for each file! -----------------------------
##/////////////////////////////////////////////////////////

files <- list.files(path="./dat_files/")
#length(files)
z_table <- NULL
for (i in 1:length(files)){
  file_name <- files[i]
  # 読み込み
  dat <- read.csv(paste0("dat_files/", file_name), header=T, row.names = 1) %>% 
    as.data.frame() %>% na.omit()
  dat <- dat %>% mutate(id = 1:nrow(dat))
  ## pref number
  pref_num <- dat$Pref %>% unique() %>% length()
  
  if(pref_num > 1 && nrow(dat) > 10){ ## more than 1 pref and more than 15 records
    
    # 面積データなら栽培面積の割合newYに変える。面積データでないならそのまま与える
    if (str_detect(file_name, "menseki")){
      dat <- dat %>%  mutate(newY = incidence/Area) # ratio
    }else{
      dat <- dat %>%  mutate(newY = incidence)
    }
    
    set.seed(100)
    # まずはシャッフル
    dat_normal <- dat %>% sample_frac(size=1)
    # 次に10分割
    spl_dat <- dat_normal %>% rsample::vfold_cv(v=10)
    #analysis(spl_dat$splits[[1]]) # 学習データを表示する
    #assessment(spl_dat$splits[[1]]) # validationデータを表示する
    
    train_list <- apply(spl_dat, 1, function(x)analysis(x$splits))
    #train_list[[1]] %>% dim() # 657 rows
    
    valid_list <- apply(spl_dat, 1, function(x)rsample::assessment(x$splits))
    
    z_rmse <- calc_10cv(train_list, valid_list)
    z_table <- rbind(z_table, c(file_name, z_rmse))
    
  }else{ # pref num ==1, or nrow(dat) =< 15
    z_table <- rbind(z_table, c(file_name, rep(NA, 21)))  
  }
  
  
}## end for 203 files

## 保存
write.csv(z_table, file="est_files/z_table.csv")

