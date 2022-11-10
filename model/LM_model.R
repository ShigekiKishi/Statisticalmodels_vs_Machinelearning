### R script for multiple linear regression model

#options(repos = "http://mirror.maffin.ad.jp/CRAN/")
library(tidyverse)
library(magrittr)
library(Metrics)
library(rsample)

factor_name <- c("Month", "longi", "lati", "airtemp", "precip")


##/////////////////////////////////////////////////////////
# 2-fold cross validation to choice explanatory variables -----------------

calc_2cv <- function(train_d, valid_d){ 


  spl_2dats <- train_d %>% rsample::vfold_cv(v=2)
  train_d1 <- analysis(spl_2dats$splits[[1]])
  train_d2 <- analysis(spl_2dats$splits[[2]])

  if(nrow(train_d1) != nrow(train_d2)){
    if(nrow(train_d1) > nrow(train_d2)){
      train_d1 <- train_d1[1:nrow(train_d2),]
    }else{
      train_d2 <- train_d2[1:nrow(train_d1),]
    }
  }

  rmse_table <- NULL

  res <- glm(newY~1, data=train_d1)
  pred_res <- predict(res, newdata=train_d2)
  rmse2 <- rmse(train_d2$newY, pred_res)
  
  res <- glm(newY~1, data=train_d2)
  pred_res <- predict(res, newdata=train_d1)  
  rmse1 <- rmse(train_d1$newY, pred_res)
  rmse_m <- mean(c(rmse1, rmse2))
  term_factor <- c(0,0,0,0,0)
  rmse_table <- rbind(rmse_table, c(term_factor, rmse_m))  
  
  for(i in 1:5){ 
    factor_mat <- combn(factor_name, i)
  
    for(j in 1:ncol(factor_mat)){
      glm_factors <- factor_mat[,j]

      if(sum(train_d1$newY)==0){ 
        rmse2<- NA
      }
      else{
        factor_form <- paste(glm_factors, collapse="+")
        eval(parse(text = paste0("res <- glm(newY~", factor_form, ", family='gaussian', data=train_d1)")))
        pred_res <- predict(res, newdata= train_d2) 
        rmse2 <- rmse(train_d2$newY, pred_res)
      }
      
      if(sum(train_d2$newY)==0){  
        rmse1 <- NA
      }
      else{
        factor_form <- paste(glm_factors, collapse="+")
        eval(parse(text = paste0("res <- glm(newY~", factor_form, ", family='gaussian', data=train_d2)")))
        pred_res <- predict(res, newdata= train_d1) 
        rmse1 <- rmse(train_d1$newY, pred_res)
      }
      
      # rmse mean
      rmse_m <- mean(c(rmse1, rmse2))
      term_factor <- ifelse(match(factor_name, glm_factors) %>% is.na(), 0, 1)
      rmse_table <- rbind(rmse_table, c(term_factor, rmse_m))
    
    }
  }
  
  #na.omit()
  min_row <- rmse_table %>% as.data.frame() %>% set_colnames(c(factor_name, "RMSE")) %>%
    na.omit() %>% slice(which.min(RMSE))
  # all NA?
  if(is.null(min_row[6])){
    term_factor <- c(0,0,0,0,0)
    rmse3 <- NA
  } else {

    used_factors <- min_row[1:5]
    valid_dat <- valid_d
    glm_factors <- factor_name[ifelse(used_factors==1, TRUE, FALSE)]
    
    
    if(length(glm_factors) == 0){
      res <- glm(newY~1, data=train_d)
      pred_res <- predict(res, newdata = valid_dat)
      rmse3 <- rmse(valid_dat$newY, pred_res)
      term_factor <- c(0,0,0,0,0)
    }else{ 
      factor_form <- paste(glm_factors, collapse="+")
      eval(parse(text = paste0("res <- glm(newY~", factor_form, ", family='gaussian', data=train_d)")))
      
      pred_res <- predict(res, newdata = valid_dat)
      rmse3 <- rmse(valid_dat$newY, pred_res)
      term_factor <- ifelse(match(factor_name, glm_factors) %>% is.na(), 0, 1)
      
    }
  }

  return(c(term_factor,rmse3))

} ## end of 2-fold cross validation function

##/////////////////////////////////////////////////////////
# z_score calculation by 10-fold cross validation -----------------


calc_10cv <- function(x, y){ 
  rmse_normal <- mapply(calc_2cv, x, y) 
  train_nlist <- lapply(x, function(a){
    newx <- a
    newx[,c(factor_name, "newY")] %<>% apply(2, sample)
    newx %<>% as.data.frame()
    return(newx)
  })
  valid_nlist <- lapply(y, function(a){
    newx <- a
    if(nrow(newx)>1){
      newx[,c(factor_name, "newY")] %<>% apply(2, sample)
      newx %<>% as.data.frame()
    } 
    return(newx)
  })
  rmse_null <-  mapply(calc_2cv, train_nlist, valid_nlist)
  z_score <- (mean(rmse_normal[6,])-mean(rmse_null[6,]))/sd(rmse_null[6,])
  
  return(c(z_score, rmse_normal[6,], rmse_null[6,]))
} 


##/////////////////////////////////////////////////////////
# applying 10-fold cv function for each file! -----------------------------

files <- list.files(path="./dat_files/")

#z_table <- NULL
set.seed(300)
rand_mat <- NULL
for(j in 1:10){
  z_table <- NULL
  rand_vec <- as.integer(runif(length(files), min=1, max=2000))
  rand_mat <- rbind(rand_mat, rand_vec)
}
rand_mat %>% head()
rownames(rand_mat) <- NULL
rand_mat %>% dim() 

for(j in 1:10){ 
  z_table <- NULL
  
  for (i in 1:length(files)){
    file_name <- files[i]
    dat <- read.csv(paste0("dat_files/", file_name), header=T, row.names = 1) %>% 
      as.data.frame() %>% na.omit()
    dat <- dat %>% mutate(id = 1:nrow(dat))
    pref_num <- dat$Pref %>% unique() %>% length()
  
    if(pref_num > 1 && nrow(dat) > 15){ 
      if (str_detect(file_name, "menseki")){
        dat <- dat %>%  mutate(newY = incidence/Area) 
      }else{
        dat <- dat %>%  mutate(newY = incidence)
      }
  
    set.seed(rand_mat[j, i])
    dat_normal <- dat %>% sample_frac(size=1)
    spl_dat <- dat_normal %>% rsample::vfold_cv(v=10)

    train_list <- apply(spl_dat, 1, function(x)analysis(x$splits))
    valid_list <- apply(spl_dat, 1, function(x)rsample::assessment(x$splits))

    z_rmse <- calc_10cv(train_list, valid_list)
    z_table <- rbind(z_table, c(file_name, z_rmse))
  
    }else{ 
      z_table <- rbind(z_table, c(file_name, rep(NA, 21)))  
    }
  
  }
  
  write.csv(z_table, file=paste0("est_files/z_table", j, ".csv"))
}


