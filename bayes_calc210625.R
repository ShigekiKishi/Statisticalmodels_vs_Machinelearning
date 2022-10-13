# ベイズ推定を使うためのRスクリプト

# ファイルの読み込み
# Data
print(file_name) # このfile_nameは外から与える

library(tidyverse)
library(rstan)
library(magrittr)
library(bayesplot)
#install.packages("Metrics", repos = "https://cran.ism.ac.jp/")
library(Metrics)
library(furrr)
library(rsample)
#plan(multicore)
plan(multiprocess) # for スパコン用
## 普通はplan(multisession) or plan(multicore)

#### Stan settings ####
iterN <- 100 #Num iterations
warmupN <- 20 # Num warmup
chainsN <- 1 # Num chins
#repN <- 100 ## 100回計算

###　おまじない
rstan_options(auto_write =T)
options(mc.cores = parallel::detectCores()) 



##################　作付面積を考慮する　#################################################

### ファイルの読み込み

dat <- read.csv(paste0("dat_files/",file_name, ".csv"), header=T, row.names = 1) %>% as.data.frame() %>% na.omit()
# 気温、降水量
airtemp_mat <- read.csv("sihou/airtemp_mat.csv", header=T, row.names = 1) %>% as.matrix()
precip_mat <- read.csv("sihou/precip_mat.csv", header=T, row.names = 1) %>% as.matrix()
# 県名を小文字にする
colnames(airtemp_mat) <- colnames(airtemp_mat) %>% tolower() # 小文字にする
colnames(precip_mat) <- colnames(precip_mat) %>% tolower() # 小文字にする

#都道府県の位置情報
p_locat <- read.csv("sihou/p_locat.csv", header=T, row.names = 1)
# 県名を小文字にする
p_locat$Pref %<>% tolower()



#県名番号
dat %<>% left_join(p_locat[,1:2], by = "Pref") # これでPref_Nが増える

#ここから別操作
#dat %>% head()
#dat$Pref_N %>% unique() %>% sort()
p_locat2 <- p_locat %>% slice(dat$Pref_N %>% unique() %>% sort()) # 出現する都道府県のみ抽出
p_locat2 %<>% mutate(Pref_N2 = 1:length(p_locat2$Pref)) # この中だけで使う県の番号付与

dat %<>% left_join(p_locat2 %>% dplyr::select(Pref, Pref_N2), by = "Pref") # Pref_N2 追加

# 気温、降水量も今回用のものを用意する
airtemp_mat2 <- airtemp_mat %>% as.data.frame() %>% dplyr::select(p_locat2$Pref) %>% as.matrix()
precip_mat2 <- precip_mat %>% as.data.frame() %>% dplyr::select(p_locat2$Pref) %>% as.matrix()


#///////////////////////////////////////////////////////////////////////
# Stan!! ------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////

set.seed(100)

#stanモデルをコンパイルしておく
model3 <- rstan::stan_model(file = "sihou/model3.stan")


# データ作成

# データが面積であるかどうかによって目的変数を変える
menseki <- file_name %>% str_detect(pattern="menseki")
if(menseki){
  dat %<>% mutate(res_val = incidence/Area) #面積で割る
}else{
  dat %<>% mutate(res_val = incidence) # そのまま
}

# スタン用データ作成
# 10 fold cross validation
set.seed(100)
# まずはシャッフル
dat_normal <- dat %>% sample_frac(size=1)
# 次に分割
spl_dat <- dat_normal %>% rsample::vfold_cv(v=10)
#analysis(spl_dat$splits[[1]]) # 学習データを表示する
#assessment(spl_dat$splits[[1]]) # validationデータを表示する

train_list <- apply(spl_dat, 1, function(x)analysis(x$splits))
# train_list[[1]] %>% dim() # 657 rows

valid_list <- apply(spl_dat, 1, function(x)rsample::assessment(x$splits))
# valid_list[[1]] %>% dim()  # 73 rows


# データランダマイズ
data_list <- NULL 
trainN <- train_list[[1]] %>% nrow()

for (i in 1:10){ # 10分割したデータを使う
  dat_piece <- rbind(train_list[[i]], valid_list[[i]]) %>% as.data.frame()
  # データの準備
  stan_data <- list(
    I = length(dat_piece$res_val),
    I_train = trainN,
    P = max(dat_piece$Pref_N2), #県の数
    M = 12, #12months
    mu_zero = 0,
    Pref = as.integer(dat_piece$Pref_N2),
    Month = as.integer(dat_piece$Month),
    Temp = airtemp_mat2,
    Precip = precip_mat2,
    Longi = as.numeric(p_locat2$longi),
    Lati = as.numeric(p_locat2$lati),
    Incidence = dat_piece$res_val # res_val
  )
  data_list %<>%  append(list(stan_data))
}


fit_list <- future_map(data_list, rstan::sampling, object = model3, iter=iterN, warmup=warmupN, chains=chainsN, verbose=F)
# ようやくいった。。。


############## 結果の要約 #########################

##rhat
rhat_summary <- fit_list %>% sapply(rhat) %>% apply(2, summary)

##extract
exfit_list <- fit_list %>% lapply(rstan::extract)

# m_muの平均値
mmu_mean <- exfit_list %>% sapply(function(x)x$m_mu %>% apply(2, mean)) %>% apply(1, mean)
##この増減は平均というか全国共通

# para_temp
paratemp_mean <- exfit_list %>% sapply(function(x)x$para_temp %>% mean()) %>% mean()
# para_precip
paraprecip_mean <- exfit_list %>% sapply(function(x)x$para_precip %>% mean()) %>% mean()
# para_longi
paralongi_mean <- exfit_list %>% sapply(function(x)x$para_longi %>% mean()) %>% mean()
# para_lati
paralati_mean <- exfit_list %>% sapply(function(x)x$para_lati %>% mean()) %>% mean()
### lambda
lambda_mean <- exfit_list %>% sapply(function(x)x$lambda %>% apply(2, mean)) %>% apply(1, mean)

# 予備的なもの
params <- exfit_list %>% sapply(function(x)sapply(x, mean)) %>% apply(1, mean)



#################################  RMSE  #################################

#pred_incidence予測値
pred_incidence <- exfit_list %>% sapply(function(x)x %$% Incidence_new %>% apply(2, mean))
#ori_incidence観測値
ori_incidence <- data_list %>% sapply(function(x)x$Incidence) #データフレームにする

#RMSE テストデータ
ori_i <- ori_incidence[(trainN+1):nrow(ori_incidence),]
pred_i <- pred_incidence[(trainN+1):nrow(pred_incidence),]
ori_pred_list <- list(ori_i, pred_i) # リストにまとめる
rmse_res <- sapply(1:ncol(ori_i), function(x)rmse(ori_pred_list[[1]][,x], ori_pred_list[[2]][,x])) # Metrics::rmse(actual, predicted)

# データフレームにまとめる
ori_res_data <- data.frame(RMSE = rmse_res, t(rhat_summary), t(params))
ori_res_data

# テストデータの観測値と予測値を保存しておく
ori_obs_test <- ori_i
ori_pred_test <- pred_i

#全体（教師＋テスト）の観測値と予測値を保存しておく
ori_obs_all <- ori_incidence
ori_pred_all <- pred_incidence

# //////////////////////////////////////////////////////////////////////////////
################################# 全体の予測値を算出する ########################

# 47都道府県×12月のマトリクスをつくって埋めていく
est_heat <- matrix(ncol=47, nrow=12) #カラのマトリクス

# 県の番号と位置情報を読み込み # pref_name_location.csv
pref_NL <- read.csv("sihou/pref_name_location.csv", header=T)
pref_NL$pref %<>% tolower()

rownames(est_heat) <- 1:12 #行：月の名前の略名
colnames(est_heat)<- pref_NL$pref #列：都道府県名


####### lambda（気温、降水量、緯度、経度の影響）のマトリクス
est_lam <- matrix(nrow=12, ncol=47) #行：月×列：県として NAを並べたもの
rownames(est_lam) <- 1:12 #行に月の名前の略名
pref_name1 <- pref_NL$pref
colnames(est_lam)<- pref_name1 #列に県名

for(i in 1:12){
  for(j in 1:47){
    est_lam[i,j]<-paratemp_mean*airtemp_mat[i,j]+paraprecip_mean*precip_mat[i,j]+paralongi_mean*p_locat[j,3]+paralati_mean*p_locat[j,4]
  }
}


###### lambdaのマトリクス（est_lam）とm_muをつかってest_heatデータをつくる
##1月　マトリクス1行目
for(i in 1:47){
  est_heat[1,i] <- mmu_mean[1]+est_lam[1,i]
}

## 2月以降　マトリクス2行目以降
for (i in 2:12){
  for(j in 1:47){
    est_heat[i, j] <- est_heat[i-1,j]+mmu_mean[i]+est_lam[i,j]
  }
}

###負の値が入ってしまった場合に0に戻す
est_heat %<>% ifelse(.<0, 0, .)

########## 予測値のデータを保存
write.csv(est_heat, paste0("est_files/PRED__", file_name, ".csv"))

exdat_normal <- list(rhat_summary, params, mmu_mean, est_lam, pred_incidence, ori_incidence)

#重いオブジェクトを削除する
rm(fit_list) 
rm(exfit_list) 


# //////////////////////////////////////////////////////////////////////////////
# 帰無モデル -------------------------------------------------------------------
# //////////////////////////////////////////////////////////////////////////////

# スタン用データ作成
# 10 fold cross validation
set.seed(100)
# まずはシャッフル
dat_null <- dat %>% sample_frac(size=1)
#dat_null
## ランダマイズ
dat_null <- dat_null %>% apply(2, sample) %>% as.data.frame()
#これだとすべて文字列になってしまうので戻す
dat_null[,1:2] <- dat_null[,1:2] %>% apply(2, as.numeric)
dat_null[,5:13] <- dat_null[,5:13] %>% apply(2, as.numeric)

#dat_null %>% as_tibble() %>% head() # OK

# 次に分割
spl_dat <- dat_null %>% rsample::vfold_cv(v=10)
#analysis(spl_dat$splits[[1]]) # 学習データを表示する
#assessment(spl_dat$splits[[1]]) # validationデータを表示する

train_list <- apply(spl_dat, 1, function(x)analysis(x$splits))
# train_list[[1]] %>% dim() # 657 rows

valid_list <- apply(spl_dat, 1, function(x)rsample::assessment(x$splits))
# valid_list[[1]] %>% dim()  # 73 rows


# データランダマイズ
data_list <- NULL 
trainN <- train_list[[1]] %>% nrow()

for (i in 1:10){ # 10分割したデータを使う
  dat_piece <- rbind(train_list[[i]], valid_list[[i]]) %>% as.data.frame()
  # データの準備
  stan_data <- list(
    I = length(dat_piece$res_val),
    I_train = trainN,
    P = max(dat_piece$Pref_N2), #県の数
    M = 12, #12months
    mu_zero = 0,
    Pref = as.integer(dat_piece$Pref_N2),
    Month = as.integer(dat_piece$Month),
    Temp = airtemp_mat2,
    Precip = precip_mat2,
    Longi = as.numeric(p_locat2$longi),
    Lati = as.numeric(p_locat2$lati),
    Incidence = dat_piece$res_val # res_val
  )
  data_list %<>%  append(list(stan_data))
}

# stan をlapply計算 
fit_list <- future_map(data_list, rstan::sampling, object = model3, iter=iterN, warmup=warmupN, chains=chainsN, verbose=F)


#結果の要約
rhat_summary <- fit_list %>% sapply(rhat) %>% apply(2, summary)

exfit_list <- fit_list %>% lapply(rstan::extract)
#parameters
params <- exfit_list %>% sapply(function(x)sapply(x, mean))
#params

#pred_incidence予測値
pred_incidence <- exfit_list %>% sapply(function(x)x %$% Incidence_new %>% apply(2, mean))

#ori_incidence観測値
ori_incidence <- data_list %>% sapply(function(x)x$Incidence) #データフレームにする

exdat_rand <- list(rhat_summary, params, pred_incidence, ori_incidence)

#RMSE テストデータ
ori_i <- ori_incidence[(trainN+1):nrow(ori_incidence),]
pred_i <- pred_incidence[(trainN+1):nrow(pred_incidence),]
ori_pred_list <- list(ori_i, pred_i) # リストにまとめる
rmse_res <- sapply(1:ncol(ori_i), function(x)rmse(ori_pred_list[[1]][,x], ori_pred_list[[2]][,x])) # Metrics::rmse(actual, predicted)

# データフレームにまとめる
null_res_data <- data.frame(RMSE = rmse_res, t(rhat_summary), t(params))

# テストデータの観測値と予測値を保存しておく
null_obs_test <- ori_i
null_pred_test <- pred_i

#全体（教師＋テスト）の観測値と予測値を保存しておく（あとから取り出せるように）
null_obs_all <- pred_incidence
null_pred_all <- ori_incidence


#重いオブジェクトを削除する
rm(fit_list) 
rm(exfit_list) 


# /////////////////////////////////////////////////
# 保存 ----------------------------------------------------------------------
# /////////////////////////////////////////////////


rmse_table <- data.frame(rmse_ori = ori_res_data$RMSE, rmse_null = null_res_data$RMSE) #面積/作付面積
# 帰無モデル平均値
null_mean <- rmse_table$rmse_null %>% mean()
# 帰無モデル標準偏差
null_sd <- rmse_table$rmse_null %>% sd()

rmse_table <- rmse_table %>% mutate(z_score = (rmse_ori - null_mean)/null_sd)

write.csv(rmse_table, paste0("est_files/RMSE__", file_name, ".csv"))



##残りのオブジェクトをまとめて保存
save.image(paste0("est_files/", file_name, ".RData"))

