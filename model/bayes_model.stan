data {
  int I; // データ数
  int I_train; //教師データ数
  int P; //出てくる県の数 47
  int M; // 12カ月
  real mu_zero; //もう０を与えることに

  int<lower=0, upper=P> Pref[I]; //県
  int<lower=0, upper=M> Month[I]; //月
  matrix[M, P] Temp; //月・県平均気温
  matrix[M, P] Precip; //月・県降水量
  vector[P] Longi; //県の緯度（県庁所在地）
  vector[P] Lati; //県の経度（同）
  vector[I] Incidence; //発病葉率
}

parameters {
  matrix<lower=0, upper=100>[M, P] mu; //県別、月別の発病葉率
  vector<lower=-100, upper=100>[M] m_mu; //月別の増減を表す効果（全国共通）
  real<lower = 0> s_a; //状態変化の分散
  real<lower = 0> s_o; //観測誤差
  
  real para_temp; //降水量が発病葉率に与える影響
  real para_precip; //同・降水量
  real para_longi; //同・緯度
  real para_lati; //同・経度

}

transformed parameters{ //lambda
  matrix[M, P] lambda; //各月・各県の影響
  for (i in 1:M){
    for (j in 1:P){
        lambda[i, j] = para_temp*Temp[i,j]+para_precip*Precip[i,j]+para_longi*Longi[j]+para_lati*Lati[j];
    }
  }
}

model {
  for(i in 1:P){
      mu[1,i] ~normal(mu_zero+m_mu[1]+lambda[1,i], s_a);//県ごとの初期値の設定
  }
  //各県の状態空間変化
  for(j in 2:M){
    for(k in 1:P){
      mu[j, k]~normal(mu[j-1,k]+m_mu[j]+lambda[j,k], s_a);//各県について状態変化する
    }
  }
  
  //観測できた発病率
  for(i in 1:I_train){
    Incidence[i]~normal(mu[Month[i], Pref[i]], s_o); //観測値
  }

}

generated quantities{ //予測値を事後サンプリング 全部
  vector[I] Incidence_new;
  for(i in 1:I){
    Incidence_new[i]=normal_rng(mu[Month[i], Pref[i]], s_o);
  }
}

