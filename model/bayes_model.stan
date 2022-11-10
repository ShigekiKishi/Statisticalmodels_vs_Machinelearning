data {
  int I; // �f�[�^��
  int I_train; //���t�f�[�^��
  int P; //�o�Ă��錧�̐� 47
  int M; // 12�J��
  real mu_zero; //�����O��^���邱�Ƃ�

  int<lower=0, upper=P> Pref[I]; //��
  int<lower=0, upper=M> Month[I]; //��
  matrix[M, P] Temp; //���E�����ϋC��
  matrix[M, P] Precip; //���E���~����
  vector[P] Longi; //���̈ܓx�i�������ݒn�j
  vector[P] Lati; //���̌o�x�i���j
  vector[I] Incidence; //���a�t��
}

parameters {
  matrix<lower=0, upper=100>[M, P] mu; //���ʁA���ʂ̔��a�t��
  vector<lower=-100, upper=100>[M] m_mu; //���ʂ̑�����\�����ʁi�S�����ʁj
  real<lower = 0> s_a; //��ԕω��̕��U
  real<lower = 0> s_o; //�ϑ��덷
  
  real para_temp; //�~���ʂ����a�t���ɗ^����e��
  real para_precip; //���E�~����
  real para_longi; //���E�ܓx
  real para_lati; //���E�o�x

}

transformed parameters{ //lambda
  matrix[M, P] lambda; //�e���E�e���̉e��
  for (i in 1:M){
    for (j in 1:P){
        lambda[i, j] = para_temp*Temp[i,j]+para_precip*Precip[i,j]+para_longi*Longi[j]+para_lati*Lati[j];
    }
  }
}

model {
  for(i in 1:P){
      mu[1,i] ~normal(mu_zero+m_mu[1]+lambda[1,i], s_a);//�����Ƃ̏����l�̐ݒ�
  }
  //�e���̏�ԋ�ԕω�
  for(j in 2:M){
    for(k in 1:P){
      mu[j, k]~normal(mu[j-1,k]+m_mu[j]+lambda[j,k], s_a);//�e���ɂ��ď�ԕω�����
    }
  }
  
  //�ϑ��ł������a��
  for(i in 1:I_train){
    Incidence[i]~normal(mu[Month[i], Pref[i]], s_o); //�ϑ��l
  }

}

generated quantities{ //�\���l������T���v�����O �S��
  vector[I] Incidence_new;
  for(i in 1:I){
    Incidence_new[i]=normal_rng(mu[Month[i], Pref[i]], s_o);
  }
}

