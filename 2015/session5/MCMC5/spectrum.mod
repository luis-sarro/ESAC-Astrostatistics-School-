data {
  int<lower=1> N;
  int nlines;
  vector[nlines] center;
  vector[N] observed;
  vector[N] lambda;
  matrix[N,N] Sigma;
}

parameters {
  real h1;
  real h2;
  real h3;
  real h4;
  real h5;
  real w1;
  real w2;
  real w3;
  real w4;
  real w5;
  real cont_const;
  real cont_slope;
}

model {
  vector[N] predicted;
  real tmp;
  
  h1 ~ uniform(0.01,1.0);
  h2 ~ uniform(0.01,1.0);
  h3 ~ uniform(0.01,1.0);
  h4 ~ uniform(0.01,1.0);
  h5 ~ uniform(0.01,1.0);
  w1 ~ uniform(1.0,5);
  w2 ~ uniform(1.0,5);
  w3 ~ uniform(1.0,5);
  w4 ~ uniform(1.0,5);
  w5 ~ uniform(1.0,5);
  cont_const ~ uniform(0,0.1);
  cont_slope ~ uniform(-0.5,0.5);
  
  for (i in 1:N)
    {
    predicted[i] <- cont_const+cont_slope*lambda[i]+
    h1*exp(normal_log(lambda[i],center[1],w1))+
    h2*exp(normal_log(lambda[i],center[2],w2))+
    h3*exp(normal_log(lambda[i],center[3],w3))+
    h4*exp(normal_log(lambda[i],center[4],w4))+
    h5*exp(normal_log(lambda[i],center[5],w5));
   }

  increment_log_prob(multi_normal_log(observed,predicted,Sigma));

}
