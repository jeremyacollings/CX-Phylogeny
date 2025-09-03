data {
  int N; // number of plants
  int FS; // number of focal species
  int CS; // number of comp species
  int sp[N]; // species of focal plants
  int comp[N]; // species of competitor plants
  int count[N]; // number of competitor plants
  int fec[N]; // fecundity of focal plants=
}

parameters { 
  vector<lower=0>[FS] lambda; // population lambda per sp
  vector<lower=0>[FS] phi; // dispersion of fec per sp
  matrix<lower=0>[FS,CS] alpha; // population alphas
}

transformed parameters {
  real mu[N];
for(n in 1:N){
  mu[n] = lambda[sp[n]]*exp(-(alpha[sp[n], comp[n]]*count[n]));
  }
}

model {
  for(s in 1:FS){
  lambda[s] ~ normal(100, 100);
  phi[s] ~ exponential(1);
  }
  for(fs in 1:FS){
    for(cs in 1:CS){
      alpha[fs, cs] ~ normal(0,.5);
    }
  }
  for(n in 1:N){
  fec[n] ~ neg_binomial_2(mu[n], phi[sp[n]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = neg_binomial_2_lpmf(fec[n] | mu[n], phi[sp[n]]);
  }
}
