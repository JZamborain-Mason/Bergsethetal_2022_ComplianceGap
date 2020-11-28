data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 int pa_tp[res];  // PA top predators in reserves
 int pa_tp2[fis];  // PA top predators in fished
}
parameters {
 real I_reserves;
 real I_fished;
}
transformed parameters {
 vector[res] mu;
 vector[fis] mu2;
//reserve component
for (i in 1:res){ 
  mu[i] = I_reserves;
 }
//fished component
for (i in 1:fis){ 
  mu2[i] = I_fished;
 }
}
model {
 //priors
 I_fished ~ normal(-5,10);
 I_reserves ~ normal(-5,10);
 //likelihoods  
 for(n in 1:res){
      pa_tp[n] ~ bernoulli_logit(mu[n]);
}
 for(n in 1:fis){
       pa_tp2[n] ~ bernoulli_logit(mu2[n]);
}
}
generated quantities {
vector[res+fis] log_lik; //log-likelihood (for loo)
for (i in 1:res) {
 log_lik[i] = bernoulli_logit_lpmf(pa_tp[i]| mu[i]);
}
for (i in 1:fis) {
 log_lik[i+res] = bernoulli_logit_lpmf(pa_tp2[i]| mu2[i]);
 }
}
