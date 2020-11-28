data {
 int<lower=1> fis;  
 int<lower=1> res;  
 vector[res] b;  
 vector[fis] b2;  
}
parameters {
 real I_reserves;
 real I_fished;
 real<lower=0> sigma_e; 
 real<lower=0> sigma_f; 
}
transformed parameters {
 vector[res] mu;
 vector[fis] mu2;
 //reserve component
for (i in 1:res){ 
  mu[i] =  I_reserves;
 }
//fished component
for (i in 1:fis){ 
  mu2[i] = I_fished;
  }
}
model {
 //priors
 sigma_e ~ cauchy(0,2.5); 
 sigma_f ~ cauchy(0,2.5);
 I_fished ~ normal(5,10);
 I_reserves ~ normal(5,10);
 //likelihoods  
 for(n in 1:res){
      b[n] ~ normal(mu[n],sigma_e);
}
 for(n in 1:fis){
      b2[n] ~ normal(mu2[n],sigma_f);
}
}
generated quantities {
vector[res+fis] log_lik; //log-likelihood (for loo)
 for (i in 1:res) {
 log_lik[i] = normal_lpdf(b[i]| mu[i], sigma_e);
}
for (i in 1:fis) {
 log_lik[i+res] = normal_lpdf(b2[i]| mu2[i], sigma_f);
 }
}
