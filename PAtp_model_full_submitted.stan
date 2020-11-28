data {
 int<lower=1> fis;  // number of observations fished
 int<lower=1> res;  // number of observations reserves
 int pa_tp[res];  // PA top predators in reserves
 int pa_tp2[fis];  // PA top predators in fished
//explanatory variables for each component
 real ag[res]; //age reserve (only for reserves)
 real si[res]; //predictor reserve size (only for reserves)
 real sa[res]; //predictor sampling area
 real op[res]; //predictor ocean productivity 
 real cs[res]; //climate stress
 real pg[res]; //predictor popgrowth 
 real gtot[res]; //predictor gravity
 real rl[res]; //predictor reef landings
 real ps[res]; //predictor pop size
 real hdi[res]; //predictor hdi
 int at[res];//predictor atoll
 int rh_b[res];//predictor backreef/lagoon habitat
 int rh_f[res];//predictor flat habitat
 int rh_c[res];//predictor crest habitat
 int cm_pc[res];//predictor census method point count
 int ds[res]; //predictor depth shallow
 int dd[res]; //predictor depth deep
 int hc[res]; //predictor highcompliance (only for reserves)
 real sa2[fis]; //predictor sampling area
 real op2[fis]; //predictor ocean productivity 
 real cs2[fis]; //climate stress
 real pg2[fis]; //predictor popgrowth 
 real gtot2[fis]; //predictor grav
 real rl2[fis]; //predictor reef landings
 real ps2[fis]; //predictor pop size
 real hdi2[fis]; //predictor hdi
 int at2[fis];//predictor atoll
 int rh_b2[fis];//predictor backreef/lagoon habitat
 int rh_f2[fis];//predictor flat habitat
 int rh_c2[fis];//predictor crest habitat
 int cm_pc2[fis];//predictor census method point count
 int cm_ds[fis]; //predictor census method distance sampling
 int ds2[fis]; //predictor depth shallow
 int dd2[fis]; //predictor depth deep
 int mr[fis]; //predictor restrictions
//random effects level 1: reef cluster(nested within larger)
 int<lower=1> RC; //number of reef clusters for reserves (groups)
 int<lower=1, upper=RC> prc[res]; //reef clusters id 
 int<lower=1> RC2; //number of reef clusters fished(groups)
 int<lower=1, upper=RC2> prc2[fis]; //reef clusters id 
//random effects level 2:larger
 int<lower=1> R; //number of data regions reserves (groups)
 int<lower=1, upper=R> pr[RC]; //region id for each reef cluster
 int<lower=1> R2; //number of data regions fished(groups)
 int<lower=1, upper=R2> pr2[RC2]; //region id for each reef cluster
}
parameters {
 vector[20] beta; //effect sizes
 real I_reserves; //intercept log-odds reserves
 real I_fished;//intercept log-odds fished
 vector[R] u; // random effects for reserves (deviation between regions)
 vector[R2] u2; // random effects for fished (deviation between regions)
 real<lower=0> sigma_u; //deviation sd 
 real<lower=0> sigma_u2; //deviation sd 
 vector[RC] urc; // random effects for reserve (deviation between reef clusters within regions)
 vector[RC2] urc2; // random effects for remote (deviation between reef clusters within regions)
 real<lower=0> sigma_urc; //deviation sd 
 real<lower=0> sigma_urc2; //deviation sd 
}
transformed parameters {
 vector[res] mu;//log odds
 vector[fis] mu2;//log odds
 //varying intercepts by reef clusters and regions
  real intercept_reserves_rc[RC];
  real intercept_reserves_r[R];
  real intercept_fished_rc[RC2];
  real intercept_fished_r[R2];
  //compute the varying intercept at the region level (non-centered)
  for(i in 1:R){
    intercept_reserves_r[i] = I_reserves + (sigma_u*u[i]);
   }
  for(i in 1:R2){
    intercept_fished_r[i] = I_fished + (sigma_u2*u2[i]);
   }
  //compute varying intercept at the reef cluster within region level
  for(i in 1:RC){
     intercept_reserves_rc[i] = intercept_reserves_r[pr[i]] + (sigma_urc*urc[i]);
  }
  for(i in 1:RC2){
     intercept_fished_rc[i] = intercept_fished_r[pr2[i]] +  (sigma_urc2*urc2[i]);
  }
//reserve component
for (i in 1:res){ 
  mu[i] = intercept_reserves_rc[prc[i]] +beta[1]*dd[i]+beta[2]*ds[i]+beta[3]*rh_c[i]+beta[4]*rh_f[i]+beta[5]*rh_b[i]+beta[6]*cm_pc[i]+beta[8]*sa[i]+ beta[9]*si[i]+ beta[10]*ag[i]+ beta[11]*at[i]+ beta[12]*op[i]+ beta[13]*cs[i]+ beta[14]*pg[i]+ beta[15]*gtot[i]+  beta[16]*hc[i]+beta[17]*rl[i]+ beta[18]*ps[i]+ beta[19]*hdi[i];
 }
//fished component
for (i in 1:fis){ 
  mu2[i] = intercept_fished_rc[prc2[i]]+beta[1]*dd2[i]+beta[2]*ds2[i]+beta[3]*rh_c2[i]+beta[4]*rh_f2[i]+beta[5]*rh_b2[i]+beta[6]*cm_pc2[i]+beta[7]*cm_ds[i]+beta[8]*sa2[i]+ beta[11]*at2[i]+ beta[12]*op2[i]+ beta[13]*cs2[i]+ beta[14]*pg2[i]+ beta[15]*gtot2[i]+ beta[17]*rl2[i]+ beta[18]*ps2[i]+ beta[19]*hdi2[i]+  beta[20]*mr[i];
  }
}
model {
 //priors
 beta[1] ~ normal (0,10); //prior slope
 beta[2] ~ normal (0,10); //prior slope
 beta[3] ~ normal (0,10); //prior slope
 beta[4] ~ normal (0,10); //prior slope
 beta[5] ~ normal (0,10); //prior slope
 beta[6] ~ normal (0,10); //prior slope
 beta[7] ~ normal (0,10); //prior slope
 beta[8] ~ normal (0,10); //prior slope
 beta[9] ~ normal (0,10); //prior slope
 beta[10] ~ normal (0,10); //prior slope
 beta[11] ~ normal (0,10); //prior slope
 beta[12] ~ normal (0,10); //prior slope
 beta[13] ~ normal (0,10); //prior slope
 beta[14] ~ normal (0,10); //prior slope
 beta[15] ~ normal (0,10); //prior slope
 beta[16] ~ normal (0,10); //prior slope
 beta[17] ~ normal (0,10); //prior slope
 beta[18] ~ normal (0,10); //prior slope
 beta[19] ~ normal (0,10); //prior slope
 beta[20] ~ normal (0,10); //prior slope
 sigma_u ~ cauchy (0,2.5); //prior sd for group varying intercept
 sigma_u2 ~ cauchy (0,2.5); //prior sd for group varying intercept
 u ~ normal(0,  1); //prior re (distribution of the varying intercept)
 u2 ~ normal(0,  1); //prior re
 sigma_urc ~ cauchy (0,2.5); //prior sd for group varying intercept
 sigma_urc2 ~ cauchy (0,2.5); //prior sd for group varying intercept
 urc ~ normal(0,  1); //prior re
 urc2 ~ normal(0,  1); //prior re
 I_fished ~ normal(-5,10); //prior in//prior intercepttercept
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
vector[res] y_rep;
vector[fis] y_rep2;
vector[res+fis] mu_fished; // log-odds for our sites if they were fished reefs
vector[res+fis] mu_lowcompliance; // log-odds for our sites if they were low compliance reserves at average reserve size and age
vector[res+fis] mu_highcompliance; //log-odds for our sites if they were high compliance reserves at average reserve size and age
vector[res+fis] p_fished; // prob observing top predators for our sites if they were fished reefs
vector[res+fis] p_lowcompliance; // prob observing top predators for our sites if they were low compliance reserves
vector[res+fis] p_highcompliance; //prob observing top predators for our sites if they were high compliance reserves
 for (i in 1:res) {
 log_lik[i] = bernoulli_logit_lpmf(pa_tp[i]| mu[i]);
 y_rep[i] =bernoulli_logit_rng(mu[i]);
 mu_fished[i]=I_fished+ beta[11]*at[i]+ beta[12]*op[i]+ beta[13]*cs[i]+ beta[14]*pg[i]+ beta[15]*gtot[i]+ beta[17]*rl[i]+ beta[18]*ps[i]+ beta[19]*hdi[i];
 p_fished[i] =bernoulli_logit_rng(mu_fished[i]);
 mu_lowcompliance[i]=I_reserves+ beta[11]*at[i]+ beta[12]*op[i]+ beta[13]*cs[i]+ beta[14]*pg[i]+ beta[15]*gtot[i]+ beta[17]*rl[i]+ beta[18]*ps[i]+ beta[19]*hdi[i];
 p_lowcompliance[i] =bernoulli_logit_rng(mu_lowcompliance[i]);
 mu_highcompliance[i]=I_reserves+beta[16]+ beta[11]*at[i]+ beta[12]*op[i]+ beta[13]*cs[i]+ beta[14]*pg[i]+ beta[15]*gtot[i]+ beta[17]*rl[i]+ beta[18]*ps[i]+ beta[19]*hdi[i];
 p_highcompliance[i] =bernoulli_logit_rng(mu_highcompliance[i]);
}
for (i in 1:fis) {
 log_lik[i+res] = bernoulli_logit_lpmf(pa_tp2[i]| mu2[i]);
 y_rep2[i] =bernoulli_logit_rng(mu2[i]);
 mu_fished[i+res]=I_fished+ beta[11]*at2[i]+ beta[12]*op2[i]+ beta[13]*cs2[i]+ beta[14]*pg2[i]+ beta[15]*gtot2[i]+ beta[17]*rl2[i]+ beta[18]*ps2[i]+ beta[19]*hdi2[i];
 p_fished[i+res] =bernoulli_logit_rng(mu_fished[i+res]);
 mu_lowcompliance[i+res]=I_reserves+ beta[11]*at2[i]+ beta[12]*op2[i]+ beta[13]*cs2[i]+ beta[14]*pg2[i]+ beta[15]*gtot2[i]+ beta[17]*rl2[i]+ beta[18]*ps2[i]+ beta[19]*hdi2[i];
 p_lowcompliance[i+res] =bernoulli_logit_rng(mu_lowcompliance[i+res]);
 mu_highcompliance[i+res]=I_reserves+beta[16]+ beta[11]*at2[i]+ beta[12]*op2[i]+ beta[13]*cs2[i]+ beta[14]*pg2[i]+ beta[15]*gtot2[i]+ beta[17]*rl2[i]+ beta[18]*ps2[i]+ beta[19]*hdi2[i];
 p_highcompliance[i+res] =bernoulli_logit_rng(mu_highcompliance[i+res]);
}
}
