#Code for: Closing the compliance gap in marine reserves
#ARC Centre of Excellence for Coral Reef Studies
#R version 3.4.2 (2017-09-28)


#clean the entire environemnt
rm(list=ls()) 

#set the working directory
setwd("c:/Users/jzamb/Dropbox/PhD/Research Worker Josh 2018/Compliance gap/R") 

##load required libraries
library(ggplot2) ##H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016
library(MuMIn) #Kamil Barton (2019) #. MuMIn: Multi-Model Inference. R package version 1.43.6.
library(ggpubr) #Alboukadel Kassambara (2018). ggpubr: 'ggplot2' Based Publication Ready Plots. 
library(rstan) #Stan Development Team (2018). RStan: the R interface to Stan.
library(bayesplot) #Jonah Gabry and Tristan Mahr (2018). bayesplot: Plotting for Bayesian Models. 
library(tidyverse) #Hadley Wickham (2017). tidyverse: Easily Install and Load the 'Tidyverse'. 
library(broom) #David Robinson and Alex Hayes (2019). broom: Convert Statistical Analysis Objects into Tidy Tibbles.
library(coda) #Martyn Plummer, Nicky Best, Kate Cowles and Karen Vines (2006). CODA: Convergence Diagnosis and Output.Analysis for MCMC, R News, vol 6, 7-11
library(data.table) #Matt Dowle and Arun Srinivasan (2019). data.table: Extension of `data.frame`
library(lme4) #Douglas Bates, Martin Maechler, Ben Bolker, Steve Walker (2015) #. Fitting Linear Mixed-Effects Models.Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
library(DHARMa)#Florian Hartig (2019). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.2.4. https://CRAN.R-project.org/package=DHARMa
library(plyr)#Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.
library(loo)#Vehtari A, Gabry J, Yao Y, Gelman A (2019). "loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models." R package version 2.1.0, <URL: https://CRAN.R-project.org/package=loo>.
library(rworldmap)#South, Andy 2011 rworldmap: A New R package for Mapping Global Data. The R Journal Vol. 3/1 : 35-43.
library(rworldxtra)#  Andy South (2012). rworldxtra: Country boundaries at high resolution.. R package version 1.01. https://CRAN.R-project.org/package=rworldxtra


###   FUNCTIONS   ###

#Standardize function for continuous variables
standardize = function(x){(x-mean(x, na.rm=T))/(2*sd(x, na.rm=T))} 

#variance inflation factors  function for mixed models
vif.mer = function (fit) {
  ## adapted from rms::vif
  v <-vcov(fit)
  nam = names(fixef(fit))
  ## exclude intercepts
  ns = sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v = v[-(1:ns), -(1:ns), drop = FALSE]
    nam = nam[-(1:ns)]
  }
  d = diag(v)^0.5
  v = diag(solve(v/(d %o% d)))
  names(v) = nam
  v
}
#correlation function
panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y, method = "pearson",use = "complete.obs"))
  txt = format(c(r, 0.123456789), digits=digits)[1]
  txt = paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor = 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r*2)
}


#upload data
alldata<-read.csv("compliancegap_data.csv",header=T) 


##Explanatory variables: relevel categorical variables and standardize continuous variables

#reef-scale covariates
alldata$DepthCategory<-relevel(alldata$DepthCategory,ref="4-10m")
alldata$CleanHabitat<-relevel(alldata$CleanHabitat,ref="Slope")
alldata$Protection<-relevel(alldata$FinalProtection4,ref="Fished")
alldata$CensusMethod<-relevel(alldata$CensusMethod,ref="Standard belt transect")
alldata$sTotal_sampling_area<-standardize(log(alldata$Total.sampling.Area))
alldata$sgrav_tot2<-standardize(log(alldata$gravtot_500km+min(alldata$gravtot_500km[alldata$gravtot_500km>0])))

#reef cluster-scale covariates
alldata$sOcean_prod<-standardize(log(alldata$Ocean_prod))
alldata$sClimate_stress<-standardize(alldata$Climate_stress)
alldata$sRegional_population_growth<-standardize(alldata$Regional_population_growth)

#nation/state-scale covariates
alldata$sReef_fish_landings_per_km2<-standardize(log(alldata$spatialcatch_tkm2+1))
alldata$sLarger_pop_size<-standardize(log(alldata$Larger_pop_size+1))
alldata$sHDI<-standardize(alldata$HDI)

#add dummy variables for factors
alldata$crest=ifelse(alldata$CleanHabitat=="Crest",1,0)
alldata$backreef=ifelse(alldata$CleanHabitat=="Lagoon_Back reef",1,0)
alldata$flat=ifelse(alldata$CleanHabitat=="Flat",1,0)
alldata$distancesampling=ifelse(alldata$CensusMethod=="Distance sampling",1,0)
alldata$pointcount=ifelse(alldata$CensusMethod=="Point intercept",1,0)
alldata$shallow=ifelse(alldata$DepthCategory=="0-4m",1,0)
alldata$deep=ifelse(alldata$DepthCategory==">10m",1,0)

#colinaerity
VIF.table<-as.data.frame(vif.mer(lmer(log(Biomass...some.families)~
                                        DepthCategory+ CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                                        log(MPAage+1)+log(NTZarea+1)+
                                        sRegional_population_growth+sOcean_prod+sClimate_stress+
                                        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger_2/SocialSIte4km), data=alldata)))
colnames(VIF.table)="VIF"
print(VIF.table)
#write.csv(VIF.table,"VIF_fullmodel_mpa.csv")
#problem of doing that way is the correlation between management and reserve age/size

VIF.table<-as.data.frame(vif.mer(lmer(log(Biomass...some.families)~
                                        DepthCategory+ CleanHabitat+Protection+CensusMethod+sTotal_sampling_area+sgrav_tot2+
                                        log(NTZarea+1)+
                                        sRegional_population_growth+sOcean_prod+sClimate_stress+
                                        sHDI+sLarger_pop_size+sReef_fish_landings_per_km2+(1|Larger_2/SocialSIte4km), data=alldata)))
colnames(VIF.table)="VIF"
print(VIF.table)
#write.csv(VIF.table,"VIF_fullmodel_mpa_nosize.csv")


## Response variables  ##

#targeted biomass and presence absence of top predators
length(alldata$Targeted.biomass[alldata$Targeted.biomass==0])/length(alldata$Targeted.biomass)
hist(alldata$PA_tp)
hist(log(alldata$Targeted.biomass+1))
alldata$tTargBiomass<-log(alldata$Targeted.biomass+min(alldata$Targeted.biomass[alldata$Targeted.biomass>0]))

#countries included
countries<-as.data.frame(unique(alldata$Larger_2))
#write.csv(countries, "countries.included.csv",row.names=F)


## Prepare data for stan  ##

#separate model subcomponents 
reserve_data=alldata[(alldata$Protection=="UnfishedHigh"|alldata$Protection=="UnfishedLow" )& alldata$MPAage >0 & alldata$NTZarea>0,]
reserve_data=droplevels(reserve_data)
fished_data=alldata[alldata$Protection=="Fished"|alldata$Protection=="Restricted",]
fished_data=droplevels(fished_data)

#add parts that are only for certain model subcomponents
fished_data$restricted=ifelse(fished_data$Protection=="Restricted", 1, 0)
reserve_data$highcompliance=ifelse(reserve_data$Protection=="UnfishedHigh", 1, 0)
reserve_data$sNTZarea=standardize(log(reserve_data$NTZarea))
reserve_data$sMPAage=standardize(log(reserve_data$MPAage))


#add indexes for reef cluster (within larger) and larger
reserve_data$indexrc=as.numeric(reserve_data$SocialSIte4km)
fished_data$indexrc=as.numeric(fished_data$SocialSIte4km)
nlevels(reserve_data$Larger_2)
nlevels(reserve_data$SocialSIte4km)
nlevels(fished_data$Larger_2)
nlevels(fished_data$SocialSIte4km)
nrow(reserve_data)
nrow(fished_data)
reserve_data_rc=ddply(reserve_data,.(SocialSIte4km),summarize, Larger_2=Larger_2[1],indexrc=indexrc[1]) 
reserve_data_rc$indexj=as.numeric(reserve_data_rc$Larger_2)
fished_data_rc=ddply(fished_data,.(SocialSIte4km),summarize, Larger_2=Larger_2[1],indexrc=indexrc[1]) 
fished_data_rc$indexj=as.numeric(fished_data_rc$Larger_2)
reserve_data=merge(reserve_data,reserve_data_rc[,c("SocialSIte4km","indexj")],by="SocialSIte4km",all.x=T)
fished_data=merge(fished_data,fished_data_rc[,c("SocialSIte4km","indexj")],by="SocialSIte4km",all.x=T)

# Sort by reefcluster id to make it easier handling in Stan
reserve_data<- arrange(reserve_data, indexrc, indexj)
fished_data<- arrange(fished_data, indexrc, indexj)

regionLookupVec <- unique(reserve_data[c("indexrc","indexj")])[,"indexj"]
regionLookupVec2 <- unique(fished_data[c("indexrc","indexj")])[,"indexj"]

##Models###########################################################################################################################
rstan_options(auto_write = T)
options(mc.cores = parallel::detectCores())
#backup_options <- options()
#options(backup_options)

##Targeted biomass.................................................................................................................
stanDat_TB<- list(res=nrow(reserve_data),b = reserve_data$tTargBiomass,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                     dd=reserve_data$deep,  ds=reserve_data$shallow,
                     pg=reserve_data$sRegional_population_growth, gtot=reserve_data$sgrav_tot2,
                     rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,
                     hdi=reserve_data$sHDI,
                     at=reserve_data$Atoll,cs=reserve_data$sClimate_stress,
                     op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                     cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area, hc=reserve_data$highcompliance,
                     fis=nrow(fished_data),b2 = fished_data$tTargBiomass,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                     dd2=fished_data$deep,  ds2=fished_data$shallow,
                     pg2=fished_data$sRegional_population_growth, gtot2=fished_data$sgrav_tot2,
                     rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,
                     hdi2=fished_data$sHDI,
                     at2=fished_data$Atoll,cs2=fished_data$sClimate_stress,
                     op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                     cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                     R=nlevels(reserve_data$Larger_2),
                     R2=nlevels(fished_data$Larger_2),
                     RC=nlevels(reserve_data$SocialSIte4km),
                     RC2=nlevels(fished_data$SocialSIte4km),
                     pr=regionLookupVec,
                     pr2=regionLookupVec2,
                     prc=reserve_data$indexrc,
                     prc2=fished_data$indexrc)

#fit null and full model
Fit_null_tb <- stan(file = "biomass_model_null.stan", data = stanDat_TB,chains = 4,iter=30000,warmup=20000,thin=10)
Fit_full_tb <- stan(file = "biomass_model_full_submitted.stan", data = stanDat_TB,chains = 4,iter=30000,warmup=20000,thin=10,control = list(adapt_delta = 0.9,stepsize = 0.005,max_treedepth = 20))
#save(Fit_full_tb,file="Fit_full_tb.RData")

#model comparison
full_loglik<- extract_log_lik(Fit_full_tb, merge_chains = F)
r_eff_full <- relative_eff(exp(full_loglik)) 
loo_full <- loo(full_loglik,r_eff =r_eff_full)


null_loglik<- extract_log_lik(Fit_null_tb, merge_chains = F)
r_eff_null <- relative_eff(exp(null_loglik)) 
loo_null <- loo(null_loglik,r_eff =r_eff_null)
comp <- loo_compare(loo_null, loo_full)
print(comp, simplify=F)
modelcomparison=as.data.frame(print(comp, simplify=F))
row.names(modelcomparison)=c("TargetedBiomass_Full","TargetedBiomass_Null")

#best-fit model posterior
Fit_full_tb_summary <- summary(Fit_full_tb,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_full_tb<- as.data.frame(Fit_full_tb_summary$summary)
output_Fit_full_tb$parameter=row.names(output_Fit_full_tb)
list_of_draws_full <- as.data.frame(Fit_full_tb) #rows are iterations and columns are parameters
tb_posterior<- rstan::extract(Fit_full_tb)

#model diagnostics
a=stan_trace(Fit_full_tb, pars=c("beta","I_reserves","I_fished"))
b=stan_rhat(Fit_full_tb)
c=stan_ess (Fit_full_tb)
diag2=ggarrange(b,c,nrow=2,ncol=1)
windows()
ggarrange(a,diag2,widths=c(2,1),nrow=1,ncol=2)+ggtitle("Targeted biomass model diagmostics")

#model fit
row.names(output_Fit_full_tb)
pred_full<- output_Fit_full_tb[861:(861+nrow(reserve_data)+nrow(fished_data)-1),1]
resid_full<- c(reserve_data$tTargBiomass, fished_data$tTargBiomass)-pred_full

a_full<- ggplot(data=NULL,aes(x=pred_full,y=resid_full))+geom_point()+theme_classic()+ggtitle("")+xlab("Fitted ")+ylab("Residuals ")
b_full<- ggplot(NULL, aes(x = resid_full)) +
  geom_histogram(colour = "white", fill = "black") +theme_classic()+ggtitle("")+ylab("Count")+xlab("Residuals ")+xlim(c(-4,4))

#posterior predictive checks
joined_sim <- rstan::extract(Fit_full_tb)
n_sims <- length(joined_sim $lp__)
y_rep_reserves <- array(NA, c(n_sims, nrow(reserve_data)))
y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

for (s in 1:n_sims){
  y_rep_reserves[s,] <- rnorm(nrow(reserve_data), joined_sim$mu[s,], joined_sim$sigma_e[s])
  y_rep_fished[s,] <- rnorm(nrow(fished_data), joined_sim$mu2[s,], joined_sim$sigma_f[s])
}
bayesplot::color_scheme_set(scheme = "gray")
a<- bayesplot::ppc_dens_overlay(reserve_data$tTargBiomass,y_rep_reserves[1:4000,])+ggtitle("Reserves ")+guides(col=F)+ labs(y="Density",x = expression ("Biomass (log("~kg/ha*"))"))
b<- bayesplot::ppc_dens_overlay(fished_data$tTargBiomass,y_rep_fished[1:4000,])+ggtitle("Fished")+ labs(y="",x = expression ("Biomass (log("~kg/ha*"))"))

resid_fig<- ggarrange(a_full,b_full,nrow=1, ncol=2,labels=c("a","b"))
ppcheckfig<- ggarrange(a,b,nrow=1,ncol=2, widths=c(1,1.2),labels=c("c","d"))
windows()
ggarrange(resid_fig,ppcheckfig,nrow=2,ncol=1)

#posterior parameters
betas<- output_Fit_full_tb[1:20,c("50%","5%","95%")]
betas$variable<- c( "Depth >10m","Depth <4m", "Crest","Flat","Backreef/lagoon","Point count","Distance sampling","Sampling area", "Reserve size","Reserve age", "Atoll", "Ocean productivity","Climate stress","Population growth","Total gravity","High compliance","Reef fish landings","Population size","HDI", "Restricted fishing")
betas$sign<- ifelse(betas$`5%`<0 & betas$`95%` <0, "negative",ifelse(betas$`5%`>0 & betas$`95%`>0, "positive", "no effect"))
betas$strength<-ifelse(betas$sign=="no effect", "open", "closed")

betasfig_tb<- ggplot(betas,aes(x=variable,y=`50%`,ymin=`5%`,ymax=`95%`))+
  geom_pointrange(size=0.7,aes(col=sign,shape=strength))+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none")+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+coord_flip()+
  ylab("Effect size")+xlab("")+ggtitle("Targeted biomass")


##Presence/absence top predators ......................................................................................................

stanDat_tp<- list(res=nrow(reserve_data),pa_tp = reserve_data$PA_tp,ag=reserve_data$sMPAage, si=reserve_data$sNTZarea,
                  dd=reserve_data$deep,  ds=reserve_data$shallow,
                  pg=reserve_data$sRegional_population_growth, gtot=reserve_data$sgrav_tot2,
                  rl=reserve_data$sReef_fish_landings_per_km2, ps=reserve_data$sLarger_pop_size,
                  hdi=reserve_data$sHDI,
                  at=reserve_data$Atoll,cs=reserve_data$sClimate_stress,
                  op=reserve_data$sOcean_prod,rh_c=reserve_data$crest,rh_b=reserve_data$backreef,rh_f=reserve_data$flat,
                  cm_pc=reserve_data$pointcount, sa=reserve_data$sTotal_sampling_area, hc=reserve_data$highcompliance,
                  fis=nrow(fished_data),pa_tp2 = fished_data$PA_tp,mr=fished_data$restricted, cm_ds=fished_data$distancesampling,
                  dd2=fished_data$deep,  ds2=fished_data$shallow,
                  pg2=fished_data$sRegional_population_growth, gtot2=fished_data$sgrav_tot2,
                  rl2=fished_data$sReef_fish_landings_per_km2, ps2=fished_data$sLarger_pop_size,
                  hdi2=fished_data$sHDI,
                  at2=fished_data$Atoll,cs2=fished_data$sClimate_stress,
                  op2=fished_data$sOcean_prod,rh_c2=fished_data$crest,rh_b2=fished_data$backreef,rh_f2=fished_data$flat,
                  cm_pc2=fished_data$pointcount, sa2=fished_data$sTotal_sampling_area,
                  R=nlevels(reserve_data$Larger_2),
                  R2=nlevels(fished_data$Larger_2),
                  RC=nlevels(reserve_data$SocialSIte4km),
                  RC2=nlevels(fished_data$SocialSIte4km),
                   pr=regionLookupVec,
                  pr2=regionLookupVec2,
                  prc=reserve_data$indexrc,
                  prc2=fished_data$indexrc)

Fit_null_tp <- stan(file = "PAtp_model_null.stan", data = stanDat_tp,chains = 4,iter=30000,warmup=20000,thin=10)
Fit_full_tp <- stan(file = "PAtp_model_full_submitted.stan", data = stanDat_tp,chains = 4,iter=30000,warmup=20000,thin=10,control = list(adapt_delta = 0.9,stepsize = 0.005,max_treedepth = 20))
#save(Fit_full_tp,file="Fit_full_tp.RData")

#model comparison
full_loglik<- extract_log_lik(Fit_full_tp, merge_chains = F)
r_eff_full <- relative_eff(exp(full_loglik)) 
loo_full <- loo(full_loglik,r_eff =r_eff_full)


null_loglik<- extract_log_lik(Fit_null_tp, merge_chains = F)
r_eff_null <- relative_eff(exp(null_loglik)) 
loo_null <- loo(null_loglik,r_eff =r_eff_null)
comp2 <- loo_compare(loo_null, loo_full)
print(comp2, simplify=F)

row.names(comp2)=c("PA_toppredators_Full","PA_toppredators_Null")
modelcomparison=rbind(modelcomparison,as.data.frame(print(comp2, simplify=F)))
#write.csv(modelcomparison,"model_comparison_loo.csv")

#bestfit model posterior
Fit_full_tp_summary <- summary(Fit_full_tp,probs=c(0.05,0.25,0.5,0.75, 0.95))
output_Fit_full_tp<- as.data.frame(Fit_full_tp_summary$summary)
list_of_draws_full <- as.data.frame(Fit_full_tp) #rows are iterations and columns are parameters
tp_posterior<- rstan::extract(Fit_full_tp)

#model diagnostics
a=stan_trace(Fit_full_tp, pars=c("beta","I_reserves","I_fished"))
b=stan_rhat(Fit_full_tp)
c=stan_ess (Fit_full_tp)
diag2=ggarrange(b,c,nrow=2,ncol=1)
windows()
ggarrange(a,diag2,widths=c(2,1),nrow=1,ncol=2)+ggtitle("Top predator model diagmostics")

#model fit
joined_sim <- rstan::extract(Fit_full_tp)
x = createDHARMa(simulatedResponse=t(joined_sim$y_rep), observedResponse=reserve_data$PA_tp)
windows()
plot(x)
a<-ggplot(NULL)+geom_histogram(aes(x=x$scaledResiduals))+xlab("Scaled residuals top predator model")+ggtitle("Reserves")
b<-ggplot(data=NULL)+
  geom_point(aes(y=x$scaledResiduals, x=x$fittedPredictedResponse),pch=21,fill="darkgrey",alpha=0.4)+ylab("Scaled residuals")+xlab("Fitted ")+ggtitle("Reserves")

x2 = createDHARMa(simulatedResponse=t(joined_sim$y_rep2), observedResponse=fished_data$PA_tp)
windows()
plot(x2)
c<-ggplot(NULL)+geom_histogram(aes(x=x2$scaledResiduals))+xlab("Scaled residuals top predator model")+ggtitle("Fished")
d<-ggplot(data=NULL)+
  geom_point(aes(y=x2$scaledResiduals, x=x2$fittedPredictedResponse),pch=21,fill="darkgrey",alpha=0.4)+ylab("Scaled residuals")+xlab("Fitted ")+ggtitle("Fished")


n_sims <- length(joined_sim $lp__)
y_rep_reserves <- array(NA, c(n_sims, nrow(reserve_data)))
y_rep_fished <- array(NA, c(n_sims, nrow(fished_data)))

for (s in 1:n_sims){
  y_rep_reserves[s,] <-joined_sim$y_rep[s,]
  y_rep_fished[s,] <- joined_sim$y_rep2[s,]
}
bayesplot::color_scheme_set(scheme = "gray")
e<- bayesplot::ppc_dens_overlay(reserve_data$PA_tp ,y_rep_reserves[0:4000,])+ggtitle("Reserves ")+guides(col=F)+  xlab("Presence/absence top predators")
f<- bayesplot::ppc_dens_overlay(fished_data$PA_tp,y_rep_fished[0:4000,])+ggtitle("Fished")+ xlab("Presence/absence top predators")
windows()
ggarrange(a,b,c,d,e,f,nrow=3,ncol=2, widths=c(1,1.2),labels=c("a","b","c","d","e","f"))

#posterior parameters
betas_tp<- output_Fit_full_tp[1:20,c("50%","5%","95%")]
betas_tp$variable<- c( "Depth >10m","Depth <4m", "Crest","Flat","Backreef/lagoon","Point count","Distance sampling","Sampling area", "Reserve size","Reserve age", "Atoll", "Ocean productivity","Climate stress","Population growth","Total gravity","High compliance","Reef fish landings","Population size","HDI", "Restricted fishing")
betas_tp$sign<- ifelse(betas_tp$`5%`<0 & betas_tp$`95%` <0, "negative",ifelse(betas_tp$`5%`>0 & betas_tp$`95%`>0, "positive", "no effect"))
betas_tp$strength<-ifelse(betas_tp$sign=="no effect", "open", "closed")

betasfig_tp<- ggplot(betas_tp,aes(x=variable,y=`50%`,ymin=`5%`,ymax=`95%`))+
  geom_pointrange(size=0.7,aes(col=sign,shape=strength))+
  scale_color_manual(values=c("no effect"="grey48","negative"="darkred","positive"="darkgreen"),
                     labels=c("Negative","No effect","Positive"), guide=F)+
  scale_shape_manual(values=c("open"=1,"closed"=19), guide=F)+
  theme(legend.position = "none",axis.text.y = element_blank())+
  geom_hline(aes(yintercept=0),colour="darkgrey",linetype="dashed")+coord_flip()+
  ylab("Log odds")+xlab("")+ggtitle("Top predator presence")
windows()
ggarrange(betasfig_tb,betasfig_tp,ncol=2,nrow=1,labels=c("a","b"),widths=c(1.5,1))


###Compliance Gap######################################################################################################################

#order datasets
fished_data$sMPAage=rep(NA,nrow(fished_data))
fished_data$sNTZarea=rep(NA,nrow(fished_data))
reserve_data$restricted=rep(NA,nrow(reserve_data))
fished_data$highcompliance=rep(NA,nrow(fished_data))
reserve_data2=reserve_data[,c("SocialSIte4km","UniqueSite","Larger_2","Site_Lat","Site_Long","Protection", "Atoll", "sTotal_sampling_area","sgrav_tot2","sOcean_prod", "sClimate_stress","sRegional_population_growth","sReef_fish_landings_per_km2","sLarger_pop_size","sHDI",  "crest","backreef","flat","distancesampling","pointcount","shallow","deep","restricted","highcompliance", "indexrc","indexj","sMPAage","sNTZarea","PA_tp","tTargBiomass","Targeted.biomass" )]
fished_data2=fished_data[,c("SocialSIte4km","UniqueSite","Larger_2","Site_Lat","Site_Long","Protection", "Atoll", "sTotal_sampling_area","sgrav_tot2","sOcean_prod", "sClimate_stress","sRegional_population_growth","sReef_fish_landings_per_km2","sLarger_pop_size","sHDI",  "crest","backreef","flat","distancesampling","pointcount","shallow","deep","restricted","highcompliance", "indexrc","indexj","sMPAage","sNTZarea","PA_tp","tTargBiomass","Targeted.biomass" )]
ordereddata=rbind(reserve_data2,fished_data2)

#map of our sites
newmap <- getMap(resolution = "high")
#jitter points 
ordereddata$Site_Lat2<- ordereddata$Site_Lat+runif(length(ordereddata$Site_Lat), min=0, max=3)
ordereddata$Site_Long2<- ordereddata$Site_Long+runif(length(ordereddata$Site_Long), min=0, max=4)
ordereddata$Site_Lat2<- ifelse(ordereddata$Site_Lat2>23.5, ordereddata$Site_Lat,ordereddata$Site_Lat2)

#map of defined protection
ggplot()+geom_polygon(data = newmap, aes(x=long, y = lat, group = group),fill = "black", color = "black", alpha=0.7)+
  coord_fixed(xlim = c(-180, 180),  ylim = c(35, -35), ratio = 1.3)+
  geom_point(data=ordereddata,aes(x=Site_Long2, y=Site_Lat2, fill = as.factor(ordereddata$Protection)),colour="black", pch=21,size=3)+
  scale_fill_manual (name="Protection",values=c( "Fished"="darkorchid1","Restricted"="goldenrod","UnfishedLow"="red", "UnfishedHigh"="turquoise1"))+geom_hline(yintercept =23.43695, lty=2)+
  geom_point(data=ordereddata[ordereddata$Protection=="UnfishedHigh",],aes(x=Site_Long2, y=Site_Lat2),fill="turquoise1",colour="black", pch=21,size=3, alpha=0.7)+
  geom_point(data=ordereddata[ordereddata$Protection=="UnfishedLow",],aes(x=Site_Long2, y=Site_Lat2),fill="red",colour="black", pch=21,size=3, alpha=0.7)+
  geom_hline(yintercept =-23.43695, lty=2)+
  xlab("")+
  ylab("")+theme( panel.background = element_rect(fill="lightcyan"),panel.border = element_rect(colour = "black", fill=NA, size=1),
  axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


#extract posterior simulations for biomass
#fished
simulated_tb_fished=tb_posterior$b_fished
ordereddata$simulated_tb_fished=exp(matrixStats::colMedians(simulated_tb_fished))-min(ordereddata$Targeted.biomass[ordereddata$Targeted.biomass>0])
#low compliance
simulated_tb_lowcompliance=tb_posterior$b_lowcompliance
ordereddata$simulated_tb_lowcompliance=exp(matrixStats::colMedians(simulated_tb_lowcompliance))-min(ordereddata$Targeted.biomass[ordereddata$Targeted.biomass>0])
#high compliance
simulated_tb_highcompliance=tb_posterior$b_highcompliance
ordereddata$simulated_tb_highcompliance=exp(matrixStats::colMedians(simulated_tb_highcompliance))-min(ordereddata$Targeted.biomass[ordereddata$Targeted.biomass>0])

#change in arithmetic units
ordereddata$diff_tb_fishedlow=(ordereddata$simulated_tb_lowcompliance-ordereddata$simulated_tb_fished)
ordereddata$diff_tb_fishedhigh=(ordereddata$simulated_tb_highcompliance-ordereddata$simulated_tb_fished)
ordereddata$diff_tb_lowhigh=(ordereddata$simulated_tb_highcompliance-ordereddata$simulated_tb_lowcompliance)

#percent increase
ordereddata$percinc_tb_fishedlow=((ordereddata$simulated_tb_lowcompliance-ordereddata$simulated_tb_fished)/ordereddata$simulated_tb_fished)*100
ordereddata$percinc_tb_fishedhigh=((ordereddata$simulated_tb_highcompliance-ordereddata$simulated_tb_fished)/ordereddata$simulated_tb_fished)*100
ordereddata$percinc_tb_lowhigh=((ordereddata$simulated_tb_highcompliance-ordereddata$simulated_tb_lowcompliance)/ordereddata$simulated_tb_lowcompliance)*100


#extract posterior simulations for top predators
#fished
simulated_tp_fished=tp_posterior$p_fished
simulated_logoddstp_fished=tp_posterior$mu_fished
ordereddata$simulated_probtp_fished=exp(matrixStats::colMedians(simulated_logoddstp_fished))/(1+exp(matrixStats::colMedians(simulated_logoddstp_fished)))
ordereddata$simulated_PAtp_fished=summary(matrixStats::colMedians(simulated_tp_fished))
#low compliance
simulated_tp_lowcompliance=tp_posterior$p_lowcompliance
simulated_logoddstp_lowcompliance=tp_posterior$mu_lowcompliance
ordereddata$simulated_probtp_lowcompliance=exp(matrixStats::colMedians(simulated_logoddstp_lowcompliance))/(1+exp(matrixStats::colMedians(simulated_logoddstp_lowcompliance)))
ordereddata$simulated_PAtp_lowcompliance=summary(matrixStats::colMedians(simulated_tp_lowcompliance))
#high compliance
simulated_tp_highcompliance=tp_posterior$p_highcompliance
simulated_logoddstp_highcompliance=tp_posterior$mu_highcompliance
ordereddata$simulated_probtp_highcompliance=exp(matrixStats::colMedians(simulated_logoddstp_highcompliance))/(1+exp(matrixStats::colMedians(simulated_logoddstp_highcompliance)))
ordereddata$simulated_PAtp_highcompliance=summary(matrixStats::colMedians(simulated_tp_highcompliance))

#change un arithmetic units
ordereddata$diff_probtp_fishedlow=(ordereddata$simulated_probtp_lowcompliance-ordereddata$simulated_probtp_fished)
ordereddata$diff_probtp_fishedhigh=(ordereddata$simulated_probtp_highcompliance-ordereddata$simulated_probtp_fished)
ordereddata$diff_probtp_lowhigh=(ordereddata$simulated_probtp_highcompliance-ordereddata$simulated_probtp_lowcompliance)

#percent increase
ordereddata$percinc_probtp_fishedlow=((ordereddata$simulated_probtp_lowcompliance-ordereddata$simulated_probtp_fished)/ordereddata$simulated_probtp_fished)*100
ordereddata$percinc_probtp_fishedhigh=((ordereddata$simulated_probtp_highcompliance-ordereddata$simulated_probtp_fished)/ordereddata$simulated_probtp_fished)*100
ordereddata$percinc_probtp_lowhigh=((ordereddata$simulated_probtp_highcompliance-ordereddata$simulated_probtp_lowcompliance)/ordereddata$simulated_probtp_lowcompliance)*100



###Create manuscript figures##############################################################################################################

#for low compliance sites
a=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="UnfishedLow",], fill="black", col="black",aes(x=diff_tb_lowhigh),alpha=0.8)+xlab("Change in targeted biomass (kg/ha)")+
 ylab("")
b=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="UnfishedLow",], fill="black", lty=2,col="black",aes(x=percinc_tb_lowhigh),alpha=0.3)+xlab("Percent increase in targeted biomass")+
 ylab("")
c=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="UnfishedLow",], fill="deeppink3", col="black",aes(x=diff_probtp_lowhigh),alpha=0.8)+xlab("Change in probability of observing top predators")+
 ylab("")
d=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="UnfishedLow",], fill="deeppink3",lty=2, col="black",aes(x=percinc_probtp_lowhigh),alpha=0.3)+xlab("Percent increase in probability of observing top predators")+
 ylab("")

windows()
finalfig=ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d"))
annotate_figure(finalfig,left="Posterior density")

median(ordereddata$percinc_tb_lowhigh[ordereddata$Protection=="UnfishedLow"])
median(ordereddata$diff_tb_lowhigh[ordereddata$Protection=="UnfishedLow"])
median(ordereddata$percinc_probtp_lowhigh[ordereddata$Protection=="UnfishedLow"])
median(ordereddata$diff_probtp_lowhigh[ordereddata$Protection=="UnfishedLow"])


#for fished sites
a=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="Fished",], fill="black", col="black",aes(x=diff_tb_lowhigh),alpha=0.8)+xlab("Change in targeted biomass (kg/ha)")+
 ylab("")
b=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="Fished",], fill="black", lty=2,col="black",aes(x=percinc_tb_lowhigh),alpha=0.3)+xlab("Percent increase in targeted biomass")+
 ylab("")
c=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="Fished",], fill="deeppink3", col="black",aes(x=diff_probtp_lowhigh),alpha=0.8)+xlab("Change in probability of observing top predators")+
 ylab("")
d=ggplot(NULL)+geom_density(data=ordereddata[ordereddata$Protection=="Fished",], fill="deeppink3",lty=2, col="black",aes(x=percinc_probtp_lowhigh),alpha=0.3)+xlab("Percent increase in probability of observing top predators")+
 ylab("")

finalfig2=ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d"))
annotate_figure(finalfig2,left="Posterior density")

median(ordereddata$percinc_tb_lowhigh[ordereddata$Protection=="Fished"])
median(ordereddata$diff_tb_lowhigh[ordereddata$Protection=="Fished"])
median(ordereddata$percinc_probtp_lowhigh[ordereddata$Protection=="Fished"])
median(ordereddata$diff_probtp_lowhigh[ordereddata$Protection=="Fished"])


#for all sites
a=ggplot(NULL)+geom_density(data=ordereddata, fill="black", col="black",aes(x=diff_tb_lowhigh),alpha=0.8)+xlab("Change in targeted biomass (kg/ha)")+
 ylab("")
b=ggplot(NULL)+geom_density(data=ordereddata, fill="black", lty=2,col="black",aes(x=percinc_tb_lowhigh),alpha=0.3)+xlab("Percent increase in targeted biomass")+
 ylab("")
c=ggplot(NULL)+geom_density(data=ordereddata, fill="deeppink3", col="black",aes(x=diff_probtp_lowhigh),alpha=0.8)+xlab("Change in probability of observing top predators")+
 ylab("")
d=ggplot(NULL)+geom_density(data=ordereddata, fill="deeppink3",lty=2, col="black",aes(x=percinc_probtp_lowhigh),alpha=0.3)+xlab("Percent increase in probability of observing top predators")+
 ylab("")


finalfig=ggarrange(a,b,c,d,nrow=2,ncol=2,labels=c("a","b","c","d"))
annotate_figure(finalfig,left="Posterior density")


#save.image(file='Compliancegap_submitted.RData')
#load('Compliancegap_submitted.RData')
