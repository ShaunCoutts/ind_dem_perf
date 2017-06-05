
## Non-spatial version of the reporductive model

require(brms)
library(rstan)

# local file locations
dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

# file locations for Zhu's Mac
#dat_loc = '/Users/shauncoutts/ind_perf_pred/Data' 
#out_ob_loc = '/Users/shauncoutts/ind_perf_pred/objects'
#code_loc = '/Users/shauncoutts/ind_perf_pred/scripts'

# ICEBERG run
# dat_loc = '/home/bo1src/ind_dem_perf/Data' 
# out_ob_loc = '/home/bo1src/ind_dem_perf/output/objects'
# code_loc = '/home/bo1src/ind_dem_perf/scripts'

# get data
setwd(dat_loc)
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactor = FALSE)
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

# start with survival
rep_dat = vr_loc_gen[!is.na(vr_loc_gen$rep), c('uID', 'uLoc', 'rep', 'year', 'height', 'X', 'Y')]
rep_dat = rep_dat[!is.na(rep_dat$height),]
rep_dat = rep_dat[rep_dat$rep <= 1, ]

# recode year and patch to act as indicies 
rep_dat$year_num = (rep_dat$year - (max(rep_dat$year))) + max(abs(rep_dat$year - (max(rep_dat$year)))) + 1    
rep_dat$gapID = sapply(strsplit(rep_dat$uID, split = ':', fixed = TRUE), FUN = function(x) x[1])
inds = unique(rep_dat$gapID)
rep_dat$gapID_num = NA
for(i in 1:length(inds)) rep_dat$gapID_num[rep_dat$gapID %in% inds[i]] = i

rep_dat$height_mod = rep_dat$height - mean(rep_dat$height)

# prepare the data for stan
dd = data.frame(y = rep_dat$rep,
  year = rep_dat$year,
  height = rep_dat$height_mod)

stan_model = "
// We recommend generating the data with the 'make_standata' function. 
functions {
} 
data { 
  int<lower=1> N;  // total number of observations
  int rep[N];  // response variable 
  vector[N] height_cen;  // height 
  // data for group-level effects of ID 1 
  int<lower=1> year[N]; 
  int<lower=1> N_years;  
  int<lower=1> M_1;
} 
transformed data { 
} 
parameters { 
  real h_ef;  // height effects 
  //vector<lower=0>[M_1] sd_1; // group-level standard deviations 
  vector[N_years] z_1[M_1];  // unscaled group-level effects
} 
transformed parameters { 
  //year effect
  vector[N_years] year_int;
  //linear predcitor
  vector[N] eta; 
  
  //year effect
  year_int = 20.0 * (z_1[1]);
  
  //linear predictor fixed effects 
  eta = height_cen * h_ef; 
  // add both random effect of year and space
  for (n in 1:N) { 
    eta[n] = eta[n] + year_int[year[n]]; 
  } 
} 
model {
  // prior specifications 
  // some divergence problems try a slightly different parameterization 
  // draw the standard deviation of year effect from a half standard normal, sd = 20 in transformed parameters block
  //sd_1[1] ~ normal(0, 1);
  z_1[1] ~ normal(0, 1);
  
  // likelihood contribution 
  rep ~ bernoulli_logit(eta); 
} 
generated quantities { 
}"              

input = make_standata(formula = y ~ 0 + height + (1|year) ,
  data=dd, family="bernoulli")

stan_dat = list(N = input$N, rep = input$Y, height_cen = dd$height, N_years = input$N_1, 
  year = input$J_1, M_1 = input$M_1)

NS_rep_stan <- stan(model_code = stan_model, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.90), 
  pars = c('h_ef', 'year_int', 'eta')) 

setwd(out_ob_loc)
save(NS_rep_stan , file = 'NS_rep_stan.Rdata')
#load('NS_rep_stan.Rdata')
