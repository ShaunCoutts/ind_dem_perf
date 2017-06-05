
## Nearest Knot version of the SPP model. With  a careful knot placement stratergy
## This means I do not have to the back interpolation step. This is computationaly
## expensive and eats RAM like fat man at a free buffet. This is pretty close to
## the vinailla version of a spatial errors model, but is much faster to run
## since the distance matrix is smaller and easier to invert or decompose, and 
## faster to build.

require(brms)
library(rstan)

# local file locations
dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
#out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

# file locations for Zhu's Mac
#dat_loc = '/Users/shauncoutts/ind_perf_pred/Data' 
#out_ob_loc = '/Users/shauncoutts/ind_perf_pred/objects'
#code_loc = '/Users/shauncoutts/ind_perf_pred/scripts'

# ICEBERG run
#dat_loc = '/home/bo1src/ind_dem_perf/Data' 
#out_ob_loc = '/home/bo1src/ind_dem_perf/output/objects'
#code_loc = '/home/bo1src/ind_dem_perf/scripts'

# read in the jags model and helper functions 
setwd(code_loc)
source('dist_neigh_setup.R')

# get data
setwd(dat_loc)
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactor = FALSE)
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

# start with survival
sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]
# neighbour data we can use all data, even the first years observed height 
neigh_dat = sur_dat

# take out data from the first year observed since nothing can be observed dead in the first year
first_year = sapply(seq_along(sur_dat$year), FUN = function(x){

  min_year_group = min(sur_dat$year[sur_dat$uID == sur_dat$uID[x]])
  return(ifelse(sur_dat$year[x] == min_year_group, FALSE, TRUE))

})
sur_dat = sur_dat[first_year, ]

# recode year and patch to act as indicies 
sur_dat$year_num = (sur_dat$year - (max(sur_dat$year))) + max(abs(sur_dat$year - (max(sur_dat$year)))) + 1    
sur_dat$gapID = sapply(strsplit(sur_dat$uID, split = ':', fixed = TRUE), FUN = function(x) x[1])
inds = unique(sur_dat$gapID)
sur_dat$gapID_num = NA
for(i in 1:length(inds)) sur_dat$gapID_num[sur_dat$gapID %in% inds[i]] = i

# height: this requires some assumptions because 1) individuals change size over the course of a year, 
# 2) we do not have the final size for indivudals that died. 3) we have indivduals that were only found when they were very large.
# Best we can do here is assue that mortality occurs just after census. 
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]
sur_dat$height_mod = sur_dat$height_prev - mean(sur_dat$height_prev)

# prepare the data for stan
dd = data.frame(y = sur_dat$sur,
  year = sur_dat$year,
  height = sur_dat$height_mod)

stan_model = "
// We recommend generating the data with the 'make_standata' function. 
functions {
} 
data { 
  int<lower=1> N;  // total number of observations
  int sur[N];  // response variable 
  vector[N] height_cen;  // population-level design matrix 
  // data for group-level effects of ID 1 
  int<lower=1> year[N]; 
  int<lower=1> N_years;  
  int<lower=1> M_1;
} 
transformed data { 
} 
parameters { 
  real h_ef;  // height effects 
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
  eta = h_ef * height_cen;
  // add both random effect of year and space
  for (n in 1:N) { 
    eta[n] = eta[n] + year_int[year[n]];
  } 
} 
model {
  // prior specifications 
  z_1[1] ~ normal(0, 1);
  
  // likelihood contribution 
  sur ~ bernoulli_logit(eta); 
} 
generated quantities { 
}"              

input = make_standata(formula = y ~ 0 + height + (1|year) ,
  data=dd, family="bernoulli")

stan_dat = list(N = input$N, sur = input$Y, height_cen = dd$height, N_years = input$N_1, year = input$J_1, 
  M_1 = input$M_1)

NS_sur_stan <- stan(model_code = stan_model, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.8), 
  pars = c('h_ef', 'year_int', 'eta')) 

setwd(out_ob_loc)
save(NS_sur_stan , file = 'NS_sur_stan.Rdata')
#load('NS_sur_stan.Rdata')
