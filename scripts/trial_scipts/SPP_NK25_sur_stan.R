## Fitting non-spatial model in stan. This includes both spatial scales (>25m and fine scale).
## Both seem to pick up some pattern so see if I can use both in the same model. But it will 
## be a minor miricle if this model is identifiable.
require(brms)
library(rstan)

# local file locations
dat_loc = '/home/bo1src/ind_dem_perf/Data' 
out_ob_loc = '/home/bo1src/ind_dem_perf/output/objects'
code_loc = '/home/bo1src/ind_dem_perf/scripts'

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
# Best we can do here is assue that mortality occurs just after census, so that we can use previous height as a predictor,
# and in the case of indivduals in their first year we say they were 0 
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]
sur_dat$height_mod = sur_dat$height_prev - mean(sur_dat$height_prev)

# create knots at 10m grid, cut out all those more than 20m from nearerst data point
knots = knot_coords(x_min = min(sur_dat$X), x_max = max(sur_dat$X), y_min = min(sur_dat$Y), 
 y_max = max(sur_dat$Y), res = 25, dat_x = sur_dat$X, dat_y = sur_dat$Y, iso_dist = 30)

dm_knots = distance_matrix(knots[,1], knots[,2])
dm_knots_obs = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knots[, 1], 
  y2 = knots[, 2])

# create knots at each data point, cluster some less than 50cm apart
# results in 575 knots. 50% of knots have another knot within 61cm and 75% have
# another knot within 0.82cm. 50% of knots have 8 or more neigbour knots within 2m
# and 75% have 5 or more neighbouring knots within 2m.
knot_cl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5) 

dm_knots_cl = distance_matrix(knot_cl[,1], knot_cl[,2])
dm_knots_obs_cl = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knot_cl[, 1], 
  y2 = knot_cl[, 2])
# get index of nearerst knot for each data point
near_knot = apply(dm_knots_obs_cl, MARGIN = 1, FUN = function(x) which(x == min(x)))

# prepare the data for stan
dd = data.frame(y = sur_dat$sur,
  year = sur_dat$year,
  height = sur_dat$height_mod)

stan_model = "// This Stan code was generated with the R package 'brms'. 
// We recommend generating the data with the 'make_standata' function. 
functions {
 
  //constructs a correlation matrix
  matrix cov_var_maker(matrix dist_mat, int mat_size, real dd_spp, real sigma_spp){
  
    matrix[mat_size, mat_size] C;
    
    for(k in 1:(mat_size - 1)){
      for(j in (k + 1):mat_size){
	
	C[k, j] = sigma_spp * exp(-(dd_spp * dist_mat[k, j])); 
	C[j, k] = C[k, j];
	
      }
    }
    //fill diagonals
    for(k in 1:mat_size){
      C[k, k] = sigma_spp;
    }
    
    return C;
    
  }

} 
data { 
  int<lower=1> N;  // total number of observations 
  int<lower=1> N_knots25;
  int<lower=1> N_knotsNN;
  matrix[N_knots25, N_knots25] knot_dist25; //distance matrix between N_knots course grid
  matrix[N_knotsNN, N_knotsNN] knot_distNN; //distance matrix between N_knots fine grid
  matrix[N, N_knots25] ob_knot25_dist; // distance between each data point and each knot on course grid
  int<lower=1, upper=N_knotsNN> nn_knot_ind[N]; // index of nearerst knot on fine grid
  
  int sur[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data for group-level effects of ID 1 
  int<lower=1> year[N]; 
  int<lower=1> N_years;  
  int<lower=1> M_1;
} 
transformed data { 
  vector[N_knots25] mu25;
  
  for(k in 1:N_knots25){
    mu25[k] = 0; //set mean of SPP course to 0
  }
} 
parameters { 
  vector[K] b;  // population-level effects 
  vector<lower=0>[M_1] sd_1; // group-level standard deviations 
  vector[N_years] z_1[M_1];  // unscaled group-level effects
  
  //Predictive Proccess parameters 
  real<lower=0.00001> sigma_spp25; // light trunction on these to make sure that C_spp > 0
  real<lower=3, upper=5000> dd_spp_inv25; // also truncate on lower bound as past a certian point making this larger does not change the function over the spatial domain 
  vector[N_knots25] spp25;
  
  real<lower=0.00001> sigma_sppNN;
  real<lower=0.1, upper=1000> dd_spp_invNN;
  vector[N_knotsNN] spp_raw;
} 
transformed parameters { 
  matrix[N_knots25, N_knots25] C_spp25;
  matrix[N_knotsNN, N_knotsNN] C_sppNN;
  real<lower=0> dd_spp25;
  real<lower=0> dd_sppNN;
  matrix[N_knots25, N_knots25] C_spp_inv25;
  matrix[N, N_knots25] C_spp_s25;
  vector[N] spp_int25;
  vector[N_knotsNN] spp_NN;
  vector[N_years] year_int; 
  vector[N] eta;
  
  // year effect
  year_int = sd_1[1] * (z_1[1]);
  
  //spatial predictive process course scale
  dd_spp25 = inv(dd_spp_inv25);
  C_spp25 = cov_var_maker(knot_dist25, N_knots25, dd_spp25, sigma_spp25);
  C_spp_inv25 = inverse_spd(C_spp25);

  // Interpolate back from the knot points to the observed data 
  for(n in 1:N){
    for(k in 1:N_knots25){
      C_spp_s25[n, k] = sigma_spp25 * exp(-(dd_spp25 * ob_knot25_dist[n, k]));
    }
  }
  
  spp_int25 = C_spp_s25 * C_spp_inv25 * spp25;  
  
  //spatial error fine scale
  dd_sppNN = inv(dd_spp_invNN);
  C_sppNN = cov_var_maker(knot_distNN, N_knotsNN, dd_sppNN, sigma_sppNN);
  
  spp_NN = cholesky_decompose(C_sppNN) * spp_raw; 
 
  // linear predictor
  eta = X * b; 
  for (n in 1:N) { 
    eta[n] = eta[n] + year_int[year[n]] + spp_int25[n] + spp_NN[nn_knot_ind[n]]; 
  } 
  
} 
model {
  // prior specifications 
  sigma_spp25 ~ cauchy(0, 5);
  dd_spp_inv25 ~ cauchy(0, 5);
  sigma_sppNN ~ cauchy(0, 5);
  dd_spp_invNN ~ cauchy(0, 5);
  sd_1 ~ student_t(3, 0, 10); 
  z_1[1] ~ normal(0, 1); 
  
  spp25 ~ multi_normal_prec(mu25, C_spp_inv25); 
  spp_raw ~ normal(0, 1);
  
  // likelihood contribution 
  sur ~ bernoulli_logit(eta); 

} 
generated quantities { 
}"              

input = make_standata(formula = y ~ 0 + intercept + height + (1|year) ,
  data=dd, family="bernoulli")
input$Z_1_1 = NULL

stan_dat = list(N = input$N, sur = sur_dat$sur, K = input$K, X = input$X, 
  N_years = input$N_1, year = input$J_1, N_knots25 = dim(dm_knots)[1], knot_dist25 = dm_knots, 
  ob_knot25_dist = dm_knots_obs, N_knotsNN = dim(dm_knots_cl)[1], knot_distNN = dm_knots_cl, 
  nn_knot_ind = near_knot, M_1 = input$M_1)

SPP_sur_stan <- stan(model_code = stan_model, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.8), include = TRUE, 
  pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv25', 'sigma_spp25', 'spp_int25',
  'dd_spp_invNN', 'sigma_sppNN', 'spp_NN','eta')) 

setwd(out_ob_loc)
save(SPP_sur_stan, file = 'SPP_25NK_sur_stan.Rdata')
