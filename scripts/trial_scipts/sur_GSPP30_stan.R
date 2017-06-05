## Fitting non-spatial model in stan
require(brms)
library(rstan)

# local file locations
dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

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

plot(sur_dat$X, sur_dat$Y, pch = 19, cex = 0.25, col = 'grey')
points(knots, pch = 3, cex= 0.25, col = 'red')

dm_knots = distance_matrix(knots[,1], knots[,2])
dm_knots_obs = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knots[, 1], 
  y2 = knots[, 2])

# prepare the data for stan
dd = data.frame(y = sur_dat$sur,
  year = sur_dat$year,
  height = sur_dat$height_mod)


stan_model2 = "// This Stan code was generated with the R package 'brms'. 
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
  int<lower=1> N_knots;
  matrix[N_knots, N_knots] knot_dist; //distance matrix between N_knots
  matrix[N, N_knots] ob_knot_dist; // distance between each data point and each knot
  int sur[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data for group-level effects of ID 1 
  int<lower=1> year[N]; 
  int<lower=1> N_years;  
  int<lower=1> M_1;
} 
transformed data { 
  vector[N_knots] mu;
  for(k in 1:N_knots){
    mu[k] = 0; //set mean of SPP to 0
  }
} 
parameters { 
  vector[K] b;  // population-level effects 
  vector<lower=0>[M_1] sd_1; // group-level standard deviations 
  vector[N_years] z_1[M_1];  // unscaled group-level effects
  
  //Predictive Proccess parameters 
  real<lower=0.00001> sigma_spp; // light trunction on these to make sure that C_spp > 0
  real<lower=0.1, upper=1000> dd_spp_inv; // also truncate on lower bound as past a certian point making this larger does not change the function over the spatial domain 
  vector[N_knots] spp;
  
} 
transformed parameters { 
  matrix[N_knots, N_knots] C_spp;
  real<lower=0> dd_spp;
  matrix[N_knots, N_knots] C_spp_inv;
  matrix[N, N_knots] C_spp_s;
  vector[N] SPP_int;
  vector[N_years] year_int; 
  vector[N] eta;
  
  // year effect
  year_int = sd_1[1] * (z_1[1]);
  
  //spatial predictive process
  dd_spp = inv(dd_spp_inv);
  C_spp = cov_var_maker(knot_dist, N_knots, dd_spp, sigma_spp);
  C_spp_inv = inverse_spd(C_spp);

  // Interpolate back from the knot points to the observed data 
  for(n in 1:N){
    for(k in 1:N_knots){
      C_spp_s[n, k] = sigma_spp * exp(-(dd_spp * ob_knot_dist[n, k]));
    }
  }
  SPP_int = C_spp_s * C_spp_inv * spp;  
  
  // linear predictor
  eta = X * b; 
  for (n in 1:N) { 
    eta[n] = eta[n] + year_int[year[n]] + SPP_int[n]; 
  } 
  
} 
model {
  // prior specifications 
  sigma_spp ~ cauchy(0, 5);
  dd_spp_inv ~ cauchy(0, 5);
  sd_1 ~ student_t(3, 0, 10); 
  z_1[1] ~ normal(0, 1); 
  
  spp ~ multi_normal_prec(mu, C_spp_inv); 
  
  // likelihood contribution 
  sur ~ bernoulli_logit(eta); 

} 
generated quantities { 
}"              

input = make_standata(formula = y ~ 0 + intercept + height + (1|year) ,
  data=dd, family="bernoulli")
input$Z_1_1 = NULL

stan_dat = list(N = length(sur_dat$sur), sur = sur_dat$sur, K = input$K, X = input$X, 
  N_years = input$N_1, year = input$J_1, N_knots = dim(dm_knots)[1], knot_dist = dm_knots, 
  ob_knot_dist = dm_knots_obs, M_1 = input$M_1)

SPP_sur_stan <- stan(model_code = stan_model2, data = stan_dat, iter = 100, warmup = 50, 
  cores = 1, chains = 1, control = list(adapt_delta = 0.8), include = TRUE, 
  pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv', 'sigma_spp', 'spp')) 

setwd(out_ob_loc)
save(SPP_sur_stan, file = 'SPP_sur_stan.Rdata')

stan_diag(SPP_sur_stan)

stan_trace(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv', 'sigma_spp', 'spp'))
pairs(SPP_sur_stan)


ext = extract(SPP_sur_stan)
is.positive.definite(ext$C_spp_inv[1, , ])
