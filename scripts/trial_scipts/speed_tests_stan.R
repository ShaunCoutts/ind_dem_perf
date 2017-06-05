## Stan speed test to try and remove a few of the bottle necks in the code with 
## some toy examples 

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
  y_max = max(sur_dat$Y), res = 25, dat_x = sur_dat$X, dat_y = sur_dat$Y, iso_dist = 700)

dm_knots = distance_matrix(knots[,'x_cord'], knots[, 'y_cord'])
dm_all = distance_matrix(jitter(sur_dat$X, amount = 0.05), jitter(sur_dat$Y, amount = 0.05))


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
  matrix[N, N] dist_mat; //distance matrix between N_knots
} 
transformed data { 
} 
parameters { 
  //Predictive Proccess parameters 
  real<lower=0.00001> sigma_spp;
  real<lower=00.1, upper=1000> dd_spp_inv;
  vector[N] spp;
  
} 
transformed parameters { 
  // spatial error
  matrix[N, N] C;
  vector[N] mu;
  real<lower=0> dd_spp;
  //vector[N] spp;
  
  //spatial error 
  dd_spp = inv(dd_spp_inv);
  C = cov_var_maker(dist_mat, N, dd_spp, sigma_spp);
  
  //spp = cholesky_decompose(C) * spp_raw; 
  
  for(k in 1:N){
    mu[k] = 0; //set mean of SPP to 0
  }
  
} 
model {
  // prior specifications 
  sigma_spp ~ cauchy(0, 5);
  dd_spp_inv ~ cauchy(0, 5);
  
  //estimated spatial effect 
  spp ~ multi_normal(mu, C); 
  //spp_raw ~ normal(0, 1);
} 
generated quantities { 
}"              

stan_dat = list(N = dim(dm_knots)[1], dist_mat = dm_knots)

SPP_sur_stan <- stan(model_code = stan_model2, data = stan_dat, iter = 200, warmup = 100, cores = 1, chains = 1, 
  control = list(adapt_delta = 0.8), pars = c('sigma_spp', 'dd_spp_inv')) 

  
stan_trace(SPP_sur_stan)
 

