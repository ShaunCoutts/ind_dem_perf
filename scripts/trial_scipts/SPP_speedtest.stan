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
  
  matrix ob_knot_C(matrix ob_knot_dm, real dd_spp, real sigma_spp){
  
    matrix[rows(ob_knot_dm), cols(ob_knot_dm)] C_knot_ob;
  
    C_knot_ob = sigma_spp * exp(-dd_spp * ob_knot_dm);
    
    return C_knot_ob;
  
  }
  
  vector interp_back(vector spp, matrix C_inv, matrix C_knot_ob, int N_obs){
  
    vector[N_obs] SPP_int;
    SPP_int = C_knot_ob * C_inv * spp;  
    return SPP_int; 
  
  }
  
  vector inter_weight_ave(vector spp, matrix knot_ob_dm, int N_obs, int N_nn, 
    real dd_spp){
    
    vector[N_obs] SPP_int;
    vector[N_nn] w;
    vector[N_nn] w_spp;
    
    for(i in 1:N_obs){
      for(j in 1:N_nn){
      
	w[j] = exp(-dd_spp * knot_ob_dm[i, j]);
	w_spp[j] = spp[j] * w[j];
      }
      SPP_int[i] = sum(w_spp) / sum(w);
    }
    
    return SPP_int;
  
  }

  vector inter_NK(vector spp, matrix knot_ob_dm, int[,] nn_ind, int N_nn, int N_obs, real dd_spp){
    
    vector[N_obs] SPP_int;
    vector[N_nn] w;
    vector[N_nn] w_spp;
    
    for(i in 1:N_obs){
      for(j in 1:N_nn){
      
	w[j] = exp(-dd_spp * knot_ob_dm[i, j]);
	w_spp[j] = spp[nn_ind[i, j]] * w[j];
      }
      SPP_int[i] = sum(w_spp) / sum(w);
    }
    
    return SPP_int;
  
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
  C_spp_s = ob_knot_C(ob_knot_dist, dd_spp, sigma_spp);
  SPP_int = interp_back(spp, C_spp_inv, C_spp_s, N);
  
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
}              
