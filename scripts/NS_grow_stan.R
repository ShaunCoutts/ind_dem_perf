
## Non-spatial growth model, where year, and potentially size can affect growth rate

require(brms)
library(rstan)
library(ggplot2)

# local file locations
dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
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
gr_dat = vr_loc_gen[!is.na(vr_loc_gen$height), c('uID', 'uLoc', 'rep', 'year', 'height', 'height_prev', 'X', 'Y')]
gr_dat = gr_dat[!is.na(gr_dat$height_prev), ]


# prepare the data for stan
dd = data.frame(y = gr_dat$height,
  year = gr_dat$year,
  height_prev = gr_dat$height_prev)

# ggplot(dd, aes(height_prev, y, group = year)) + geom_point() + facet_wrap(~year)
  
stan_mod = "
// We recommend generating the data with the 'make_standata' function. 
functions {
} 
data { 
  int<lower=1> N;  // total number of observations
  real height[N];  // response variable 
  vector[N] height_prev;  
  // data for year effect on growth rate
  int<lower=1> year[N];
  int<lower=1> N_years;  
} 
transformed data {
  //might log things in here
} 
parameters { 
  real<lower=0> sigma; // standard deviaiton on height
  vector[N_years] b0; // year specifc intercept
  vector[N_years] gr_rate; // year specifc slope
} 
transformed parameters { 
  //linear predcitor
  vector[N] mu; 
  
  //linear predictor fixed effects
  for (n in 1:N) { 
    mu[n] = b0[year[n]] + height_prev[n] * gr_rate[year[n]];
  } 
} 
model {
  // prior specifications 
  sigma ~ cauchy(0, 20); 
  
  // likelihood contribution 
  height ~ normal(mu, sigma); 
} 
generated quantities { 
}"         

input = make_standata(formula = y ~ 0 + height_prev + (1|year) ,
  data=dd)

stan_dat = list(N = input$N, height = input$Y, height_prev = dd$height_prev,  
  N_years = input$N_1, year = input$J_1)

NS_gr_stan <- stan(model_code = stan_mod, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.90), pars = c('b0', 'gr_rate', 'sigma', 'mu')) 

setwd(out_ob_loc)
save(NS_gr_stan , file = 'NS_gr_stan.Rdata')
#load('NS_gr_stan.Rdata')
