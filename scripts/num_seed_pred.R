# seed per fruit prediction 
# get data (in spss)

library(foreign)
library(dplyr)
library(tidyr)
library(ggplot2)
require(brms)
library(rstan)

# local file locations
dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

setwd(dat_loc)
seeds98 = read.spss('seeds98.sav', to.data.frame = TRUE) 
seeds96 = read.spss('SEED96.SAV', to.data.frame = TRUE) 
seeds95 = read.spss('SEED95.SAV', to.data.frame = TRUE) 
burn_dat = read.csv('hcdem_ABS_2015_checked.csv')

## labels
## TSF = 
## SITE = study site (mine is site 103)
## PLANT = plant ID
## TOTAL = total number of seeds 
## DEVELOPE = It has developed but may not be good) 
## ABORTED = seed tiny, white, not viable
## BROWN = good seed (this is the one I want)
## NOZERO = just markes the 0's

# can only use the 95 data as that is the only one with height, seems like there are
# multipule sites and often two measures per plant, will need some sort of heiracical model maybe, 
# or just average over plants, since only two obervations per plant just average is going to work better
# take the coloumns I need
seed95_useful = select(seeds95, SITE, PLANT, HEIGHT, BROWN)
seed95_useful = mutate(seed95_useful, ID = paste0(PLANT, ':', SITE))
plant_group = group_by(seed95_useful, ID)
seed95_site = summarise(plant_group, n_obs = n(), site = unique(SITE), height = mean(HEIGHT),
  n_seed = mean(BROWN))
# ggplot(seed95_site, aes(height, n_seed)) + geom_point() + facet_grid(.~ site)
  

stan_model = "
// We recommend generating the data with the 'make_standata' function. 
functions {
} 
data { 
  int<lower=1> N;  // total number of observations
  vector[N] seed_num;  // response variable 
  vector[N] height;  // height 
  // data for group-level effects of ID 1 
  int<lower=1> site[N]; 
  int<lower=1> N_sites;  
} 
transformed data {
  vector[N] ln_seed;
  ln_seed = log(seed_num + 1);
} 
parameters {
  real glob_int;
  real h_ef;  // height effects
  real<lower=0> sigma; //global sd
  real<lower=0> site_sd; // group-level standard deviations 
  vector[N_sites] site_raw;  // unscaled group-level effects
} 
transformed parameters { 
  //site effect
  vector[N_sites] site_int;
  //linear predcitor
  vector[N] mu; 
  
  //year effect
  site_int = site_sd * site_raw;
  
  //linear predictor fixed effects 
  mu = glob_int + height * h_ef; 
  // add both random effect of year and space
  for (n in 1:N) { 
    mu[n] = mu[n] + site_int[site[n]]; 
  } 
} 
model{
  // prior specifications 
  site_raw ~ normal(0, 1);
  
  // likelihood contribution 
  ln_seed ~ normal(mu, sigma); 
} 
generated quantities { 
}"              

input = make_standata(formula = n_seed ~ height + (1|site),
  data = seed95_site, family = "normal")

stan_dat = list(N = input$N, seed_num = input$Y, height = input$X[, 2], N_sites = input$N_1, 
  site = input$J_1)

seed_num_stan <- stan(model_code = stan_model, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 3, chains = 3, control = list(adapt_delta = 0.90), 
  pars = c('glob_int', 'h_ef', 'site_int', 'site_sd', 'sigma', 'mu')) 

setwd(out_ob_loc)
save(seed_num_stan , file = 'seed_num_stan.Rdata')
#load('seed_num_stan.Rdata')

## try a model with better data to estimate number of fruits (assume seed production is proportional to fruit production.
# first get a clean set of data and put into long form, only take years from 1994-2003 so that most indivudals have repeated 
# meaures for them so individual level effects can be fitted 
fruit_ht_dat = select(burn_dat, year, burn_yr, site, gap, tag, ht94, ht95, ht96, ht97, ht98, ht99, ht00,
  ht01, ht02, ht03)
fruit_rep_dat = select(burn_dat, year, site, gap, tag, rep94, rep95, rep96, rep97, rep98, rep99, rep00, rep01, rep02,
  rep03)
  
ht_long = gather(fruit_ht_dat, ht_lab, height, ht94:ht03)
rep_long = gather(fruit_rep_dat, rep_lab, fruit_num, rep94:rep03)

# drop NA's
ht_long = ht_long[ht_long$height != '#NULL!', ]
rep_long = rep_long[rep_long$fruit_num != '#NULL!', ]

# make the year labels numeric 
ht_long$m_year = sapply(ht_long$ht_lab, FUN = function(x){
  a = as.numeric(strsplit(x, split = 'ht')[[1]][2])
  if(a >= 50) return(1900 + a) else return(2000 + a)
})

rep_long$m_year = sapply(rep_long$rep_lab, FUN = function(x){
  a = as.numeric(strsplit(x, split = 'rep')[[1]][2])
  if(a >= 50) return(1900 + a) else return(2000 + a)
})

# add some columns for ID 
ht_long = mutate(ht_long, ID = paste0(site, ':', gap, ':', tag),
  join_ID = paste0(m_year, ':', site, ':', gap, ':', tag), 
  height = as.numeric(height))
  
rep_long = mutate(rep_long, ID = paste0(site, ':', gap, ':', tag),
  join_ID = paste0(m_year, ':', site, ':', gap, ':', tag),
  fruit_num = as.numeric(fruit_num))

# merg the data frames 
fruit_dat = inner_join(ht_long, rep_long, by = 'join_ID')

fruit_dat = select(fruit_dat, burn_yr = burn_yr, site = site.x, ID = ID.x,
  join_ID = join_ID, m_year = m_year.x, height = height, fruit_num = fruit_num)  
# only take the rows with fruit number greater than 1, as these are counts of fruits 
fruit_dat = filter(fruit_dat, fruit_num >= 2)
# get time since fire for each observation
fruit_dat = mutate(fruit_dat, time_since_fire = m_year - burn_yr, site = as.character(site))

# explore the data a bit
plot(fruit_dat[, c('m_year', 'height', 'time_since_fire', 'fruit_num')])

site_group = group_by(fruit_dat, ID)
fruit_ID_sum = summarise(site_group, n_obs = n())
hist(fruit_ID_sum$n_obs)

# the majority of my indivudals only have one obervation, making it 
# very hard to fit an indivdual random effect, next most common number of 
# observations is 2. So need to have a site level random effect only
# and realise that not all my obervations are independent.	
# I also want a differnt slope estimated for time since fire for each site 
# since times since fire can vary widely between sites, so some of those lines 
# might be very flat


fruit_model = "
// We recommend generating the data with the 'make_standata' function. 
functions {
} 
data { 
  int<lower=1> N;  // total number of observations
  vector[N] fruit_num;  // response variable 
  vector[N] height;  // height 
  vector[N] TSF; // time since fire
  
  // data for group-level effects of site 
  int<lower=1> site[N]; 
  int<lower=1> N_sites;  
} 
transformed data {
  vector[N] ln_fruit;
  ln_fruit = log(fruit_num);
} 
parameters {
  real glob_int;
  real h_ef;  // height effects
  real glob_TSF_ef; // global effect of time since fire
  
  real<lower=0> sigma; //global sd
  real<lower=0> site_int_sd; // group-level standard deviation on intercept 
  real<lower=0> site_TSF_sd; //group-level standard deviation on slope
  vector[N_sites] site_int_raw;  // unscaled group-level intercept
  vector[N_sites] site_TSF_raw; // unscaled group level slope on TSF
  
} 
transformed parameters { 
  //site effect
  vector[N_sites] site_int;
  vector[N_sites] site_TSF;
  //linear predcitor
  vector[N] mu; 
  
  //site effect
  site_int = site_int_sd * site_int_raw;
  site_TSF = site_TSF_sd * site_TSF_raw + glob_TSF_ef;
  
  //linear predictor fixed effects 
  mu = glob_int + height * h_ef; 
  
  // add random effect of 1+TSF|site
  for (n in 1:N) { 
    mu[n] = mu[n] + site_int[site[n]] + TSF[n] * site_TSF[site[n]]; 
  } 
} 
model{
  // prior specifications 
  site_int_raw ~ normal(0, 1);
  site_TSF_raw ~ normal(0, 1);
  
  // likelihood contribution 
  ln_fruit ~ normal(mu, sigma); 
} 
generated quantities { 
}"              

input = make_standata(formula = fruit_num ~ height + time_since_fire + (time_since_fire|site),
  data = fruit_dat, family = "normal")

stan_dat = list(N = input$N, fruit_num = input$Y, height = input$X[, 2], TSF = input$X[, 3], 
  N_sites = input$N_1, site = input$J_1)

fruit_num_stan <- stan(model_code = fruit_model, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.90), 
  pars = c('glob_int', 'h_ef', 'glob_TSF_ef', 'site_TSF_raw', 'site_TSF_sd', 'site_int_raw', 
  'site_int_sd', 'site_int', 'site_TSF', 'sigma', 'mu')) 

traceplot(fruit_num_stan, pars = c('glob_int', 'h_ef', 'glob_TSF_ef', 'site_TSF_raw', 'site_TSF_sd', 'site_int_raw', 
  'site_int_sd', 'site_int', 'site_TSF', 'sigma'))  
  
setwd(out_ob_loc)
save(fruit_num_stan , file = 'fruit_num_stan.Rdata')
#load('fruit_num_stan.Rdata')

# fitted vs predicted plot

str(fruit_num_stan)
lp_means = apply(extract(fruit_num_stan, pars = 'mu')$mu, MARGIN = 2, FUN = mean)
plot(log(fruit_dat$fruit_num), lp_means)
abline(0, 1)

# try and model time since fire in a different way, maybe non-linear effect
# so that I can try and fit a single relationship across sites for the effect of time since fire
# first step is to plot TSF against ln(fruit numer) for each site
ggplot(fruit_dat, aes(time_since_fire, log(fruit_num))) + geom_point() + facet_grid(site ~.)
ggplot(fruit_dat, aes(jitter(time_since_fire), log(fruit_num))) + geom_point() 
# seems like a very weak realtionship if it exists 

ggplot(fruit_dat, aes(height, log(fruit_num))) + geom_point() + facet_grid(site ~.)
# strong relationship, try a model with just height 
######################################################################################################################################
fruit_model = "
// We recommend generating the data with the 'make_standata' function. 
functions {
} 
data { 
  int<lower=1> N;  // total number of observations
  vector[N] fruit_num;  // response variable 
  vector[N] height;  // height 
  
  // data for group-level effects of site 
  int<lower=1> site[N]; 
  int<lower=1> N_sites;  
} 
transformed data {
  vector[N] ln_fruit;
  vector[N] height_cen;
  
  ln_fruit = log(fruit_num);
  height_cen = height - mean(height);
} 
parameters {
  real h_ef;  // height effects
  real<lower=0> sigma; //global sd
  vector[N_sites] site_int;  // unscaled group-level intercept
  
} 
transformed parameters { 
  //linear predcitor
  vector[N] mu; 
  
  for (n in 1:N) { 
    mu[n] = site_int[site[n]] + height_cen[n] * h_ef; 
  } 
} 
model{
  
  // likelihood contribution 
  ln_fruit ~ normal(mu, sigma); 
} 
generated quantities { 
}"              

input = make_standata(formula = fruit_num ~ height + (1|site),
  data = fruit_dat, family = "normal")

stan_dat = list(N = input$N, fruit_num = input$Y, height = input$X[, 2], 
  N_sites = input$N_1, site = input$J_1)

fruit_num_stan <- stan(model_code = fruit_model, data = stan_dat, iter = 3000, warmup = 1000, 
  cores = 4, chains = 4, control = list(adapt_delta = 0.90), 
  pars = c('h_ef', 'site_int', 'sigma', 'mu')) 

traceplot(fruit_num_stan, pars = c('h_ef', 'site_int', 'sigma'))  
pairs(fruit_num_stan, pars = c('h_ef', 'site_int', 'sigma'))  
  
setwd(out_ob_loc)
save(fruit_num_stan , file = 'fruit_num_stan.Rdata')
#load('fruit_num_stan.Rdata')

# fitted vs predicted plot

str(fruit_num_stan)
lp_means = apply(extract(fruit_num_stan, pars = 'mu')$mu, MARGIN = 2, FUN = mean)
plot(log(fruit_dat$fruit_num), lp_means)
abline(0, 1)
# residual plot
plot(lp_means, lp_means - log(fruit_dat$fruit_num))

# model over estimates fruit production of very small plants, but there are very few of these plants in the model
# but there are very few of these small plants in the data, since mostly plants that small will not produce fruit
# and residual plot looks really good