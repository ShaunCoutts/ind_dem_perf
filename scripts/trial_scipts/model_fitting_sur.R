# model fitting script, data set up and call to jags
library(R2jags)
library(coda)

# local file locations
dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

# Zhu computer file locations 
dat_loc = '/Users/shauncoutts/ind_perf_pred/Data' 
out_plot_loc = '/Users/shauncoutts/ind_perf_pred/plots'
out_ob_loc = '/Users/shauncoutts/ind_perf_pred/objects'
code_loc = '/Users/shauncoutts/ind_perf_pred/scripts'

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


# read in the jags model and helper functions 
setwd(code_loc)
source('sur_model.R')
source('dist_neigh_setup.R')

jags_data <- list(sur = sur_dat$sur, year = sur_dat$year_num, height = sur_dat$height_mod,
  num_years = max(sur_dat$year_num), num_obs = length(sur_dat$sur), gapID = sur_dat$gapID_num,
  num_gaps = max(sur_dat$gapID_num))

#intial values
inital <- function() list(b_height = runif(1, -20, 20), 
  int_year = runif(jags_data$num_years, -20, 20),
  gap_int = runif(jags_data$num_gaps, -20, 20)) 

# Parameters monitored
monitored<- c('b_height', 'int_year', 'gap_int', 'sigma_gap', 'z_lp')

#mcmc params
nIter <- 500
nThin <- 1
nBurnIn <- 100
nChains <- 3

#call to JAGS
sur_mod <- jags(jags_data, inital, monitored, "sur_model.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

setwd(out_ob_loc)
save(sur_mod, file = 'sur_mod.Rdata')
#load('sur_mod.Rdata')

jags_obj_mc <- as.mcmc(sur_mod)
plot(jags_obj_mc, ask = TRUE)
sur_mod$BUGSoutput$summary[1:100,]


## TODO: IT IS WORKING WITH ONLY THE YEAR INTERCEPT (SO ALL PARAMETERS ARE IDENTIFIABLE) AND HEIGHT. CANNOT HAVE AN INDIVIDUAL 
## INTERCEPT SINCE BY DEFINITION EACH INDIVDUAL ONLY DIES ONCE.
## ADD A DENSITY KERNEL EFFECT B_DEN * SUM(HEIGHT[K] * DIST[K] * EXP(-W)) TO MODEL
## ADD SPATIAL EFFECT TO MODEL

# get a simple model of 
dist_height = dist_height(xLoc = sur_dat$X, yLoc = sur_dat$Y, height = sur_dat$height, year = sur_dat$year)

jags_data <- list(sur = sur_dat$sur, year = sur_dat$year_num, height = sur_dat$height_mod,
  num_years = max(sur_dat$year_num), num_obs = length(sur_dat$sur), dists = dist_height$dist,
  num_cont = dist_height$num_contemp, height_mat = dist_height$height_mat)

#intial values
inital <- function() list(b_height = runif(1, -2, 2), int_year = runif(jags_data$num_years, -2, 2)) 

# Parameters monitored
monitored<- c('b_height', 'b_den', 'dd_den', 'int_year', 'dist_height')

#mcmc params
nIter <- 1000
nThin <- 1
nBurnIn <- 50
nChains <- 3

#call to JAGS
set.seed(123)
sur_dd_mod <- jags(jags_data, inital, monitored, "sur_model_den.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

jags_obj_mc <- as.mcmc(sur_dd_mod)
plot(jags_obj_mc, ask = TRUE)
sur_dd_mod$BUGSoutput$summary[1:10,]

# Simple version works, try one with an interaction with height, as density might only matter for small plants
#The DIC is lower for the density model

## Spatial predictive proccess
# also try a spatial error version

dist_mat = distance_matrix(sur_dat$X, sur_dat$Y)

knots = knot_coords(x_min = min(sur_dat$X), x_max = max(sur_dat$X), y_min = min(sur_dat$Y), 
  y_max = max(sur_dat$Y), res = 10, dat_x = sur_dat$X, dat_y = sur_dat$Y, iso_dist = 20)

#look at knots
par(mar = c(4, 4, 0, 0))
plot(x = sur_dat$X, y = sur_dat$Y, bty = 'n', pch = 19, col = grey(0.7), tck = 0.015, xlab = 'x location (m)',
  xaxt = 'n', ylab = 'y location (m)', yaxt = 'n', xlim = c(-20, 160), ylim = c(-350, 60), cex = 0.5)
axis(side = 1, at = c(-20, 20, 60, 100, 140), tck = 0.015)
axis(side = 2, at = seq(-350, 50, 100), tck = 0.015)
points(knots, pch = 3, col = 'red', cex = 0.25)

dm_knots = distance_matrix(knots[,'x_cord'], knots[, 'y_cord'])
dm_knots_obs = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knots[,'x_cord'], y2 = knots[, 'y_cord'])
dm_obs = distance_matrix(sur_dat$X, sur_dat$Y)

# fit the spatial predictive process model
jags_data <- list(sur = sur_dat$sur, year = sur_dat$year_num, num_years = max(sur_dat$year_num), 
  num_obs = length(sur_dat$sur), height = sur_dat$height_mod, num_knots = dim(dm_knots)[1], s_knot_dists = dm_knots, 
  s_knot_ob_dist = dm_knots_obs)

#intial values
inital <- function() list(int_year = runif(jags_data$num_years, -2, 2)) 

# Parameters monitored
monitored<- c('int_year', 'dd_spp', 'z_lp', 'sigma_spp')

#mcmc params
nIter <- 1000
nThin <- 1
nBurnIn <- 1
nChains <- 3

#call to JAGS
setwd(code_loc)
sur_spp_mod <- jags(jags_data, inital, monitored, "sur_model_GSPP.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

setwd(out_ob_loc)
save(sur_spp_mod, file = 'sur_spp10_mod.Rdata')
#load('sur_spp10_mod.Rdata')

jags_obj_mc <- as.mcmc(sur_spp_mod)
plot(jags_obj_mc, ask = TRUE)
sur_dd_mod$BUGSoutput$summary[1:10,]

# converges very quickly, but finds no spatial auto-correlation in survival, but I know there is 
# a spatial effect at a fairly small scale, about the 1m resolution, which is a bit 
# beyond the GPP computationally try a spatially lagged model

## spatially auto-correlation model
# fit the non-spatial model, then fit a spatial lag model to the residual 
sac_sur_dists = dist_res_neigh(xLoc = sur_dat$X, yLoc = sur_dat$Y, year = sur_dat$year, win_r = 10)
 ## This works if win_r >15 

# fit the spatial predictive process model
jags_data <- list(sur1 = sur_dat$sur, sur2 = sur_dat$sur, year = sur_dat$year_num, num_years = max(sur_dat$year_num), 
  num_obs = length(sur_dat$sur), num_gaps = max(sur_dat$gapID_num), height = sur_dat$height_mod, 
  num_neigh = sac_sur_dists$num_contemp, inds_neigh = sac_sur_dists$ind_mat, dist_neigh = sac_sur_dists$dist, 
  gapID = sur_dat$gapID_num)
  
inital <- function() list(b_height = runif(1, -20, 20), 
  int_year = runif(jags_data$num_years, -20, 20))
  #gap_int = runif(jags_data$num_gaps, -20, 20)) 

# Parameters monitored
monitored<- c('b_height', 'b_sac', 'dd_sac', 'int_year', 'gap_int', 'sigma_gap', 'z_lp')

#mcmc params
nIter <- 1000
nThin <- 1
nBurnIn <- 100
nChains <- 3

#call to JAGS
sur_sac_mod <- jags(jags_data, inital, monitored, "sur_model_sac.jags", n.chains = nChains, 
  n.thin = nThin, n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

setwd(out_ob_loc)
save(sur_sac_mod, file = 'sur_sac_mod.Rdata')
#load('sur_sac_mod.Rdata')

jags_obj_mc <- as.mcmc(sur_sac_mod)
plot(jags_obj_mc, ask = TRUE)
sur_sl_mod$BUGSoutput$summary[1:10,]

## Sort of works, but confounds the year and the space effect, so need to take the survival prob across all num_years
## try a window effect on the SPP model 
knots = knot_coords(x_min = min(sur_dat$X), x_max = max(sur_dat$X), y_min = min(sur_dat$Y), 
 y_max = max(sur_dat$Y), res = 0.75, dat_x = sur_dat$X, dat_y = sur_dat$Y, iso_dist = 3)
 
dm_knots = distance_matrix(knots[,'x_cord'], knots[, 'y_cord'])
dm_knots_obs = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knots[,'x_cord'], y2 = knots[, 'y_cord'])

dm_knots_sparse = distmat_2_sparse(dm_knots, win_r = 5.5)
dm_knots_obs_sparse = distmat_2_sparse2(dm_knots_obs, win_r = 5.5)

jags_data <- list(sur = sur_dat$sur, year = sur_dat$year_num, num_years = max(sur_dat$year_num), 
  num_obs = length(sur_dat$sur), num_knots = dim(dm_knots)[1], knot_zero_ent = dim(dm_knots_sparse$zero_ind)[2], 
  k_zero_ind = dm_knots_sparse$zero_ind, knot_ent = dim(dm_knots_sparse$ent_ind)[2], 
  k_non_zero_ind = dm_knots_sparse$ent_ind, knot_dists = dm_knots_sparse$dist, 
  obs_knot_zero = dim(dm_knots_obs_sparse$zero_ind)[2], ob_zero_ind = dm_knots_obs_sparse$zero_ind, 
  obs_knots = dim(dm_knots_obs_sparse$ent_ind)[2], ob_non_zero_ind = dm_knots_obs_sparse$ent_ind, 
  s_knot_ob_dist = dm_knots_obs_sparse$dist, win_r = 5.5)
  
#intial values
inital <- function() list(int_year = runif(jags_data$num_years, -2, 2)) 

# Parameters monitored
monitored<- c('int_year', 'dd_spp', 'lp', 'sigma_spp')

#mcmc params
nIter <- 10
nThin <- 1
nBurnIn <- 1
nChains <- 2

#call to JAGS
sur_spp_win_mod <- jags(jags_data, inital, monitored, "sur_model_SPP_window.jags", n.chains = nChains, n.thin = nThin, 
  n.iter = nIter, n.burnin = nBurnIn, working.directory = getwd())

jags_obj_mc <- as.mcmc(sur_spp_win_mod)
plot(jags_obj_mc, ask = TRUE)
sur_dd_mod$BUGSoutput$summary[1:10,]

## This has many porblems inverting the knot distance matrix, possibly due to the grid 
## being made irregular by the window. Also it was very slow, took about 50 minuets to finally fail.

## So far my best plan to to use the GSPP to show that there is not much happening 
## at larger scales (say at a knot resolution of 10m, this should run in a resonable time).
## Then use a spatial lag proccess with a 10m window (to cut down the number of distance
## calulations) to pick up the fine scale spatial effects. Need to be careful gnerating the 
## predictor set and the same indivudals are entered multipule times, so those that survive 
## longer will get up weighted in prediction as they have more entries. To give every neighbour 
## equal weight need one estimated Y per individual
