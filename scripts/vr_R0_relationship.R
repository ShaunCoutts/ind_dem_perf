# script to try and untangle how each vital rate affects R0

library(dplyr)
library(tidyr)
library(colorspace)
library(rstan)
library(arrayhelpers)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)
library(lme4)
library(lmtest)
library(optimx)
library(MuMIn)

dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

setwd(dat_loc)
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
setwd(code_loc)
source('dist_neigh_setup.R')
source('model_perf_output_helper.R')
source('demographic_perf.R')

# First step is to see how the spatial efect for each vital rate correlates with each other to get 
# a sense of what is going on. 
# this will take some data massaging to work out which knot point is in which gap, and match up knot
# locations between the vital rates, since each vital rate has a slightly differnt data set, and so 
# different knot locations

# get the data sets used to fit the model
setwd(dat_loc)
sur_dat = read.csv('sur_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
rep_dat = read.csv('rep_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
gr_dat = read.csv('gr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
setwd(out_ob_loc)
ob_name = load('R_E_samp_df.Rdata')

setwd(out_ob_loc)
# survival
ob_name = load('NS_sur_stan.Rdata')
NS_sur = NS_sur_stan
ob_name = load('SPP_NK_sur_stan.Rdata') 
SPP_sur = SPP_NK_sur_stan
# reproduction
ob_name = load('SPP_NS_rep_stan.Rdata')
NS_rep = SPP_NS_rep_stan
ob_name = load('SPP_NK_rep_stan_ap.Rdata') 
SPP_rep = SPP_NK_rep_stan
# growth
ob_name = load('NS_gr_stan.Rdata')
NS_gr = NS_gr_stan
ob_name = load('SPP_NK_gr_stan.Rdata') 
SPP_gr = SPP_gr_stan
# fruit production
ob_name = load('fruit_num_stan.Rdata')
fruit_stan = fruit_num_stan

# pull the spp from the models
sur_spp = apply(extract_flat(SPP_sur, pars = c('spp')), MARGIN = 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))
sur_kl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5) 

rep_spp = apply(extract_flat(SPP_rep, pars = c('spp')), MARGIN = 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))
rep_kl = knot_coords2(dat_x = rep_dat$X, dat_y = rep_dat$Y, min_dist = 0.5) 

gr_spp = apply(extract_flat(SPP_gr, pars = c('spp')), MARGIN = 2, FUN = quantile, probs = c(0.025, 0.5, 0.975))
gr_kl = knot_coords2(dat_x = gr_dat$X, dat_y = gr_dat$Y, min_dist = 0.5) 

# find the knots common to all vital rates, and ensure they get the right knot index  
sur_kl_df = data.frame(X = round(sur_kl[, 1], 3), Y = round(sur_kl[, 2], 3))
sur_kl_df = mutate(sur_kl_df, locID = paste0(X, ':', Y))
sur_kl_df$rowID = 1:length(sur_kl_df$X)

rep_kl_df = data.frame(X = rep_kl[, 1], Y = rep_kl[, 2])
rep_kl_df = mutate(rep_kl_df, locID = paste0(X, ':', Y))
rep_kl_df$rowID = 1:length(rep_kl_df$X)

gr_kl_df = data.frame(X = gr_kl[, 1], Y = gr_kl[, 2])
gr_kl_df = mutate(gr_kl_df, locID = paste0(X, ':', Y))
gr_kl_df$rowID = 1:length(gr_kl_df$X)


# associate a gap ID to each knot location
sur_dat$gapID = sapply(sur_dat$uID, FUN = function(x) strsplit(x, split = ':')[[1]][1])
rep_dat$gapID = sapply(rep_dat$uID, FUN = function(x) strsplit(x, split = ':')[[1]][1])
gr_dat$gapID = sapply(gr_dat$uID, FUN = function(x) strsplit(x, split = ':')[[1]][1])

sur_kl_df$gapID = NA
for(i in 1:length(sur_kl_df$X)){

  gap = unique(sur_dat$gapID[find_near(sur_kl_df$X[i], sur_kl_df$Y[i], sur_dat$X, sur_dat$Y, dist = 0.25)])
  sur_kl_df$gapID[i] = ifelse(length(gap) > 0, gap, NA) 

}

rep_kl_df$gapID = NA
for(i in 1:length(rep_kl_df$X)){

  gap = unique(rep_dat$gapID[find_near(rep_kl_df$X[i], rep_kl_df$Y[i], rep_dat$X, rep_dat$Y, dist = 0.25)])
  rep_kl_df$gapID[i] = ifelse(length(gap) > 0, gap, NA) 

}

gr_kl_df$gapID = NA
for(i in 1:length(gr_kl_df$X)){

  gap = unique(gr_dat$gapID[find_near(gr_kl_df$X[i], gr_kl_df$Y[i], gr_dat$X, gr_dat$Y, dist = 0.25)])
  gr_kl_df$gapID[i] = ifelse(length(gap) > 0, gap, NA) 

}


# build a common set of knot locations by taking vital rate with least number of 
# knots (gr_kl_df), then finding all the knot locatons in the other 2 that are within
# 25 cm of the knot locations in gr_kl, and add those indexes to a key
kl_common = list()
count = 1
for(i in 1:length(gr_kl_df$X)){

  sur_inds = find_near(gr_kl_df$X[i], gr_kl_df$Y[i], sur_kl_df$X, sur_kl_df$Y, 0.25)
  rep_inds = find_near(gr_kl_df$X[i], gr_kl_df$Y[i], rep_kl_df$X, rep_kl_df$Y, 0.25)
  
  if(length(sur_inds) > 0 & length(rep_inds) > 0 ){
  
    kl_common[[count]] = list(gr_ind = i, sur_ind = sur_inds, rep_ind = rep_inds)
    count = count + 1
    
  }

}

# make pair wise knot location index key for survial and reproduction
kl_sur_rep_key = list()
count = 1
for(i in 1:length(sur_kl_df$X)){

  rep_inds = find_near(sur_kl_df$X[i], sur_kl_df$Y[i], rep_kl_df$X, rep_kl_df$Y, 0.25)
  
  if(length(rep_inds) > 0){
  
    kl_sur_rep_key[[count]] = list(sur_ind = i, rep_ind = rep_inds)
    count = count + 1
    
  }

}

# make data frame with knot locations, gap ID and spatial effect for survival and reproduction
sur_rep_kl_df = sur_kl_df 
sur_rep_kl_df$sur_medGPP = NA
sur_rep_kl_df$sur_lqGPP = NA
sur_rep_kl_df$sur_uqGPP = NA
sur_rep_kl_df$rep_medGPP = NA
sur_rep_kl_df$rep_lqGPP = NA
sur_rep_kl_df$rep_uqGPP = NA

for(i in 1:length(kl_sur_rep_key)){
  
  sur_rep_kl_df[i, c('sur_lqGPP', 'sur_medGPP', 'sur_uqGPP')] = sur_spp[, kl_sur_rep_key[[i]]$sur_ind] 
  sur_rep_kl_df[i, c('rep_lqGPP', 'rep_medGPP', 'rep_uqGPP')] = rep_spp[, kl_sur_rep_key[[i]]$rep_ind] 
  
}

# make data frame with knot locations, gap ID and spatial effect for growth survival and reproduction, 
# and estimate of fruit_prod
vr_kl_df = data.frame(row_ID = 1:length(kl_common), gapID = NA, X = NA, Y = NA, 
  sur_medGPP = NA, sur_lqGPP = NA, sur_uqGPP = NA, rep_medGPP = NA, 
  rep_lqGPP = NA, rep_uqGPP = NA, gr_medGPP = NA, gr_lqGPP = NA,
  gr_uqGPP = NA, med_E_Rf = NA, lq_E_Rf = NA, uq_E_Rf = NA)

for(i in 1:length(kl_common)){
  
  vr_kl_df$gapID[i] = gr_kl_df$gapID[kl_common[[i]]$gr_ind]
  vr_kl_df$X[i] = gr_kl_df$X[kl_common[[i]]$gr_ind]
  vr_kl_df$Y[i] = gr_kl_df$Y[kl_common[[i]]$gr_ind]
  vr_kl_df[i, c('sur_lqGPP', 'sur_medGPP', 'sur_uqGPP')] = sur_spp[, kl_common[[i]]$sur_ind] 
  vr_kl_df[i, c('rep_lqGPP', 'rep_medGPP', 'rep_uqGPP')] = rep_spp[, kl_common[[i]]$rep_ind] 
  vr_kl_df[i, c('gr_lqGPP', 'gr_medGPP', 'gr_uqGPP')] = gr_spp[, kl_common[[i]]$gr_ind] 
  vr_kl_df[i, c('med_E_Rf', 'lq_E_Rf', 'uq_E_Rf')] = fn_df[i, c('median', 'lq', 'uq')]
}

####################################################################################################
# now I have the data all set up I can start plotting
plot(vr_kl_df[, c('sur_medGPP', 'rep_medGPP', 'gr_medGPP', 'med_E_Rf')])

# look at each vital rate correlated with each other, also look within each patch 
ggplot(vr_kl_df, aes(sur_medGPP, rep_medGPP, colour = gapID)) + geom_point() +
    facet_wrap(~ gapID, nrow = 6)

# make a reduced data set with only those data points in patches shared by more than three locations
gap_groups = group_by(vr_kl_df, gapID)
gap_counts = summarize(gap_groups, count = n(), cor_sur_rep = cor(sur_medGPP, rep_medGPP), 
  cor_sur_gr = cor(sur_medGPP, gr_medGPP), cor_rep_gr = cor(rep_medGPP, gr_medGPP))

gaps_n4 = gap_counts$gapID[gap_counts$count > 2]

vr_kl_redu = filter(vr_kl_df, gapID %in% gaps_n4)
vr_kl_n4 = filter(gap_counts, gapID %in% gaps_n4)

setwd(out_plot_loc)
pdf(file = 'gap_vr_correlation.pdf', width = 15, height = 15)
  
  ggplot(vr_kl_redu, aes(sur_medGPP, rep_medGPP, colour = gapID)) + geom_point() +
    annotate('text', label = round(cor(vr_kl_redu$sur_medGPP, vr_kl_redu$rep_medGPP), 2), x = -2.5, y = -2.5)
  ggplot(vr_kl_redu, aes(sur_medGPP, rep_medGPP)) + geom_point() +
    geom_text(data = vr_kl_n4, aes(label = round(cor_sur_rep, 2), x = 1, y = 1)) +
    facet_wrap(~ gapID, nrow = 6)

  ggplot(vr_kl_redu, aes(sur_medGPP, gr_medGPP, colour = gapID)) + geom_point() +
    annotate('text', label = round(cor(vr_kl_redu$sur_medGPP, vr_kl_redu$gr_medGPP), 2), x = -2.5, y = -2.5)
  ggplot(vr_kl_redu, aes(sur_medGPP, gr_medGPP)) + geom_point() +
    geom_text(data = vr_kl_n4, aes(label = round(cor_sur_gr, 2), x = 1, y = 1)) +
    facet_wrap(~ gapID, nrow = 6)
    
  ggplot(vr_kl_redu, aes(rep_medGPP, gr_medGPP, colour = gapID)) + geom_point() + 
    annotate('text', label = round(cor(vr_kl_redu$rep_medGPP, vr_kl_redu$gr_medGPP), 2), x = -2.5, y = -2.5)
  ggplot(vr_kl_redu, aes(rep_medGPP, gr_medGPP)) + geom_point() +
    geom_text(data = vr_kl_n4, aes(label = round(cor_rep_gr, 2), x = 1, y = 1)) +
    facet_wrap(~ gapID, nrow = 6)

dev.off()

# plot these same vital rates against fruit production
ggplot(vr_kl_redu, aes(sur_medGPP, med_E_Rf)) + geom_point() +
  facet_wrap(~ gapID, ncol = 4)

ggplot(vr_kl_redu, aes(rep_medGPP, med_E_Rf)) + geom_point() +
  facet_wrap(~ gapID, ncol = 4)

ggplot(vr_kl_redu, aes(gr_medGPP, med_E_Rf)) + geom_point() +
  facet_wrap(~ gapID, ncol = 4)

# try and see how much variance in R0 is explained by the spatial effect of each vitial rate 
# (only part of the R0 equation that can change).

Rf_full = lmer(med_E_Rf ~ sur_medGPP + rep_medGPP + gr_medGPP + (sur_medGPP + rep_medGPP + gr_medGPP | gapID), 
  data = vr_kl_redu, REML = TRUE)
  
# think about how to work out how much variance explained with just pairs, might not be doable
# try taking patches with 3 locations see if that helps
  
# try take out the random effects
Rf_rand_int = lmer(med_E_Rf ~ sur_medGPP + rep_medGPP + gr_medGPP + (1 | gapID), 
  data = vr_kl_redu, REML = TRUE)
  
lrtest(Rf_full, Rf_rand_int)  
# cannot take out random effect
summary(Rf_full)
ranef(Rf_full)
# take out each vr to see how much R squared 

## two tasks, first calculate R^2 for full model and then models without each vital rate to assess how much spatial effect of each vr 
## has on R0. Second redo but with data generated under the null model that in each location vital rates are independent (that is no trade-offs).
## If a vital rate explains very little of the variation in the observed, and more under the null, then it is evidence that in the observed
## data there are spatial tradeoffs, becuase if all vr were independent they would be less correlated 

# I need to re-implement a bit of stuff so a lot of prep shit first to get the right data and some constants from the data

setwd(dat_loc)
burn_dat = read.csv('hcdem_ABS_2015_checked.csv')

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
 
mean_height = mean(fruit_dat$height) 
site_num = 2

# get the rep data to get mean height from this data for the centering
# get rep data
setwd(dat_loc)
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactor = FALSE)
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

# find the mean height for centering
rep_dat = vr_loc_gen[!is.na(vr_loc_gen$rep), c('uID', 'uLoc', 'rep', 'year', 'height', 'X', 'Y')]
rep_dat = rep_dat[!is.na(rep_dat$height),]
rep_dat = rep_dat[rep_dat$rep <= 1, ]

rep_mean_height = mean(rep_dat$height)

# find the mean height for centering survival
sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y')]
# neighbour data we can use all data, even the first years observed height 
neigh_dat = sur_dat

# take out data from the first year observed since nothing can be observed dead in the first year
first_year = sapply(seq_along(sur_dat$year), FUN = function(x){

  min_year_group = min(sur_dat$year[sur_dat$uID == sur_dat$uID[x]])
  return(ifelse(sur_dat$year[x] == min_year_group, FALSE, TRUE))

})
sur_dat = sur_dat[first_year, ]
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]

sur_mean_height = mean(sur_dat$height_prev)


# first step create matrix of observed R0 from the stan objects
# I will have to re-implement some fo these calculations to do this 

dz = 1
Z = seq(0, max(rep_dat$height) * 1.3, dz)
num_samps = 500

# get intial distribution
gr_dat_firt = filter(gr_dat, year == 2003)
z0_den = density(gr_dat_firt$height_prev)
# look at log normal distribution
z0_mean = mean(log(gr_dat_firt$height_prev))
z0_sd = sd(log(gr_dat_firt$height_prev))

RJJ_space = matrix(NA, nrow = length(kl_common), ncol = num_samps)
sur_space = matrix(NA, nrow = length(kl_common), ncol = num_samps)
rep_space = matrix(NA, nrow = length(kl_common), ncol = num_samps)
gr_space = matrix(NA, nrow = length(kl_common), ncol = num_samps)

for(i in 1:length(kl_common)){
  
  print(i)
  out_ob = R_E_samp(SPP_sur, SPP_gr, SPP_rep, fruit_stan, Z, dz, z0_mean, z0_sd,
    mean_height, rep_mean_height, sur_mean_height, i, kl_common, num_samps)

  RJJ_space[i, ] = out_ob$RJJ
  sur_space[i, ] = out_ob$sur_spp
  rep_space[i, ] = out_ob$rep_spp
  gr_space[i, ] = out_ob$gr_spp
  
}

RJJ_space_list = list(RJJ = RJJ_space, sur_spp = sur_space, rep_spp = rep_space, gr_spp = gr_space)

# this takes ages to produce so save for future used
setwd(out_ob_loc)
save(RJJ_space_list, file = 'R0_space_dist_list.Rdata')
#space_ob_name = load('R0_null_dist_list.Rdata')


# randomization over parameter uncertianty under the null distribution where vital rates are not associated with a location 
num_rands = num_samps * length(kl_common)

RJJ_null = R_E_null(SPP_sur, SPP_gr, SPP_rep, fruit_stan, Z, dz, z0_mean, z0_sd,
    mean_height, rep_mean_height, sur_mean_height, kl_common, num_rands, num_samps)

# this takes ages to produce so save for future used
setwd(out_ob_loc)
save(RJJ_null, file = 'R0_null_dist_list.Rdata')
#null_ob_name = load('R0_null_dist_list.Rdata')


# first make a simple histagram to just see how the two distribtuions look over all 
# turn the dists of Rjj over all space and under the null dist into a data frame 
hist_df = data.frame(model = rep(c('spatially structured', 'null model'), each = length(RJJ_space_list$RJJ)), 
  RJJ = c(as.numeric(RJJ_space_list$RJJ), as.numeric(RJJ_null$RJJ)))

setwd(out_plot_loc)
pdf(file = 'RJJ_hist_space_V_null.pdf')
  
  ggplot(hist_df, aes(RJJ, color = model, fill = model)) + 
    geom_density(alpha = 0.2) + xlim(0, 750) +
    theme(legend.position = c(0.8, 0.5))
    
dev.off()


# now make a plot for the different R^2 
# first set up the matricies so the rows only refer to locations in gaps that have 3 or more locations 
gaps_n3 = gap_counts$gapID[gap_counts$count > 2]
gaps_n3_test = vr_kl_df$gapID %in% gaps_n3

RJJ_space_n3 = list(RJJ = RJJ_space_list$RJJ[gaps_n3_test, ], sur_spp = RJJ_space_list$sur_spp[gaps_n3_test, ],
  rep_spp = RJJ_space_list$rep_spp[gaps_n3_test, ], gr_spp = RJJ_space_list$gr_spp[gaps_n3_test, ])

RJJ_null_n3 = list(RJJ = RJJ_null2$RJJ[gaps_n3_test, ], sur_spp = RJJ_null2$sur_spp[gaps_n3_test, ],
  rep_spp = RJJ_null2$rep_spp[gaps_n3_test, ], gr_spp = RJJ_null2$gr_spp[gaps_n3_test, ])

# set up the different models 
form_full = 'RJJ ~ sur_spp + rep_spp + gr_spp + (sur_spp + rep_spp + gr_spp | gapID)'
form_sur_rep = 'RJJ ~ sur_spp + rep_spp + (sur_spp + rep_spp | gapID)'
form_sur_gr = 'RJJ ~ sur_spp + gr_spp + (sur_spp + gr_spp | gapID)'
form_rep_gr = 'RJJ ~ rep_spp + gr_spp + (rep_spp + gr_spp | gapID)'

# fit model to each sample in posterior and put in list  
full_mod_list = lmer_fitter(form_full, RJJ_space_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(full_mod_list, file = 'RJJ_model_list_space.Rdata')

full_mod_list_null = lmer_fitter(form_full, RJJ_null_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(full_mod_list_null, file = 'RJJ_model_list_null.Rdata')

sur_rep_mod_list = lmer_fitter(form_sur_rep, RJJ_space_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(sur_rep_mod_list, file = 'RJJ_sur_rep_model_list_space.Rdata')

sur_rep_mod_list_null = lmer_fitter(form_sur_rep, RJJ_null_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(sur_rep_mod_list_null, file = 'RJJ_sur_rep_model_list_null.Rdata')

sur_gr_mod_list = lmer_fitter(form_sur_gr, RJJ_space_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(sur_gr_mod_list, file = 'RJJ_sur_gr_model_list_space.Rdata')

sur_gr_mod_list_null = lmer_fitter(form_sur_gr, RJJ_null_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(sur_gr_mod_list_null, file = 'RJJ_sur_gr_model_list_null.Rdata')

rep_gr_mod_list = lmer_fitter(form_rep_gr, RJJ_space_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(rep_gr_mod_list, file = 'RJJ_rep_gr_model_list_space.Rdata')

rep_gr_mod_list_null = lmer_fitter(form_rep_gr, RJJ_null_n3, vr_kl_redu$gapID, control = lmerControl(optimizer = 'optimx', optCtrl = list(method = 'nlminb', maxit = 1000)))
setwd(out_ob_loc)
save(rep_gr_mod_list_null, file = 'RJJ_rep_gr_model_list_null.Rdata')

# load all these back in to save time 
setwd(out_ob_loc)
load('RJJ_model_list_space.Rdata')
load('RJJ_model_list_null.Rdata')
load('RJJ_sur_rep_model_list_space.Rdata')
load('RJJ_sur_rep_model_list_null.Rdata')
load('RJJ_sur_gr_model_list_space.Rdata')
load('RJJ_sur_gr_model_list_null.Rdata')
load('RJJ_rep_gr_model_list_space.Rdata')
load('RJJ_rep_gr_model_list_null.Rdata')

# make a violin plot of the R^2 for each model to see how much each of the fixed effects are explaining 
Rsq_full_space = sapply(full_mod_list, FUN = r.squaredGLMM)
Rsq_full_null = sapply(full_mod_list_null, FUN = r.squaredGLMM)

Rsq_sur_rep_space = sapply(sur_rep_mod_list, FUN = r.squaredGLMM)
Rsq_sur_rep_null = sapply(sur_rep_mod_list_null, FUN = r.squaredGLMM)

Rsq_sur_gr_space = sapply(sur_gr_mod_list, FUN = r.squaredGLMM)
Rsq_sur_gr_null = sapply(sur_gr_mod_list_null, FUN = r.squaredGLMM)

Rsq_rep_gr_space = sapply(rep_gr_mod_list, FUN = r.squaredGLMM)
Rsq_rep_gr_null = sapply(rep_gr_mod_list_null, FUN = r.squaredGLMM)

Rsq_mar_df = data.frame(struct = rep(c('spatial', 'null'), each = num_samps * 4), 
  model = rep(rep(c('sur + rep + gro', 'sur + rep', 'sur + gro', 'rep + gro'), each = num_samps), times = 2),
  Rsq = c(Rsq_full_space[1, ],  Rsq_sur_rep_space[1, ], Rsq_sur_gr_space[1, ],   Rsq_rep_gr_space[1, ], 
  Rsq_full_null[1, ],  Rsq_sur_rep_null[1, ], Rsq_sur_gr_null[1, ],   Rsq_rep_gr_null[1, ]))

Rsq_con_df = data.frame(struct = rep(c('spatial', 'null'), each = num_samps * 4), 
  model = rep(rep(c('sur + rep + gro', 'sur + rep', 'sur + gro', 'rep + gro'), each = num_samps), times = 2),
  Rsq = c(Rsq_full_space[2, ],  Rsq_sur_rep_space[2, ], Rsq_sur_gr_space[2, ],   Rsq_rep_gr_space[2, ], 
  Rsq_full_null[2, ],  Rsq_sur_rep_null[2, ], Rsq_sur_gr_null[1, ],   Rsq_rep_gr_null[2, ]))

  
# make the violin plots 
ggplot(Rsq_mar_df, aes(model, Rsq, colour = struct, fill = struct)) +
  geom_violin(alpha = 0.3)
  
ggplot(Rsq_con_df, aes(model, Rsq, colour = struct, fill = struct)) +
  geom_violin(alpha = 0.3)
  
  
# try fix up the null list so I can play around with it. Regenerate it tonight for real
RJJ_null2 = list(RJJ = matrix(RJJ_null$RJJ, ncol = num_samps), 
  sur_spp = matrix(RJJ_null$sur_spp, ncol = num_samps),
  rep_spp = matrix(RJJ_null$rep_spp, ncol = num_samps),
  gr_spp = matrix(RJJ_null$gr_spp, ncol = num_samps))

