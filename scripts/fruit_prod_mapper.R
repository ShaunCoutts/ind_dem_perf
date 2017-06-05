# Make a map of study period fruit production spatial structure 
# comapire it to fruit production predicted by a non-spatial model

library(dplyr)
library(tidyr)
library(colorspace)
library(rstan)
library(arrayhelpers)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)

dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

# get some helper functions 
setwd(code_loc)
source('demographic_perf.R')
source('dist_neigh_setup.R')

# set up a few constants from the data
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

# get the locations of the knots
setwd(dat_loc)
sur_dat = read.csv('sur_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
rep_dat = read.csv('rep_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
gr_dat = read.csv('gr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)

sur_kl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5) 
rep_kl = knot_coords2(dat_x = rep_dat$X, dat_y = rep_dat$Y, min_dist = 0.5) 
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


inner_join(sur_kl_df, rep_kl_df, by  = 'locID')

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

# get the intail height distribution from the 2002 height data
gr_dat_firt = filter(gr_dat, year == 2003)
z0_den = density(gr_dat_firt$height_prev)
# look at log normal distribution
z0_mean = mean(log(gr_dat_firt$height_prev))
z0_sd = sd(log(gr_dat_firt$height_prev))
z0 = rlnorm(1000, z0_mean, z0_sd)
plot(z0_den)
hist(z0, add = TRUE, prob = TRUE, breaks = seq(0, 100, 2))
# log normal looks pretty good, use that

# pull in the stan objets to sample the postieros from 
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

# hold the results
num_samps = 1000
RJJ_space = matrix(NA, nrow = length(kl_common), ncol = num_samps)

# make the distributions for the spatial model
for(i in 1:length(kl_common)){
  
  print(i)
  RJJ_space[i, ] = R_samp_space(SPP_sur, SPP_gr, SPP_rep, fruit_stan, num_samps, 
    z0_mean, z0_sd, mean_height, rep_mean_height, sur_mean_height, i, kl_common)

}

# make the distribution for the non-spatial model
RJJ_NS = R_samp_nospace(NS_sur, NS_gr, NS_rep, fruit_stan, num_samps, 
  z0_mean, z0_sd, mean_height, rep_mean_height, sur_mean_height)

# find difference between the spatial and non-spatial model
fruit_X = 0:20000
CDF_NS = ecdf(RJJ_NS)
pred_NS = CDF_NS(fruit_X)

desty_shift = numeric(length(kl_common))

for(i in 1:length(kl_common)){

  CDF = ecdf(RJJ_space[i, ])
  desty_shift[i] = sum(pred_NS - CDF(fruit_X)) * (1 / max(fruit_X)) 

}
hist(desty_shift)

# data frame of locations and the shift in predicted fruit number between 
# spatial and non-spatial model
loc_inds = sapply(kl_common, FUN = function(x) x$gr_ind)

R_sim_df = data.frame(X = gr_kl_df$X[loc_inds], Y = gr_kl_df$Y[loc_inds], den_shift = desty_shift)

R_sim_map = ggplot(R_sim_df, aes(X, Y)) + geom_point(aes(color = den_shift), size = 0.5) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = "right") + 
  labs(title = bquote(''*R[j0]^J*''), colour = 'shift') + scale_colour_gradient2(space="Lab") +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m') 

setwd(out_plot_loc)
pdf('life_fruit_sim_SPP_vs_NS.pdf', width = 10, height = 10)
  
  R_sim_map
  
dev.off()

# try plot the actual numbers of fruits rather than the shift, also put the upper and lower of each dist on the 
# side 
fruit_nums = apply(RJJ_space, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975))

# turn to a data frame
fn_df = data.frame(X = gr_kl_df$X[loc_inds], Y = gr_kl_df$Y[loc_inds], median = fruit_nums[2, ],
  lq = fruit_nums[1, ], uq = fruit_nums[3, ])

# set up my three ggplots
R_sim_med = ggplot(fn_df, aes(X, Y)) + geom_point(aes(color = median), size = 1) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0),
    plot.title = element_text(hjust = 0)) + 
  labs(title = bquote('a) median '*R[j0]^J*''), colour = '') + scale_colour_gradient2(space="Lab") +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m')
  
R_sim_lq = ggplot(fn_df, aes(X, Y)) + geom_point(aes(color = lq), size = 1) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0),
    plot.title = element_text(hjust = 0)) + 
  labs(title = bquote('c) lower 95% CI of '*R[j0]^J*''), colour = '') + scale_colour_gradient2(space="Lab")

R_sim_uq = ggplot(fn_df, aes(X, Y)) + geom_point(aes(color = uq), size = 1) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0),
    plot.title = element_text(hjust = 0)) + 
  labs(title = bquote('b) upper 95% CI of '*R[j0]^J*''), colour = '') + scale_colour_gradient2(space="Lab") +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m') 


setwd(out_plot_loc)
pdf('life_fruit_num_sim_map.pdf', width = 15, height = 10)
  
  ggdraw() + 
    draw_plot(R_sim_med, 0, 0, 0.65, 1) + 
    draw_plot(R_sim_lq, 0.65, 0, 0.35, 0.5) +
    draw_plot(R_sim_uq, 0.65, 0.5, 0.35, 0.5)

dev.off()

#############################################################################################
## Now try and make the same plot with the expectations under an IPM framework

#height domaine 
dz = 1
Z = seq(0, max(rep_dat$height) * 1.3, dz)
num_samps = 500

RJJ_space = matrix(NA, nrow = length(kl_common), ncol = num_samps)

for(i in 1:length(kl_common)){
  
  print(i)
  RJJ_space[i, ] = R_E_samp(SPP_sur, SPP_gr, SPP_rep, fruit_stan, Z, dz, z0_mean, z0_sd,
    mean_height, rep_mean_height, sur_mean_height, i, kl_common, num_samps)

}


# try plot the actual numbers of fruits rather than the shift, also put the upper and lower of each dist on the 
# side 
fruit_nums = apply(RJJ_space, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.5, 0.975))

# turn to a data frame
loc_inds = sapply(kl_common, FUN = function(x) x$gr_ind)
fn_df = data.frame(X = gr_kl_df$X[loc_inds], Y = gr_kl_df$Y[loc_inds], median = fruit_nums[2, ],
  lq = fruit_nums[1, ], uq = fruit_nums[3, ])

# set up my three ggplots
R_sim_med = ggplot(fn_df, aes(X, Y)) + geom_point(aes(color = median), size = 1) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0),
    plot.title = element_text(hjust = 0)) + 
  labs(title = bquote('a) median '*R[j0]^J*''), colour = '') + scale_colour_gradient2(space="Lab") +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m')
  
R_sim_lq = ggplot(fn_df, aes(X, Y)) + geom_point(aes(color = lq), size = 1) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0),
    plot.title = element_text(hjust = 0)) + 
  labs(title = bquote('c) lower 95% CI of '*R[j0]^J*''), colour = '') + scale_colour_gradient2(space="Lab")

R_sim_uq = ggplot(fn_df, aes(X, Y)) + geom_point(aes(color = uq), size = 1) + xlim(-20, 140) + ylim(-200, 50) + 
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0),
    plot.title = element_text(hjust = 0)) + 
  labs(title = bquote('b) upper 95% CI of '*R[j0]^J*''), colour = '') + scale_colour_gradient2(space="Lab") +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m') 


setwd(out_plot_loc)
pdf('life_fruit_num_E_map.pdf', width = 15, height = 10)
  
  ggdraw() + 
    draw_plot(R_sim_med, 0, 0, 0.65, 1) + 
    draw_plot(R_sim_lq, 0.65, 0, 0.35, 0.5) +
    draw_plot(R_sim_uq, 0.65, 0.5, 0.35, 0.5)

dev.off()

# takes a long time to generate the distribution of expected fruit number so save it for future use
setwd(out_ob_loc )
save(fn_df, file = 'R_E_samp_df.Rdata')
