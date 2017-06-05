# plotting script for individual performance over space

library(plyr)
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

setwd(dat_loc)
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
setwd(code_loc)
source('dist_neigh_setup.R')
source('model_perf_output_helper.R')

#need location data so drop the obervations that do not have a location
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X),]

#aggregate the data for plotting 
plot_dat = ddply(vr_loc_gen, .(uID), summarise, X = unique(X), Y = unique(Y),
  sur = min(sur), height = max(height, na.rm = TRUE), rep = sum(rep, na.rm = TRUE),
  num_years = unique(num_years))

# start with survival
sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]

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

# create knots at each data point, cluster some less than 50cm apart
# results in 575 knots. 50% of knots have another knot within 61cm and 75% have
# another knot within 0.82cm. 50% of knots have 8 or more neigbour knots within 2m
# and 75% have 5 or more neighbouring knots within 2m.
knot_cl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5) 

plot_dat_sur = ddply(sur_dat, .(uID), summarise, X = unique(X), Y = unique(Y),
  sur = min(sur))
  
# plot the vital rates in space 
setwd(out_loc)
pdf('spatial_vital_rates.pdf', height = 7.5, width = 11.25)
  par(mfrow = c(1, 3))
  # start with survival 
  plot(x = plot_dat_sur$X[plot_dat_sur$sur == 1], y = plot_dat_sur$Y[plot_dat_sur$sur == 1], 
    xlim = c(-20, 140), ylim = c(-345, 46), xlab = 'location (m)', ylab = 'loaction (m)',
    pch = 19, , bty = 'n', tck = 0.015, col = grey(0.7), main = 'survival')
  points(x = plot_dat_sur$X[plot_dat$sur == 0], y = plot_dat_sur$Y[plot_dat$sur == 0],
    pch = 19)
  points(knot_cl, col = 'red', pch = 3, cex = 0.25)
mtext('a)', side = 3, adj = 0)	
    
  #height, color by max height 
  col_height = sequential_hcl(max(round(plot_dat$height)), h = 260, c. = c(0, 100), l = c(10, 50))
  cols_plotted = col_height[round(plot_dat$height)]
  plot(x = plot_dat$X, y = plot_dat$Y, xlim = c(-20, 140), ylim = c(-345, 46), 
    xlab = 'location (m)', ylab = '', pch = 19, , bty = 'n', tck = 0.015, 
    col = cols_plotted, main = 'max height')
mtext('b)', side = 3, adj = 0)	

  #reproduction, color by number of years reproductive 
  col_height = sequential_hcl(max(plot_dat$rep) + 1, h = 260, c. = c(0, 100), l = c(10, 50))
  cols_plotted = col_height[plot_dat$rep + 1]
  plot(x = plot_dat$X, y = plot_dat$Y, xlim = c(-20, 140), ylim = c(-345, 46), 
    xlab = 'location (m)', ylab = '', pch = 19, , bty = 'n', tck = 0.015, 
    col = cols_plotted, main = 'reproduction')
mtext('c)', side = 3, adj = 0)	
dev.off()    
    
    
hist(plot_dat$num_years)  # not many indviduals with only one uear of data, and lots with 3, 4, and 5

# need to find a proxy for survival probability. A good candidate might be max height observed. It can be 
# observed for each indivudual. To see if it is a good proxy take the all the idividuals that were observed >10cm in the 
# first year and where their death was observed.

setwd(dat_loc)
vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)

year1_height = vr_loc_gen[vr_loc_gen$age == 0, ]
hist(year1_height$height, breaks = seq(0, 50, 1))

first_height = tapply(seq_along(vr_loc_gen$uID), INDEX = vr_loc_gen$uID, FUN = function(x){ 
  heights = vr_loc_gen$height[x]
  ages = which(vr_loc_gen$age[x] == 0)
  return(c(ifelse(length(ages) > 0, heights[ages], NA), min(vr_loc_gen$sur[x])))
})

suit = sapply(first_height, FUN = function(x) x[1] < 10 & x[2] == 0)
suit = ifelse(is.na(suit), FALSE, suit)
uIDs = names(suit)
uIDs = uIDs[suit]

vr_full_lc = vr_loc_gen[vr_loc_gen$uID %in% uIDs,]
plot_dat = ddply(vr_full_lc, .(uID), summarise, max_age = max(age, na.rm = TRUE), 
  max_height = max(height, na.rm = TRUE))

plot(x = plot_dat$max_age, y = log(plot_dat$max_height))  
cor(plot_dat$max_age, plot_dat$max_height) # 0.825  
h_mod = lm(log(max_height) ~ max_age + I(max_age ^ 2), data = plot_dat)
summary(h_mod)
#plot(h_mod)

new_dat = data.frame(max_age = seq(min(plot_dat$max_age), max(plot_dat$max_age), 0.1)) 
pred_age = predict(h_mod, newdata = new_dat, interval = 'prediction', level = 0.95) 

setwd(out_loc)
pdf('age_vs_height.pdf', height = 8, width = 8)
  plot(x = plot_dat$max_age, y = log(plot_dat$max_height), type = 'n', bty = 'n', xlab = 'max age observed (years)', 
    ylab = 'ln(max height observed (cm))', ylim = c(-1, 5))
  polygon(x = c(new_dat$max_age, rev(new_dat$max_age)), y = c(pred_age[,2], rev(pred_age[,3])), 
    col = 'lightgoldenrod', border = NA)
  points(x = jitter(plot_dat$max_age, amount = 0.1), y = log(plot_dat$max_height), pch = 19, cex = 0.5)  
  lines(x = new_dat$max_age, y = pred_age[,1])
dev.off()
##NOTE: used prediction interval not confidence interval

# high correlation, (R^2 = 0.79 for log model with quadratic term). Lots of spread but max height should act as 
# okay porxy for age reached, it will at least discrimnate between first year, second year (to an extent) 
# and then 3 and 4 together (line flattens out).

## observed versus predicted for the 4 survival models (no space), 25m grid, fine grid, combined 25m grid and fine grid.
# fitted versus predicted for 4 survival models 
setwd(out_ob_loc)
ob_name = load('NS_sur_stan.Rdata')
NS_sur = NS_sur_stan
ob_name = load('SPP_25_sur_stan.Rdata') 
SPP25_sur = SPP_sur_stan
ob_name = load('SPP_NK_sur_stan.Rdata') 
SPPfine_sur = SPP_NK_sur_stan
ob_name = load('SPP_25NK_sur_stan.Rdata') 
SPP_25fine_sur = SPP_sur_stan 

## split out the parameter samples and linear predictor samples
num_warmup = NS_sur@sim$warmup
num_samp = NS_sur@sim$iter
NS_sur_param = extract_flat(NS_sur, pars = c('b', 'sd_1', 'year_int'))
NS_sur_pred = extract_flat(NS_sur, pars = c('eta'))
SPP25_sur_pred = extract_flat(SPP25_sur, pars = c('eta'))
SPPfine_sur_pred = extract_flat(SPPfine_sur, pars = c('eta'))
SPP25fine_sur_pred = extract_flat(SPP_25fine_sur, pars = c('eta'))

non_space = obs_v_pred_bin(obs = sur_dat$sur, lp_samp = NS_sur_pred, num_resamp = 1000, jitter_am = 0.2,
  labels = labs(list(title = 'non-spatial', x= 'predicted', y = 'observed')))
GPP25 = obs_v_pred_bin(obs = sur_dat$sur, lp_samp = SPP25_sur_pred, num_resamp = 1000, jitter_am = 0.2,
  labels = labs(list(title = 'GPP 25m res', x= 'predicted', y = '')))
GPPfine = obs_v_pred_bin(obs = sur_dat$sur, lp_samp = SPPfine_sur_pred, num_resamp = 1000, jitter_am = 0.2,
  labels = labs(list(title = 'GPP fine scale', x= 'predicted', y = '')))
GPP25fine = obs_v_pred_bin(obs = sur_dat$sur, lp_samp = SPP25fine_sur_pred, num_resamp = 1000, jitter_am = 0.2,
  labels = labs(list(title = 'GPP 25m res and fine scale', x= 'predicted', y = '')))

ROC_sur_NS = ROC_curve(obs = sur_dat$sur, pred = NS_sur_pred, num_resamp = 200, th_res = 0.025)
ROC_NS = ROC_plot(ROC_sur_NS, labels = labs(list(title = '', x= 'False positive rate', y = 'True positive rate'))) 
AUC_NS = round(quantile(auc_dist(sur_dat$sur, NS_sur_pred, 1000), probs = c(0.5, 0.025, 0.975)), 3)
ROC_NS = ROC_NS + annotate('text', x = 0.5, y = 0.25, label = paste0('AUC: ',  AUC_NS[1], ' (', AUC_NS[2], '-', AUC_NS[3], ')'))

ROC_sur_GPP25 = ROC_curve(obs = sur_dat$sur, pred = SPP25_sur_pred, num_resamp = 200, th_res = 0.025)
ROC_GPP25 = ROC_plot(ROC_sur_GPP25, labels = labs(list(title = '', x= 'False positive rate', y = 'True positive rate'))) 
AUC_GPP25 = round(quantile(auc_dist(sur_dat$sur, SPP25_sur_pred, 1000), probs = c(0.5, 0.025, 0.975)), 3)
ROC_GPP25 = ROC_GPP25 + annotate('text', x = 0.5, y = 0.25, 
  label = paste0('AUC: ',  AUC_GPP25[1], ' (', AUC_GPP25[2], '-', AUC_GPP25[3], ')'))
  
ROC_sur_GPPfine = ROC_curve(obs = sur_dat$sur, pred = SPPfine_sur_pred, num_resamp = 200, th_res = 0.025)
ROC_GPPfine = ROC_plot(ROC_sur_GPPfine, labels = labs(list(title = '', x= 'False positive rate', y = 'True positive rate'))) 
AUC_GPPfine = round(quantile(auc_dist(sur_dat$sur, SPPfine_sur_pred, 1000), probs = c(0.5, 0.025, 0.975)), 3)
ROC_GPPfine = ROC_GPPfine + annotate('text', x = 0.5, y = 0.25, 
  label = paste0('AUC: ',  AUC_GPPfine[1], ' (', AUC_GPPfine[2], '-', AUC_GPPfine[3], ')'))
  
ROC_sur_GPP25fine = ROC_curve(obs = sur_dat$sur, pred = SPP25fine_sur_pred, num_resamp = 200, th_res = 0.025)
ROC_GPP25fine = ROC_plot(ROC_sur_GPP25fine, labels = labs(list(title = '', x= 'False positive rate', y = 'True positive rate'))) 
AUC_GPP25fine = round(quantile(auc_dist(sur_dat$sur, SPP25fine_sur_pred, 1000), probs = c(0.5, 0.025, 0.975)), 3)
ROC_GPP25fine = ROC_GPP25fine + annotate('text', x = 0.5, y = 0.25, 
  label = paste0('AUC: ',  AUC_GPP25fine[1], ' (', AUC_GPP25fine[2], '-', AUC_GPP25fine[3], ')'))
  
grid.arrange(non_space, GPP25, GPPfine, GPP25fine, ROC_NS, ROC_GPP25, ROC_GPPfine, ROC_GPP25fine, 
  ncol = 4, nrow = 2)  

# Make violin plot of year effects on prob scale. 
NS_sur_param = data.frame(extract_flat(NS_sur, pars = c('b', 'sd_1', 'year_int')))
GPP25_sur_param = data.frame(extract_flat(SPP25_sur, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv', 'sigma_spp')))
GPPfine_sur_param = data.frame(extract_flat(SPPfine_sur, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv', 'sigma_spp')))
GPP25fine_sur_param = data.frame(extract_flat(SPP_25fine_sur, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv25', 'sigma_spp25', 
  'dd_spp_invNN', 'sigma_sppNN')))

NS_sur_year = gather(NS_sur_param, year, effect_size, year_int.1.:year_int.5.)
df_year = data.frame(model = rep('non-spatial', dim(NS_sur_year)[1]), year = NS_sur_year$year, 
  effect_size = NS_sur_year$effect_size)
GPP25_sur_year = gather(GPP25_sur_param, year, effect_size, year_int.1.:year_int.5.)
df_year = rbind(df_year, data.frame(model = rep('GPP 25m res', dim(GPP25_sur_year)[1]), year = GPP25_sur_year$year, 
  effect_size = GPP25_sur_year$effect_size))
GPPfine_sur_year = gather(GPPfine_sur_param, year, effect_size, year_int.1.:year_int.5.)
df_year = rbind(df_year, data.frame(model = rep('GPP fine scale', dim(GPPfine_sur_year)[1]), year = GPPfine_sur_year$year, 
  effect_size = GPPfine_sur_year$effect_size))
#GPP25fine_sur_year = gather(GPP25fine_sur_param, year, effect_size, year_int.1.:year_int.5.)
#df_year = rbind(df_year, data.frame(model = rep('GPP 25m res + GPP fine scale', dim(GPP25fine_sur_year)[1]), year = GPP25fine_sur_year$year, 
#  effect_size = GPP25fine_sur_year$effect_size))

year_plt = ggplot(df_year, aes(year, effect_size)) + 
  geom_violin(fill = hcl(h = 180), color = hcl(h = 180, c = 70)) + 
  geom_abline(intercept = 0, slope = 0) + theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9))) +
  scale_x_discrete(labels = c('2003', '2004', '2005', '2006', '2007')) +
  facet_grid(.~model, scales="free", space="free") 

year_plt


# plot the decay curves
decay_plotter = function(x, dd_inv){
  dd = 1 / dd_inv
  return(exp(-dd * x))
}

# data for decay plot
x_ax = seq(0, 50, 1)
num_resamp = 200
rs_row = sample.int(dim(GPP25_sur_param)[1], num_resamp)
GPP25_sur_dd_rs = GPP25_sur_param$dd_spp_inv[rs_row]
GPPfine_sur_dd_rs = GPPfine_sur_param$dd_spp_inv[rs_row]

df_dd_plot = data.frame(model = rep(c('GPP 25m res', 'GPP fine scale'), each = length(x_ax) * num_resamp), 
  rep_lab = rep(rep(as.character(1:length(GPP25_sur_param$dd_spp_inv)), each = length(x_ax)), times = 2), 
  dd_inv = c(rep(GPP25_sur_dd_rs, each = length(x_ax)), rep(GPPfine_sur_dd_rs, each = length(x_ax))), 
  dist = rep(x_ax, times = 2 * num_resamp))
df_dd_plot$cor_curve = decay_plotter(df_dd_plot$dist, df_dd_plot$dd_inv)  


dd_plt = ggplot(df_dd_plot, aes(dist, cor_curve, group = rep_lab)) + 
  geom_line(color = hcl(h = 180, c = 70)) + 
  theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9))) +
  facet_grid(model~., scales="free", space="free") 

dd_plt

# try a color mat version 
GPP25_dd_quant = quantile(GPP25_sur_param$dd_spp_inv, probs = c(0.025, 0.5, 0.975))
GPPfine_dd_quant = quantile(GPPfine_sur_param$dd_spp_inv, probs = c(0.025, 0.5, 0.975))

x_ax = seq(0, 50, 1)
df_dd_quant = data.frame(model = rep(c('GPP 25m res', 'GPP fine scale'), each = length(x_ax) * 2),
  dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)), 
  quant_dd = c(rep(GPP25_dd_quant[1], length(x_ax)), rep(GPP25_dd_quant[3], length(x_ax)), 
    rep(GPPfine_dd_quant[1], length(x_ax)), rep(GPPfine_dd_quant[3], length(x_ax))), 
  med_dd = c(rep(GPP25_dd_quant[2], length(x_ax) * 2), rep(GPPfine_dd_quant[2], length(x_ax) * 2))) 

df_dd_quant$quant95 = decay_plotter(df_dd_quant$dist, df_dd_quant$quant_dd)  
df_dd_quant$quant50 = decay_plotter(df_dd_quant$dist, df_dd_quant$med_dd)  

dd_quant_plt = ggplot(df_dd_quant, aes(dist, quant95)) + 
  geom_polygon(aes(fill = hcl(h = 180), color = hcl(h = 180)), alpha = 0.5) +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_line(colour = grey(1)),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none") +
  geom_line(aes(dist, quant50)) + labs(x = 'distance', y = 'correlation') +
  facet_grid(model~., scales="free", space="free") 

dd_quant_plt

############################################################################################################################
## plot to look at performance of models and likelihood, mapped also by year 
# get the data sets used to fit the model
setwd(dat_loc)
sur_dat = read.csv('sur_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
rep_dat = read.csv('rep_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)
gr_dat = read.csv('gr_loc_gen_postburn.csv', header = TRUE, stringsAsFactors = FALSE)

# build up obs vs predicted for spatial and non-spatial model for survival rep and growth
# get the predictions that we need
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

NS_sur_pred = extract_flat(NS_sur, pars = c('eta'))
SPP_sur_pred = extract_flat(SPP_sur, pars = c('eta'))

NS_rep_pred = extract_flat(NS_rep, pars = c('eta'))
SPP_rep_pred = extract_flat(SPP_rep, pars = c('eta'))

NS_gr_pred = extract_flat(NS_gr, pars = c('mu'))
SPP_gr_pred = extract_flat(SPP_gr, pars = c('mu'))
NS_gr_sigma = extract_flat(NS_gr, pars = c('sigma'))
SPP_gr_sigma = extract_flat(SPP_gr, pars = c('sigma'))

pred_obs_sur = data.frame(model = rep(c('non-spatial', 'SPP'), each = length(sur_dat$sur)), 
  ob = jitter(rep(sur_dat$sur, times = 2), amount = 0.1),
  pred = c(as.vector(apply(NS_sur_pred, MARGIN = 2, FUN = function(x) plogis(mean(x)))),
    as.vector(apply(SPP_sur_pred, MARGIN = 2, FUN = function(x) plogis(mean(x))))))
  
p <- ggplot(pred_obs_sur, aes(ob, pred))
pvo_sur <- p + geom_point() + 
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none", plot.title = element_text(size = 20), axis.title = element_text(size = 20),
    strip.background = element_blank(), strip.text = element_blank()) + 
  annotate('text', x = 0.5, y = 0.80, label = paste0('logLik ', 
    c(round(sum(logLik_bin(post_mat = NS_sur_pred, obs = sur_dat$sur)), 2), 
    round(sum(logLik_bin(post_mat = SPP_sur_pred, obs = sur_dat$sur)), 2))), size = 5) +
  labs(title = 'survival', x = 'observed', y = 'predicted') + ylim(0, 1) +
  facet_grid(model~.) 

pred_obs_rep = data.frame(model = rep(c('non-spatial', 'SPP'), each = length(rep_dat$rep)), 
  ob = jitter(rep(rep_dat$rep, times = 2), amount = 0.1),
  pred = c(as.vector(apply(NS_rep_pred, MARGIN = 2, FUN = function(x) plogis(mean(x)))),
    as.vector(apply(SPP_rep_pred, MARGIN = 2, FUN = function(x) plogis(mean(x))))))
  
p <- ggplot(pred_obs_rep, aes(ob, pred))
pvo_rep <- p + geom_point() + 
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none", plot.title = element_text(size = 20), axis.title = element_text(size = 20),
    strip.background = element_blank(), strip.text = element_blank(),
    axis.title.y = element_blank()) + 
  annotate('text', x = 0.5, y = 0.80, label = paste0('logLik ', 
    c(round(sum(logLik_bin(post_mat = NS_rep_pred, obs = rep_dat$rep)), 2), 
    round(sum(logLik_bin(post_mat = SPP_rep_pred, obs = rep_dat$rep)), 2))), size = 5) +
  labs(title = 'reproduction', x = 'observed', y = 'predicted') + ylim(0, 1) +
  facet_grid(model~.) 

pred_obs_gr = data.frame(model = rep(c('non-spatial', 'SPP'), each = length(gr_dat$height)), 
  ob = jitter(rep(gr_dat$height, times = 2), amount = 0.1),
  pred = c(as.vector(apply(NS_gr_pred, MARGIN = 2, FUN = function(x) mean(x))),
    as.vector(apply(SPP_gr_pred, MARGIN = 2, FUN = function(x) mean(x)))))
  
p <- ggplot(pred_obs_gr, aes(ob, pred))
pvo_gr <- p + geom_point() + 
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.9)),
    legend.position = "none", plot.title = element_text(size = 20), axis.title = element_text(size = 20),
    axis.title.y = element_blank(), strip.text = element_text(size = 20)) + 
  annotate('text', x = 60, y = 10, label = paste0('logLik ', 
    c(round(sum(logLik_cont(post_mu = NS_gr_pred, obs = gr_dat$height, post_sd = NS_gr_sigma)), 2), 
    round(sum(logLik_cont(post_mu = SPP_gr_pred, obs = gr_dat$height, post_sd = SPP_gr_sigma)), 2))), size = 5) +
  labs(title = 'growth', x = 'observed', y = 'predicted') + geom_abline(intercept = 0, slope = 1) + ylim(0, 65) +
  facet_grid(model~.) 

setwd(out_plot_loc)
pdf('prev_vs_obs.pdf', height = 10, width = 15)
  grid.arrange(pvo_sur, pvo_rep, pvo_gr, nrow = 1, ncol = 3)
dev.off()


# map of likliehood for each point for each year for both non-spatial and spatial models
# survival first
sur_pred_map_df = data.frame(X = sur_dat$X, Y = sur_dat$Y, year = sur_dat$year, 
  pred = apply(plogis(SPP_sur_pred), MARGIN = 2, FUN = mean))

sur_pmap = ggplot(sur_pred_map_df, aes(X, Y)) + geom_point(aes(color = pred), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) + 
  labs(color = 'prob surv') + 
  annotate('segment', x = 0, xend = 50, y =-250, yend = -250, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 25, y =-240, label = c('50m', '', '', '', '')) + 
  facet_grid(.~year, scales="free", space="free") 
  
sur_lik_map_df = data.frame(X = sur_dat$X, Y = sur_dat$Y, year = sur_dat$year, 
  logLik = logLik_bin(post_mat = SPP_sur_pred, obs = sur_dat$sur))

sur_lmap = ggplot(sur_lik_map_df, aes(X, Y)) + geom_point(aes(colour = exp(logLik)), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) + 
  labs(color = 'likliehood') + facet_grid(.~year, scales="free", space="free") 

setwd(out_plot_loc)
pdf('sur_prev_lik_map.pdf', height = 8, width = 13)
  grid.arrange(sur_pmap, sur_lmap, ncol = 1, nrow = 2)  
dev.off()

# make a zomed in verison see how patches look closer in 
sur_pmap_zoom = sur_pmap + xlim(30,  45) + ylim(-55, -20) + 
  annotate('segment', x = 31, xend = 36, y =-21, yend = -21, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 33.5, y =-20, label = c('5m', '', '', '', '')) 
  
sur_lmap_zoom = sur_lmap + xlim(30,  45) + ylim(-55, -20)  

setwd(out_plot_loc)
pdf('sur_prev_lik_zoom_map.pdf', height = 8, width = 13)
  grid.arrange(sur_pmap_zoom, sur_lmap_zoom, ncol = 1, nrow = 2)  
dev.off()

# reproduction
rep_pred_map_df = data.frame(X = rep_dat$X, Y = rep_dat$Y, year = rep_dat$year, 
  pred = apply(plogis(SPP_rep_pred), MARGIN = 2, FUN = mean))

rep_pmap = ggplot(rep_pred_map_df, aes(X, Y)) + geom_point(aes(color = pred), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) + 
  labs(color = 'prob repo') + 
  annotate('segment', x = 0, xend = 50, y =-250, yend = -250, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 25, y =-240, label = c('50m', '', '', '', '', '')) + 
  facet_grid(.~year, scales="free", space="free") 
  
rep_lik_map_df = data.frame(X = rep_dat$X, Y = rep_dat$Y, year = rep_dat$year, 
  logLik = logLik_bin(post_mat = SPP_rep_pred, obs = rep_dat$rep))

rep_lmap = ggplot(rep_lik_map_df, aes(X, Y)) + geom_point(aes(colour = exp(logLik)), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) + 
  labs(color = 'likliehood') + facet_grid(.~year, scales="free", space="free") 

setwd(out_plot_loc)
pdf('rep_pred_lik_map.pdf', height = 8, width = 13)
  grid.arrange(rep_pmap, rep_lmap, ncol = 1, nrow = 2)  
dev.off()

# make a zomed in verison see how patches look closer in 
rep_pmap_zoom = rep_pmap + xlim(30,  45) + ylim(-55, -20) + 
  annotate('segment', x = 31, xend = 36, y =-21, yend = -21, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 33.5, y =-20, label = c('5m', '', '', '', '', '')) 
  
rep_lmap_zoom = rep_lmap + xlim(30,  45) + ylim(-55, -20)  

setwd(out_plot_loc)
pdf('rep_pred_lik_zoom_map.pdf', height = 8, width = 13)
  grid.arrange(rep_pmap_zoom, rep_lmap_zoom, ncol = 1, nrow = 2)  
dev.off()

# growth
gr_pred_map_df = data.frame(X = gr_dat$X, Y = gr_dat$Y, year = gr_dat$year, 
  pred = apply(SPP_gr_pred, MARGIN = 2, FUN = mean))

gr_pmap = ggplot(gr_pred_map_df, aes(X, Y)) + geom_point(aes(color = pred), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) + 
  labs(color = 'height (cm)') + 
  annotate('segment', x = 0, xend = 50, y =-250, yend = -250, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 25, y =-240, label = c('50m', '', '', '', '')) + 
  facet_grid(.~year, scales="free", space="free") 
  
gr_lik_map_df = data.frame(X = gr_dat$X, Y = gr_dat$Y, year = gr_dat$year, 
  logLik = logLik_cont(post_mu = SPP_gr_pred, obs = gr_dat$height, post_sd = SPP_gr_sigma))

gr_lmap = ggplot(gr_lik_map_df, aes(X, Y)) + geom_point(aes(colour = exp(logLik)), size = 0.5) + xlim(-20, 140) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.key.width = unit(1.5, 'cm'), legend.key.height = unit(1.5, 'cm')) + 
  labs(color = 'likliehood') + facet_grid(.~year, scales="free", space="free") 

setwd(out_plot_loc)
pdf('gr_pred_lik_map.pdf', height = 8, width = 13)
  grid.arrange(gr_pmap, gr_lmap, ncol = 1, nrow = 2)  
dev.off()

# make a zomed in verison see how patches look closer in 
gr_pmap_zoom = gr_pmap + xlim(30,  45) + ylim(-55, -20) + 
  annotate('segment', x = 31, xend = 36, y =-21, yend = -21, color = c(grey(0),grey(0.95),grey(0.95),grey(0.95),grey(0.95))) +
  annotate('text', x = 33.5, y =-20, label = c('5m', '', '', '', '')) 
  
gr_lmap_zoom = gr_lmap + xlim(30,  45) + ylim(-55, -20)  

setwd(out_plot_loc)
pdf('gr_pred_lik_zoom_map.pdf', height = 8, width = 13)
  grid.arrange(gr_pmap_zoom, gr_lmap_zoom, ncol = 1, nrow = 2)  
dev.off()

## todo mapp the spp term for each vital rate, year invariant so only need 3 plots, all need a differnt z scale so cant facet  
# spatial effect plot
# pull the spp from the models
sur_spp = apply(extract_flat(SPP_sur, pars = c('spp')), MARGIN = 2, FUN = mean)
sur_kl = knot_coords2(dat_x = sur_dat$X, dat_y = sur_dat$Y, min_dist = 0.5) 

rep_spp = apply(extract_flat(SPP_rep, pars = c('spp')), MARGIN = 2, FUN = mean)
rep_kl = knot_coords2(dat_x = rep_dat$X, dat_y = rep_dat$Y, min_dist = 0.5) 

gr_spp = apply(extract_flat(SPP_gr, pars = c('spp')), MARGIN = 2, FUN = mean)
gr_kl = knot_coords2(dat_x = gr_dat$X, dat_y = gr_dat$Y, min_dist = 0.5) 

sur_dd = extract_flat(SPP_sur, pars = c('dd_spp_inv'))
rep_dd = extract_flat(SPP_rep, pars = c('dd_spp_inv'))
gr_dd = extract_flat(SPP_gr, pars = c('dd_spp_inv'))


# spp on survial big map
sur_spp_map_df = data.frame(X = sur_kl[,1], Y = sur_kl[,2], spp = sur_spp)

sur_spp_map = ggplot(sur_spp_map_df, aes(X, Y)) + geom_point(aes(color = spp), size = 0.5) + xlim(-20, 140) + ylim(-200, 50) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = "none") + 
  labs(title = ' survival') + scale_colour_gradient2(space="Lab") +
  annotate('segment', x = 0, xend = 50, y =-150, yend = -150, color = grey(0)) +
  annotate('text', x = 25, y =-140, label = '50 m') 

# reproduction big map
rep_spp_map_df = data.frame(X = rep_kl[,1], Y = rep_kl[,2], spp = rep_spp)

rep_spp_map = ggplot(rep_spp_map_df, aes(X, Y)) + geom_point(aes(color = spp), size = 0.5) + xlim(-20, 140) + ylim(-200, 50) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = "none") + 
  labs(title = 'reproduction') + scale_colour_gradient2(space="Lab") + 
  annotate("rect", xmin = 12, xmax = 43, ymin = -58, ymax = 13, colour = grey(0), alpha = 0.0) + 
  annotate("text", x = 12, y = 20, label = 'panels d,e,f', hjust = 0) 
  
# growth big map
gr_spp_map_df = data.frame(X = gr_kl[,1], Y = gr_kl[,2], spp = gr_spp)

gr_spp_map = ggplot(gr_spp_map_df, aes(X, Y)) + geom_point(aes(color = spp), size = 0.5) + xlim(-20, 140) + ylim(-200, 50) +
  theme(panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)), 
    axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(),
    axis.title.y = element_blank(),axis.line.x = element_blank(), axis.line.y = element_blank(), 
    axis.text.y = element_blank(), legend.position = c(0.95, 0.05), legend.justification = c(1,0)) + 
  labs(title = 'growth', color = 'GPP') + scale_colour_gradient2(space = "Lab", breaks = c(-4, 0, 4), labels = c('', '0', ''))
  
# zoomed in areas
sur_spp_zoom_map = sur_spp_map  + xlim(15,  40) + ylim(-55, 10) + labs(title = '') +
  annotate('segment', x = 30, xend = 35, y = 5, yend = 5, color = grey(0)) +
  annotate('text', x = 32.5, y = 2, label = '5 m') 
rep_spp_zoom_map = rep_spp_map  + xlim(15,  40) + ylim(-55, 10) + labs(title = '') 
gr_spp_zoom_map = gr_spp_map  + xlim(15,  40) + ylim(-55, 10) + labs(title = '') + theme(legend.position = 'none') 


# show the distance decay curves 
# plot the decay curves
decay_plotter = function(x, dd_inv){
  dd = 1 / dd_inv
  return(exp(-dd * x))
}

sur_dd_quant = quantile(sur_dd, probs = c(0.025, 0.5, 0.975))

x_ax = seq(0, 40, 1)
sur_dd_df = data.frame(dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)), 
  quant_dd = c(rep(sur_dd_quant[1], length(x_ax)), rep(sur_dd_quant[3], length(x_ax))), 
  med_dd = rep(sur_dd_quant[2], length(x_ax) * 2)) 

sur_dd_df$quant95 = decay_plotter(sur_dd_df$dist, sur_dd_df$quant_dd)  
sur_dd_df$quant50 = decay_plotter(sur_dd_df$dist, sur_dd_df$med_dd)  

sur_dd_plt = ggplot(sur_dd_df, aes(dist, quant95)) + 
  geom_polygon(aes(fill = hcl(h = 180), color = hcl(h = 180)), alpha = 0.5) +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_line(colour = grey(1)),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    legend.position = "none") +
  geom_line(aes(dist, quant50)) + labs(x = 'distance (m)', y = 'correlation') 

rep_dd_quant = quantile(rep_dd, probs = c(0.025, 0.5, 0.975))
rep_dd_df = data.frame(dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)), 
  quant_dd = c(rep(rep_dd_quant[1], length(x_ax)), rep(rep_dd_quant[3], length(x_ax))), 
  med_dd = rep(rep_dd_quant[2], length(x_ax) * 2)) 

rep_dd_df$quant95 = decay_plotter(rep_dd_df$dist, rep_dd_df$quant_dd)  
rep_dd_df$quant50 = decay_plotter(rep_dd_df$dist, rep_dd_df$med_dd)  

rep_dd_plt = ggplot(rep_dd_df, aes(dist, quant95)) + 
  geom_polygon(aes(fill = hcl(h = 180), color = hcl(h = 180)), alpha = 0.5) +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_line(colour = grey(1)),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    legend.position = "none") +
  geom_line(aes(dist, quant50)) + labs(x = 'distance (m)', y = '') 

gr_dd_quant = quantile(gr_dd, probs = c(0.025, 0.5, 0.975))
gr_dd_df = data.frame(dist = c(x_ax, rev(x_ax), x_ax, rev(x_ax)), 
  quant_dd = c(rep(gr_dd_quant[1], length(x_ax)), rep(gr_dd_quant[3], length(x_ax))), 
  med_dd = rep(gr_dd_quant[2], length(x_ax) * 2)) 

gr_dd_df$quant95 = decay_plotter(gr_dd_df$dist, gr_dd_df$quant_dd)  
gr_dd_df$quant50 = decay_plotter(gr_dd_df$dist, gr_dd_df$med_dd)  

gr_dd_plt = ggplot(gr_dd_df, aes(dist, quant95)) + 
  geom_polygon(aes(fill = hcl(h = 180), color = hcl(h = 180)), alpha = 0.5) +
  theme(panel.grid.major.y = element_line(colour = grey(1)), panel.grid.major.x = element_line(colour = grey(1)),
    panel.grid.minor = element_blank(), panel.background = element_rect(fill = grey(0.95)),
    legend.position = "none") +
  geom_line(aes(dist, quant50)) + labs(x = 'distance (m)', y = '') 

setwd(out_plot_loc)
pdf('SPP_plot.pdf', height = 15, width = 15)
  grid.arrange(sur_spp_map, rep_spp_map, gr_spp_map, sur_spp_zoom_map, 
    rep_spp_zoom_map, gr_spp_zoom_map, sur_dd_plt, rep_dd_plt, gr_dd_plt, nrow = 3, ncol = 3)
  grid.text(paste0(letters[1:9], ')'), x = rep(c(0.01, 0.34, 0.67), times = 3), y = rep(c(1, 0.66, 0.33), each = 3), vjust = 1, hjust = 1)  
dev.off()


## parameter plots
# make violine plots of the coeffiecents 
NS_sur_param = data.frame(extract_flat(NS_sur, pars = c('h_ef', 'year_int')))
SPP_sur_param = data.frame(extract_flat(SPP_sur, pars = c('h_ef', 'year_int', 'sigma_spp')))
NS_rep_param = data.frame(extract_flat(NS_rep, pars = c('h_ef', 'year_int')))
SPP_rep_param = data.frame(extract_flat(SPP_rep, pars = c('h_ef', 'year_int', 'sigma_spp')))

NS_sur_long = gather(NS_sur_param, param, effect_size, h_ef:year_int.5.)
#make the slope height effect of an individual of mean height to put on same scale as intercepts 
NS_sur_long[NS_sur_long$param == 'h_ef', 'effect_size'] = NS_sur_long[NS_sur_long$param == 'h_ef', 'effect_size'] * mean(sur_dat$height_prev) 
NS_sur_long$model = 'non-spatial'
NS_sur_long$vrate = 'survival'
#only keep the 99% central part of the distribution
NS_sur_long = NS_sur_long[unlist(tapply(NS_sur_long$effect_size, INDEX = NS_sur_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

SPP_sur_long = gather(SPP_sur_param, param, effect_size, h_ef:sigma_spp)
#make the slope height effect of an individual of mean height to put on same scale as intercepts 
SPP_sur_long[SPP_sur_long$param == 'h_ef', 'effect_size'] = SPP_sur_long[SPP_sur_long$param == 'h_ef', 'effect_size'] * mean(sur_dat$height_prev) 
SPP_sur_long$model = 'SPP'
SPP_sur_long$vrate = 'survival'
#only keep the 99% central part of the distribution
SPP_sur_long = SPP_sur_long[unlist(tapply(SPP_sur_long$effect_size, INDEX = SPP_sur_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

NS_rep_long = gather(NS_rep_param, param, effect_size, h_ef:year_int.6.)
#make the slope height effect of an individual of mean height to put on same scale as intercepts 
NS_rep_long[NS_rep_long$param == 'h_ef', 'effect_size'] = NS_rep_long[NS_rep_long$param == 'h_ef', 'effect_size'] * mean(rep_dat$height) 
# recode years so they start at year 0 for reproduction
for(i in 1:6){
  NS_rep_long[NS_rep_long$param == paste0('year_int.', i, '.'), 'param'] = paste0('year_int.', i - 1, '.')
}
NS_rep_long$model = 'non-spatial'
NS_rep_long$vrate = 'reproduction'
#only keep the 99% central part of the distribution
NS_rep_long = NS_rep_long[unlist(tapply(NS_rep_long$effect_size, INDEX = NS_rep_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

SPP_rep_long = gather(SPP_rep_param, param, effect_size, h_ef:sigma_spp)
#make the slope height effect of an individual of mean height to put on same scale as intercepts 
SPP_rep_long[SPP_rep_long$param == 'h_ef', 'effect_size'] = SPP_rep_long[SPP_rep_long$param == 'h_ef', 'effect_size'] * mean(rep_dat$height) 
# recode years so they start at year 0 for reproduction
for(i in 1:6){
  SPP_rep_long[SPP_rep_long$param == paste0('year_int.', i, '.'), 'param'] = paste0('year_int.', i - 1, '.')
}
SPP_rep_long$model = 'SPP'
SPP_rep_long$vrate = 'reproduction'
SPP_rep_long = SPP_rep_long[unlist(tapply(SPP_rep_long$effect_size, INDEX = SPP_rep_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

par_long = rbind(NS_sur_long, SPP_sur_long, NS_rep_long, SPP_rep_long)

par_sr_plt = ggplot(par_long, aes(param, effect_size)) + geom_abline(intercept = 0, slope = 0) + 
  geom_violin(fill = hcl(h = 180), color = hcl(h = 180, c = 70)) + 
  theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9)),
    strip.background = element_blank(),
    strip.text.y = element_blank(), 
    strip.text.x = element_text(size = 20)) +
  scale_x_discrete(labels = c('height', 'spp_sd', 'y_2002', 'y_2003', 'y_2004', 'y_2005', 'y_2006', 'y_2007')) +
  facet_grid(model ~ vrate, space="free") 
  
# setwd(out_plot_loc)
# pdf('pars_sur_rep_plot.pdf', height = 15, width = 15)
#   par_sr_plt
# dev.off()
 
# make growth seperate since there is a different set of predictors and the effect sizes are on a different scale
NS_gr_param = data.frame(extract_flat(NS_gr, pars = c('b0', 'gr_rate', 'sigma')))
SPP_gr_param = data.frame(extract_flat(SPP_gr, pars = c('b0', 'gr_rate', 'sigma_spp', 'sigma')))
#make the slope height effect of an individual of mean height to put on same scale as intercepts 
NS_gr_param[, 6:10] = NS_gr_param[, 6:10] * mean(gr_dat$height_prev) 
SPP_gr_param[, 6:10] = SPP_gr_param[, 6:10] * mean(gr_dat$height_prev) 

NS_gr_long = gather(NS_gr_param, param, effect_size, b0.1.:sigma)
NS_gr_long$model = 'non-spatial'
NS_gr_long$vrate = 'growth'
NS_gr_long = NS_gr_long[unlist(tapply(NS_gr_long$effect_size, INDEX = NS_gr_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

SPP_gr_long = gather(SPP_gr_param, param, effect_size, b0.1.:sigma)
SPP_gr_long$model = 'SPP'
SPP_gr_long$vrate = 'growth'
SPP_gr_long = SPP_gr_long[unlist(tapply(SPP_gr_long$effect_size, INDEX = SPP_gr_long$param, FUN = function(x){
  quant99 = quantile(x, probs = c(0.01, 0.99))
  return(x > quant99[1] & x < quant99[2])
})), ]

par_long = rbind(NS_gr_long, SPP_gr_long)

par_gr_plt = ggplot(par_long, aes(param, effect_size)) + geom_abline(intercept = 0, slope = 0) + 
  geom_violin(fill = hcl(h = 180), color = hcl(h = 180, c = 70)) + 
  theme(
    panel.grid.major.y = element_line(colour = grey(1)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = grey(0.9)),
    strip.background = element_blank(),
    axis.title.y=element_blank(),
    strip.text.y = element_text(size = 20),
    strip.text.x = element_text(size = 20)) +
  scale_x_discrete(labels = c(sprintf('y_%02d', 3:7), sprintf('gr_%02d', 3:7), 'sigma', 'spp_sd')) +
  facet_grid(model ~ vrate, space="free") 

# setwd(out_plot_loc)
# pdf('pars_gr_plot.pdf', height = 15, width = 10)
#   par_gr_plt
# dev.off()
 

setwd(out_plot_loc)
pdf('pars_vr_plot.pdf', height = 12, width = 21)
  grid.arrange(par_sr_plt, par_gr_plt, nrow = 1, widths = c(2, 1.2))
dev.off()
