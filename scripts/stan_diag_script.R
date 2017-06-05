## script to run diagnostics on stan objects
library(rstan)
#install.packages('extrafont', repos = "https://cloud.r-project.org")
library(extrafont)
#font_import()
#loadfonts()

hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
} 

dat_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/Data' 
out_plot_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/plots'
out_ob_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/output/objects'
code_loc = '/home/shauncoutts/Dropbox/projects/ind_perform_correlation/scripts'

## survival diagnosis
setwd(out_ob_loc)
ob_name = load('NS_sur_stan.Rdata')
setwd(out_plot_loc)
pdf('NS_sur_diagnosis.pdf', width = 20, height = 20)
  stan_trace(NS_sur_stan, inc_warmup = FALSE, pars = c('h_ef', 'year_int'))
  stan_trace(NS_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(NS_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(NS_sur_stan, pars = c('h_ef', 'year_int'))
  stan_ess(NS_sur_stan, pars = c('h_ef', 'year_int'), binwidth = 0.05)
  hist_treedepth(NS_sur_stan)
  pairs(NS_sur_stan, pars = c('h_ef', 'year_int'))
dev.off()
# looks perfect

setwd(out_ob_loc)
ob_name = load('SPP_NK_sur_stan.Rdata')
setwd(out_plot_loc)
pdf('SPP_NK_sur_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_NK_sur_stan, inc_warmup = FALSE, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = FALSE, pars = paste0('spp[', 1:100, ']'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = FALSE, pars = paste0('spp[', 400:500, ']'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP_NK_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP_NK_sur_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP_NK_sur_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'), binwidth = 0.05)
  hist_treedepth(SPP_NK_sur_stan)
  pairs(SPP_NK_sur_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
dev.off()
# all looks really nice and converged, effective sample size is a bit small, but at least 1600 for all parameters (from 8000 chains). 
# some correlation between year effects, model is only moderatly identifiable. 

setwd(out_ob_loc)
ob_name = load('SPP_25_sur_stan.Rdata')
setwd(out_plot_loc)
pdf('SPP_25_sur_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_sur_stan, inc_warmup = FALSE, pars = c('b', 'sd_1', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP_sur_stan, inc_warmup = FALSE, pars = c('spp'))
  stan_trace(SPP_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  hist_treedepth(SPP_sur_stan)
  pairs(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'sigma_spp', 'dd_spp_inv'))
dev.off()
# all looks really nice and converged, sigma_spp and dd_spp_inv have pretty low efefctive sample sizes, but still over 400 and 800 respectivley
stan_ess(SPP_sur_stan, pars = c('b', 'sd_1'), binwidth = 0.025, xlim = c(0, 1))
stan_ess(SPP_sur_stan, pars = c('year_int'))
stan_ess(SPP_sur_stan, pars = c('sigma_spp', 'dd_spp_inv'))


setwd(out_ob_loc)
ob_name = load('SPP_25NK_sur_stan.Rdata')
setwd(out_plot_loc)
pdf('SPP_25NK_sur_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_sur_stan, inc_warmup = FALSE, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv25', 'sigma_spp25', 
    'dd_spp_invNN', 'sigma_sppNN'))
  stan_trace(SPP_sur_stan, inc_warmup = FALSE, pars = c('spp_int25'))
  stan_trace(SPP_sur_stan, inc_warmup = FALSE, pars = paste0('spp_NN[', 1:100, ']'))
  stan_trace(SPP_sur_stan, inc_warmup = FALSE, pars = paste0('spp_NN[', 400:500, ']'))
  stan_trace(SPP_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP_sur_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'sigma_spp25', 'dd_spp_inv25', 'spp_int25', 'dd_spp_invNN', 'sigma_sppNN', 'spp_NN'))
  stan_ess(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'sigma_spp25', 'dd_spp_inv25', 'dd_spp_invNN', 'sigma_sppNN'), binwidth = 0.01)
  hist_treedepth(SPP_sur_stan)
  pairs(SPP_sur_stan, pars = c('b', 'sd_1', 'year_int', 'dd_spp_inv25', 'sigma_spp25', 
    'dd_spp_invNN', 'sigma_sppNN'))
dev.off()
# all looks really nice and converged, sigma_spp25 and dd_spp_inv25 have very low efefctive sample sizes, 100-300
stan_ess(SPP_sur_stan, pars = c('b', 'sd_1'), binwidth = 0.025)
stan_ess(SPP_sur_stan, pars = c('year_int'))
stan_ess(SPP_sur_stan, pars = c('sigma_spp25', 'dd_spp_inv25'))


## reproduction diagnosis setwd(out_ob_loc)
setwd(out_ob_loc)
ob_name = load('SPP_NS_rep_stan.Rdata')
setwd(out_plot_loc)
pdf('NS_rep_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_NS_rep_stan, inc_warmup = FALSE, pars = c('h_ef', 'year_int'))
  stan_trace(SPP_NS_rep_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP_NS_rep_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP_NS_rep_stan, pars = c('h_ef', 'year_int'))
  stan_ess(SPP_NS_rep_stan, pars = c('h_ef', 'year_int'))
  hist_treedepth(SPP_NS_rep_stan)
  pairs(SPP_NS_rep_stan, pars = c('h_ef', 'year_int'))
dev.off()
# perfect, new parameterisation works much better

setwd(out_ob_loc)
ob_name = load('SPP_25_rep_stan.Rdata')
setwd(out_plot_loc)
pdf('SPP_25_rep_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP25_rep_stan, inc_warmup = FALSE, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP25_rep_stan, inc_warmup = FALSE, pars = c('spp'))
  stan_trace(SPP25_rep_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP25_rep_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP25_rep_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP25_rep_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  hist_treedepth(SPP25_rep_stan)
  pairs(SPP25_rep_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
dev.off()
# a handful of divergent samples on dd_spp_inv, with big spikes, so try re-run with higher adapt_delta
# also lots of auto-correlation in the spatial effect chains, but still 400 effective samples, and all 
# well converged
# some divergence problems
# Trying an alternative parameterisation that should work better with more cleanly identifiable parameters 
# for the year effect

setwd(out_ob_loc)
ob_name = load('SPP_NK_rep_stan_ap.Rdata')
setwd(out_plot_loc)
pdf('SPP_NK_rep_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_NK_rep_stan, inc_warmup = FALSE, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP_NK_rep_stan, inc_warmup = FALSE, pars = c(paste0('spp[', 1:100, ']')))
  stan_trace(SPP_NK_rep_stan, inc_warmup = FALSE, pars = c(paste0('spp[', 400:500, ']')))
  stan_trace(SPP_NK_rep_stan, inc_warmup = TRUE, pars = paste0('eta[', 1:100, ']'))
  stan_trace(SPP_NK_rep_stan, inc_warmup = TRUE, pars = paste0('eta[', 1000:1100, ']'))
  stan_rhat(SPP_NK_rep_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP_NK_rep_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
  hist_treedepth(SPP_NK_rep_stan)
  pairs(SPP_NK_rep_stan, pars = c('h_ef', 'year_int', 'sigma_spp', 'dd_spp_inv'))
dev.off()

# This one looks very nice, some with autocorrelated chains but more than 1600 effective samples
# some correlation in samples between parameters but is better than all other parameterisations I have tried 
# Could be the model is only weakly identifiable

## growth models 
# non-spatial model 
setwd(out_ob_loc)
ob_name = load('NS_gr_stan.Rdata')
setwd(out_plot_loc)
pdf('NS_gr_diagnosis.pdf', width = 20, height = 20)
  stan_trace(NS_gr_stan, inc_warmup = FALSE, pars = c('b0', 'gr_rate', 'sigma'))
  stan_trace(NS_gr_stan, inc_warmup = FALSE, pars = paste0('mu[', 1:100, ']'))
  stan_rhat(NS_gr_stan, pars = c('b0', 'gr_rate', 'sigma'))
  stan_ess(NS_gr_stan, pars = c('b0', 'gr_rate', 'sigma'))
  hist_treedepth(NS_gr_stan)
  pairs(NS_gr_stan, pars = c('b0', 'gr_rate', 'sigma'))
dev.off()
# All looks very good 

setwd(out_ob_loc)
ob_name = load('SPP_NK_gr_stan.Rdata')
setwd(out_plot_loc)
pdf('SPP_NK_gr_diagnosis.pdf', width = 20, height = 20)
  stan_trace(SPP_gr_stan, inc_warmup = FALSE, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv'))
  stan_trace(SPP_gr_stan, inc_warmup = FALSE, pars = c(paste0('spp[', 1:50, ']')))
  stan_trace(SPP_gr_stan, inc_warmup = FALSE, pars = c(paste0('spp[', 400:450, ']')))
  stan_trace(SPP_gr_stan, inc_warmup = TRUE, pars = paste0('mu[', 1:50, ']'))
  stan_trace(SPP_gr_stan, inc_warmup = TRUE, pars = paste0('mu[', 1000:1050, ']'))
  stan_rhat(SPP_gr_stan, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv', 'spp'))
  stan_ess(SPP_gr_stan, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv'))
  hist_treedepth(SPP_gr_stan)
  pairs(SPP_gr_stan, pars = c('b0', 'gr_rate', 'sigma', 'sigma_spp', 'dd_spp_inv'))
dev.off()
# all looks good, eveything converged, some auto-correlation but all chains have at least 2000 effective samples 


## seed number prediction
setwd(out_ob_loc)
ob_name = load('seed_num_stan.Rdata')
setwd(out_plot_loc)
pdf('seed_num_diagnosis.pdf', width = 20, height = 20)
  stan_trace(seed_num_stan, inc_warmup = FALSE, pars = c('glob_int', 'h_ef', 'site_int'))
  stan_rhat(seed_num_stan, pars = c('glob_int', 'h_ef', 'site_int'))
  stan_ess(seed_num_stan, pars = c('glob_int', 'h_ef', 'site_int'))
  hist_treedepth(seed_num_stan)
  pairs(seed_num_stan, pars = c('glob_int', 'h_ef', 'site_int'))
dev.off()
# all looks very good


setwd(out_ob_loc)
ob_name = load('fruit_num_stan.Rdata')
setwd(out_plot_loc)
pdf('fruit_num_diagnosis.pdf', width = 20, height = 20)
  stan_trace(fruit_num_stan, inc_warmup = FALSE, pars = c('h_ef', 'site_int', 'sigma'))
  stan_rhat(fruit_num_stan, pars = c('h_ef', 'site_int', 'sigma'))
  stan_ess(fruit_num_stan, pars = c('h_ef', 'site_int', 'sigma'))
  stan_ess(fruit_num_stan, pars = c('site_int'))
  hist_treedepth(fruit_num_stan)
  pairs(fruit_num_stan, pars = c('h_ef', 'site_int', 'sigma'))
dev.off()
# all looks very good, 


  

