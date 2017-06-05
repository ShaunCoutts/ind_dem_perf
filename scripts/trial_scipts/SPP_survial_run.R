## Spatial prediction proccess for survial, find spatial effect at course scale 
library(R2jags)

# file locations for Zhu's Mac
dat_loc = '/Users/shauncoutts/ind_perf_pred/Data' 
out_ob_loc = '/Users/shauncoutts/ind_perf_pred/objects'
code_loc = '/Users/shauncoutts/ind_perf_pred/scripts'

# DATA PREP
setwd(dat_loc)

vr_loc_gen = read.csv('vr_loc_gen_postburn.csv', header = TRUE, stringsAsFactor = FALSE)
vr_loc_gen = vr_loc_gen[!is.na(vr_loc_gen$X), ]

# start with survival
sur_dat = vr_loc_gen[!is.na(vr_loc_gen$sur), c('uID', 'uLoc', 'year','sur', 'height', 'height_prev',
  'X', 'Y', 'MNR', 'MDH1', 'MDH3', 'X6PGD', 'IDH')]
# take out data from the first year observed since nothing can be observed dead in the first year
first_year = sapply(seq_along(sur_dat$year), FUN = function(x){

  min_year_group = min(sur_dat$year[sur_dat$uID == sur_dat$uID[x]])
  return(ifelse(sur_dat$year[x] == min_year_group, FALSE, TRUE))

})
sur_dat = sur_dat[first_year, ]

sur_dat$height_mod = NA

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
# assumed order: new plants emerege -> observation -> mortality -> growth -> observation
sur_dat = sur_dat[!is.na(sur_dat$height_prev), ]
sur_dat$height_mod = sur_dat$height_prev - mean(sur_dat$height_prev)

# read in the jags model and helper functions 
setwd(code_loc)
source('sur_model.R')
source('dist_neigh_setup.R')

# create knots at 10m grid, cut out all those more than 20m from nearerst data point
knots = knot_coords(x_min = min(sur_dat$X), x_max = max(sur_dat$X), y_min = min(sur_dat$Y), 
  y_max = max(sur_dat$Y), res = 10, dat_x = sur_dat$X, dat_y = sur_dat$Y, iso_dist = 20)

dm_knots = distance_matrix(knots[,'x_cord'], knots[, 'y_cord'])
dm_knots_obs = distance_matrix2(x1 = sur_dat$X, y1 = sur_dat$Y, x2 = knots[,'x_cord'], 
  y2 = knots[, 'y_cord'])

jags_data <- list(sur = sur_dat$sur, year = sur_dat$year_num, num_years = max(sur_dat$year_num), 
  num_obs = length(sur_dat$sur), height = sur_dat$height_mod, num_knots = dim(dm_knots)[1], 
  s_knot_dists = dm_knots, s_knot_ob_dist = dm_knots_obs)
  
inital <- function() list(b_height = runif(1, -20, 20), 
  int_year = runif(jags_data$num_years, -20, 20))
#  gap_int = runif(jags_data$num_gaps, -20, 20)) 

# Parameters monitored
monitored<- c('int_year', 'b_height', 'dd_spp', 'z_lp', 'sigma_spp')

#call to JAGS
setwd(code_loc)
sur_spp_mod <- jags.parallel(data = jags_data, inits = inital, parameters.to.save = monitored, 
  "sur_model_GSPP.jags", n.chains = 3, n.thin = 50, n.iter = 50000, n.burnin = 10000, 
  working.directory = getwd(), export_obj_names = c('jags_data', 'inital', 'monitored'))

setwd(out_ob_loc)
save(sur_spp_mod, file = 'sur_spp10NG2_mod.Rdata')
