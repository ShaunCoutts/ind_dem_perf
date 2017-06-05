# jags model for survival 

sink("sur_model.jags")
  cat('
    model{
      
      # PRIORS 
      # year effect
       for(i in 1:num_years){
	 int_year[i] ~ dnorm(0, 0.0001) 
       }
      
      # slope on height
      b_height ~ dnorm(0, 0.0001)
      
      # random intercept for each gap
      for(g in 1:num_gaps){
	gap_int[g] ~ dnorm(0, tua_gap)
      }
      tua_gap ~ dgamma(0.001, 0.001)
      sigma_gap = 1 / tua_gap
      
      for(n in 1:num_obs){
      
	sur[n] ~ dbern(p[n])
	p[n] <- max(0.000001, min(0.99999999, p_raw[n]))
	p_raw[n] <- 1 / (1 + exp(-z_lp[n]))	
	z_lp[n] <- int_year[year[n]] + b_height * height[n] + gap_int[gapID[n]] 
      
      }

    }',fill = TRUE)
sink()

# model with density dependence, very simple model that combines distance and size
sink("sur_model_den.jags")
  cat('
    model{
      
      # PRIORS 
      # year effect
       for(i in 1:num_years){
	 int_year[i] ~ dnorm(0, 0.0001) 
       }
      
      # slopes
      b_height ~ dnorm(0, 0.0001)
      b_den ~ dnorm(0, 0.0001)
      
      # DISTANCE EFFECTS
      dd_den ~ dunif(0.0000001, 5)
      
      # find the summed exp(-dd_den * dist) * height
      # var dist_weight[num_obs, max_num_cont]
      for(n in 1:num_obs){
	for(col in 1:num_cont[n]){
	
	  dist_weight[n, col] <- exp(-dd_den * dists[n, col])
	
	}
      
	dist_height[n] <- sum(dist_weight[n, 1:num_cont[n]] * height_mat[n, 1:num_cont[n]])  
      
      }
      
      #LIKLIEHOOD
      for(n in 1:num_obs){
      
	sur[n] ~ dbern(p[n])
	p[n] <- max(0.000001, min(0.99999999, p_raw[n]))
	p_raw[n] <- 1 / (1 + exp(-z_lp[n]))
	z_lp[n] <- int_year[year[n]] + b_height * height[n] + b_den * dist_height[n] 
      
      }

    }',fill = TRUE)
sink()

# spatial guassian proccess model following Viana et al. 2013

sink("sur_model_GSPP.jags")
  cat('
    model{
      
      # PRIORS 
      # year effect
      for(i in 1:num_years){
	int_year[i] ~ dnorm(0, 0.0001) 
      }
      
      # slopes
      b_height ~ dnorm(0, 0.0001)
       
      # DISTANCE EFFECTS
      dd_spp ~ dunif(0.00001, 5)
      sigma_spp ~ dunif(0.0001, 10)
      #sigma_spp = 1 / sigma_spp_inv
      #sigma_spp_inv ~ dgamma(0.001, 0.001)


      # spatial predictive proccess 
      spp ~ dmnorm(mu_spp, C_spp_inv)
      C_spp_inv = inverse(C_spp)
      for(k in 1:num_knots){
 	
 	mu_spp[k] = 0
 	C_spp[k, k] = sigma_spp
 	  
 	for(j in 1:(k - 1)){
 	  
 	  C_spp[k, j] = sigma_spp * exp(-(dd_spp * s_knot_dists[k, j])) 
 	  C_spp[j, k] = C_spp[k, j]
 	  
 	}
      }
      # Interpolate back from the knot points to the observed data 
      for(n in 1:num_obs){
	for(k in 1:num_knots){
	  C_spp_s[n, k] = sigma_spp * exp(-(dd_spp * s_knot_ob_dist[n, k]))
	}
      }
      spp_interp = C_spp_s %*% C_spp_inv %*% spp
 
#       # random intercept for each gap
#       for(g in 1:num_gaps){
# 	gap_int[g] ~ dnorm(0, tua_gap)
#       }
#       tua_gap ~ dgamma(0.001, 0.001)
#       sigma_gap = 1 / tua_gap
      
      #LIKLIEHOOD
      for(n in 1:num_obs){
      
	sur[n] ~ dbern(p[n])
	p[n] <- max(0.000001, min(0.99999999, p_raw[n]))
	p_raw[n] <- 1 / (1 + exp(-z_lp[n]))
	z_lp[n] <- int_year[year[n]] + b_height * height[n] + spp_interp[n] #gap_int[gapID[n]]   
      
      }

    }',fill = TRUE)
sink()

# The above is very slow and very memory hungry, so createa window version that only does the interperlation 
# from knots that are with some window of the target observation (say 20m). This will also allow the use of 
# sparser matricies, so may be lighter on memory. At high knot resolution (0.5m or higher) there are more
# knots than there are data points and it is better to just cut out the interperlation step and 
# get the covariance matrix for each data point. Note that care is needed in constucting the 
# sparse matrix. the array of indicies cannot include diagonal, and to make the number of computation 
# manageable only the lower triangle is used (then reflected when constructing the co-variance matrix)
sink("sur_model_SPP_window.jags")
  cat('
    model{
      
      ## PRIORS 
      # year effect
      for(i in 1:num_years){
	int_year[i] ~ dnorm(0, 0.0001) 
      }
#       
#       # slopes
#       b_height ~ dnorm(0, 0.0001)
#       b_den ~ dnorm(0, 0.0001)
#       
#       ## DISTANCE EFFECTS
#       dd_den ~ dunif(0.0000001, 5)
      dd_spp ~ dunif(0.00001, 5)
      sigma_spp = 1 / sigma_spp_inv
      sigma_spp_inv ~ dgamma(0.001, 0.001)

#       # find the summed exp(-dd_den * dist) * height
#       # var dist_weight[num_obs, max_num_cont]
#       for(n in 1:num_obs){
# 	for(col in 1:num_cont[n]){
# 	
# 	  dist_weight[n, col] <- exp(-dd_den * dists[n, col])
# 	
# 	}
#       
# 	dist_height[n] <- sum(dist_weight[n, 1:num_cont[n]] * height_mat[n, 1:num_cont[n]])  
#       
#       }

      ## SPP (spatial predictive proccess) 
      spp ~ dmnorm(mu_spp, C_spp_inv)
      C_spp_inv = inverse(C_spp)
      
      # fill the diagonal
      for(k in 1:num_knots){
      
	mu_spp[k] = 0
	C_spp[k, k] = sigma_spp
	
      }
      
      # fill in the off-diag 0 entries 
      # assumes if distance greater than win_r then dist = win_r, and so covar = sigma * exp(-dd_spp * win_r))
      out_win_covar = sigma_spp * exp(-dd_spp * win_r)
      for(j in 1:knot_zero_ent){
      
	C_spp[k_zero_ind[1, j] , k_zero_ind[2, j]] = out_win_covar
	C_spp[k_zero_ind[2, j] , k_zero_ind[1, j]] = out_win_covar
      
      }
	
      # fill non-zero elements (those within win_r of observation k)
      for(j in 1:knot_ent){
      
	C_spp[k_non_zero_ind[1, j], k_non_zero_ind[2, j]] = sigma_spp * exp(-dd_spp * knot_dists[j]) 
	C_spp[k_non_zero_ind[2, j], k_non_zero_ind[1, j]] = C_spp[k_non_zero_ind[1, j], k_non_zero_ind[2, j]]
      
      }
 
      ## Interpolate back from the knot points to the observed data 
      # fill zero elements 
      for(j in 1:obs_knot_zero){
      
	C_spp_s[ob_zero_ind[1, j], ob_zero_ind[2, j]] = out_win_covar
	
      }
      
      # fill non-zero elements
      for(j in 1:obs_knots){
      
	C_spp_s[ob_non_zero_ind[1, j], ob_non_zero_ind[2, j]] = sigma_spp * exp(-(dd_spp * s_knot_ob_dist[j]))
	
      }
      
      spp_interp = C_spp_s %*% C_spp_inv %*% spp
      
      #LIKLIEHOOD
      for(n in 1:num_obs){
      
	sur[n] ~ dbern(p[n])
	p[n] <- max(0.000001, min(0.99999999, p_raw[n]))
	p_raw[n] <- 1 / (1 + exp(-lp[n]))
	lp[n] <- int_year[year[n]] + spp_interp[n] #+ b_height * height[n] + b_den * dist_height[n] 
      
      }

    }',fill = TRUE)
sink()


# spatial auto-correlation model model
sink("sur_model_sac.jags")
  cat('
    model{
      
      # PRIORS 
      # year effect
      for(i in 1:num_years){
	int_year[i] ~ dnorm(0, 0.0001) 
      }
       
      # slopes
      b_sac ~ dnorm(0, 0.0001)
      b_height ~ dnorm(0, 0.0001)
      
      # random intercept for each gap
#       for(g in 1:num_gaps){
# 	gap_int[g] ~ dnorm(0, tua_gap)
#       }
#       tua_gap ~ dgamma(0.001, 0.001)
#       sigma_gap = 1 / tua_gap
      
      # fit the non-spatial model to get the residuals 
      for(n in 1:num_obs){
      
	r_ns[n] = sur1[n] - p_ns[n] 
	p_ns[n] = 1 / (1 + exp(-z_lp_ns[n]))
	z_lp_ns[n] = int_year[year[n]] + b_height * height[n]# + gap_int[gapID[n]]  
      
      }
      
      # DISTANCE EFFECTS
      # fit a model to the residuals
      dd_sac ~ dunif(0.0001, 5)

      # spatial predictor term 
      for(n in 1:num_obs){
	for(j in 1:num_neigh[n]){
	
	  w[n, j] = exp(-dd_sac * dist_neigh[n, j])
	  sac_raw[n, j] = r_ns[inds_neigh[n, j]] * w[n, j]
	
	}
	
	sac[n] = sum(sac_raw[n, 1:num_neigh[n]]) / sum(w[n, 1:num_neigh[n]]) 
      }
      
      #LIKLIEHOOD
      for(n in 1:num_obs){
      
	sur2[n] ~ dbern(p[n])
	p[n] = max(0.000001, min(0.99999999, p_raw[n]))
	p_raw[n] = 1 / (1 + exp(-z_lp[n]))
	z_lp[n] = z_lp_ns[n] + b_sac * sac[n] 
      
      }

    }',fill = TRUE)
sink()
