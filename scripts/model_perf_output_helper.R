# set of functions to work out diganostics and model performance measures 

library(colorspace)
library(pROC)

# for binary reponse use ROC curve, pred is a length(obs) by N_samples matrix giving the uncertianity 
# around the prediction for observation. Number samples gives the number of times the prediction chains 
# are sampled to capture the uncertianty. True Positive and False Positive is calculated at each threshold
# level to generate the curve.
ROC_curve = function(obs, pred, num_resamp, th_res){
  
  thresholds = seq(0, 1, th_res)
  num_MCMC = dim(pred)[1]
  num_obs = dim(pred)[2]
  num_pos = sum(obs)
  num_neg = sum(!obs)
  
  rs_row = sample.int(dim(pred)[1], num_resamp)
  pred_obs = data.frame(ob = rep(obs, each = num_resamp), 
    pred_lp = as.vector(apply(pred, MARGIN = 2, FUN = function(x) x[rs_row])))
   
  # do this as gg plot density plot object 
   
  ROC = data.frame(ts_ID =rep(1:num_samp, each  = length(thresholds)), true_pos = NA, false_pos = NA,
    thres = NA)  
  count = 1
  #array(NA, dim = c(num_samp, length(thresholds), 3))
  row_inds = sample.int(num_MCMC, num_samp)
  for(i in 1:num_samp){
    #draw a set of random predictions from the chains
    pred_samp = pred[row_inds[i], ]
    # put the linear predictor on a probability scale 
    pred_prob = 1 / (1 + exp(-pred_samp))
    
    for(k in seq_along(thresholds)){
      # find all predicted positives for each threshold 
      pred_pos = pred_prob >= thresholds[k]
      ROC$true_pos[count] = sum(obs[pred_pos]) / num_pos 
      ROC$false_pos[count] = sum(!obs[pred_pos]) / num_neg 
      ROC$thres[count] = thresholds[k]
      count = count + 1
    }
  }
  return(ROC)
}

## plot the ROC curve from the array returned by ROC_curve
ROC_plot = function(ROC_df, labels = labs(list(title = 'title', x = 'x_lab', y = 'y_lab')),
  AUC_text){
 
  lin_col = hcl(h = 100, c = 35, l = 85, alpha = 1)
  plt = ggplot(ROC_df, aes(x = false_pos, y = true_pos, group = ts_ID, color = lin_col)) +
    geom_line() + theme(legend.position = "none") + xlim(0, 1) + ylim(0, 1) + labels +
    geom_abline(intercept = 0, slope = 1)
    
  return(plt)

}

#calculate the auc for the ROC curve, take a matrix for prediction, one set of predictions for each sample.
# pred is a num_samples matrix by nrow(obs)  
auc_dist = function(obs, pred, samp_num){
  
  auc_hold = numeric(samp_num)
  rs_row = sample.int(dim(pred)[1], samp_num)
  for(i in 1:samp_num) auc_hold[i] = auc(obs, pred[rs_row[i], ])
  
  return(auc_hold)
  
}

## observed vs predicted for binary response, use 2Dhist to give density and jitter to visulise uncertianty
obs_v_pred_bin = function(obs, lp_samp, num_resamp, jitter_am, labels = labs(list(title = 'title', x = 'x_lab', y = 'y_lab'))){
  
  obs_jit = jitter(obs, amount = jitter_am)
  # now need to expand each obs so there are num_resamp predictions for each obs, the -1 dummy data is to trick the image function
  # into coloring the whole mat
  rs_row = sample.int(dim(lp_samp)[1], num_resamp)
  
  pred_obs = data.frame(ob = rep(obs_jit, each = num_resamp), 
    pred_lp = as.vector(apply(lp_samp, MARGIN = 2, FUN = function(x) x[rs_row])))
  pred_obs$p_pred = 1 / (1 + exp(-pred_obs$pred_lp))
  
  p <- ggplot(pred_obs, aes(p_pred, ob))
  h3 <- p + stat_bin2d() + theme(legend.position = "none") + xlim(0, 1) + ylim(0, 1) + labels
  h3
  
}

## wrapper on rstan extract to flatten but not permute the output of extract
extract_flat = function(rstan_ob, ...){
  
  ext = extract(rstan_ob, permuted = FALSE, inc_warmup = FALSE, ...)
  ext_flat = sapply(1:dim(ext)[3], FUN = function(x) as.numeric(ext[, , x]))
  colnames(ext_flat) <- colnames(ext[1, , ])
  
  return(ext_flat)
}

## calculate the logLik from a set of postieriers of the linear predictor of a model
# post_mat is a num_samp by num_obs matrix of the postierier on the logit scale (plogis() to turn to prob)
logLik_bin = function(post_mat, obs){

  # turn the obs 1 into 0's and 0's into 1
  obs_flip = 1 - obs
  
  pred_prob = plogis(post_mat)
  
  return(log(apply(abs(obs_flip - t(pred_prob)), MARGIN = 1, FUN = mean)))

}

logLik_cont = function(post_mu, post_sd, obs){
  
  lik = numeric(length(obs))
  
  for(i in seq_along(obs)){
    lik[i] = mean(dnorm(obs[i], post_mu[,i], post_sd)) 
  }
  
  return(log(lik))

}
