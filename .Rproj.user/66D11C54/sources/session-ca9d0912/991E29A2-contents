####functions for the simulation
simu_ESDGENlm_alphaTest = function(parallel_idx = 1,
                                   num_sig = 10,
                                   k = 10,
                                   alpha = 0.05,
                                   total_audio = 50,
                                   npatient=80,
                                   vec_lam)
{
  require(glmnet)
  vec_res = rep(NA,2)
  names(vec_res)=c("GENlm","Lasso")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 66.95,
      medium.mean        = 66.95,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig,
      seed_idx=parallel_idx
    )
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(total_audio + 4)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(total_audio + 4)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(total_audio + 4), 5:(total_audio + 4)] * num_sigma
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam[1:k])
  vec_res[1] = length(vec_idx)
  # browser()
  # obj_lasso = cv.glmnet(
  #   x = as.matrix(ll_dat$df.XY[, 1:(4+total_audio)]),
  #   y = ll_dat$df.XY$Y,
  #   lambda = NULL,
  #   family = "gaussian",
  #   penalty.factor = c(rep(0, 4), rep(1, total_audio))
  # )
  # vec_outlier = which((coef(obj_lasso)[6:(total_audio+5)]!=0))
  # vec_res[2] = length(vec_outlier)
  
  return(vec_res)
}

simu_ESDGENlm_alphaTest_mean = function(parallel_idx = 1,
                                   num_sig = 10,
                                   k = 10,
                                   alpha = 0.05,
                                   total_audio = 50,
                                   npatient=80,
                                   vec_lam)
{
  # browser()
  require(glmnet)
  vec_res = rep(NA,2)
  names(vec_res)=c("GENlm","Lasso")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 66.95,
      medium.mean        = 66.95,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig,
      seed_idx=parallel_idx
    )
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(total_audio + 4)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(total_audio + 4)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(total_audio + 4), 5:(total_audio + 4)] * num_sigma
  ll_Ri = Ri_ESDChi_mean(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam[1:k])
  vec_res[1] = length(vec_idx)
  # browser()
  # obj_lasso = cv.glmnet(
  #   x = as.matrix(ll_dat$df.XY[, 1:(4+total_audio)]),
  #   y = ll_dat$df.XY$Y,
  #   lambda = NULL,
  #   family = "gaussian",
  #   penalty.factor = c(rep(0, 4), rep(1, total_audio))
  # )
  # vec_outlier = which((coef(obj_lasso)[6:(total_audio+5)]!=0))
  # vec_res[2] = length(vec_outlier)
  
  return(vec_res)
}


###model misspecify, missing age2
simu_ESDGENlm_alphaTest_miscov = function(parallel_idx = 1,
                                   num_sig = 10,
                                   k = 10,
                                   alpha = 0.05,
                                   total_audio = 50,
                                   npatient=80,
                                   vec_lam,
                                   num_cov = 3)
{
  require(glmnet)
  vec_res = rep(NA,2)
  names(vec_res)=c("GENlm","Lasso")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 66.95,
      medium.mean        = 66.95,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig,
      seed_idx=parallel_idx
    )
  ll_dat$coef = ll_dat$coef[-2]
  ll_dat$df.XY = ll_dat$df.XY[,-2]
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[(num_cov+1):(total_audio + num_cov)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(total_audio + num_cov)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[
    (num_cov+1):(total_audio + num_cov), (num_cov+1):(total_audio + num_cov)] * num_sigma
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam[1:k])
  vec_res[1] = length(vec_idx)

  return(vec_res)
}
############################################################################
##functions for power Test
############################################################################
simu_ESDGENlm_PowerTest = function(parallel_idx = 1,
                                      num_sig = 10,
                                      k = 10,
                                      alpha = 0.05,
                                      total_audio = 50,
                                      npatient=80,
                                      vec_lam)
{
  require(glmnet)
  warning("Critical value specific version! use with caution!")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig, seed_idx=parallel_idx
    )
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(4 + total_audio)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(4 + total_audio)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(4 + total_audio), 5:(4 + total_audio)] * num_sigma
  
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam)

  ll_res=list()
  if (length(vec_idx) == 0)
  {
    ll_res$GENlm=c(0)
  } else{
    ll_res$GENlm = ll_Ri$Ic[[max(vec_idx)]]
  }
  
  # obj_lasso = cv.glmnet(
  #   x = as.matrix(ll_dat$df.XY[, 1:(4+total_audio)]),
  #   y = ll_dat$df.XY$Y,
  #   lambda = NULL,
  #   family = "gaussian",
  #   penalty.factor = c(rep(0, 4), rep(1, total_audio))
  # )
  # vec_outlier = which((coef(obj_lasso)[6:(total_audio+5)]!=0))
  # if(length(vec_outlier)==0) {
  #   ll_res$Lasso=c(0)
  # } else{
  #   ll_res$Lasso = vec_outlier
  # }
  ll_res
}

simu_ESDGENlm_PowerTest_mean = function(parallel_idx = 1,
                                   num_sig = 10,
                                   k = 10,
                                   alpha = 0.05,
                                   total_audio = 50,
                                   npatient=80,
                                   vec_lam)
{
  require(glmnet)
  warning("Critical value specific version! use with caution!")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig, seed_idx=parallel_idx
    )
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(4 + total_audio)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(4 + total_audio)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(4 + total_audio), 5:(4 + total_audio)] * num_sigma
  
  ll_Ri = Ri_ESDChi_mean(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam)
  
  ll_res=list()
  if (length(vec_idx) == 0)
  {
    ll_res$GENlm=c(0)
  } else{
    ll_res$GENlm = ll_Ri$Ic[[max(vec_idx)]]
  }
  
  # obj_lasso = cv.glmnet(
  #   x = as.matrix(ll_dat$df.XY[, 1:(4+total_audio)]),
  #   y = ll_dat$df.XY$Y,
  #   lambda = NULL,
  #   family = "gaussian",
  #   penalty.factor = c(rep(0, 4), rep(1, total_audio))
  # )
  # vec_outlier = which((coef(obj_lasso)[6:(total_audio+5)]!=0))
  # if(length(vec_outlier)==0) {
  #   ll_res$Lasso=c(0)
  # } else{
  #   ll_res$Lasso = vec_outlier
  # }
  ll_res
}

simu_ESDGENlm_PowerTest_miscov = function(parallel_idx = 1,
                                   num_sig = 10,
                                   k = 10,
                                   alpha = 0.05,
                                   total_audio = 50,
                                   npatient=80,
                                   vec_lam,
                                   num_cov = 3)
{
  require(glmnet)
  warning("Critical value specific version! use with caution!")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig, seed_idx=parallel_idx
    )
  ll_dat$coef = ll_dat$coef[-2]
  ll_dat$df.XY = ll_dat$df.XY[,-2]
  
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[(num_cov+1):(num_cov + total_audio)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(num_cov + total_audio)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[
    (num_cov+1):(num_cov + total_audio), (num_cov+1):(num_cov + total_audio)] * num_sigma
  
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam)
  
  ll_res=list()
  if (length(vec_idx) == 0)
  {
    ll_res$GENlm=c(0)
  } else{
    ll_res$GENlm = ll_Ri$Ic[[max(vec_idx)]]
  }

  ll_res
}

################################################################
#multi outcome tests
################################################################
simu_ESDGENlm_alphaTest_multi = function(parallel_idx = 1,
                                         num_sig = 10,
                                         num_rho = 0.3,
                                         k = 10,
                                         alpha = 0.05,
                                         total_audio = 50,
                                         npatient=80,
                                         vec_lam,
                                         var_type="rob")
{
  ####gene data with same betas
  ll_dat <-
    gene.data.gee(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 66.95,
      medium.mean        = 66.95,
      medium_number = 3,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig,
      rho = num_rho,
      seed_idx = parallel_idx
    )
  temp = GEEfit(
    my.data = ll_dat$df.XY,
    Y.name = "Y",
    X.name = c("age", "age.sq", "goodhear", "littletrouble"),
    Test_reviewer = colnames(ll_dat$df.XY)[-c(1:6)],
    id = "id",
    corstr = "exchangeable"
  )
  vec_betahat = temp$coefficients[5:(total_audio + 4)]
  if(var_type=="rob") {
    mat_cov_beta = temp$robust.variance[5:(total_audio + 4), 5:(total_audio + 4)]
  } else {
    mat_cov_beta = temp$naive.variance[5:(total_audio + 4), 5:(total_audio + 4)]
  }
  rm(temp)
  # browser()
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam[1:k])
  
  vec_res=rep(NA,2)
  names(vec_res)=c("GENlm","Lasso")
  vec_res[1] = length(vec_idx)
}





simu_ESDGENlm_PowerTest_multi = function(parallel_idx = 1,
                                         num_sig = 10,
                                         num_rho = 0.3,
                                         k = 10,
                                         alpha = 0.05,
                                         total_audio = 50,
                                         npatient=80,
                                         vec_lam,
                                         var_type="rob")
{
  ####gene data with same betas
  ll_dat <-
    gene.data.gee(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig,
      rho = num_rho
    )
  temp = GEEfit(
    my.data = ll_dat$df.XY,
    Y.name = "Y",
    X.name = c("age", "age.sq", "goodhear", "littletrouble"),
    Test_reviewer = colnames(ll_dat$df.XY)[-c(1:6)],
    id = "id",
    corstr = "exchangeable"
  )
  vec_betahat = temp$coefficients[5:(total_audio + 4)]
  if(var_type=="rob") {
    mat_cov_beta = temp$robust.variance[5:(total_audio + 4), 5:(total_audio + 4)]
  } else {
    mat_cov_beta = temp$naive.variance[5:(total_audio + 4), 5:(total_audio + 4)]
  }
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = npatient, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > vec_lam)
  if (length(vec_idx) == 0)
  {
    return(0)
  } else{
    return(ll_Ri$Ic[[max(vec_idx)]])
  }
}


# simu_ESDGENlm_alphaTest_multi_v2 = function(parallel_idx = 1,
#                                          num_sig = 10,
#                                          num_rho = 0.3,
#                                          k = 20,
#                                          alpha = 0.05,
#                                          total_audio = 100,
#                                          var_type="rob")
# {
#   ####gene data with same betas
#   ll_dat <-
#     gene.data.gee(
#       num_audio          = total_audio,
#       patient = 40,
#       non_zero_audio     = 5,
#       abn.audio.mean = 66.95,
#       medium.mean        = 66.95,
#       medium_number = 3,
#       nor.audio.mean     = 66.95,
#       age.mean = 56.56,
#       age.std            = 4.36,
#       goodhear.p = 0.4358,
#       littletrouble.p    = 0.2519,
#       coef_age = -2.73,
#       coef_age.sq        = 0.03,
#       coef_goodhear = 0.03,
#       coef_littletrouble = 3.32,
#       RMSE = num_sig,
#       rho = num_rho,
#       seed_idx=parallel_idx
#     )
#   temp = GEEfit(
#     my.data = ll_dat$df.XY,
#     Y.name = "Y",
#     X.name = c("age", "age.sq", "goodhear", "littletrouble"),
#     Test_reviewer = colnames(ll_dat$df.XY)[-c(1:6)],
#     id = "id",
#     corstr = "exchangeable"
#   )
#   vec_betahat = temp$coefficients[5:(total_audio + 4)]
#   if(var_type=="rob") {
#     mat_cov_beta = temp$robust.variance[5:(total_audio + 4), 5:(total_audio + 4)]
#   } else {
#     mat_cov_beta = temp$naive.variance[5:(total_audio + 4), 5:(total_audio + 4)]
#   }
#   rm(temp)
#   ll_ord = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
#   vec_lambda = critival_value_GENlm_mvnorm_v2(vec_betahat, mat_cov_beta, ll_ord, k, alpha, F)
#   ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
#   
#   vec_idx = which(ll_Ri$Ri[1:k] > vec_lambda)
#   length(vec_idx)
# }
# 
# simu_ESDGENlm_PowerTest = function(parallel_idx = 1,
#                                    num_sig = 10,
#                                    k = 20,
#                                    alpha = 0.05)
# {
#   ####gene data with same betas
#   ll_dat <-
#     gene.data(
#       num_audio          = 100,
#       patient = 40,
#       non_zero_audio     = 5,
#       abn.audio.mean = 75.10,
#       medium.mean        = 70.10,
#       medium_number = 5,
#       nor.audio.mean     = 66.95,
#       age.mean = 56.56,
#       age.std            = 4.36,
#       goodhear.p = 0.4358,
#       littletrouble.p    = 0.2519,
#       coef_age = -2.73,
#       coef_age.sq        = 0.03,
#       coef_goodhear = 0.03,
#       coef_littletrouble = 3.32,
#       RMSE = num_sig
#     )
#   obj_lm = lm(Y ~ 0 + .,
#               data = ll_dat$df.XY)
#   vec_betahat = obj_lm$coefficients[5:104]
#   num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
#   mat_X = as.matrix(ll_dat$df.XY[, 1:104])
#   mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:104, 5:104] * num_sigma
#   
#   ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = 100, k = k)
#   vec_lam = critival_value_GENlm_mvnorm(vec_betahat,
#                                         mat_X,
#                                         ll_Ri$Ic[[k]],
#                                         k = k,
#                                         alpha = alpha)
#   vec_idx = which(ll_Ri$Ri[1:k] > vec_lam)
#   if (length(vec_idx) == 0)
#   {
#     return(0)
#   } else{
#     return(ll_Ri$Ic[[max(vec_idx)]])
#   }
# }
# 

# 
# simu_ESDGENlm_alphaTest_v2 = function(parallel_idx = 1,
#                                       num_sig = 10,
#                                       k = 10,
#                                       alpha = 0.05,
#                                       total_audio = 50,
#                                       npatient=80,
#                                       vec_lam)
# {
#   warning("Critical value specific version! use with caution!")
#   ####gene data with same betas
#   ll_dat <-
#     gene.data(
#       num_audio          = total_audio,
#       patient = npatient,
#       non_zero_audio     = 5,
#       abn.audio.mean = 66.95,
#       medium.mean        = 66.95,
#       medium_number = 5,
#       nor.audio.mean     = 66.95,
#       age.mean = 56.56,
#       age.std            = 4.36,
#       goodhear.p = 0.4358,
#       littletrouble.p    = 0.2519,
#       coef_age = -2.73,
#       coef_age.sq        = 0.03,
#       coef_goodhear = 0.03,
#       coef_littletrouble = 3.32,
#       RMSE = num_sig
#     )
#   obj_lm = lm(Y ~ 0 + .,
#               data = ll_dat$df.XY)
#   vec_betahat = obj_lm$coefficients[5:(total_audio + 4)]
#   num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
#   mat_X = as.matrix(ll_dat$df.XY[, 1:(total_audio + 4)])
#   mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(total_audio + 4), 5:(total_audio + 4)] * num_sigma
#   # browser()
#   ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
#   # vec_lam = critival_value_GENlm_mvnorm(vec_betahat, mat_X, ll_Ri$Ic[[k]],k = k, alpha = alpha)
#   vec_idx = which(ll_Ri$Ri[1:k] > vec_lam)
#   length(vec_idx)
# }



#################################################################
#unified lambda
#################################################################
fun_uni_lambda_Ri = function(
    parallel_idx = 1,
    mat_beta_simu,
    mat_sigma,
    k,
    M) 
{
  require(mvtnorm)
  vec_beta = mat_beta_simu[parallel_idx, ]
  mat_omega_cur = mat_sigma
  vec_Ri_simu = rep(NA, k)
  for (k0 in 0:(k - 1))
  {
    vec_Ri = Ri_ESDChi_Hl(vec_beta, mat_omega_cur, M = M - k0)
    int_out = sample(M - k0, 1)
    ##record the value
    vec_Ri_simu[k0+1] = max(vec_Ri)
    ##update the mat_omega matrix and beta
    mat_omega_cur = mat_omega_cur[-int_out, ]
    mat_omega_cur = mat_omega_cur[, -int_out]
    vec_beta = vec_beta[-int_out]
  }
  vec_Ri_simu
}

fun_uni_lambda_Ri_v2 = function(
    parallel_idx = 1,
    mat_beta_simu,
    mat_sigma,
    k,
    M
    )
{
  require(mvtnorm)
  vec_beta = mat_beta_simu[parallel_idx, ]
  mat_omega_cur = mat_sigma
  ll_Ri = Ri_ESDChi(vec_beta, mat_omega_cur, M = M, k = k)
  ll_Ri$Ri[1:k]
}

fun_exact_lambda_Ri = function(
    parallel_idx = 1,
    vec_ord,
    mat_beta_simu,
    mat_sigma,
    k,
    M
    )
{
  ll_Rmat = list()
  require(mvtnorm)
  vec_beta_simu = mat_beta_simu[parallel_idx,]
  for (k1 in (k-1):0) {
    if(k1>0)
    {
      vec_keep = setdiff(1:M,vec_ord[1:k1])
    } else {
      vec_keep = 1:M
    }
    
    submat_sigma = mat_sigma[vec_keep,vec_keep]
    ll_Ri = Ri_ESDChi(
      vec_beta_simu[vec_keep],
      submat_sigma,
      M = M-k1,
      k = k-k1)
    ll_Rmat[[k-k1]] = ll_Ri$Ri[1:(k-k1)]
  }
  ll_Rmat
}

f_exact_critical_bsearch = function(
    obj_Ri,alf,tol = 1e-3)
{
  k = length(obj_Ri[[1]])
  mat_lambda = NULL
  n_simus = length(obj_Ri)
  for (k1 in (k-1):0) {
    mat_Ri_simu = matrix(ncol = k-k1, nrow = n_simus)
    for(ii in 1:n_simus) {
      mat_Ri_simu[ii,] = obj_Ri[[ii]][[k-k1]]
    }
    ##reset binary search parameter
    lft = 0
    rgt = max(mat_Ri_simu)
    diff = 99
    while (diff > tol) {
      mid = (lft + rgt) / 2
      cuttst = cbind(rep(mid, n_simus), mat_lambda)
      idxmat = (mat_Ri_simu > cuttst)
      f.mid = mean(apply(idxmat, 1, any))
      if (f.mid <= alf) {
        rgt = mid
      }
      if (f.mid > alf) {
        lft = mid
      }
      diff = rgt - lft
    }
    mat_lambda = cuttst
  }
  mat_lambda[1,]
}

f_critical_bsearch = function(
    mat_Ri,alpha)
{
  Lb = 0
  Ub = max(mat_Ri)
  cur = (Lb + Ub) / 2
  n_bsearch = 100
  ct = 0
  while (ct < n_bsearch)
  {
    num_er = type1error(cur, mat_Ri)
    if (abs(num_er - alpha) < 1e-4)
    {
      break
    }
    if (num_er < alpha)
    {
      Ub = cur
    } else
    {
      Lb = cur
    }
    cur = (Lb + Ub) / 2
    ct = ct + 1
  }
  num_lambda = cur
  vec_lambda = apply(mat_Ri, 2, quantile, probs = 1-alpha)
  list(num_lambda,vec_lambda)
}

fun_uni_lambda_typeI = function(
    parallel_idx = 1,
    lambda,
    k = 20,
    M = 100,
    npatient,
    num_noise = num_sig)
{
  ll_dat <-
    gene.data(
      num_audio          = M,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 66.95,
      medium.mean        = 66.95,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_noise,
      seed_idx = parallel_idx
    )
  
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(4+M)]
  vec_order = order(obj_lm$coefficients[5:(4+M)], decreasing = T)
  
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(4+M)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(4+M), 5:(4+M)] * num_sigma
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = total_audio, k = k)
  any(ll_Ri$Ri[1:k] > lambda)
}

fun_uni_lambda_power = function(
    parallel_idx = 1,
    lambda,
    k = 20,
    M = 100,
    npatient,
    num_noise = num_sig)
{
  ##power
  ll_dat <-
    gene.data(
      num_audio          = M,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_noise,
      seed_idx = parallel_idx
    )
  
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(4 + M)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(4 + M)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(4 + M), 5:(4 + M)] * num_sigma
  
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = M, k = k)
  vec_idx = which(ll_Ri$Ri[1:k] > lambda)
  
  vec_out_idx = c()
  if (length(vec_idx) == 0)
  {
    vec_out_idx = c(0)
  } else{
    vec_out_idx = ll_Ri$Ic[[max(vec_idx)]]
  }
  
  vec_out_idx = ll_Ri$Ic[[max(vec_idx)]]
  
  num_trp = sum(vec_out_idx %in% 1:10) / 10
  num_tnr = length(setdiff(11:100, vec_out_idx)) / 90
  num_nrej = length(vec_out_idx)
  c(num_trp,num_tnr, num_nrej)
}

##multi simulation
fun_multi_lambda_typeI = function(
    parallel_idx = 1,
    lambda,
    k = 20,
    M = 100,
    num_rho,
    npatient,
    num_noise,
    var_type)
{
  ll_dat <-
    gene.data.gee(
      num_audio          = M,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 66.95,
      medium.mean        = 66.95,
      medium_number = 3,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_noise,
      rho = num_rho,
      seed_idx = parallel_idx
    )
  temp = GEEfit(
    my.data = ll_dat$df.XY,
    Y.name = "Y",
    X.name = c("age", "age.sq", "goodhear", "littletrouble"),
    Test_reviewer = colnames(ll_dat$df.XY)[-c(1:6)],
    id = "id",
    corstr = "exchangeable"
  )
  vec_betahat = temp$coefficients[5:(M + 4)]
  if(var_type=="rob") {
    mat_cov_beta = temp$robust.variance[5:(M + 4), 5:(M + 4)]
  } else {
    mat_cov_beta = temp$naive.variance[5:(M + 4), 5:(M + 4)]
  }
  
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = M, k = k)
  # browser()
  any(ll_Ri$Ri[1:k] > lambda)
}

fun_multi_lambda_power = function(
    parallel_idx = 1,
    lambda,
    k,
    M,
    num_rho,
    npatient,
    num_noise,
    var_type)
{
  ##power
  ll_dat <-
    gene.data.gee(
      num_audio          = M,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_noise,
      rho = num_rho,
      seed_idx = parallel_idx
    )
  temp = GEEfit(
    my.data = ll_dat$df.XY,
    Y.name = "Y",
    X.name = c("age", "age.sq", "goodhear", "littletrouble"),
    Test_reviewer = colnames(ll_dat$df.XY)[-c(1:6)],
    id = "id",
    corstr = "exchangeable"
  )
  vec_betahat = temp$coefficients[5:(M + 4)]
  if(var_type=="rob") {
    mat_cov_beta = temp$robust.variance[5:(M + 4), 5:(M + 4)]
  } else {
    mat_cov_beta = temp$naive.variance[5:(M + 4), 5:(M + 4)]
  }
  ll_Ri = Ri_ESDChi(vec_betahat, mat_cov_beta, M = npatient, k = k)
  
  vec_idx = which(ll_Ri$Ri[1:k] > lambda)
  
  vec_out_idx = c()
  if (length(vec_idx) == 0)
  {
    vec_out_idx = c(0)
  } else{
    vec_out_idx = ll_Ri$Ic[[max(vec_idx)]]
  }
  
  vec_out_idx = ll_Ri$Ic[[max(vec_idx)]]
  
  num_trp = sum(vec_out_idx %in% 1:10) / 10
  num_tnr = length(setdiff(11:100, vec_out_idx)) / 90
  num_nrej = length(vec_out_idx)
  c(num_trp,num_tnr, num_nrej)
}

pyDataGen = function(
  parallel_idx = 1,
  num_sig = 10,
  k = 10,
  total_audio = 50,
  npatient=80
)
{
  require(glmnet)
  warning("Critical value specific version! use with caution!")
  ####gene data with same betas
  ll_dat <-
    gene.data(
      num_audio          = total_audio,
      patient = npatient,
      non_zero_audio     = 5,
      abn.audio.mean = 75.10,
      medium.mean        = 70.10,
      medium_number = 5,
      nor.audio.mean     = 66.95,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = num_sig, seed_idx=parallel_idx
    )
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  vec_betahat = obj_lm$coefficients[5:(4 + total_audio)]
  num_sigma = (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
  mat_X = as.matrix(ll_dat$df.XY[, 1:(4 + total_audio)])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(4 + total_audio), 5:(4 + total_audio)] * num_sigma
  return(as.numeric(vec_betahat))
}