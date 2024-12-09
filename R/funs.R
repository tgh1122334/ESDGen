#' Simulating data generation, single outcome.
#'
#' Generate the Data which is similar with the CHEARS AAA data.
#' The predictors are \emph{age}, \emph{age}^2, \emph{hearing status}
#' (3 levels, good, a little trouble, and bad). \emph{age} is generated from a
#' normal distribution, mean is \emph{age.mean} and
#'  the standard deviation is \emph{age.std}.
#' \emph{Hearing status} is generated from multinormial distribution, with
#' probabilities \emph{goodhear.p}, and \emph{littletrouble.p}.
#'
#' @param num_audio number of audiologists.
#' @param patient number of patients evaluated by each audiologist.
#' @param non_zero_audio number of significant outlier audiologists.
#' @param abn.audio.mean effect of the significant outlier audiologists
#' @param medium.mean effect of the intermediate outlier audiologists
#' @param medium_number number of intermediate outlier audiologists.
#' @param age.mean mean of simulated covariate age for the patients.
#' @param age.std std of simulated covariate age for the patients.
#' @param goodhear.p probability that hearing status is good.
#' @param littletrouble.p probability that hearing status is having a little trouble.
#' @param coef_age coefficient of age.
#' @param coef_age.sq coefficient of age^2.
#' @param coef_goodhear coefficient of hearing status good.
#' @param coef_littletrouble coefficient of hearing status having a little trouble.
#' @param RMSE standard deviation of random effects.
#' @param seed_idx random seed for random effects.
#' @import glmnet
#' @export
#' @return A dataframe, with number of columns=\emph{4+num_audio},
#' rows=\emph{patient}*\emph{num_audio}.
#' @examples
#' add(1, 1) #add comments
#' add(10, 1)
gene.data <- function(num_audio,
                      patient,
                      non_zero_audio,
                      abn.audio.mean,
                      medium.mean,
                      medium_number ,
                      nor.audio.mean,
                      age.mean,
                      age.std,
                      goodhear.p,
                      littletrouble.p,
                      coef_age,
                      coef_age.sq,
                      coef_goodhear,
                      coef_littletrouble,
                      RMSE,
                      seed_idx = 1) {
  # set.seed(1)
  age            <-
    round(rnorm(num_audio * patient, age.mean, age.std), 0)
  age.sq         <- age ^ 2
  goodhear       <- rbinom(num_audio * patient, 1, goodhear.p)
  littletrouble  <- rbinom(num_audio * patient, 1, littletrouble.p)
  audiologist <-
    matrix(0, ncol = num_audio, nrow = patient * num_audio)
  colnames(audiologist) <- paste0("audiologist_", 1:num_audio)
  for (i in 1:num_audio) {
    index                  <- (patient * (i - 1) + 1):(patient * i)
    audiologist[, i][index] <- 1
  }
  XX             <-
    cbind(age, age.sq, goodhear, littletrouble, audiologist)

  #### give coefficient values to age, self_report and non-zero audiologists

  coef_audio     <- rep(NA, num_audio)
  index.sig      <- 1:(non_zero_audio + medium_number)
  sig.audio      <- paste0("audiologist_", index.sig)
  nosig.audio    <-
    paste0("audiologist_", c(1:num_audio)[-index.sig])

  abn.audio.coef <- rep(abn.audio.mean, non_zero_audio)

  coef_audio[index.sig]   <-
    c(abn.audio.coef, rep(medium.mean, medium_number))
  coef_audio[-index.sig]  <-
    rep(nor.audio.mean, (num_audio - non_zero_audio - medium_number))

  coef_all   <- c(coef_age,
                  coef_age.sq,
                  coef_goodhear,
                  coef_littletrouble,
                  coef_audio)
  # set.seed(seed_idx)
  rm(.Random.seed, envir = globalenv())
  #### generating outcome variable Y, with random noise N(0,1)
  YY         <-
    XX %*% coef_all + rnorm(num_audio * patient, 0, RMSE)
  df.XY      <- data.frame(XX, Y = YY)

  return(list(df.XY = df.XY,
              coef  = coef_all))
}

mygene.data <- function(num_audio,
                        patient,
                        vec_audio_beta,
                        age.mean,
                        age.std,
                        goodhear.p,
                        littletrouble.p,
                        coef_age,
                        coef_age.sq,
                        coef_goodhear,
                        coef_littletrouble,
                        RMSE,
                        seed_idx = 1) {
  # set.seed(1)
  age            <-
    round(rnorm(num_audio * patient, age.mean, age.std), 0)
  age.sq         <- age ^ 2
  goodhear       <- rbinom(num_audio * patient, 1, goodhear.p)
  littletrouble  <- rbinom(num_audio * patient, 1, littletrouble.p)
  audiologist <-
    matrix(0, ncol = num_audio, nrow = patient * num_audio)
  colnames(audiologist) <- paste0("audiologist_", 1:num_audio)
  for (i in 1:num_audio) {
    index                  <- (patient * (i - 1) + 1):(patient * i)
    audiologist[, i][index] <- 1
  }
  XX             <-
    cbind(age, age.sq, goodhear, littletrouble, audiologist)
  coef_all   <- c(coef_age,
                  coef_age.sq,
                  coef_goodhear,
                  coef_littletrouble,
                  vec_audio_beta)
  # set.seed(seed_idx)
  rm(.Random.seed, envir = globalenv())
  #### generating outcome variable Y, with random noise N(0,1)
  YY         <-
    XX %*% coef_all + rnorm(num_audio * patient, 0, RMSE)
  df.XY      <- data.frame(XX, Y = YY)

  return(list(df.XY = df.XY,
              coef  = coef_all))
}

gene.data.gee <- function(num_audio,
                          patient,
                          non_zero_audio,
                          abn.audio.mean,
                          medium.mean,
                          medium_number,
                          nor.audio.mean,
                          age.mean,
                          age.std,
                          goodhear.p,
                          littletrouble.p,
                          coef_age,
                          coef_age.sq,
                          coef_goodhear,
                          coef_littletrouble,
                          RMSE,
                          rho,
                          seed_idx = 1) {
  require(MASS)
  require(gee)
  # set.seed(1)
  age            <-
    round(rnorm(num_audio * patient, age.mean, age.std), 0)
  age.sq         <- age ^ 2
  goodhear       <- rbinom(num_audio * patient, 1, goodhear.p)
  littletrouble  <- rbinom(num_audio * patient, 1, littletrouble.p)
  audiologist <-
    matrix(0, ncol = num_audio, nrow = patient * num_audio)
  colnames(audiologist) <- paste0("audiologist_", 1:num_audio)
  for (i in 1:num_audio) {
    index                  <- (patient * (i - 1) + 1):(patient * i)
    audiologist[, i][index] <- 1
  }
  XX             <-
    cbind(age, age.sq, goodhear, littletrouble, audiologist)
  # set.seed(seed_idx)
  rm(.Random.seed, envir = globalenv())
  #### give coefficient values to age, self_report and non-zero audiologists

  coef_audio     <- rep(NA, num_audio)
  index.sig      <- 1:(non_zero_audio + medium_number)
  sig.audio      <- paste0("audiologist_", index.sig)
  nosig.audio    <-
    paste0("audiologist_", c(1:num_audio)[-index.sig])

  abn.audio.coef <- rep(abn.audio.mean, non_zero_audio)

  coef_audio[index.sig]   <-
    c(abn.audio.coef, rep(medium.mean, medium_number))
  coef_audio[-index.sig]  <-
    rep(nor.audio.mean, (num_audio - non_zero_audio - medium_number))

  coef_all   <- c(coef_age,
                  coef_age.sq,
                  coef_goodhear,
                  coef_littletrouble,
                  coef_audio)

  #### generating outcome variable Y, with random noise N(0,1)
  mu         <- XX %*% coef_all
  XX         <- cbind(XX, id = 1:(num_audio * patient))
  mu         <- rep(mu, each = 2)
  Y          <- c()
  for (pair in 1:(length(mu) / 2)) {
    index      <- (2 * (pair - 1) + 1):(2 * pair)
    YY         <-
      mvrnorm(
        n = ,
        mu = mu[index],
        Sigma = matrix(
          c(RMSE ^ 2, rho * RMSE ^ 2,
            rho * RMSE ^
              2, RMSE ^ 2),
          byrow = TRUE,
          nrow = 2,
          ncol = 2
        )
      )
    id         <- pair
    YY         <- c(YY, pair)
    Y          <- rbind(Y, YY)
  }
  Y.long       <- data.frame(Y  = c(Y[, 1], Y[, 2]),
                             id = c(Y[, 3], Y[, 3]))
  Y.long       <- Y.long[order(Y.long$id),]

  df.XY        <- merge(Y.long, XX, by.x = "id")

  return(list(df.XY = df.XY,
              coef  = coef_all))
}

###return R_i i=1:k, k is the number of pre_set outliers
ESD_Ri = function(X, k)
{
  n = length(X)
  vec_R = rep(NA, k)
  vec_R[1] = max(abs(X - mean(X))) / sd(X)
  vec_XI = rep(NA, k)
  ll_I = list()
  ll_I[[1]] = 1:n

  for (i in 2:(k + 1))
  {
    ###renew I_i
    num_mu_It = mean(X[ll_I[[i - 1]]])
    int_t = which.max((abs(X[ll_I[[i - 1]]] - num_mu_It)))
    ll_I[[i]] = setdiff(ll_I[[i - 1]], ll_I[[i - 1]][int_t])
    vec_R[i] = max(abs(X[ll_I[[i]]] - mean(X[ll_I[[i]]]))) / sd(X[ll_I[[i]]])
    vec_XI[i] = int_t
  }
  list(R = vec_R, Xi = vec_XI, I = ll_I)
}
dat_gen = function(ab_beta = 75.1,
                   mab_beta = 70.10,
                   nor_beta = 66.95,
                   sig_sq = 10)
{
  ll_dat <-
    gene.data(
      num_audio          = 100,
      patient = 40,
      non_zero_audio     = 5,
      abn.audio.mean = ab_beta,
      # abn.audio.mean = 68,
      medium.mean        = mab_beta,
      # medium.mean        = 67,
      medium_number = 5,
      nor.audio.mean     = nor_beta,
      age.mean = 56.56,
      age.std            = 4.36,
      goodhear.p = 0.4358,
      littletrouble.p    = 0.2519,
      coef_age = -2.73,
      coef_age.sq        = 0.03,
      coef_goodhear = 0.03,
      coef_littletrouble = 3.32,
      RMSE = sig_sq
    )


  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)
  summary(obj_lm)
  num_sigma = (sum(obj_lm$residuals ^ 2) / 3896)
  mat_X = as.matrix(ll_dat$df.XY[, 1:104])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:104, 5:104] * num_sigma

  vec_betahat = obj_lm$coefficients[5:104]
  list(
    beta_true = ll_dat$coef,
    beta_pre = vec_betahat,
    conv = mat_cov_beta,
    sigma = num_sigma,
    X = mat_X
  )
}
###determing the percentage point lambda(beta) and alpha
### two side version
lambda_l_2side = function(alpha, k, n)
{
  vec_l = 0:(k - 1)
  vec_lambda_lp1 = rep(NA, k)
  for (l in 0:(k - 1))
  {
    num_p = 1 - ((alpha / 2) / (n - l))
    vec_lambda_lp1[l + 1] =
      qt(p = num_p, df = n - l - 2) * (n - l - 1) /
      sqrt((n - l - 2 + qt(p = num_p, df = n - l - 2) ^ 2) * (n - l))
  }
  vec_lambda_lp1
}
# print(lambda_l_2side(alpha = 0.05, k = 10, n = 50))
fun_simu_discontinued = function(alpha = 0.05, k = 20)
{
  ll_dat <-
    gene.data(
      num_audio          = 100,
      patient = 40,
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
      RMSE = 5
    )
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat$df.XY)

  mat_X = as.matrix(ll_dat$df.XY[, 1:104])
  mat_cov_beta = solve(t(mat_X) %*% mat_X)
  #R'R = x
  mat_cov_beta_sqr = chol(mat_cov_beta)
  ###verify the diag
  vec_beta_trans = t(solve(mat_cov_beta_sqr)) %*% obj_lm$coefficients
  vec_beta_trans = vec_beta_trans[5:104]
  # browser()
  obj_Ri = ESD_Ri(vec_beta_trans, k)
  vec_lambda = lambda_l_2side(alpha = alpha, k = k, n = 100)
  vec_idx = which(obj_Ri$R[2:(k + 1)] > vec_lambda)
  if (length(vec_idx) == 0)
  {
    list(nOutliers = 0, IdxOutliers = NA)
  }
  ##return the last continuous true
  total_out = 0
  # browser()
  for (i in 1:length(vec_idx))
  {
    if (vec_idx[i] == i)
    {
      total_out = i
    } else{
      break
    }
  }
  list(nOutliers = total_out,
       IdxOutliers = setdiff(1:100, obj_Ri$I[[total_out + 1]]))
}

parallel_simu = function(parallel_idx = 1,
                         alpha = 0.05,
                         k = 20)
{
  fun_simu_discontinued(alpha = alpha, k = k)
}


mean_L = function(M){
  matL = matrix(-1 / M, ncol = M, nrow = M)
  diag(matL) = 1.0 + diag(matL)
  matL
}

mean_L_j = function(M,j){
  matL = matrix(-1 / M, ncol = M, nrow = M)
  diag(matL) = 1.0 + diag(matL)
  matL[j,]
}

truncated_L = function(M, Mdelta)
{
  matL = matrix(0.0, ncol = M, nrow = M)
  diag(matL) = 1.0
  set_A = c(1:Mdelta, (M - Mdelta + 1):M)
  for (l in setdiff(1:M, set_A))
  {
    matL[, l] = matL[, l] - 1 / (M - 2 * Mdelta)
  }
  matL
}

trun_mean_v1 = function(xx, M, Mdelta)
{
  xx = sort(xx)
  mean(xx[-c(1:Mdelta, (M - Mdelta + 1):M)])
}

trun_mean_v2 = function(xx, j, M, Mdelta)
{
  vec_w = truncated_L(M, Mdelta)[j,]
  vec_w[j] = vec_w[j] - 1
  # print(vec_w)
  vec_rank = rank(xx)
  sum(xx * vec_w[vec_rank])
}

truncated_L_j = function(vec_rank, j, M, Mdelta)
{
  vec_w = truncated_L(M, Mdelta)[1,]
  vec_w[1] = vec_w[1] - 1
  vec_w = vec_w[vec_rank]
  vec_w[j] = vec_w[j] + 1
  vec_w
}

CHi_sq_stat = function(vecL, vec_beta, mat_cov_beta)
{
  vec_L_beta = t(vecL) %*% vec_beta
  # t(vec_L_beta) %*% solve(t(vecL) %*% mat_cov_beta %*% vecL) %*% vec_L_beta
  vec_L_beta*vec_L_beta/(t(vecL) %*% mat_cov_beta %*% vecL)
  # browser()
}

Ri_ESDChi_Hl = function(vec_betahat,
                        mat_cov_beta,
                        M = length(vec_betahat))
{
  k = 1
  Mdelta = floor(M * 0.2)
  vec_R = rep(NA, k)
  ###calculate the mean and cov of first Rj:R0
  vec_temp = rep(1:length(vec_betahat))
  for (j in 1:length(vec_temp))
  {
    vec_L = truncated_L_j(rank(vec_betahat), j, M = M, Mdelta = Mdelta)
    vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                              vec_beta = vec_betahat,
                              mat_cov_beta = mat_cov_beta)
  }
  return(vec_temp)
}

Ri_ESDChi_1 = function(
    vec_betahat,
    mat_cov_beta,
    M = length(vec_betahat),
    k = 20,
    trim_proc = 0.2
) {
  M = length(vec_betahat)
  Mdelta = floor(M * 0.2)
  ###calculate the mean and cov of first Rj:R0
  vec_temp = rep(1:length(vec_betahat))
  for (j in 1:length(vec_temp))
  {
    vec_L = truncated_L_j(rank(vec_betahat), j, M = M, Mdelta = Mdelta)
    vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                              vec_beta = vec_betahat,
                              mat_cov_beta = mat_cov_beta)
  }

  num_R = max(vec_temp)
  num_idx = which.max(vec_temp)
  c(num_R, num_idx)
}

Ri_ESDChi = function(vec_betahat,
                     mat_cov_beta,
                     M = length(vec_betahat),
                     k = 20)
{
  M = length(vec_betahat)
  Mdelta = floor(M * 0.2)
  vec_R = rep(NA, k)
  ###calculate the mean and cov of first Rj:R0
  vec_temp = rep(1:length(vec_betahat))
  for (j in 1:length(vec_temp))
  {
    vec_L = truncated_L_j(rank(vec_betahat), j, M = M, Mdelta = Mdelta)
    vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                              vec_beta = vec_betahat,
                              mat_cov_beta = mat_cov_beta)
    # browser()
  }

  ###deleteing process
  vec_R[1] = max(vec_temp)
  t = 1
  ll_I = list()
  ll_Ic = list()
  ll_Ic[[1]] = numeric(0)
  ll_I[[1]] = 1:M ### start at including all betas
  for (t in 2:(k + 1))
  {
    ###update the cov mat
    mat_cov_beta_temp = mat_cov_beta[ll_I[[t - 1]], ll_I[[t - 1]]]
    ###update the beta vec, rank vec
    vec_beta_temp = vec_betahat[ll_I[[t - 1]]]
    vec_rank = rank(vec_beta_temp)
    ###create a temp vec for potential Ri
    vec_temp = rep(NA, M - t + 2)
    ### fill in the vec
    for (j in 1:length(vec_temp))
    {
      vec_L = truncated_L_j(vec_rank, j, M - t + 2, Mdelta)
      vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                                vec_beta = vec_beta_temp,
                                mat_cov_beta = mat_cov_beta_temp)
    }
    ###find the largest Ri
    int_t = which.max(vec_temp)
    ###thinning the Ii
    ll_I[[t]] = setdiff(ll_I[[t - 1]], ll_I[[t - 1]][int_t])
    # browser()
    ll_Ic[[t]] = append(ll_Ic[[t - 1]], ll_I[[t - 1]][int_t])
    ###store the Ri on thinned Ii
    vec_temp = rep(NA, M - t + 1)
    mat_cov_beta_temp = mat_cov_beta[ll_I[[t]], ll_I[[t]]]
    ###update the beta vec, rank vec
    vec_beta_temp = vec_betahat[ll_I[[t]]]
    vec_rank = rank(vec_beta_temp)
    for (j in 1:length(vec_temp))
    {
      vec_L = truncated_L_j(vec_rank, j, M - t + 1, Mdelta)
      vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                                vec_beta = vec_beta_temp,
                                mat_cov_beta = mat_cov_beta_temp)
    }
    vec_R[t] = max(vec_temp)
  }
  # browser()
  list(Ri = vec_R, I = ll_I[-1], Ic = ll_Ic[2:(k + 1)])
}

Ri_ESDChi_mean = function(vec_betahat,
                     mat_cov_beta,
                     M = length(vec_betahat),
                     k = 20)
{
  M = length(vec_betahat)
  Mdelta = floor(M * 0.2)
  vec_R = rep(NA, k)
  ###calculate the mean and cov of first Rj:R0
  vec_temp = rep(1:length(vec_betahat))
  for (j in 1:length(vec_temp))
  {
    vec_L = mean_L_j(M = M, j = j)
    vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                              vec_beta = vec_betahat,
                              mat_cov_beta = mat_cov_beta)
    # browser()
  }

  ###deleteing process
  vec_R[1] = max(vec_temp)
  t = 1
  ll_I = list()
  ll_Ic = list()
  ll_Ic[[1]] = numeric(0)
  ll_I[[1]] = 1:M ### start at including all betas
  for (t in 2:(k + 1))
  {
    ###update the cov mat
    mat_cov_beta_temp = mat_cov_beta[ll_I[[t - 1]], ll_I[[t - 1]]]
    ###update the beta vec, rank vec
    vec_beta_temp = vec_betahat[ll_I[[t - 1]]]
    vec_rank = rank(vec_beta_temp)
    ###create a temp vec for potential Ri
    vec_temp = rep(NA, M - t + 2)
    ### fill in the vec
    for (j in 1:length(vec_temp))
    {
      vec_L = mean_L_j(M = length(vec_temp), j = j)
      vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                                vec_beta = vec_beta_temp,
                                mat_cov_beta = mat_cov_beta_temp)
    }
    ###find the largest Ri
    int_t = which.max(vec_temp)
    ###thinning the Ii
    ll_I[[t]] = setdiff(ll_I[[t - 1]], ll_I[[t - 1]][int_t])
    # browser()
    ll_Ic[[t]] = append(ll_Ic[[t - 1]], ll_I[[t - 1]][int_t])
    ###store the Ri on thinned Ii
    vec_temp = rep(NA, M - t + 1)
    mat_cov_beta_temp = mat_cov_beta[ll_I[[t]], ll_I[[t]]]
    ###update the beta vec, rank vec
    vec_beta_temp = vec_betahat[ll_I[[t]]]
    vec_rank = rank(vec_beta_temp)
    for (j in 1:length(vec_temp))
    {
      vec_L = mean_L_j(M = length(vec_temp), j = j)
      vec_temp[j] = CHi_sq_stat(vecL = vec_L,
                                vec_beta = vec_beta_temp,
                                mat_cov_beta = mat_cov_beta_temp)
    }
    vec_R[t] = max(vec_temp)
  }
  # browser()
  list(Ri = vec_R, I = ll_I[-1], Ic = ll_Ic[2:(k + 1)])
}

simu_ChiESD_sample = function(mat_beta_conv,
                              nsample = 1e3,
                              k = 20,
                              mu)
{
  require(mvtnorm)
  beta_sample = rmvnorm(mean = mu, n = nsample, sigma = mat_beta_conv)
  mat_Ri_sample = matrix(ncol = k + 1, nrow = nsample)
  for (i in 1:nsample)
  {
    temp_Ri = Ri_ESDChi(beta_sample[i, ], mat_beta_conv, k = k)
    # browser()
    mat_Ri_sample[i,] = temp_Ri$Ri
    cat("\r", i, "/", nsample)
  }
  mat_Ri_sample
}

simu_ChiESD_sample_parallel = function(mat_beta_conv,
                                       nsample = 1e3,
                                       k = 20,
                                       mu,
                                       ncpus = 4)
{
  # browser()

  require(mvtnorm)
  require(snowfall)
  try(sfStop())
  sfInit(cpus = ncpus, parallel = T)
  sfSource("funs.R")
  sfExportAll()
  beta_sample = rmvnorm(mean = mu, n = nsample, sigma = mat_beta_conv)
  fun_parallel = function(i, beta_sample, mat_cov_beta, M, k)
  {
    Ri_ESDChi(beta_sample[i, ],
              mat_beta_conv,
              M = dim(beta_sample)[2],
              k = k)$Ri
  }
  temp = sfClusterApplyLB(
    1:nsample,
    fun_parallel,
    beta_sample = beta_sample,
    mat_cov_beta = mat_beta_conv,
    M = dim(beta_sample)[2],
    k = k
  )
  mat_Ri_sample = matrix(unlist(temp), ncol = k + 1, byrow = TRUE)
  sfStop()
  mat_Ri_sample
}

critival_value_ESD = function(vec_b,
                              mat_Ri_sample,
                              nsample = dim(mat_Ri_sample)[1],
                              k = dim(mat_Ri_sample)[2])
{
  vec_alpha = rep(NA, length(vec_b))
  mat_lambda_i = matrix(ncol = length(vec_b), nrow = k)
  colnames(mat_lambda_i) = vec_b
  for (j in 1:length(vec_b))
  {
    num_beta = vec_b[j]
    vec_lambda_i = rep(NA, k)
    for (i in 1:k)
    {
      vec_lambda_i[i] = quantile(mat_Ri_sample[, i], 1 - num_beta)
    }
    mat_lambda_i[, j] = vec_lambda_i
    ###calculate alpha
    vec_alp_idx = integer(0)
    for (i in 1:k)
    {
      vec_alp_idx = union(vec_alp_idx, which(mat_Ri_sample[, i] > vec_lambda_i[i]))
    }
    vec_alpha[j] = length(vec_alp_idx) / nsample
  }
  ###return lambda_i and alpha
  list(lambda = mat_lambda_i,
       alpha = vec_alpha,
       b = vec_b)
}

critival_value_GEN = function(vec_b,
                              mat_Ri_sample,
                              nsample = dim(mat_Ri_sample)[1],
                              k = dim(mat_Ri_sample)[2])
{
  vec_alpha = rep(NA, length(vec_b))
  mat_lambda_i = matrix(ncol = length(vec_b), nrow = k + 1)
  colnames(mat_lambda_i) = vec_b
  for (j in 1:length(vec_b))
  {
    num_beta = vec_b[j]
    vec_lambda_i = rep(NA, k + 1)
    for (i in 1:k)
    {
      vec_lambda_i[i] = quantile(mat_Ri_sample[, i], 1 - num_beta)
    }
    mat_lambda_i[, j] = vec_lambda_i
    ###calculate alpha
    vec_alp_idx = integer(0)
    for (i in 1:k)
    {
      vec_alp_idx = union(vec_alp_idx, which(mat_Ri_sample[, i] > vec_lambda_i[i]))
    }
    vec_alpha[j] = length(vec_alp_idx) / nsample
  }
  ###return lambda_i and alpha
  list(lambda = mat_lambda_i,
       alpha = vec_alpha,
       b = vec_b)
}

Out_hypothesis = function(vec_Ri, ll_I, mat_lambda_i, vec_alpha, k, M)
{
  #dim(mat_lambda_i) == c(k, length(alpha))
  mat_outs = matrix(0.0, ncol = length(vec_alpha), nrow = k)
  colnames(mat_outs) = vec_alpha
  for (i in 1:length(vec_alpha))
  {
    idx_l_max = max(which(vec_Ri[1:k] > mat_lambda_i[, i]))
    # browser()
    if (is.infinite(idx_l_max))
    {
      next
    }
    vec_I = ll_I[[idx_l_max]]
    vec_outs = setdiff(1:M, vec_I)
    mat_outs[1:length(vec_outs), i] = vec_outs
  }
  mat_outs
}

parallel_metric = function(parallel_idx = 1,
                           sig_sq = 10,
                           k = 15,
                           obj_crit)
{
  ll_para = dat_gen(sig_sq = sig_sq)
  ll_Ri = Ri_ESDChi(ll_para$beta_pre, ll_para$conv, k = k)
  mat_outs = Out_hypothesis(ll_Ri$Ri,
                            ll_Ri$I,
                            obj_crit$lambda,
                            obj_crit$alpha,
                            k = k,
                            M = 100)
  vec_tpr = vec_fdr = rep(NA, length(obj_crit$alpha))
  for (i in 1:length(obj_crit$alpha))
  {
    vec_tpr[i] = sum(mat_outs[, i] %in% 1:10) / 10
    vec_fdr[i] = sum(mat_outs[, i] %in% 11:100) / sum(mat_outs[, i] != 0)
    if (sum(mat_outs[, i] != 0) == 0)
    {
      vec_fdr[i] = 0
    }
  }
  data.frame(TPR = vec_tpr,
             FDR = vec_fdr,
             alpha = obj_crit$alpha)
}

parallel_verify_T1 = function(parallel_idx = 1,
                              sig_sq = 10,
                              k = 15,
                              obj_crit)
{
  ll_para = dat_gen(
    ab_beta = 66.95,
    mab_beta = 66.95,
    nor_beta = 66.95,
    sig_sq = sig_sq
  )
  ll_Ri = Ri_ESDChi(ll_para$beta_pre, ll_para$conv, k = k)
  mat_outs = Out_hypothesis(ll_Ri$Ri,
                            ll_Ri$I,
                            obj_crit$lambda,
                            obj_crit$alpha,
                            k = k,
                            M = 100)

  vec_reg_num = rep(NA, length(obj_crit$alpha))
  names(vec_reg_num) = obj_crit$alpha
  for (i in 1:length(vec_reg_num))
  {
    vec_reg_num[i] = sum(mat_outs[, i] != 0)
  }
  vec_reg_num
}

lambda_l = function(l, alpha, n)
{
  num_p = 1 - ((alpha / 2) / (n - l))
  qt(p = num_p, df = n - l - 2) * (n - l - 1) /
    sqrt((n - l - 2 + qt(p = num_p, df = n - l - 2) ^ 2) * (n - l))
}

check_alpha = function(k, n, alpha, mat_Ri)
{
  vec_idx = c()
  for (i in 0:(k - 1))
  {
    vec_idx = union(vec_idx, which(mat_Ri[, i + 1] > lambda_l(i, alpha, n)))
  }
  length(vec_idx) / dim(mat_Ri)[1]
}

GEN_Ri_simu = function(niter = 1e4,
                       n = 50,
                       k = 10)
{
  require(parallel)
  cl = makePSOCKcluster(8)
  clusterEvalQ(cl, {
    source("funs.R")
  })
  fun_para = function(para_idx = 1, X, k)
  {
    ll = ESD_Ri(X[para_idx,], k)
    ll$R
  }
  mat_nsample = matrix(rnorm(n * niter), ncol = n)
  temp = clusterApply(cl, 1:niter, fun_para, X = mat_nsample, k = k)
  mat_Ri = matrix(unlist(temp), ncol = k + 1, byrow = T)
  stopCluster(cl)
  mat_Ri
}

###############begin GENlm approach
###critical value obtained by direct simulation
critival_value_GENlm_simu = function(vec_beta_hat, mat_X, vec_order, k, alpha = 0.05)
{
  M = length(vec_beta_hat)
  Mdelta = floor(M * 0.2)
  ###determine the order by quad form
  # vec_temp = rep(NA,length(vec_beta_hat))
  # for(j in 1:length(vec_temp))
  # {
  #   vec_L = truncated_L_j(rank(vec_beta_hat),j, M = M, Mdelta = Mdelta)
  #   vec_temp[j] = CHi_sq_stat(
  #     vecL = vec_L,
  #     vec_beta = vec_beta_hat,
  #     mat_cov_beta = solve(t(mat_X) %*% mat_X)[5:(4+M),5:(4+M)]
  #   )
  # }
  # vec_order = order(vec_temp, decreasing = T)
  vec_lambda = rep(NA, k)
  # browser()
  for (l in 0:(k - 1))
  {
    if (l > 0)
    {
      ###delete some columns
      vec_delete = c()
      for (j in 1:l)
      {
        ##append each column
        vec_delete = c(vec_delete, which(mat_X[, 4 + vec_order[1:l]] == 1))
      }
      ##remove the rows
      mat_Xl = mat_X[-vec_delete,]
      ##remove the columns
      mat_Xl = mat_Xl[,-(vec_order[1:l] + 4)]
    } else{
      mat_Xl = mat_X
    }
    # browser()
    mat_conv_beta = solve(t(mat_Xl) %*% mat_Xl)[5:(M + 4 - l), 5:(M + 4 - l)]
    # print(dim(mat_conv_beta))
    mat_Ri = simu_ChiESD_sample_parallel(
      mat_conv_beta,
      1000,
      k = 1,
      mu = rep(0, M - l),
      ncpus = 8
    )
    vec_lambda[l + 1] = quantile(mat_Ri[, 1], probs = 1 - alpha)
    cat(l)
  }
  vec_lambda
}

###critical value obtained by chisq extreme value approximation
critival_value_GENlm_chisq = function(M, k, alpha = 0.05)
{
  vec_lambda = rep(NA, k)
  int_nsamples = 1e3
  for (l in 0:(k - 1))
  {
    temp = matrix(rchisq(n = (M - l) * int_nsamples, df = 1),
                  ncol = M - l,
                  nrow = int_nsamples)
    vec_lambda[l + 1] = quantile(apply(temp, 1, max), 1 - alpha)
  }
  vec_lambda
}

###critical value by integration mvtnorm
critival_value_GENlm_mvnorm = function(vec_betahat,
                                       mat_X,
                                       vec_order,
                                       k,
                                       alpha = 0.05,
                                       return_cov = F,
                                       num_covar = 4)
{
  require(mvtnorm)
  M = length(vec_betahat)
  Mdelta = floor(M * 0.2)

  vec_lambda = rep(NA, k)
  for (l in 0:(k - 1))
  {
    if (l > 0)
    {
      ###delete some columns
      vec_delete = c()
      for (j in 1:l)
      {
        ##append each column
        vec_delete = c(vec_delete, which(mat_X[, num_covar + vec_order[1:l]] == 1))
      }
      ##remove the rows
      mat_Xl = mat_X[-vec_delete,]
      ##remove the columns
      mat_Xl = mat_Xl[,-(vec_order[1:l] + num_covar)]
      vec_beta_sub = vec_betahat[-vec_order[1:l]]
    } else{
      mat_Xl = mat_X
      vec_beta_sub = vec_betahat
    }
    mat_conv_beta = solve(
      t(mat_Xl) %*% mat_Xl)[
        (num_covar+1):(M + num_covar - l),
        (num_covar+1):(M + num_covar - l)]
    # browser()
    mat_A = matrix(ncol = M - l, nrow = M - l)
    vec_rank = rank(vec_beta_sub)
    for (i in 1:length(vec_beta_sub))
    {
      vec_L = truncated_L_j(vec_rank, i, M - l, Mdelta)
      mat_A[i,] = vec_L / sqrt(c((t(vec_L) %*% mat_conv_beta %*% vec_L)))
    }
    mat_cov_betatilta = mat_A %*% mat_conv_beta %*% t(mat_A)
    vec_lambda[l + 1] = qmvnorm(1 - alpha, tail = "both.tails", sigma = mat_cov_betatilta)$quantile ^
      2
    cat("\r", l, "out of ", k - 1, "finished")
    if (l == 0)
    {
      mat_cov_betatilta0 = mat_cov_betatilta
    }
  }
  if (return_cov)
  {
    list(lambda = vec_lambda, cov = mat_cov_betatilta0)
  } else
  {
    vec_lambda
  }
}

####use for multiple measurements
critival_value_GENlm_mvnorm_v2 = function(vec_betahat,
                                          beta_cov,
                                          obj_ord,
                                          k,
                                          alpha = 0.05,
                                          return_cov = F)
{
  require(mvtnorm)
  M = length(vec_betahat)
  Mdelta = floor(M * 0.2)
  vec_lambda = rep(NA, k)
  for (l in 0:(k - 1))
  {
    if (l == 0)
    {
      mat_conv_beta = beta_cov
      vec_beta_sub = vec_betahat
    } else{
      mat_conv_beta = beta_cov[obj_ord$I[[l]], obj_ord$I[[l]]]
      vec_beta_sub = vec_betahat[obj_ord$I[[l]]]
    }
    # browser()
    mat_A = matrix(ncol = M - l, nrow = M - l)
    vec_rank = rank(vec_beta_sub)
    # vec_rank = 1:(M - l)
    for (i in 1:length(vec_beta_sub))
    {
      vec_L = truncated_L_j(vec_rank, i, M - l, Mdelta)
      mat_A[i,] = vec_L * sqrt(c(solve(t(vec_L) %*% mat_conv_beta %*% vec_L)))
    }
    # browser()
    mat_cov_betatilta = mat_A %*% mat_conv_beta %*% t(mat_A)
    vec_lambda[l + 1] = qmvnorm(1 - alpha, tail = "both.tails", sigma = mat_cov_betatilta)$quantile ^
      2
    cat("\r", l, "out of ", k - 1, "finished")
    if (l == 0)
    {
      mat_cov_betatilta0 = mat_cov_betatilta
    }
  }
  if (return_cov)
  {
    list(lambda = vec_lambda, cov = mat_cov_betatilta0)
  } else
  {
    vec_lambda
  }
}

##use simulation for critical values
critival_value_GENlm_simu_intercept = function(M,
                                               k,
                                               omega,
                                               alf,
                                               trimfac,
                                               numreps,
                                               tol,
                                               BIG,
                                               cutmax,
                                               bs_max = 100,
                                               fac) {
  # M = number of units
  # k = max number of outlier units
  # omega = covariance matrix
  # alf = desired Type I error level
  # trimfac = percent of trimming on each side
  # numreps = number of simulation replications
  # tolerance for critical value search
  # BIG = outlier value
  # cutmax = upper bound for critical values
  # fac = factor for critical value search

  #PRELIMINARIES
  inputs = cbind(M, k, alf, trimfac, numreps, tol, BIG, cutmax, fac)
  cat('\n')
  print(noquote('inputs'))
  print(inputs)
  cat('\n')
  Rmat = matrix(0, numreps, k)
  cutmat = NULL

  #GENERATE SIMULATED BETAS
  betmat = rmvnorm(numreps, mean = rep(0, M), sigma = omega)

  #loop over k1 values
  # browser()
  for (k1 in (k - 1):0) {
    diff = 99
    ##parallel version
    require(parallel)
    cl = makeCluster(8)
    clusterExport(
      cl,
      varlist = list(
          "Ri_ESDChi",
          "truncated_L_j",
          "CHi_sq_stat",
          "truncated_L"))
    para_ful = function(ii, k1,
                        matbeta,
                        mat_cov_beta,
                        M,
                        k) {
      Ri_ESDChi(matbeta[ii, (k1 + 1):M],
                mat_cov_beta[(k1 + 1):M, (k1 + 1):M],
                M = M-k1, k = k)$Ri[1]
    }
    row_Rmat = clusterApply(
      cl,
      1:numreps,
      para_ful,
      k1 = k1,
      matbeta = betmat,
      mat_cov_beta = omega,
      M = M-k1,
      k = k
    )
    Rmat[, (k1 + 1)] = unlist(row_Rmat)
    stopCluster(cl)
    lft = 0
    rgt = max(Rmat[, (k1 + 1):k])
    bs_int = 0
    #bisection search for critical value
    while ((diff > tol) & (bs_int <= bs_max)) {
      mid = (lft + rgt) / 2
      cuttst = cbind(rep(mid, numreps), cutmat)
      idxmat = (Rmat[, (k1 + 1):k] > cuttst)

      f.mid =  mean(apply(idxmat, 1, any))
      if (f.mid <= alf) {
        rgt = mid
      }
      if (f.mid > alf) {
        lft = mid
      }
      diff = rgt - lft
      bs_int = bs_int + 1
      # cat(f.mid)
    }
    cutcur = mid
    cutmat = cuttst
  }

  #return result
  cutvec = cutmat[1, ]
  # browser()
  return(cutvec)

}


critival_value_GENlm_simu_intercept_v2 = function(
  M,
  k,
  omega,
  alf,
  trimfac,
  numreps,
  tol,
  bs_max = 100) {
  # M = number of units
  # k = max number of outlier units
  # omega = covariance matrix
  # alf = desired Type I error level
  # trimfac = percent of trimming on each side
  # numreps = number of simulation replications
  # tolerance for critical value search
  # BIG = outlier value
  # cutmax = upper bound for critical values
  # fac = factor for critical value search

  #PRELIMINARIES
  inputs = cbind(M, k, alf, trimfac, numreps, tol, BIG, cutmax, fac)
  cat('\n')
  print(noquote('inputs'))
  print(inputs)
  cat('\n')
  Rmat = matrix(0, numreps, k)
  cutmat = NULL

  #GENERATE SIMULATED BETAS
  betmat = rmvnorm(numreps, mean = rep(0, M), sigma = omega)

  #loop over k1 values
  # browser()
  for (k1 in (k - 1):0) {
    diff = 99
    ##parallel version
    require(parallel)
    cl = makeCluster(8)
    clusterExport(
      cl,
      varlist = list(
        "Ri_ESDChi",
        "truncated_L_j",
        "CHi_sq_stat",
        "truncated_L"))
    para_ful = function(ii, k1,
                        matbeta,
                        mat_cov_beta,
                        M,
                        k) {
      Ri_ESDChi(matbeta[ii, (k1 + 1):M],
                mat_cov_beta[(k1 + 1):M, (k1 + 1):M],
                M = M-k1, k = k)$Ri[1]
    }
    row_Rmat = clusterApply(
      cl,
      1:numreps,
      para_ful,
      k1 = k1,
      matbeta = betmat,
      mat_cov_beta = omega,
      M = M-k1,
      k = k
    )
    Rmat[, (k1 + 1)] = unlist(row_Rmat)
    stopCluster(cl)
    lft = 0
    rgt = max(Rmat[, (k1 + 1):k])
    bs_int = 0
    #bisection search for critical value
    while ((diff > tol) & (bs_int <= bs_max)) {
      mid = (lft + rgt) / 2
      cuttst = cbind(rep(mid, numreps), cutmat)
      idxmat = (Rmat[, (k1 + 1):k] > cuttst)

      f.mid =  mean(apply(idxmat, 1, any))
      if (f.mid <= alf) {
        rgt = mid
      }
      if (f.mid > alf) {
        lft = mid
      }
      diff = rgt - lft
      bs_int = bs_int + 1
      # cat(f.mid)
    }
    cutcur = mid
    cutmat = cuttst
  }

  #return result
  cutvec = cutmat[1, ]
  # browser()
  return(cutvec)

}
###calculate the test statistic
Ri_ESDChi_GENlm = function(ll_dat, k = 20)
{
  M = dim(ll_dat)[2] - 4 - 1
  Mdelta = floor(M * 0.2)
  ###fit a full model and determine the order
  obj_lm = lm(Y ~ 0 + .,
              data = ll_dat)
  vec_betahat = obj_lm$coefficients[5:(M + 4)]
  vec_order = order(vec_betahat, decreasing = T)
  vec_Ri = rep(NA, k)
  for (l in 0:(k - 1))
  {
    if (l > 0)
    {
      vec_delete = c()
      for (j in 1:l)
      {
        ##append each column
        vec_delete = c(vec_delete, which(ll_dat[, 4 + which(vec_order == j)] == 1))
      }
      ##remove the rows
      ll_datl = ll_dat[-vec_delete,]
      ##remove the columns
      ll_datl = ll_datl[,-(which(vec_order <= l) + 4)]
    } else{
      ll_datl = ll_dat
    }
    mat_X = as.matrix(ll_datl[, 1:(M - l + 4)])
    ###fit the lm with remaining covariates
    obj_lm = lm(Y ~ 0 + .,
                data = ll_datl)
    mat_conv = solve(t(mat_X) %*% mat_X)[5:(4 + M - l), 5:(4 + M - l)] *
      (sum(obj_lm$residuals ^ 2) / obj_lm$df.residual)
    vec_betahat = obj_lm$coefficients[5:(M - l + 4)]
    # vec_temp = rep(NA, M - l)
    # for(j in 1:length(vec_temp))
    # {
    #   vec_L = truncated_L_j(rank(vec_betahat),j, M = M - l, Mdelta = Mdelta)
    #   vec_temp[j] = CHi_sq_stat(
    #     vecL = vec_L,
    #     vec_beta = vec_betahat,
    #     mat_cov_beta = mat_conv
    #   )
    # }
    # vec_Ri[l + 1] = max(vec_temp)
    vec_Ri[l + 1] = max(Ri_ESDChi_Hl(vec_betahat, mat_conv))
  }
  # browser()
  return(vec_Ri)
}

GEEfit <- function(my.data,
                   Y.name,
                   X.name,
                   Test_reviewer,
                   id,
                   corstr,
                   trunc.prop = NA) {
  require(gee)
  #### First stage no intercept linear regression
  regressors      <- c(X.name, Test_reviewer)
  X.length        <- length(X.name)
  Tester.length   <- length(Test_reviewer)
  Full.length     <- length(regressors)
  regressor.part  <- paste0(regressors, collapse = "+")
  formula.full    <-
    as.formula(paste0(Y.name, "~ -1 + ", regressor.part))
  gee.fit         <-
    gee(
      formula = formula.full,
      id = id,
      corstr = corstr,
      data = my.data,
      silent = T
    )
  gee.fit
}

GEEfit_v2 <- function(my.data,
                      Y.name,
                      X.name,
                      Test_reviewer,
                      id,
                      corstr,
                      trunc.prop = NA) {
  require(geepack)
  #### First stage no intercept linear regression
  regressors      <- c(X.name, Test_reviewer)
  X.length        <- length(X.name)
  Tester.length   <- length(Test_reviewer)
  Full.length     <- length(regressors)
  regressor.part  <- paste0(regressors, collapse = "+")
  formula.full    <-
    as.formula(paste0(Y.name, "~ -1 + ", regressor.part))
  gee.fit         <-
    geeglm(
      formula = formula.full,
      id = id,
      corstr = corstr,
      data = my.data
    )
  gee.fit
}

###############################################################
##used for unified lambda
type1error = function(num_lambda, mat_Ri) {
  mat_diff = mat_Ri - num_lambda
  sum(apply(mat_diff, 1, function(x)
    any(x > 0))) / dim(mat_Ri)[1]
}
