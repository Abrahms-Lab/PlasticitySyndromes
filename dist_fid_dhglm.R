# Code for Johansson et al. 2024 -- distance/fidelity DHGLM

# 0. Set up
# 1. Model definition
# 2. Data prep and model fitting
# 3. Extact model results

#############################
# 0. Set up
#############################

# load packages
library(rstan)
library(mgcv)
library(MASS)
library(Matrix)
library(tidyverse)
library(matrixStats)
library(lubridate)
library(bayestestR)
library(deming)

# load data
dist.fid.df <- read.csv("dist.fid.dhglm.df")

#############################
# 1. Model definition
#############################

dist.fid.dhglm <- as.character("
  data {
    int N_obs;                   // number of observations
    int N_peng;                  // number of penguins
    real N_peng2;
    int K;                       // number of fixed effects
    int pengID[N_obs];           // penguin id vector
    int N_yr;                    // Number of years in dataset
    int yrID[N_obs];             // year id vector
    real sst_f[N_obs];           // sst predictor fidelity
    real sst_d[N_obs];		       // sst predictor distance
    int sex_f[N_obs];            // fixed effect for fideltiy
    int sex_d[N_obs];            // fixed effect for distance
    real y[N_obs];               // response values
    int mm[N_obs];               // which model component (0=trip distance, 1=trip fidelity)
  }
  
  parameters {
    vector[K] beta;                         // Fixed effects
    vector<lower=0>[4] sd_peng;             // sd for random bird effects 
    vector<lower=0>[2] sd_resid;            // residual sigma
    vector<lower=0>[2] sd_yr;               // sd for random yr effects
    
    // Random year effects
    matrix[2, N_yr] z_yr;                   // matrix of random yr effects (z-scores)
    cholesky_factor_corr[2] cor_yr;         // cholesky corr matrix for rand yr effs
    
    // Random penguin effects
    matrix[4, N_peng] z_peng;               // matrix of int/slope deviates (z-scores) per bird
    cholesky_factor_corr[4] cor_peng;       // cholesky corr matrix for rand bird effs
  }
  
  transformed parameters {
    matrix[4, N_peng] d_peng;                                 // Rand bird effs per bird
    matrix[2, N_yr] d_yr;                                     // Random yr effects
    d_yr = diag_pre_multiply(sd_yr, cor_yr) * z_yr;
    d_peng = diag_pre_multiply(sd_peng, cor_peng) * z_peng;   // Converts z-scores to scaled deviates

  }
  
  model {
    vector[N_obs] mu;                           // vector of expected values
    
    // priors
    beta ~ normal(0, 100);                      // Prior for fixed effects
    
    // Random effects (sigma)
    sd_resid ~ normal(0, 100);                  // Prior for residual SD
    sd_peng ~ normal(0, 100);                   // Prior for among bird SD
    sd_yr ~ normal(0, 100);                     // Prior for among yr SD
    
    // Correlation matrices
    cor_peng ~ lkj_corr_cholesky(1);            // Prior for correlation matrix
    cor_yr ~ lkj_corr_cholesky(1);              // Prior for correlation matrix 
    
    // Random effects expressed as z-scores
    to_vector(z_peng) ~ normal(0, 1);           // Prior for penguin z-score deviates
    to_vector(z_yr) ~ normal(0, 1);             // Prior for year z-score deviates
    
    // likelihood
    for(i in 1:N_obs) {
      if(mm[i] == 0) {
        // Sampling statement for distance model
        mu[i] = beta[1] + beta[5]*sex_d[i] + beta[2]*sst_d[i] + d_yr[1,yrID[i]] + d_peng[1,pengID[i]] + d_peng[2,pengID[i]]*sst_d[i];
        y[i] ~ normal(mu[i],sd_resid[1]);
      }
      else {
        // Sampling statement for fidelity model
        mu[i] = beta[3] + beta[6]*sex_f[i] + beta[4]*sst_f[i] + d_yr[2,yrID[i]] + d_peng[3,pengID[i]] + d_peng[4,pengID[i]]*sst_f[i]; 
        y[i] ~ normal(mu[i],sd_resid[2]);
      }
    }
  }
  
  generated quantities {
    matrix[4,4] rho_peng;                 // rho = corr matrix 
    matrix[2,2] rho_yr;
    matrix[4,N_peng] re_peng;             // Penguin Random effects, including popn effect

    rho_peng = multiply_lower_tri_self_transpose(cor_peng);     // Converts the chol corr matrix to corr matrix
    rho_yr = multiply_lower_tri_self_transpose(cor_yr);         // Converts the chol corr matrix to corr matrix
    re_peng[1,] = d_peng[1,] + beta[1];                         // Adds population level estimate 
    re_peng[2,] = d_peng[2,] + beta[2];
    re_peng[3,] = d_peng[3,] + beta[3];
    re_peng[4,] = d_peng[4,] + beta[4];
  }
")

#############################
# 2. Data prep and model fitting
#############################
obs_dist_data <- dist.fid.df %>% filter(behavior=='td')
obs_fid_data <- dist.fid.df %>% filter(behavior=='sf')

# load data into variables
mm <- c(rep(0,nrow(obs_dist_data)), rep(1,nrow(obs_fid_data)))
y <- c(scale(obs_dist_data$response), scale(obs_fid_data$response))
sst_d <- c(obs_dist_data$zscore_sst, rep(0,nrow(obs_fid_data)))
sst_f <- c(rep(0,nrow(obs_dist_data)), (obs_fid_data$zscore_sst)) 
yr.orig <- c(obs_dist_data$center_year, obs_fid_data$center_year)
yrID <- as.numeric(as.factor(yr.orig))
N_yr <- length(unique(yrID))
pengID.orig <- c(obs_dist_data$penguinseq,obs_fid_data$penguinseq)
pengID <- as.numeric(as.factor(pengID.orig))
pengID <- as.numeric(factor(pengID, levels = unique(pengID)))
center_yr_d <- c(obs_dist_data$center_year, rep(0,nrow(obs_fid_data)))
center_year_f <- c(rep(0,nrow(obs_dist_data)), (obs_fid_data$center_year))
sex_d <- c((obs_dist_data$sex=="M"), rep(0,nrow(obs_fid_data)))
sex_f <- c(rep(F,nrow(obs_dist_data)), obs_fid_data$sex=="TRUE")
K <- 6
N_peng <- max(pengID)
N_obs <- length(y)

# combine to a list
mod.in <- list(mm=mm, y=y, sst_f=sst_f, sst_d=sst_d, yrID=yrID, N_yr=N_yr, pengID=pengID, K=K, N_peng=N_peng, N_peng2=N_peng,
               N_obs=N_obs, sex_d=sex_d, sex_f=sex_f)

# Model run values
nchain = 4              # Number of MCMC chains
nwarm = 4000            # Number of warm-up iterations
niter = 8000            # Total number of iterations
nperm <- niter-nwarm    # Number of estimation iterations

# fit model
bimod.data.real <- stan(model_code = dist.fid.dhglm, 
                        data = mod.in, 
                        chains = nchain, warmup = nwarm, iter = niter,
                        cores = 4, verbose = FALSE, refresh=50)

#############################
# 3. Extact model results
#############################

# set parameters
b1 <- "dist"
b2 <- "fid"
params_vec <- c(1:14,16:18,21:22,26) 
fixed_efs <- c("day_d","sex_f")
obs_b1_data <- obs_dist_data
obs_b2_data <- obs_fid_data

# Extract MCMC estimates for key parameters
mcmc.pars <- rstan::extract(bimod.data.real, permuted=FALSE, inc_warmup=FALSE, 
                            pars=c("beta", "sd_peng","sd_yr","sd_resid", "rho_peng"))

# format into a single 2-d data frame
npar <- dim(mcmc.pars)[3]
params <- matrix(0, nrow=nchain*nperm, ncol=npar)

for(i in 1:npar){
  tmp1 <- mcmc.pars[,,i]
  parfill <- tmp1[,1]
  for(j in 2:nchain){
    parfill <- c(parfill, tmp1[,j])
  }
  params[,i] <- parfill
}
params <- data.frame(params)

head(params)

# several columns derived from correlation matrices are redundant
params <- params[,params_vec]

# name columns
names(params) <- c("b1_Int","b1_Slope","b2_Int","b2_Slope",fixed_efs,
                   "sdPeng.b1_Int","sdPeng.b1_Slope","sdPeng.b2_Int","sdPeng.b2_Slope",
                   "sdYr.lay","sdYr.trip",
                   "sdRes.lay","sdRes.trip",
                   "cor.b1_Int.b1_Slope","cor.b1_Int.b2_Int","cor.b1_Int.b2_Slope",
                   "cor.b1_Slope.b2_Int","cor.b1_Slope.b2_Slope",
                   "cor.b2_Int.b2_Slope")

# print metrics for all parameters
for(i in 1:ncol(params)){
  post <- params[,i]
  print(colnames(params)[i])
  print(round(median(post),digits=3))
  print(p_direction(post))
  cat("\n")
  print(ci(post, method = "HDI", ci = 0.89))
  cat("\n")
}

# Extract bird-specific effects
re.peng.est <- rstan::extract(bimod.data.real, permuted=FALSE, inc_warmup=FALSE, 
                              pars=c("re_peng"))

local.re.peng.est <- re.peng.est
bird.effs <- matrix(0, nrow=dim(local.re.peng.est)[3], ncol=nchain*nperm)

for(i in 1:dim(local.re.peng.est)[3]){
  tmp <- local.re.peng.est[,,i]   
  bird.effs[i,] <- c(tmp[,1],tmp[,2],tmp[,3],tmp[,4])
}

# inverse transformation of slope estimates
all_slopes <- bird.effs[c(seq(from=2,to=nrow(bird.effs),by=2)),]
all_slopes_v2 <- all_slopes
for(i in 1:nrow(all_slopes_v2)){
  if(mean(all_slopes_v2[i,])<0){
    all_slopes_v2[i,] <- (-1) * all_slopes_v2[i,]
  }
}
bird.effs[c(seq(from=2,to=nrow(bird.effs),by=2)),] <- all_slopes_v2

cor_v2 <-c()
for(i in 1:ncol(all_slopes_v2)){
  cor_v2 <- c(cor_v2,cor(all_slopes_v2[seq(from=1,to=nrow(all_slopes_v2),by=2),i],
                         all_slopes_v2[seq(from=2,to=nrow(all_slopes_v2),by=2),i]))
}

# reduce this to mean and 95% CI
bird.REs <- data.frame(mod.ID=rep(1:N_peng, each=4),
                       response=rep(c("b1_Int","b1_Slope","b2_Int","b2_Slope"), N_peng),
                       stan.mean=rowMeans(bird.effs), 
                       stan.sd=rowSds(bird.effs),
                       stan.lwr=rowQuantiles(bird.effs, probs=0.025),
                       stan.upr=rowQuantiles(bird.effs, probs=0.975))

# identify penguin ID match ups
pengID.lookup <- distinct(data.frame(ID=pengID.orig,mod.ID=pengID))

# Merge with look up table
bird.REs <- merge(x=bird.REs, y=pengID.lookup, by="mod.ID", all.x=TRUE)

# Identify whether the corresponding data was present
bird.REs$has.data.b1 <- bird.REs$ID %in% obs_b1_data$penguinseq
bird.REs$has.data.b2 <- bird.REs$ID %in% obs_b2_data$penguinseq

# reformat data
tmp.b1I <- bird.REs[bird.REs$response=="b1_Int",c("ID","mod.ID","stan.mean","stan.sd","stan.lwr","stan.upr",
                                                  "has.data.b1")]
names(tmp.b1I) <- c("ID","mod.ID","b1_Int.mean","b1_Int.sd","b1_Int.lwr","b1_Int.upr","has.data.b1")

tmp.b1S <- bird.REs[bird.REs$response=="b1_Slope",c("stan.mean","stan.sd","stan.lwr","stan.upr")]
names(tmp.b1S) <- c("b1_Slope.mean","b1_Slope.sd","b1_Slope.lwr","b1_Slope.upr")

tmp.b2I <- bird.REs[bird.REs$response=="b2_Int",c("stan.mean","stan.sd","stan.lwr","stan.upr",
                                                  "has.data.b2")]
names(tmp.b2I) <- c("b2_Int.mean","b2_Int.sd","b2_Int.lwr","b2_Int.upr","has.data.b2")

tmp.b2S <- bird.REs[bird.REs$response=="b2_Slope",c("stan.mean","stan.sd","stan.lwr","stan.upr")]
names(tmp.b2S) <- c("b2_Slope.mean","b2_Slope.sd","b2_Slope.lwr","b2_Slope.upr")

bird.Rand <- cbind(tmp.b1I,tmp.b1S,tmp.b2I,tmp.b2S)

bird.Rand.abs_v2 <- bird.Rand

full.Rand <- bird.Rand

# Change STAN estimates when information wasn't present to NA
full.Rand[full.Rand$has.data.b1==F,
          c("b1_Int.mean","b1_Int.sd","b1_Int.lwr","b1_Int.upr","b1_Slope.mean","b1_Slope.sd","b1_Slope.lwr","b1_Slope.upr")] <- NA
full.Rand[full.Rand$has.data.b2==F,
          c("b2_Int.mean","b2_Int.sd","b2_Int.lwr","b2_Int.upr","b2_Slope.mean","b2_Slope.sd","b2_Slope.lwr","b2_Slope.upr")] <- NA

full.Rand.abs_v2 <- full.Rand

# metrics of transformed slope correlation estimates
# median
median(cor_v2)
# HDI
ci(cor_v2, method = "HDI", ci = 0.89)
# probability of direction
p_direction(cor_v2)

#############################
# 4. DHGLM figures 
#############################

## plot results for each behavior- highlight individuals 
temp.b1 <- full.Rand.abs_v2 %>% 
  filter(has.data.b1) %>%
  dplyr::select(ID,b1_Int.mean,b1_Slope.mean) %>% 
  distinct()

temp.b2 <- full.Rand.abs_v2 %>% 
  filter(has.data.b2) %>%
  dplyr::select(ID,b2_Int.mean,b2_Slope.mean) %>% 
  distinct()

l.d.low_index.b1 <- which(temp.b1$ID==46877)
l.d.high_index.b1 <- which(temp.b1$ID==46039)

d.f.low_index.b1 <- which(temp.b1$ID==59168)
d.f.high_index.b1 <- which(temp.b1$ID==51522)

d.f.low_index.b2 <- which(temp.b2$ID==59168)
d.f.high_index.b2 <- which(temp.b2$ID==51522)

# foraging trip distaince
p <- ggplot() + geom_point(aes(
  c(min(obs_b1_data$zscore_sst),max(obs_b1_data$zscore_sst)),
  c(min(temp.b1$b1_Int.mean),max(temp.b1$b1_Int.mean)+.7)),
  alpha=0)

for(i in 1:nrow(temp.b1)){
  p <- p +
    geom_abline(intercept = temp.b1$b1_Int.mean[i],
                slope = temp.b1$b1_Slope.mean[i], alpha=.05)
}
p + 
  geom_abline(intercept = temp.b1$b1_Int.mean[l.d.low_index.b1],
              slope = temp.b1$b1_Slope.mean[l.d.low_index.b1], color="#ccbb44",lwd=1.5) + 
  geom_abline(intercept = temp.b1$b1_Int.mean[l.d.high_index.b1],
              slope = temp.b1$b1_Slope.mean[l.d.high_index.b1], color="#228833",lwd=1.5) + 
  
  geom_abline(intercept = temp.b1$b1_Int.mean[d.f.low_index.b1],
              slope = temp.b1$b1_Slope.mean[d.f.low_index.b1], color="#66CCEE",lwd=1.5) + 
  geom_abline(intercept = temp.b1$b1_Int.mean[d.f.high_index.b1],
              slope = temp.b1$b1_Slope.mean[d.f.high_index.b1], color="#4477AA",lwd=1.5) + 
  
  ylab("Foraging Distance") +
  xlab("SST") + 
  theme_classic()

# foraging site selection
p <- ggplot() + geom_point(aes(
  c(min(obs_b2_data$zscore_sst),max(obs_b2_data$zscore_sst)),
  c(min(temp.b2$b2_Int.mean)-.5,max(temp.b2$b2_Int.mean)+.5)),
  alpha=0)

for(i in 1:nrow(temp.b2)){
  p <- p +
    geom_abline(intercept = temp.b2$b2_Int.mean[i],
                slope = temp.b2$b2_Slope.mean[i], alpha=.05)
}
p + 
  geom_abline(intercept = temp.b2$b2_Int.mean[d.f.low_index.b2],
              slope = temp.b2$b2_Slope.mean[d.f.low_index.b2], color="#66CCEE",lwd=1.5) + 
  geom_abline(intercept = temp.b2$b2_Int.mean[d.f.high_index.b2],
              slope = temp.b2$b2_Slope.mean[d.f.high_index.b2], color="#4477AA",lwd=1.5) + 
  ylab("Foraging Site Selection") +
  xlab("SST") + 
  theme_classic()

## b1 plasticity vs b2 plasticity with deming best-fit line
dem_df <- bird.Rand.abs_v2 %>% filter(has.data.b1&has.data.b2)

dem_fit <- deming(b2_Slope.mean ~ b1_Slope.mean, data = dem_df,
                  xstd = b1_Slope.sd, ystd = b2_Slope.sd)

low <- dem_df[which(dem_df$ID==59168),]
high <- dem_df[which(dem_df$ID==51522),]

ggplot(dem_df, aes(x=b1_Slope.mean,
                   y=b2_Slope.mean)) +
  geom_point(alpha=0.3,stroke=0,size=2) + 
  geom_point(data=low,aes(b1_Slope.mean,b2_Slope.mean),alpha=1,size=3,color="#66CCEE") + 
  geom_point(data=high,aes(b1_Slope.mean,b2_Slope.mean),alpha=1,size=3,color="#4477AA") + 
  geom_abline(intercept=dem_fit$coefficients[[1]],
              slope=dem_fit$coefficients[[2]]) + 
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Mate Swtiching Plasticity") + xlab("Laying Date Plasticity")
