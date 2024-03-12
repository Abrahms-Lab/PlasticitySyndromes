# Code for Johansson et al. 2024 -- lay date/mate switch DHGLM

# 0. Set up
# 1. Model definition
# 2. Data prep and model fitting
# 3. Extact model results
# 4. Reproductive success modeling 
# 5. DHGLM figures 

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
library(DescTools)
library(ggeffects)
library(ggpubr)
library(interactions)
library(sandwich)
library(deming)
library(broom.mixed)
library(lme4)
library(MuMIn)

# load data
lay.mate.df <- read.csv("lay.mate.dhglm.df")
LRS.df <- read_csv("LRS.df") 
annualRS.df <- read_csv("annualRS.df") 

#############################
# 1. Model definition
#############################

lay.mate.dhglm <- as.character("
  data {
    int N_obs;                    // number of observations
    int N_peng;                   // number of penguins
    real N_peng2;                 // number of penguins
    int K;                        // number of fixed effects
    int pengID[N_obs];            // penguin id vector
    int N_yr;                     // Number of years in dataset
    int yrID[N_obs];              // year id vector
    real sam[N_obs];          // sam predictor for lay
    real rain[N_obs];         // sam predictor for mate
    real center_year_lay[N_obs];  // fixed effect for lay
    int prev_fledge[N_obs];       // fixed effect for mate
    int ss2[N_obs];               // fixed effect for mate
    int y[N_obs];                 // response values
    int mm[N_obs];                // which model component (0=lay date, 1=mate choice)
  }
  
  parameters {
    vector[K] beta;                 // Fixed effects
    vector<lower=0>[4] sd_peng;     // sd for random bird effects 
    vector<lower=0>[2] sd_resid;    // residual sigma
    vector<lower=0>[2] sd_yr;       // sd for random yr effects
    
    matrix[2, N_yr] z_yr;                  // matrix of random yr effects (z-scores)
    cholesky_factor_corr[2] cor_yr;        // cholesky corr matrix for rand yr effs


    // Random penguin effects
    matrix[4, N_peng] z_peng;               // matrix of int/slope deviates (z-scores) per bird
    cholesky_factor_corr[4] cor_peng;       // Cholesky corr matrix for rand bird effs
  }
  
  transformed parameters {
    matrix[4, N_peng] d_peng;                                 // Rand bird effs per bird
    matrix[2, N_yr] d_yr;                                     // Random yr effects

    d_yr = diag_pre_multiply(sd_yr, cor_yr) * z_yr;
    d_peng = diag_pre_multiply(sd_peng, cor_peng) * z_peng;   // Converts z-scores to scaled deviates
    
}
  
  model {
    vector[N_obs] mu;                        // vector of expected values

    // priors
    beta ~ normal(0, 100);                    // Prior for fixed effects
    
    // Random effects (sigma)
    sd_resid ~ normal(0, 100);                // Prior for residual SD
    sd_peng ~ normal(0, 100);                 // Prior for among bird SD
    sd_yr ~ normal(0, 100);                   // Prior for among yr SD
    
    // Correlation matrices
    cor_peng ~ lkj_corr_cholesky(1);          // Prior for correlation matrix
    cor_yr ~ lkj_corr_cholesky(1);            // Prior for correlation matrix

    // Random effects expressed as z-scores
    to_vector(z_peng) ~ normal(0, 1);         // Prior for penguin z-score deviates
    to_vector(z_yr) ~ normal(0, 1);           // Prior for year z-score deviates


    // likelihood
    for(i in 1:N_obs) {
      if(mm[i] == 0) {
        // Sampling statement for lay date model
        mu[i] = beta[1] + beta[5]*center_year_lay[i] + beta[2]*sam[i] + d_yr[1,yrID[i]] + d_peng[1,pengID[i]] + d_peng[2,pengID[i]]*sam[i];
        y[i] ~ normal(mu[i],sd_resid[1]);
      }
      else {
        // Sampling statement for mate choice model
        mu[i] = beta[3] + beta[4]*rain[i] + beta[6]*ss2[i] + beta[7]*prev_fledge[i] + d_yr[2,yrID[i]] + d_peng[3, pengID[i]] + d_peng[4, pengID[i]]*rain[i];
        y[i] ~ bernoulli_logit(mu[i]);
      }
    }
  }
  
  generated quantities {
    matrix[4,4] rho_peng;               // rho = corr matrix 
    matrix[2,2] rho_yr;
    matrix[4,N_peng] re_peng;           // Penguin Random effects, including popn effect

    rho_peng = multiply_lower_tri_self_transpose(cor_peng);   // Converts the chol corr matrix to corr matrix
    rho_yr = multiply_lower_tri_self_transpose(cor_yr);       // Converts the chol corr matrix to corr matrix
    re_peng[1,] = d_peng[1,] + beta[1];                       // Adds population level estimate
    re_peng[2,] = d_peng[2,] + beta[2];
    re_peng[3,] = d_peng[3,] + beta[3];
    re_peng[4,] = d_peng[4,] + beta[4];
  }
")

#############################
# 2. Data prep and model fitting
#############################
obs_lay_data <- lay.mate.df %>% filter(behavior=='ld')
obs_mate_data <- lay.mate.df %>% filter(behavior=='ms')

mm <- c(rep(0,nrow(obs_lay_data)), rep(1,nrow(obs_mate_data)))
y <- c(obs_lay_data$response, obs_mate_data$response)
sam <- c(obs_lay_data$zscore_SAM, rep(0,nrow(obs_mate_data)))
rain <- c(rep(0,nrow(obs_lay_data)), (obs_mate_data$log_rain_per_day))
yr.orig <- c(obs_lay_data$center_year, obs_mate_data$center_year)
yrID <- as.numeric(as.factor(yr.orig))
N_yr <- length(unique(yrID))
pengID.orig <- c(obs_lay_data$penguinseq,obs_mate_data$penguinseq)
pengID <- as.numeric(as.factor(pengID.orig))
pengID <- as.numeric(factor(pengID, levels = unique(pengID)))
center_year_lay <- c(obs_lay_data$center_year, rep(0,nrow(obs_mate_data)))
ss2 <- c(rep(0,nrow(obs_lay_data)), (obs_mate_data$StaySwitch2)) 
prev_fledge <- c(rep(0,nrow(obs_lay_data)), (obs_mate_data$prev_fledge)) 
K <- 7
N_peng <- max(pengID)
N_obs <- length(y)

# combine to a list
mod.in <- list(mm=mm, y=y, rain=rain, sam=sam, yrID=yrID, N_yr=N_yr, pengID=pengID, K=K, N_peng=N_peng, N_peng2=N_peng,
               N_obs=N_obs, center_year_lay=center_year_lay, ss2=ss2, prev_fledge=prev_fledge)

# Model run values
nchain = 4             # Number of MCMC chains
nwarm = 2000           # Number of warm-up iterations
niter = 4000          # Total number of iterations
nperm <- niter-nwarm   # Number of estimation iterations

# fit model
bimod.data.real <- stan(model_code = m_lay.sam.Ryr.Rbird.Rsam_mate.Ryr.Rbird.Rrain, 
                        data = mod.in, 
                        chains = nchain, warmup = nwarm, iter = niter,
                        cores = 4, verbose = FALSE, refresh=50)

#############################
# 3. Extact model results
#############################

# set parameters
params_vec <- c(1:15,17:19,22:23,27)
b1 <- "lay"
b2 <- "mate"
fixed_efs <- c( "center_year_lay","ss2","prev_fledge")
obs_b1_data<-obs_lay_data
obs_b2_data<-obs_mate_data

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

# several columns derived from correlation matrices are redundant
params <- params[,params_vec]
# name columns
names(params) <- c("b1_Int","b1_Slope","b2_Int","b2_Slope",fixed_efs,
                   "sdPeng.b1_Int","sdPeng.b1_Slope","sdPeng.b2_Int","sdPeng.b2_Slope",
                   "sdYr.lay","sdYr.mate",
                   "sdRes.lay","sdRes.mate",
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
# 4. Reproductive success modeling 
#############################

## long-term reproductive success model
# set data frame
lrs_model_df <- 
  left_join(bird.Rand.abs_v2,LRS.df,by = c("ID"="penguinseq")) %>% 
  filter(has.data.b1,has.data.b2) %>%
  rowwise() %>%
  mutate(joint_slope = CartToPol(b1_Slope.mean,b2_Slope.mean)$r) # joint plasticity metric from polar conversion

# fit model
lrs_model <- lm(f_yr ~ joint_slope, 
                weights = yrs,
                data=lrs_model_df)
summary(lrs_model)

# plot model
newx = seq(-1,4,by = 0.05)
conf_interval <- predict(lrs_model, newdata=data.frame(joint_slope=newx), interval="confidence",
                         level = 0.95)
ggplot() +
  geom_point(data=lrs_model_df,aes(joint_slope,f_yr),alpha=.5,stroke=0) +
  geom_ribbon(aes(newx, ymin=conf_interval[,2], ymax=conf_interval[,3]),alpha=.3) + 
  geom_abline(intercept = coef(lrs_model)[1], slope = coef(lrs_model)[2]) + 
  xlim(c(.25,2.75)) + ylab("Long-Term Reproductive Success") + xlab("Plasticity") +
  theme_classic() + theme(legend.position = "none") 

## reproductive success vs SAM
# set data frame
sam_rs_model_df <- left_join(annualRS.df,bird.Rand.abs_v2 %>% 
                     dplyr::select(ID,
                                   b1_Int.mean,b2_Int.mean,
                                   b1_Slope.mean,b2_Slope.mean,b1_Slope.sd,
                                   b2_Slope.sd,has.data.b1,has.data.b2),
                   by=c("penguinseq"="ID")) %>% 
  filter(has.data.b1,has.data.b2) %>%
  mutate(NFledged_01 = as.integer(NFledged>0))

# fit model with random effect for ID
SAM_rs_glmer <- glmer(NFledged_01 ~ SAM_quant + (1|penguinseq),
                      data=sam_rs_model_df,
                      family=binomial(link = "logit"),
                      control=glmerControl(optimizer = "bobyqa"))
summary(SAM_rs_glmer)

# fit model without random effect for plotting
SAM_RS_glm <- glm(NFledged_01 ~ SAM_quant,
                  data=sam_rs_model_df,
                  family=binomial(link = "logit"))
summary(SAM_RS_glm)

# plot model
preds <- data.frame(SAM_quant = seq(-3, 3.5, length.out = 100))
preds$pred <- predict(SAM_RS_glm, preds, type = "response")
preds$upper <- predict(SAM_RS_glm, preds, type = "response", se.fit = TRUE)$fit + 1.96 * predict(SAM_RS_glm, preds, type = "response", se.fit = TRUE)$se.fit
preds$lower <- predict(SAM_RS_glm, preds, type = "response", se.fit = TRUE)$fit - 1.96 * predict(SAM_RS_glm, preds, type = "response", se.fit = TRUE)$se.fit
ggplot(sam_rs_model_df, aes(x = SAM_quant, y = NFledged_01)) +
  geom_point(alpha = 0.2) +
  geom_line(data = preds, aes(x = SAM_quant, y = pred), color = "black", inherit.aes = FALSE) +
  geom_ribbon(data = preds, aes(x = SAM_quant, ymin = lower, ymax = upper), alpha = 0.2, inherit.aes = FALSE) +
  theme_classic()

## interaction between plasticity and SAM
sam_rs_model_df$SAM_qual5 <- factor(sam_rs_model_df$SAM_qual5,levels = c("norm","high","low","other"))

# set data frame
b1_b2_df <- sam_rs_model_df %>%
  filter(has.data.b1&has.data.b2) %>%
  mutate(NFledged_01 = as.integer(NFledged>0),
         b1_b2_Slope = CartToPol(b1_Slope.mean,b2_Slope.mean)$r)

# fit interaction model
RS.fit_b1_b2 <- glmer(NFledged_01 ~ b1_b2_Slope + SAM_qual5 + b1_b2_Slope*SAM_qual5 +  (1|penguinseq),
                      data=b1_b2_df,
                      family=binomial(link = "logit"),
                      control=glmerControl(optimizer = "bobyqa"))
summary(RS.fit_b1_b2)

# plot model
interact_plot(RS.fit_b1_b2, 
              pred = b1_b2_Slope, 
              modx = SAM_qual5, 
              modxvals = c("low"),
              interval = T,
              plot.points = T, colors = "#009988") +
  theme_classic() + ylab("Successful Fledging") + xlab("Plasticity")+ theme(legend.position = "none") +
  xlim(.4,2.5)
interact_plot(RS.fit_b1_b2, 
              pred = b1_b2_Slope, 
              modx = SAM_qual5, 
              modxvals = c("norm"),
              interval = T,
              plot.points = T, colors = "#33bbee") +
  theme_classic() + ylab("Successful Fledging") + xlab("Plasticity")+ theme(legend.position = "none")+
  xlim(.4,2.5)
interact_plot(RS.fit_b1_b2, 
              pred = b1_b2_Slope, 
              modx = SAM_qual5, 
              modxvals = c("high"),
              interval = T,
              plot.points = T, colors = "#ee3377") +
  theme_classic() + ylab("Successful Fledging") + xlab("Plasticity") + theme(legend.position = "none")+
  xlim(.4,2.5)

#############################
# 5. DHGLM figures 
#############################

## plot results for each behavior- highlight individuals 
temp.b1 <- sam_rs_model_df %>% 
  filter(has.data.b1) %>%
  dplyr::select(penguinseq,
         b1_Int.mean,b1_Slope.mean) %>% 
  distinct()

temp.b2 <- sam_rs_model_df %>% 
  filter(has.data.b2) %>%
  dplyr::select(penguinseq,
         b2_Int.mean,b2_Slope.mean) %>% 
  distinct()

l.d.low_index.b1 <- which(temp.b1$penguinseq==46877) # individual to highlight
l.d.high_index.b1 <- which(temp.b1$penguinseq==46039) # individual to highlight

l.m.low_index.b1 <- which(temp.b1$penguinseq==62482) # individual to highlight
l.m.high_index.b1 <- which(temp.b1$penguinseq==9820) # individual to highlight

l.m.low_index.b2 <- which(temp.b2$penguinseq==62482) # individual to highlight
l.m.high_index.b2 <- which(temp.b2$penguinseq==9820) # individual to highlight


# lay day
p <- ggplot() + geom_point(aes(
  c(min(obs_b1_data$zscore_SAM),max(obs_b1_data$zscore_SAM)),
  c(min(temp.b1$b1_Int.mean),max(temp.b1$b1_Int.mean)+1)),
  alpha=0)

for(i in 1:nrow(temp.b1)){
  p <- p +
    geom_abline(intercept = temp.b1$b1_Int.mean[i],
                slope = temp.b1$b1_Slope.mean[i], alpha=.03)
}
p + 
  geom_abline(intercept = temp.b1$b1_Int.mean[l.m.low_index.b1],
              slope = temp.b1$b1_Slope.mean[l.m.low_index.b1], color="#EE6677",lwd=1.5) + 
  geom_abline(intercept = temp.b1$b1_Int.mean[l.m.high_index.b1],
              slope = temp.b1$b1_Slope.mean[l.m.high_index.b1], color="#AA3377",lwd=1.5) + 
  
  geom_abline(intercept = temp.b1$b1_Int.mean[l.d.low_index.b1],
              slope = temp.b1$b1_Slope.mean[l.d.low_index.b1], color="#ccbb44",lwd=1.5) + 
  geom_abline(intercept = temp.b1$b1_Int.mean[l.d.high_index.b1],
              slope = temp.b1$b1_Slope.mean[l.d.high_index.b1], color="#228833",lwd=1.5) + 
  
  ylab("Lay Date Timing") +
  xlab("SAM") + 
  theme_classic()


# mate switching
p <- ggplot() + geom_point(aes(
  c(min(obs_b2_data$log_rain),max(obs_b2_data$log_rain)),
  c(min(temp.b2$b2_Int.mean)-1,max(temp.b2$b2_Int.mean)+1)),
  alpha=0)

for(i in 1:nrow(temp.b2)){
  p <- p +
    geom_abline(intercept = temp.b2$b2_Int.mean[i],
                slope = temp.b2$b2_Slope.mean[i], alpha=.03)
}
p + 
  geom_abline(intercept = temp.b2$b2_Int.mean[l.m.low_index.b2],
              slope = temp.b2$b2_Slope.mean[l.m.low_index.b2], color="#EE6677",lwd=1.5) + 
  geom_abline(intercept = temp.b2$b2_Int.mean[l.m.high_index.b2],
              slope = temp.b2$b2_Slope.mean[l.m.high_index.b2], color="#AA3377",lwd=1.5) + 
  ylab("Mate Switch Probability") +
  xlab("Rainfall") + 
  theme_classic()

## b1 plasticity vs b2 plasticity with deming best-fit line
dem_df <- bird.Rand.abs_v2 %>% filter(has.data.b1&has.data.b2)

dem_fit <- deming(b2_Slope.mean ~ b1_Slope.mean, data = dem_df,
                  xstd = b1_Slope.sd, ystd = b2_Slope.sd)

low <- dem_df[which(dem_df$ID==62482),]
high <- dem_df[which(dem_df$ID==9820),]

ggplot(dem_df, aes(x=b1_Slope.mean,
                 y=b2_Slope.mean)) +
  geom_point(alpha=0.3,stroke=0,size=2) + 
  geom_point(data=low,aes(b1_Slope.mean,b2_Slope.mean),alpha=1,size=3,color="#EE6677") + 
  geom_point(data=high,aes(b1_Slope.mean,b2_Slope.mean),alpha=1,size=3,color="#AA3377") + 
  geom_abline(intercept=dem_fit$coefficients[[1]],
              slope=dem_fit$coefficients[[2]]) + 
  theme_classic() +
  theme(legend.position = "none") +
  ylab("Mate Swtiching Plasticity") + xlab("Laying Date Plasticity")

