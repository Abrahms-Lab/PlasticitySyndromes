# Code for Johansson et al. 2024 -- Single-behavior models

# 0. Set up
# 1. Egg laying date
# 2. Mate switching
# 3. Foraging trip distance
# 4. Foraging site fidelity 

#############################
# 0. Set up
#############################

# Import libraries 
library(tidyverse)
library(lme4)
library(MuMIn)

# Load data
lay.date.data <- read_csv("lay.date.data.csv")
mate.switch.data <- read_csv("mate.switch.data.csv")
trip.dist.data <- read_csv("trip.dist.data.csv")
trip.fid.data <- read_csv("trip.fid.data.csv")

#############################
# 1. Egg laying date
#############################

# random effects structure - model fitting
model.list <- list()
model.list[[1]] <- lm(LayDay ~ zscore_SAM + center_year, data=lay.date.data)
model.list[[2]] <- lmer(LayDay ~ zscore_SAM + center_year + (1|center_year), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[3]] <- lmer(LayDay ~ zscore_SAM + center_year + (zscore_SAM|penguinseq), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[4]] <- lmer(LayDay ~ zscore_SAM + center_year + (1|center_year) + (zscore_SAM|penguinseq), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))

# random effects structure - model selection
MuMIn::model.sel(model.list)

# fixed effects structure - model fitting 
model.list <- list()
model.list[[1]] <- lmer(LayDay ~ 1 + (1|center_year) + (zscore_SAM|penguinseq), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[2]] <- lmer(LayDay ~ zscore_SAM + (1|center_year) + (zscore_SAM|penguinseq), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[3]] <- lmer(LayDay ~ center_year + (1|center_year) + (zscore_SAM|penguinseq), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[4]] <- lmer(LayDay ~ zscore_SAM + center_year + (1|center_year) + (zscore_SAM|penguinseq), data=lay.date.data,control = lmerControl(optimizer = "bobyqa"))

# fixed effects structure - model selection
MuMIn::model.sel(model.list)


#############################
# 2. Mate switching
#############################

# random effects structure - model fitting
model.list <- list()
model.list[[1]] <- glm(StaySwitch ~ center_year + log_rain_per_day + StaySwitch2 + prev_fledge,
                       data = mate.switch.data, family = binomial(link = "logit"))
model.list[[2]] <- glmer(StaySwitch ~ center_year + log_rain_per_day + StaySwitch2 + prev_fledge + (1|center_year),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[3]] <- glmer(StaySwitch ~ center_year + log_rain_per_day + StaySwitch2 + prev_fledge + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[4]] <- glmer(StaySwitch ~ center_year + log_rain_per_day + StaySwitch2 + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))

# random effects structure - model selection
MuMIn::model.sel(model.list)

# fixed effects structure - model fitting
model.list <- list()
model.list[[1]] <- glmer(StaySwitch ~ 1 + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[2]] <- glmer(StaySwitch ~ log_rain_per_day + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[3]] <- glmer(StaySwitch ~ StaySwitch2 + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[4]] <- glmer(StaySwitch ~ prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[5]] <- glmer(StaySwitch ~ log_rain_per_day + StaySwitch2  + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[6]] <- glmer(StaySwitch ~ log_rain_per_day + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[7]] <- glmer(StaySwitch ~ StaySwitch2 + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[8]] <- glmer(StaySwitch ~ log_rain_per_day + StaySwitch2 + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))

model.list[[9]] <- glmer(StaySwitch ~ center_year  + (1|center_year) + (log_rain_per_day|penguinseq),
                         data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[10]] <- glmer(StaySwitch ~ center_year + (1|center_year) + log_rain_per_day + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[11]] <- glmer(StaySwitch ~ center_year + StaySwitch2 + (1|center_year) + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[12]] <- glmer(StaySwitch ~ center_year + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[13]] <- glmer(StaySwitch ~ center_year + log_rain_per_day + StaySwitch2  + (1|center_year) + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[14]] <- glmer(StaySwitch ~ center_year + log_rain_per_day + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[15]] <- glmer(StaySwitch ~ center_year + StaySwitch2 + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))
model.list[[16]] <- glmer(StaySwitch ~ center_year + log_rain_per_day + StaySwitch2 + prev_fledge + (1|center_year) + (log_rain_per_day|penguinseq),
                          data = mate.switch.data, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"))

# fixed effects structure - model selection
MuMIn::model.sel(model.list)


#############################
# 3. Foraging trip distance
#############################

# random effects structure - model fitting
model.list <- list()
model.list[[1]] <- lm(MaxTripDistance ~ zscore_climate + center_year + sex, data=trip.dist.data)
model.list[[2]] <- lmer(MaxTripDistance ~ zscore_climate + center_year + sex + (1|center_year), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[3]] <- lmer(MaxTripDistance ~ zscore_climate + center_year + sex + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[4]] <- lmer(MaxTripDistance ~ zscore_climate + center_year + sex + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))


# random effects structure - model selection
MuMIn::model.sel(model.list)

# fixed effects structure - model fitting 
model.list <- list()
model.list[[1]] <- lmer(MaxTripDistance ~ 1 + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[2]] <- lmer(MaxTripDistance ~ zscore_climate + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[3]] <- lmer(MaxTripDistance ~ center_year + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[4]] <- lmer(MaxTripDistance ~ sex + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[5]] <- lmer(MaxTripDistance ~ sex + center_year + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[6]] <- lmer(MaxTripDistance ~ zscore_climate + center_year + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[7]] <- lmer(MaxTripDistance ~ zscore_climate + sex + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[8]] <- lmer(MaxTripDistance ~ zscore_climate + sex + center_year + (1|center_year) + (zscore_climate|penguinseq), data=trip.dist.data,control = lmerControl(optimizer = "bobyqa"))

# fixed effects structure - model selection
MuMIn::model.sel(model.list)


#############################
# 4. Foraging site fidelity 
#############################

# random effects structure - model fitting
model.list <- list()
model.list[[1]] <- lm(DistBetweenPrev ~ zscore_climate + center_year + Sex, data=trip.fid.data)
model.list[[2]] <- lmer(DistBetweenPrev ~ zscore_climate + center_year + Sex + (1|center_year), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[3]] <- lmer(DistBetweenPrev ~ zscore_climate + center_year + Sex + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[4]] <- lmer(DistBetweenPrev ~ zscore_climate + center_year + Sex + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))


# random effects structure - model selection
MuMIn::model.sel(model.list)

# fixed effects structure - model fitting 
model.list <- list()
model.list[[1]] <- lmer(DistBetweenPrev ~ 1 + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[2]] <- lmer(DistBetweenPrev ~ zscore_climate + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[3]] <- lmer(DistBetweenPrev ~ center_year + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[4]] <- lmer(DistBetweenPrev ~ Sex + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[5]] <- lmer(DistBetweenPrev ~ Sex + center_year + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[6]] <- lmer(DistBetweenPrev ~ zscore_climate + center_year + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[7]] <- lmer(DistBetweenPrev ~ zscore_climate + Sex + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))
model.list[[8]] <- lmer(DistBetweenPrev ~ zscore_climate + Sex + center_year + (1|center_year) + (zscore_climate|PenguinSeq), data=trip.fid.data,control = lmerControl(optimizer = "bobyqa"))

# fixed effects structure - model selection
MuMIn::model.sel(model.list)
