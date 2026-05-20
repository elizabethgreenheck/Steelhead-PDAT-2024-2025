####Steelhead PDAT Mortality Project^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Last Updated 2026-05-20 by EMG

#This code accompanies the model used in the following publication:

#Quantifying component mortality estimates of out-migrating juvenile steelhead (Oncorhynchus mykiss) 
#Authors: E. M. Greenheck1, C. Michel2, B. Lehman2, L. Takata2, N. Demetras2, T. R. Nelson1
#1 George Mason University, Department of Environmental Science and Policy
#2 University of California Santa Cruz in affiliation with NOAA-NMFS-SWFSC
#Canadian Journal of Fisheries and Aquatic Sciences

#Kéry and Schaub 2012 (Bayesian population analysis using WinBUGS: a hierarchical perspective) was used to build these models

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Script information

#The data needed for this script is available on GitHub (four files: 'Obs_matrix_all_new_20250711.csv','temp_it.csv','turb_it.csv','discharge_it.csv')
#This script indexes a non-staggered observation matrix (all fish enter at time=1)
#The data is censored once a fish emigrates (fish become a "0" once they pass gate 17)
#This uses 8 days condensed and models values for each day (8t/8p; real days 1-8 are modeled + days 9-52 as one informative "final state")
#The data includes release (release (p1), state at end of each day (p2-p9), final state (p10))
#There are 3 observation states and 3 model states (3o, 3s)
#detection probability is modeled for each period (t)
#three candidate covariates are modeled to estimate component mortality -  but the covariates are not used to model the final state

#Rather than creating six different scripts for each candidate model, the candidate models are toggled below (Lines 203-216).
#Each jagsfile is saved for each candidate model to create FigureS4 (whisker plot)
#The user needs only to updated the monitored parameters, and the prediction calculations (Lines 254-259) are reserved only for the best fit model.
#The code below calculates the best fit models - the plotting code only plots significant relationships from the best fit model.

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ------------------------------------------------- # Parameters
# t: Period (Periods)
# i: individual (nFish)

# Pr: predation mortality probability
# M: unknown mortality probability
# Z: overall mortality probability
# S: survival probability
# p: individual level detection probability (p)

# a0: intercept for predation mortality
# a1: slope for predation mortality
# b0: intercept for unknown mortality
# b1: slope for unknown mortality
# a2/a3/b2/b3: slopes for additive and interactive model effects (as needed)
# error: residual error
# ------------------------------------------------- # Calculated values:
# Pr.length.pred: Predicted predation mortality using slope for length effects
# Pr.disch.pred: Predicted predation mortality using slope for discharge effects
# M.length.pred: Predicted unknown mortality using slope for length effects
# M.disch.pred: Predicted unknown mortality using slope for discharge effects
# ------------------------------------------------- # States (S):
# 1 = alive
# 2 = predated (PDAT triggered)
# 3 = unknown mortality
# ------------------------------------------------- # Observations (O):
# 1 = detected alive (PDAT not triggered) (A/E in capture history)
# 2 = detected predation (PDAT triggered) (P in capture history)
# 3 = unobserved or stationary tag (U/D in capture history)
# -------------------------------------------------

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Required packages & setting the working directory
##clear environment and set working directory
rm(list = ls())
setwd("yourwd")

###Install/Load needed packages
packs <- c('purrr','coda','R2OpenBUGS','rjags','gdata','devtools','ggh4x',
           'jagsUI', 'postpack','dplyr','ggplot2','tidyr','tidyverse','bayesplot','MCMCvis')
#install.packages(packs)
devtools::install_github("bstaton1/postpack", force = TRUE) ##great package by Ben Staton for MCMC objects
lapply(packs, library, character.only = TRUE)

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Reading in the capture history - available on Dryad for download
dat <- read_csv('Obs_matrix_all_new_20250711.csv')    

#Fish that died after the 8-day cutoff but eventually died
#they need to be considered "alive" during the first 8- days
dat$`8`[dat$FishID == 724] <- 1
dat$`8`[dat$FishID == 1293] <- 1
dat$`8`[dat$FishID == 1560] <- 1 #this fish became a 4 (observed dead) after the 8-days

##Defining the number of days to be included in the model, this does not include release or final state
ndays = 8

#Identifying column names to retain
periods = c(0:ndays,"FinalState") #day 0 (release) and the additional columns defined by ndays + the final state
ivars = c("FishID","Release_Group_binned","Length","Weight") #These are the column names of our individual covariates

#Pulling the columns that we define in "periods" from our CH matrix
CH <- as.matrix(dat[,colnames(dat) %in% periods],header=F)

#Censoring out the fish when they've emigrated from the system
CH <- t(apply(CH, 1, function(row) {
  max_val <- max(row[row < 5], na.rm = TRUE)
  if (max_val %in% c(3)) {
    first_max <- which(row == max_val)[1]
    row[row == max_val] <- 0
    row[first_max] <- max_val
  }
  return(row)
}))

#Recoding states for a 3-state model:
CH[CH == 3] <- 1 #recoding emigration as survival
CH[CH == 4] <- 3 #recoding observed natural mortality (stationary tags) as natural mortality/unobserved
CH[CH == 5] <- 3 #recoding unobserved state 5 as unobserved state 3

#We need the environmental data to be in the same order so it's indexed properly
##Your scaling is incorrect here, likely why you are getting weird results
ivars <- as.matrix(dat[,colnames(dat) %in% ivars])
temp_it <- read_csv("temp_it.csv")
temp_it <- fill(temp_it,everything())
temp_it <- as.data.frame(ivars) %>% left_join(temp_it)
temp_it_ <- temp_it

#For this model, we're doing daily mortality over a 10-day period; 
#For the last "day", it's going to be the average of days 9-ncol(temp_it)
temp_it_$Temp_F <- rowMeans(temp_it_[,c(13:54)], na.rm = TRUE)
temp_it_ <- temp_it_[,c(0:8,"Temp_F")] 

temp_it_ <- as.matrix(temp_it_)

#Turbidity data
turb_it <- read_csv("turb_it.csv")
turb_it <- fill(turb_it,everything())
turb_it <- as.data.frame(ivars) %>% left_join(turb_it)
turb_it_ <- turb_it

#For this model, we're doing daily mortality over a 10-day period; 
#For the last "day", it's going to be the average of days 9-ncol(temp_it)
turb_it$Turb_NTU <- rowMeans(turb_it_[,c(13:55)], na.rm = TRUE)
turb_it_ <- turb_it[,c(0:8,"Turb_NTU")] 

turb_it_ <- as.matrix(turb_it_)

#Discharge data
disch_it <- read_csv("discharge_it.csv")
disch_it <- fill(disch_it,everything())
disch_it <- as.data.frame(ivars) %>% left_join(disch_it)
disch_it_ <- disch_it

#For this model, we're doing daily mortality over a 10-day period; 
#For the last "day", it's going to be the average of days 9-ncol(temp_it)
disch_it_$Disch_m3s <- rowMeans(disch_it_[,c(13:55)], na.rm = TRUE)
disch_it_ <- disch_it_[,c(0:8,"Disch_m3s")] 

disch_it_ <- as.matrix(disch_it_)

#Removing values when a fish was censored from temp calculations by multiplying the CH with 1s and 0s to the covariate matrix
#This will remove days when a fish isn't in the study area and, subsequently, their environmemtal coviariates
CH_0s <- CH
CH_0s[CH_0s > 0] <- 1

#For temperature
temp_it_ <- temp_it_*CH_0s
temp_it_[temp_it_ == 0] <- NA #filling in NAs in place of 0 so the 0s aren't calculated as part of rowmeans

#For discharge
disch_it_ <- disch_it_*CH_0s
disch_it_[disch_it_ == 0] <- NA #filling in NAs in place of 0 so the 0s aren't calculated as part of rowmeans

#For turbidity
turb_it_ <- turb_it_*CH_0s
turb_it_[turb_it_ == 0] <- NA #filling in NAs in place of 0 so the 0s aren't calculated as part of rowmeans

#Taking the mean values of discharge & temp and pulling the length column
disch_i_ <- rowMeans(disch_it_, na.rm = TRUE)
temp_i_ <- rowMeans(temp_it_, na.rm = TRUE)
turb_i_ <- rowMeans(turb_it_, na.rm = TRUE)
length_i_ <- ivars[,3]

#Scaling the variables
scale_i <- function(raw){
  scale_i <- scale(raw)
  scale_i_mod <- scale_i[1:length(raw)]
  cntr <- attr(scale_i,'scaled:center')
  scl <- attr(scale_i,'scaled:scale')
  predict = seq(min(raw),max(raw),length = length(raw))
  s_values <- data.frame(raw.var = raw, 
                         scl.var = scale_i_mod, 
                         pred.var = (predict[1:length(predict)]-cntr)/scl,
                         scale_param = scl,
                         cntr_param = cntr)
  return(s_values)
}

#Applying the scale function to each variable
temppredict_c <- scale_i(temp_i_)
dischpredict_c <- scale_i(disch_i_)
lengthpredict_c <- scale_i(length_i_)
turbpredict_c <- scale_i(turb_i_)

#removing things out of the environment
rm(disch_it_,disch_it,ivars,temp_it,temp_it_,turb_it_,disch_i_,length_i_,turb_i_,ndays,packs,periods,temp_i_,scale_i,CH_0s)

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##model code for JAGS####
mod = function() {
  # Calculating parameters to be used in state-transition and observation matrices
  # Each candidate model is toggled on/off; this code just selects the best model (Model 5)
  for (i in 1:nFish){  
    for (t in 1:(Periods-2)){
      #log(Pr[i,t]) <- a0 + a1*length[i] #Model 1
      #log(M[i,t]) <- b0 + b1*length[i] #Model 1
      #log(Pr[i,t]) <- a0 + a1*temp[i] #Model 2
      #log(M[i,t]) <- b0 + b1*temp[i] #Model 2
      #log(Pr[i,t]) <- a0 + a1*disch[i] #Model 3
      #log(M[i,t]) <- b0 + b1*disch[i] #Model 3
      #log(Pr[i,t]) <- a0 + a1*turb[i] #Model 4
      #log(M[i,t]) <- b0 + b1*turb[i] #Model 4
      log(Pr[i,t]) <- a0 + a1*length[i] + a2*disch[i]  #Model 5
      log(M[i,t]) <- b0 + b1*length[i] #Model 5 
      #log(Pr[i,t]) <- a0 + a1*length[i] + a2*disch[i] + a3*length[i]*disch[i] #Model 6
      #log(M[i,t]) <- b0 + b1*length[i] #Model 5 
      
      Z[i,t] <- Pr[i,t] + M[i,t]
      S[i,t] <- exp(-Z[i,t])
    } #t
  }# i
  
  ##priors/constraints for last period that encompasses day 9 until last detection
  ##We cannot measure covariate effects here, but can use this to inform true states
  ##for unobserved fish during days 1 - 8
  for (i in 1:nFish){  
    M[i,Periods-1] ~ dunif(0,1)
    Pr[i,Periods-1] ~ dunif(0,1) 
    Z[i,Periods-1] <- Pr[i,Periods-1] + M[i,Periods-1]
    S[i,Periods-1] <- exp(-Z[i,Periods-1])
  }# i
  
  #Priors for the linear models; all are uninformative and all can be listed even if not included
  a0 ~ dunif(-20,2)
  a1 ~ dnorm(0,0.01) 
  a2 ~ dnorm(0,0.01)
  a3 ~ dnorm(0,0.01)
  
  b0 ~ dunif(-20,2)
  b1 ~ dnorm(0,0.01)
  b2 ~ dnorm(0,0.01) 
  b3 ~ dnorm(0,0.01) 
  
  #Generating predictive estimates using our model intercepts (a0/b0) and slopes (a1/b1)
  #For model 5 - Pr was predicted by length & discharge so we predict length holding discharge at its mean (mean = 0 because it's scaled)
  #And then we predict discharge at the mean length (again, mean = 0)
  #For M, only length was predictive
  
  #When testing what model is best, omit this portion (only calculate for best model)
  for(j in 1:length(lengthpredict)){
    log(Pr.length.pred[j]) <- a0 + a1*lengthpredict[j] + a2*meandisch
    log(Pr.disch.pred[j]) <- a0 + a1*meanlength + a2*dischpredict[j]
    log(M.length.pred[j]) <- b0 + b1*lengthpredict[j]
  } #j
  
  #detection probability 
  p[1] <- 1 # Model conditioned on first capture, estimate separate p for remaining periods
  for (t in 2:(Periods)){
    p[t] ~ dunif(0, 1)
  }
  
  # Define state-transition and observation matrices
  for (i in 1:nFish){
    # Define probabilities of State (t+1) given State (t). First index is state at time t, next is
    #state at t+1
    for (t in first[i]:(last[i]-1)){
      ps[1,i,t,1] <- S[i,t]                     ##alive fish stays alive in system
      ps[1,i,t,2] <- Pr[i,t]*(1-S[i,t])/Z[i,t]  ##alive fish is predated (PDAT trigger)
      ps[1,i,t,3] <- M[i,t] *(1-S[i,t])/Z[i,t]  ##alive fish dies of unknown causes
      ps[2,i,t,1] <- 0                          ##predated fish becomes alive (cannot happen)
      ps[2,i,t,2] <- 1                          ##predated fish remains predated 
      ps[2,i,t,3] <- 0                          ##predated fish dies of unknown causes (cannot happen)
      ps[3,i,t,1] <- 0                          ##fish dead from unknown causes becomes alive (cannot happen)
      ps[3,i,t,2] <- 0                          ##fish dead from unknown causes becomes predated (cannot happen)
      ps[3,i,t,3] <- 1                          ##fish dead from unknown causes remains dead from unknown causes
    } #t
    for (t in first[i]:(last[i])){
      # Define probabilities of Observed (t) given State (t). First index is state, last index is
      #observed
      po[1,i,t,1] <- p[t]   ##alive fish is detected alive in system
      po[1,i,t,2] <- 0        ##alive fish is detected as predated (cannot happen)
      po[1,i,t,3] <- 1-p[t] ##alive fish is not detected
      po[2,i,t,1] <- 0        ##predated fish is detected alive (cannot happen)
      po[2,i,t,2] <- 1        ##predated fish is detected as predated; terminal state (once a 2, stays a 2)
      po[2,i,t,3] <- 0        ##predated fish is detected as an unknown mortality (cannot happen)
      po[3,i,t,1] <- 0        ##unobserved fish is detected alive (cannot happen) 
      po[3,i,t,2] <- 0        ##unobserved fish is detected as predated (cannot happen) 
      po[3,i,t,3] <- 1        ##unobserved fish remains unobserved
    } #t
  } #i
  
  # Likelihood
  for (i in 1:nFish){
    z[i,first[i]] <- 1 # Individuals have known status (alive) at first occasion in study
    for (t in (first[i]+1):last[i]){
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,]) # State process: draw State (t) given State (t-1)
    } #t
    for (t in first[i]:last[i]){
      y[i,t] ~ dcat(po[z[i,t], i, t,]) # Observation process: draw Observed (t) given State (t)
    } #t
  } #i
} #model close

# write model to a text file
model.file = "model.txt"
write.model(mod, model.file)

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Get data for model to run 
y=CH                ##Capture History from Above
nFish=dim(CH)[1]    ##Number of individuals obtained from row count
Periods=dim(CH)[2]  ##Number of time-frames obtained from column count
first <- numeric()        ##first observation of each fish
for (i in 1:dim(y)[1]) {
  first[i] <-min(which(y[i,]!=0))}
last <-numeric()          ##last observation of each fish
for (i in 1:dim(y)[1]) {
  last[i] <-max(which(y[i,]!=0))}
f <- first
l <- last

#pulling the scaled variable to model
temp = as.numeric(temppredict_c$scl.var)
length = as.numeric(lengthpredict_c$scl.var)
disch = as.numeric(dischpredict_c$scl.var)
turb = as.numeric(turbpredict_c$scl.var)

#calculating the mean for each variable (for model predictions)
meanlength = mean(length) #should be 0
meantemp= mean(temp) #0
meandisch = mean(disch)  #0
meanturb = mean(turb)  #0

#grabbing predictive variables (from scaling function)
lengthpredict = as.numeric(lengthpredict_c$pred.var)
temppredict = as.numeric(temppredict_c$pred.var)
dischpredict = as.numeric(dischpredict_c$pred.var)
turbpredict = as.numeric(turbpredict_c$pred.var)

jags.data<-list(y=y,nFish=nFish,first=first,last=last,Periods=Periods,
                length=length,temp = temp,disch=disch,turb=turb,
                lengthpredict=lengthpredict,temppredict=temppredict,dischpredict=dischpredict,turbpredict=turbpredict,
                meanlength=meanlength,meantemp=meantemp,meandisch=meandisch,meanturb=meanturb)

#initial values for 3S model
ms.init.z <- function(y, f) {
  y.iv <- y
  y.iv[y.iv==3] <- NA #changing this to our unobserved state
  
  for(i in 1:nrow(y.iv)){
    if(max(y.iv[i,],na.rm=TRUE)==1){
      y.iv[i,(f[i]+1):l[i]] <- 1}# Not detected dead or dispersed so initialize as alive
    
    if(max(y.iv[i,],na.rm=TRUE)==2){
      m <- min(which(y.iv[i,]==2))
      y.iv[i,f[i]:(m-1)] <- 1 # Not detected dead or dispersed so initialize as alive
      y.iv[i,m:l[i]] <- 2} # Not detected dead or dispersed so initialize as alive
    
  }
  y.iv[y.iv == 0] <- NA #Then, filling 0s as NAs
  for (i in 1:dim(y.iv)[1]){y.iv[i,1:f[i]] <- NA}
  return(y.iv)
} #Final initial function contains 1s and 2s


jags.inits <- function(){list(z=ms.init.z(y,f))}

##parameters to monitor during model run 
params <-c('a0','b0','a1','b1','a2','b2','a3','b3',
           "Pr",
           'Pr.length.pred',
           'Pr.disch.pred',
           'M.length.pred')

##run model in JAGS

#Plot using autojags
set.seed(121)
jagsfit <- autojags(data=jags.data, inits=jags.inits, parameters.to.save=params, model.file,
                    n.chains=3, n.adapt=NULL, iter.increment=3000, n.burnin=100000, n.thin=1,
                    save.all.iter=FALSE, modules=c('glm'), parallel=TRUE,
                    DIC=TRUE, store.data=FALSE, codaOnly=FALSE,
                    bugs.format=FALSE, Rhat.limit=1.01, max.iter=500000, verbose=TRUE)

jagsfit

#Each model iteration is saved as a .rds so they can be pulled for analyses later 
#saveRDS(jagsfit,"8t_8p_3o_3s_length_log_noerror.rds") #Model 1
#saveRDS(jagsfit,"8t_8p_3o_3s_temp_log_noerror.rds") #Model 2
#saveRDS(jagsfit,"8t_8p_3o_3s_disch_log_noerror.rds") #Model 3
#saveRDS(jagsfit,"8t_8p_3o_3s_turb_log_noerror.rds") #Model 4
#saveRDS(jagsfit,"8t_8p_3o_3s_Prlength,disch_Mlength_log_noerror.rds") #Model 5
#saveRDS(jagsfit,"8t_8p_3o_3s_Prlength:disch_Mlength_log_noerror.rds") #Model 6
saveRDS(jagsfit,"8t_8p_3o_3s_Prlength,disch_Mlength_log_predicted_noerror.rds") #Model 5 - with predictions

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Figure plotting - if starting here you'll need to read in the environmental data too
jagsfit <- readRDS("8t_8p_3o_3s_Prlength,disch_Mlength_log_predicted_noerror.rds")
jagsfit <- `8t_8p_3o_3s_Prlength,disch_Mlength_log_predicted_noerror`

#Extracting model predictions
predict_extract <- function(jagsfit){
  param_pattern <- "^(Pr|M)"
  params <- grep(param_pattern, names(jagsfit$mean), value = TRUE)
  
  param_list <- lapply(params, function(p) {
    data.frame(
      mean = jagsfit$mean[[p]],
      sd = jagsfit$sd[[p]],
      q2.5 = jagsfit$q2.5[[p]],
      q25 = jagsfit$q25[[p]],
      q50 = jagsfit$q50[[p]],
      q75 = jagsfit$q75[[p]],
      q97.5 = jagsfit$q97.5[[p]],
      param = p
    )
  })
  do.call(rbind, param_list)
}

predict <- predict_extract(jagsfit) 
predict <- predict %>%
  filter(!grepl(".D",param))

#Joining the model predictions to the actual predicted values from our initial scaling
param=rep(unique(predict$param),each=300)
predictvalues <- rbind(lengthpredict_c,dischpredict_c,lengthpredict_c)
predictvalues <- cbind(param,predictvalues)
predict <- cbind(predict,predictvalues)

predict <- predict[,c(1:9,10:14)]

#Some data tidying to get things in the proper format for plotting
predict_ <- predict %>%
  mutate(Var = case_when(grepl("length",param) ~ "'Fork Length (mm)'",
                         grepl("temp",param) ~ "Temperature (°C)",
                         grepl("disch",param) ~ "Discharge~(m^3~s^{-1})",
  )) %>%
  mutate(param = case_when(grepl("Pr",param) ~ "Pr",
                           grepl("M",param) ~ "M"))

length_ <- predict_ %>%
  filter(Var == "'Fork Length (mm)'")
length_breaks = c(-3:2)
length_labels = round((length_breaks*unique(length_$scale_param)) + unique(length_$cntr_param))
temp_ <- predict_ %>%
  filter(Var == "Temperature (°C)")
temp_breaks = c(-2:2)
temp_labels = round((temp_breaks*unique(temp_$scale_param)) + unique(temp_$cntr_param))
disch_ <- predict_ %>%
  filter(Var == "Discharge~(m^3~s^{-1})")
disch_breaks = c(-2:2)
disch_labels = round((disch_breaks*unique(disch_$scale_param)) + unique(disch_$cntr_param))

# Creating Figure 4
jagsplot_covariates <- predict_ %>%
  filter(param %in% c("M","Pr")) %>%
  mutate(Var = factor(Var, levels = c("'Fork Length (mm)'", "Discharge~(m^3~s^{-1})")),
         param = factor(param,levels=c("Pr","M"))) %>%
  ggplot(aes(x=pred.var))+
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5, fill = param),
            alpha = 0.2,
              linetype = 0)+
  geom_line(aes(y=q50,color=param))+
  guides(fill =  guide_legend(position = "inside"),
         col =    guide_legend(position = "inside"))+
  scale_color_manual(name="Parameter",values=c("#785EF0","#FD005D"),
                     labels=c(bquote('Fish predation (' * italic(Pr)^I * ')'),
                              bquote('Unknown mortality (' * italic(M)^I * ')')))+
  scale_fill_manual(name="Parameter",values=c("#785EF0","#FD005D"),
                    labels=c(bquote('Fish predation (' * italic(Pr)^I * ')'),
                             bquote('Unknown mortality (' * italic(M)^I * ')')))+
  xlab(element_blank())+
  ylab("Mortality Estimates")+
  facet_wrap(~Var,scales="free",strip.position = "bottom",labeller = label_parsed)+
  facetted_pos_scales(
    x = list(
      Var == "Discharge~(m^3~s^{-1})" ~ scale_x_continuous(breaks = disch_breaks, labels = disch_labels,expand=c(0.01,0.01)),
      Var == "'Fork Length (mm)'" ~ scale_x_continuous(breaks = length_breaks, labels = length_labels,expand=c(0.01,0.01))
    )
  )+
  scale_y_continuous(expand=c(0,0))+
  theme_linedraw()+
  theme(text=element_text(size=18,color="black"),legend.position.inside = c(0.85,0.9),
        strip.background = element_rect(color="white",fill="white"),strip.text = element_text(color="black"),
        strip.placement = "outside",panel.grid.minor = element_blank())
jagsplot_covariates

# Creating Table 5 & Figure S4
#These require that all candidate models are read in (the saved .rds files)
setwd("yourwd")
file_paths <- list.files(pattern = "\\.rds$", full.names = TRUE)
file_paths <- file_paths[grepl("log", file_paths)]
file_paths <- file_paths[!grepl("pred", file_paths)]

#I just want my best model slopes + temp on its own + the interaction on its own
# Extract file names without extensions
file_names <- tools::file_path_sans_ext(basename(file_paths))

# Read files into a list
data_list <- lapply(file_paths, readRDS)

# Assign file names to list elements
names(data_list) <- file_names

#Function to extract slopes
slopes_extract <- function(jagsfit) {
  # Get all parameters starting with 'a' or 'b' followed by a number
  param_pattern <- "^[ab][0-9]+$"
  params <- grep(param_pattern, names(jagsfit$mean), value = TRUE)
  
  # Create a list of data frames for each valid parameter
  param_list <- lapply(params, function(p) {
    data.frame(
      mean = jagsfit$mean[[p]],
      sd = jagsfit$sd[[p]],
      q2.5 = jagsfit$q2.5[[p]],
      q25 = jagsfit$q25[[p]],
      q50 = jagsfit$q50[[p]],
      q75 = jagsfit$q75[[p]],
      q97.5 = jagsfit$q97.5[[p]],
      overlap0 = jagsfit$overlap0[[p]],
      param = p
    )
  })
  
  # Combine all data frames and return
  do.call(rbind, param_list)
}

whiskers <- lapply(data_list,slopes_extract)
whiskers <- do.call("rbind",whiskers)
whiskers <- rownames_to_column(whiskers)

#getting rid of this long naming string 
patterns <- c("8t_8p_3o_3s","_log_noerror")
pattern_regex <- paste(patterns, collapse = "|")

# Remove all patterns in one mutate
whiskers <- whiskers %>%
  mutate(
    rowname_clean = str_replace_all(rowname, pattern_regex, ""),
  ) 

#I accidentally monitored "a2" & "a3" for non-additive single effects models- removing those
whiskers <- whiskers %>%
  mutate(remove = case_when(grepl("_disch.",rowname_clean) & grepl("2|3",param) ~ "y",
                            grepl("_temp.",rowname_clean) & grepl("2|3",param) ~ "y",
                            grepl("_length.",rowname_clean) & grepl("2|3",param) ~ "y",
                            grepl("_turb.",rowname_clean) & grepl("2|3",param) ~ "y",
                            grepl(",",rowname_clean) & grepl("3",param) ~ "y",
                            TRUE ~ "n")) %>%
  filter(remove == "n") %>%
  select(!remove)
                          
#unscale the link-scale slopes (first)
#Then transform them to the log-scale. 
#It's a little convoluted
Table5 <- whiskers %>%
  dplyr::select(!c(rowname,mean,sd,q25,q75)) %>%
  mutate(col= case_when(grepl("0",param) ~ "b0",
                        grepl("1",param) ~ "b1",
                        grepl("2",param) ~ "b2",
                        grepl("3",param) ~ "b3"),
         Model = case_when(grepl("_disch.",rowname_clean) ~ "D",
                           grepl("_length,disch.",rowname_clean) ~ "FL + D",
                           grepl("_length:disch.",rowname_clean) ~ "FL * D",
                           grepl("_length.",rowname_clean) ~ "FL",
                           grepl("_temp.",rowname_clean) ~ "T",
                           grepl("_turb.",rowname_clean) ~ "N",
                           grepl("Prlength,disch",rowname_clean) & grepl("a",param) ~ "Pr FL + D",
                           grepl("Prlength,disch",rowname_clean) & grepl("b",param) ~ "M FL",
                           grepl("Prlength:disch",rowname_clean) & grepl("a",param) ~ "Pr FL * D",
                           grepl("Prlength:disch",rowname_clean) & grepl("b",param) ~ "M FL",
         )) %>%
  mutate(Model = case_when(grepl("P|M",Model) ~ Model,
                           grepl("a",param) ~ paste0("Pr ",Model),
                           grepl("b",param) ~ paste0("M ",Model)),
         scale_param = case_when(grepl("length,disch.",rowname_clean) & grepl("1",col) ~ unique(lengthpredict_c$scale_param),
                                 grepl("length,disch.",rowname_clean) & grepl("2",col) ~ unique(dischpredict_c$scale_param),
                                 grepl("length:disch.",rowname_clean) & grepl("1",col) ~ unique(lengthpredict_c$scale_param),
                                 grepl("length:disch.",rowname_clean) & grepl("2",col) ~ unique(dischpredict_c$scale_param),
                                 grepl("length:disch.",rowname_clean) & grepl("3",col) ~ unique(lengthpredict_c$scale_param)*unique(dischpredict_c$scale_param),
                                 grepl("T",Model) & grepl("1",col) ~ unique(temppredict_c$scale_param),
                                 grepl("N",Model) & grepl("1",col) ~ unique(turbpredict_c$scale_param),
                                 grepl("FL",Model) & grepl("1",col) ~ unique(lengthpredict_c$scale_param),
                                 grepl("D",Model) & grepl("1",col) ~ unique(dischpredict_c$scale_param),
         )) %>%
  arrange(Model) %>%
  mutate(q2.5_unscaled = round(q2.5/scale_param),
         q50_unscaled = round(q50/scale_param),
         q97.5_unscaled = round(q97.5/scale_param)) %>%
  mutate(q2.5_raw = round(exp(q2.5_unscaled)),
         q50_raw = round(exp(q50_unscaled)),
         q97.5_raw = round(exp(q97.5_unscaled))) %>%
  mutate(link_scaled = paste0( sprintf("%.3f", q50), " [",sprintf("%.3f", q2.5), " - ",sprintf("%.3f", q97.5), "]"),
         link = paste0( sprintf("%.3f", q50_unscaled), " [",sprintf("%.3f", q2.5_unscaled), " - ",sprintf("%.3f", q97.5_unscaled), "]"),
         response = paste0( sprintf("%.3f", q50_raw), " [",sprintf("%.3f", q2.5_raw), " - ",sprintf("%.3f", q97.5_raw), "]")) %>%
  dplyr::select(!c(q2.5,q2.5_unscaled,q2.5_raw,q50,q50_unscaled,q50_raw,q97.5,q97.5_unscaled,q97.5_raw))
Table5$response[Table5$response == "NA [NA - NA]"] <- NA
Table5$link[Table5$link == "NA [NA - NA]"] <- NA

#Cleaning up the rownames
Table5$rowname_clean <- str_replace_all(Table5$rowname_clean, "\\..*", "")
Table5_ <- Table5 %>%
  mutate(overlap0=as.character(overlap0)) %>%
  pivot_longer(
    cols = c(link, response,overlap0),
    names_to = "type",
    values_to = "value"
  ) %>%
  dplyr::select(rowname_clean, Model, type, col, value) %>%
  pivot_wider(
    names_from = col,
    values_from = value
  ) %>%
  select(!b0) %>%
  mutate(BestMod = case_when(grepl("M",Model) & b1 == FALSE & is.na(b2) ~ 1,
                             grepl("Pr",Model) & b1 == FALSE & b2 == FALSE & is.na(b3) ~ 1)) %>%
  group_by(Model, rowname_clean) %>%
  fill(BestMod,.direction="downup") %>%
  mutate(BestMod = cumsum(BestMod))

Table5_$BestMod[is.na(Table5_$BestMod)] <- 0

Table5_ <- Table5_ %>%
  group_by(BestMod,Model,rowname_clean) %>%
  slice_head() %>%
  filter(type != "overlap0") %>%
  mutate(remove = case_when(rowname_clean == "_Prlength,disch_Mlength" & Model == "M FL" ~ 1,
                            rowname_clean == "_Prlength:disch_Mlength" & Model == "M FL" ~ 1,
                            TRUE ~ 0
                            )) %>%
  filter(remove == 0) %>%
  ungroup() %>%
  select(!c(remove,BestMod,rowname_clean)) %>%
  mutate(Model = factor(Model,levels=c("Pr T","Pr N","Pr D","Pr FL","Pr FL + D","Pr FL * D",
                                       "M T","M N","M D","M FL"))) %>%
  arrange(Model) 
  
Table5_ %>% write_csv("Table5.csv")
#Table 5 doesn't print exactly as it's written, some decisions had to be made with what to report

#Creating the whisker plot of all models (Figure S4)
whiskers <- whiskers %>%
  filter(param != "b0", param != "a0") %>%
  mutate(keep = case_when(grepl(",",rowname_clean) ~ "y",
                          grepl("length:disch",rowname_clean) & grepl("a3",param) ~ "y",
                          grepl("temp",rowname_clean) ~ "y",
                          grepl("turb",rowname_clean) ~ "y",
                          TRUE ~ "n"))

whiskers_ <- whiskers %>%
  filter(keep == "y") %>%
  mutate(Var = case_when(grepl(":",rowname_clean)~ "ForkLength:Discharge",
                         grepl("length",rowname_clean) & grepl("1",param) ~ "Fork Length (mm)",
                         grepl("temp",rowname_clean)~ "Temperature (°C)",
                         grepl("turb",rowname_clean)~ "Turbidity (NTU)",
                         grepl("disch",rowname_clean) & grepl("2",param) ~ "Discharge (m3s)"),
         Parameter = case_when(grepl("a",param) ~ "Pr",
                               grepl("b",param) ~ "M")) %>%
  dplyr::select(!c(rowname,rowname_clean)) %>%
  mutate(order = case_when(Var == "Fork Length (mm)" ~ 3,
                           Var == "Temperature (°C)" ~ 1,
                           Var == "Discharge (m3s)" ~ 4,
                           Var == "Turbidity (NTU)" ~ 2,
                           Var == "ForkLength:Discharge" ~ 5)) %>%
  mutate(order = case_when(grepl("a",param) ~ order,
                           TRUE ~ order+5))%>%
  mutate( 
    xlab = case_when(
      Var == "Fork Length (mm)" & Parameter == "M"  ~ 'Fork*Length~"(mm)"~"-"~italic(M)^I',
      Var == "Fork Length (mm)" & Parameter == "Pr" ~ 'Fork*Length~"(mm)"~"-"~italic(Pr)^I',
      Var == "ForkLength:Discharge" & Parameter == "M"  ~ 'ForkLength*":"*Discharge~"-"~italic(M)^I',
      Var == "ForkLength:Discharge" & Parameter == "Pr" ~ 'ForkLength*":"*Discharge~"-"~italic(Pr)^I',
      Var == "Temperature (°C)" & Parameter == "M"  ~ 'Temperature~"(°C)"~"-"~italic(M)^I',
      Var == "Temperature (°C)" & Parameter == "Pr" ~ 'Temperature~"(°C)"~"-"~italic(Pr)^I',
      Var == "Turbidity (NTU)" & Parameter == "M"  ~ 'Turbidity~"(NTU)"~"-"~italic(M)^I',
      Var == "Turbidity (NTU)" & Parameter == "Pr" ~ 'Turbidity~"(NTU)"~"-"~italic(Pr)^I',
      Var == "Discharge (m3s)" & Parameter == "M"  ~ 'Discharge~(m^3*s^{-1})~"-"~italic(M)^I',
      Var == "Discharge (m3s)" & Parameter == "Pr" ~ 'Discharge~(m^3*s^{-1})~"-"~italic(Pr)^I',
      TRUE ~ NA_character_
    )
  )%>%
  mutate(xlab = factor(xlab, levels = unique(xlab[order(-order)]))) %>%
  filter(!grepl("0",param)) %>%
  mutate(facet = case_when(order == 9 ~ "1",
                           TRUE ~ "0"))

whiskerplot1 <- whiskers_ %>%
  group_by(order) %>%
  slice_head() %>%
  filter(order != 9) %>%
  ggplot(aes(x=xlab))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_errorbar(aes(ymin=q2.5,ymax=q97.5),color="black",width=0,linewidth=0.5)+
  geom_errorbar(aes(ymin=q25,ymax=q75),color="#333333",width=0,linewidth=2)+
  geom_point(aes(y=q50,fill=overlap0),color="black",pch=21,size=3)+
  scale_fill_manual(name="Significant",values=c("black","white"),
                    labels=c("Yes","No"))+
  coord_flip()+
  xlab(element_blank())+
  ylab(element_blank())+
  theme_linedraw()+
  theme(text=element_text(size=18),legend.position = "none",
        plot.margin = unit(c(0.5, 1.5, 0.5, 0.5), "lines"))+
  scale_x_discrete(labels = function(x) parse(text = x))+
  scale_y_continuous(limits=c(-1.5,1.5))
whiskerplot1

whiskerplot2 <- whiskers_ %>%
  group_by(order) %>%
  slice_head() %>%
  filter(order == 9) %>%
  ggplot(aes(x=xlab))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_errorbar(aes(ymin=q2.5,ymax=q97.5),color="black",width=0,linewidth=0.5)+
  geom_errorbar(aes(ymin=q25,ymax=q75),color="#333333",width=0,linewidth=2)+
  geom_point(aes(y=q50),fill="white",color="black",pch=21,size=3)+
  scale_fill_manual(name="Significant",values=c("black","white"),
                    labels=c("Yes","No"))+
  coord_flip()+
  xlab(element_blank())+
  ylab("Estimates (95% CrI)")+
  theme_linedraw()+
  theme(text=element_text(size=18),legend.position = "none",
        plot.margin = unit(c(0.5, 1.5, 0.5, 0.5), "lines"))+
  scale_x_discrete(labels = function(x) parse(text = x))+
  scale_y_continuous(limits=c(-20,20),expand=c(0,0))
whiskerplot2

library(cowplot)
FigureS4 <- plot_grid(
  whiskerplot1, whiskerplot2,
  ncol = 1,
  align = "v",
  rel_heights = c(4, 1)
)
FigureS4 
ggsave(
  "FigureS4.png",
  FigureS4,
  width  = 12,
  height = 9,
  units  = "in",
  dpi    = 300
)

