####Steelhead PDAT Mortality Project^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Last Updated 2025-09-17 by EMG

#This code accompanies the model used in the following publication:

#Quantifying component mortality estimates of out-migrating juvenile steelhead (Oncorhynchus mykiss) 
#Authors: E. M. Greenheck1, C. Michel2, B. Lehman2, L. Takata2, N. Demetras2, T. R. Nelson1
#1 George Mason University, Department of Environmental Science and Policy
#2 University of California Santa Cruz in affiliation with NOAA-NMFS-SWFSC
#[in progress]

#Kéry and Schaub 2012 (Bayesian population analysis using WinBUGS: a hierarchical perspective) was used to build these models

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#Script information

#The data needed for this script is available on Dryad (one file: 'Obs_matrix_all_new_20250711.csv')
#MSM_Non_C_8t_8p_3o_3s_pg_FINAL
#This script indexes a non-staggered observation matrix (all fish enter at time=1; Non)
#The data is censored once a fish emigrates (fish become a "0" once they pass gate 17; C)
#This uses 8 days condensed and models values for each day (8t/8p; real days 1-8 are modeled + days 9-52 as one informative "final state")
#The data includes release (release (p1), state at end of each day (p2-p9), final state (p10))
#There are 3 observation states and 3 model states (3o, 3s)
#The models estimate parameters for each group (u) and detection probability is modeled per individual and group for each time (pg)

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# ------------------------------------------------- # Parameters
# t: Period (Periods)
# i: individual (nFish)
# u: group (g)

# Pr: predation mortality probability
# M: unknown mortality probability
# Z: overall mortality probability
# S: survival probability
# p: individual level detection probability (p)

# ------------------------------------------------- # Calculated values:
# Pr_timeD: discrete group-level 8-day predation mortality estimate
# Pr_D: discrete group-level overall predation mortality estimate
# Pr_alltimeD: mean discrete 8-day predation mortality estimates
# Pr_allD: mean discrete overall predation mortality estimates
# M_timeD: discrete group-level 8-day unknown mortality estimate
# M_D: discrete group-level overall unknown mortality estimate
# M_alltimeD: mean discrete 8-day unknown mortality estimates
# M_allD: mean discrete overall predation unknown estimates
# S_timeD: discrete group-level 8-day survival estimate
# S_D: discrete group-level overall survival estimate
# S_alltimeD: mean discrete 8-day survival estimates
# S_allD: mean discrete overall survival estimates
# p.all: calculated overall detection probability for all groups

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
packs <- c('coda','R2OpenBUGS','rjags','gdata','devtools','jagsUI', 'dplyr','ggplot2','tidyr')
install.packages(packs)
devtools::install_github("bstaton1/postpack", force = TRUE) ##great package by Ben Staton for MCMC objects
packs <- c('purrr','coda','R2OpenBUGS','rjags','gdata','devtools','ggh4x',
           'jagsUI', 'postpack','dplyr','ggplot2','tidyr','tidyverse','bayesplot','MCMCvis')
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

#Pulling the group information
group=dat$Release_Group_binned

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##model code for JAGS####
mod = function() {
  # Calculating parameters to be used in state-transition and observation matrices
  # Group-level estimates for each individual at Period t
  for (i in 1:nFish){  
    for (t in 1:(Periods-1)){
      Pr[i,t] <- eta.phi.Pr[group[i],t] 
      M[i,t] <- eta.phi.M[group[i],t]
      Z[i,t] <- Pr[i,t]+M[i,t]
      S[i,t] <- exp(-Z[i,t])
    }#t
  }#i
  
  # Group-level estimates for each individual at Period t
  for(u in 1:g){
    for (t in 1:(Periods-1)){
      eta.phi.Pr[u,t] <- exp(ln.eta.phi.Pr[u,t])
      ln.eta.phi.Pr[u,t] ~ dunif(-10,0) #uninformative prior for Pr
      eta.phi.M[u,t] <- exp(ln.eta.phi.M[u,t])
      ln.eta.phi.M[u,t] ~ dunif(-10,0) #uninformative prior for M
    } #t

    ##Derived sums of discrete mortality for each release month
    #See definitions on L26-61
    Pr_timeD[u] <- 1-exp(-sum(eta.phi.Pr[u,1:(Periods-2)]))
    Pr_D[u] <- 1-exp(-sum(eta.phi.Pr[u,1:(Periods-1)]))
    M_timeD[u] <- 1-exp(-sum(eta.phi.M[u,1:(Periods-2)]))
    M_D[u] <- 1-exp(-sum(eta.phi.M[u,1:(Periods-1)]))
    S_timeD[u] <- 1 - sum(Pr_timeD[u],M_timeD[u])
    S_D[u] <- 1 - sum(Pr_D[u],M_D[u])
    
  }#u
    
    # Detection probability (p)
    for(i in 1:nFish){
      p[i,1] <- 1  # first time period p=1
    } #i
    
    for(t in 2:Periods){
      for(i in 1:nFish){
        p[i,t] <- p.g[group[i],t] ##varies by group and time
      } #i
      for(u in 1:g){
        p.g[u,t] ~ dunif(0,1) #uninformative prior
      } #u
      p.all[t] <- mean(p.g[1:3,t]) ##calculating overall detection probability for all groups
    } #t
  
  #Derived sums of discrete mortality for all release months
  #See definitions on L26-61
  Pr_allD <- mean(Pr_D) 
  Pr_alltimeD <- mean(Pr_timeD)
  M_allD <- mean(M_D) 
  M_alltimeD <- mean(M_timeD)
  S_allD <- mean(S_D) 
  S_alltimeD <- mean(S_timeD)
  
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
      po[1,i,t,1] <- p[i,t]   ##alive fish is detected alive in system
      po[1,i,t,2] <- 0        ##alive fish is detected as predated (cannot happen)
      po[1,i,t,3] <- 1-p[i,t] ##alive fish is not detected
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
group = group
g = as.numeric(length(unique(group)))
first <- numeric()        ##first observation of each fish
for (i in 1:dim(y)[1]) {
  first[i] <-min(which(y[i,]!=0))}
last <-numeric()          ##last observation of each fish
for (i in 1:dim(y)[1]) {
  last[i] <-max(which(y[i,]!=0))}
f <- first
l <- last

jags.data<-list(y=y,nFish=nFish,first=first,last=last,
                Periods=Periods, group = group,g = g)

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

##may need to specify initals for lnE and lnM, but uncertain how to do that.
jags.inits <- function(){list(ln.eta.phi.Pr = matrix(nrow = g,ncol = Periods - 1,
                                                     data = runif(g*(Periods-1),-10,0)),
                              ln.eta.phi.M = matrix(nrow = g,ncol = Periods - 1,
                                                    data = runif(g*(Periods-1),-10,0)),
                              z=ms.init.z(y,f))}


##parameters to monitor during model run 
params <-c('Pr_timeD',"Pr_D","M_timeD","M_D","S_D","S_timeD",
           "Pr_allD","Pr_alltimeD","M_allD","M_alltimeD","S_allD","S_alltimeD",
           "p.g","p.all")

##run model in JAGS
jagsfit <- autojags(data=jags.data, inits=jags.inits, parameters.to.save=params, model.file,
                    n.chains=3, n.adapt=NULL, iter.increment=3000, n.burnin=100000, n.thin=1,
                    save.all.iter=FALSE, modules=c('glm'), parallel=TRUE, n.cores = 4,
                    DIC=TRUE, store.data=FALSE, codaOnly=FALSE,
                    bugs.format=FALSE, Rhat.limit=1.01, max.iter=500000, verbose=TRUE)

##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
##Plotting and pulling model estimates
#Function to extract parameter estimates
jags_extract <- function(jagsfit){
  param_pattern <- "^(Pr|M|S)"
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
jagsfits_ <- jags_extract(jagsfit)

#Adding back release group data
RG = c(rep(c(1,2,3),times=6),rep("all",each=6))
jagsfits_ <- cbind(jagsfits_,RG)

#Cleaning up the df
jagsfits_ <- jagsfits_ %>% 
  mutate(parameter = case_when(grepl("M",param) ~ "M",
                               grepl("Pr",param) ~ "Pr",
                               grepl("S",param) ~ "S"),
         parameter_type = case_when(grepl("time",param) ~ "timeD",
                                    grepl("_D|allD",param) ~ "D")) 

#Coding variables as factors for plotting
jagsfits_$parameter <- factor(jagsfits_$parameter, levels=c("S","Pr","M"),labels=c("S","Pr","M"))
jagsfits_$RG <- factor(jagsfits_$RG, levels=c("1","2","3","all"),labels=c("March","April","May","All"))
jagsfits_$parameter_type <- factor(jagsfits_$parameter_type, levels=c("timeD","D"),labels=c("8-day","Overall"))

#Plotting model estimates (Figure 3)
est.plt <- jagsfits_ %>%
  ggplot(aes(x = parameter, fill=parameter))+
  geom_boxplot(aes(ymin = q2.5, lower = q25, middle = q50, upper = q75, ymax = q97.5),
               stat="identity")+
  facet_grid(cols=vars(RG),rows=vars(parameter_type))+
  scale_fill_manual(name="State",values=c("#FFB000","#785EF0","#FD005D"),
                    labels=c(
                      bquote('Survival (' * italic(S) * ')'),
                      bquote('Fish predation (' * italic(Pr) * ')'),
                      bquote('Unknown mortality (' * italic(M) * ')')
                    ))+
  theme_linedraw()+
  scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01))+
  theme(text=element_text(size=18))+
  ylab('Estimates (95% CrI)')+
  xlab(element_blank())+
  scale_x_discrete(
    labels = c(expression(italic(S)), expression(italic(Pr)), expression(italic(M)))
  )+
  theme(strip.background = element_rect(color="white",fill="white"),strip.text=element_text(color="black"))
est.plt

#Creating Table 5
jagsfit_dat_ <- round(jagsfits_[,3:7]*100,digits=1)
jagsfit_dat <- cbind(jagsfits_[,8:11],jagsfit_dat_)
jagsfit_dat <- jagsfit_dat[,c(2:5,7,9)]
Table5 <- jagsfit_dat %>%
  mutate(Estimate = paste0(q50," (",q2.5,"-",q97.5,")"),
         header = paste0(parameter,".",parameter_type)) %>%
  dplyr::select(!c(q2.5,q50,q97.5,parameter_type,parameter)) %>%
  pivot_wider(id_cols=RG,values_from=Estimate,names_from=header) %>%
  dplyr::select(RG,`S.8-day`,S.Overall,`Pr.8-day`,Pr.Overall,`M.8-day`,
                M.Overall)

#Pulling detection probability data from the jagsfit object
p.g.all <- list()
for(i in 2:10){
  p.g <- data.frame(q2.5 = jagsfit$q2.5$p.g[,i], q25 = jagsfit$q25$p.g[,i], 
                    q50 = jagsfit$q50$p.g[,i], q75 = jagsfit$q75$p.g[,i],
                    q97.5 = jagsfit$q97.5$p.g[,i]) %>%
    mutate(P = i) 
  RG <- c(1,2,3) 
  p.g <- cbind(p.g,RG)
  p.g.all[[i]] <- p.g
}
p.g.all <- do.call("rbind",p.g.all)

p.all <- data.frame(q2.5 = jagsfit$q2.5$p.all, q25 = jagsfit$q25$p.all, 
                    q50 = jagsfit$q50$p.all, q75 = jagsfit$q75$p.all,
                    q97.5 = jagsfit$q97.5$p.all) 
p.all <- p.all %>%
  mutate(P=rownames(p.all)) %>%
  mutate(RG = "All") %>%
  drop_na()

p.all <- rbind(p.g.all,p.all)
#Filling in 1s as P[t=1]
fill <- c(1,1,1,1)
p.1 <- tibble(q2.5=fill,q25=fill,q50=fill,q75=fill,q97.5=fill,P=fill,RG=c(1,2,3,"All"))
p.all <- rbind(p.1,p.all)

#Factoring parameters for plotting
p.all$RG <- factor(p.all$RG,levels=c(1,2,3,"All"),labels=c("March","April","May","All"))
p.all$P <- factor(p.all$P,levels=c(1:10))

#Plotting detection probability (Figure S2)
p.plt <- p.all %>%
  ggplot(aes(x = P))+
  geom_boxplot(aes(ymin = q2.5, lower = q25, middle = q50, upper = q75, ymax = q97.5),
               stat="identity",fill="grey")+
  facet_grid(cols=vars(RG))+
  theme_linedraw()+
  scale_y_continuous(limits=c(0,1),expand=c(0.01,0.01))+
  theme(text=element_text(size=18))+
  ylab('Estimates (95% CrI)')+
  xlab("p[t]")+
  theme(strip.background = element_rect(color="white",fill="white"),strip.text=element_text(color="black"))
p.plt
setwd("~/Desktop/")
ggsave("Group_ps.png",height=9,width=12,unit="in",dpi=600)
