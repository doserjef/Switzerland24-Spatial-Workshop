library(spOccupancy)
library(tidyverse)


# Simulate data ----------------------------------------------------------------

# set parameters
psi <- 0.4
p <- 0.7

# simulate a dataset
set.seed(3)
simdata <- simOcc(J.x=50, J.y=20,       # 50*20 sites distributed on a grid with 1000 total sites
                  n.rep=rep(10, 50*20), # 10 visits to each site
                  beta=c(qlogis(psi),   # occupancy intercept
                         -0.4),         # coefficient of occupancy covariate, e.g. forest density
                  alpha=c(qlogis(p),    # detection intercept
                          0.6))         # coefficient of detection covariate, e.g. temperature
str(simdata)
simdata$y[1:10,1:5]
simdata$X[1:10,2]
simdata$X.p[1:10,1:5,2]

# subset detection data according to sampling design
nvisits <- 4
nsites <- 500

# select sites randomly
index <- sample(1:1000, nsites)

# bundle data for spOccupancy
data.list <- list(y = simdata$y[index, 1:nvisits],
                  occ.covs = data.frame(forest=simdata$X[index,2]),
                  det.covs = list(temperature=simdata$X.p[index,1:nvisits,2]))
str(data.list)

# fit model
out <- PGOcc(occ.formula = ~ forest,
             det.formula = ~ temperature,
             data=data.list,
             n.samples=2000, n.thin=2, n.burn=500, n.chains=3, n.report=500)

# extract estimates
str(out)
head(out$beta.samples)
mean(out$beta.samples[,1])
sd(out$beta.samples[,1])
summary(out)



# Run multiple simulations with different parameter settings -------------------


# define parameters
nsims <- 10
psi <- c(0.25, 0.5, 0.75)
p <- c(0.25, 0.5, 0.75)
total.surveys <- 2000
nvisits <- 2:10

# prepare array to hold estimated standard error of occupancy estimate
results <- array(data=NA, dim=c(length(psi), length(p), length(nvisits), nsims),
                 dimnames=list(paste('psi',psi,sep='_'),
                               paste('p',p,sep='_'),
                               paste('K',nvisits,sep='_'),
                               NULL))
dim(results)
dimnames(results)

# run simulation
# approximate run time: 1 hour
for(i in 1:length(psi)){
  for(j in 1:length(p)){
    for(s in 1:nsims){
      cat(paste('\n*** Doing simulation ', s, ' with psi=', psi[i], ' and p=', p[j], ' ***', sep=''))
      
      # simulate data (same for all sampling designs)
      simdata <- simOcc(J.x=50, J.y=20, n.rep=rep(10, 50*20),
                        beta=c(qlogis(psi[i]), # occupancy intercept
                               -0.4),
                        alpha=c(qlogis(p[j]),  # detection intercept
                                0.6))
      
      for(k in 1:length(nvisits)){
        
        # sample from the simulated data according to the sampling design
        nsites <- round(total.surveys/nvisits[k],0)
        index <- sample(1:1000, nsites)
        data.list <- list(y=simdata$y[index, 1:nvisits[k]],
                          occ.covs=data.frame(forest=simdata$X[index,2]),
                          det.covs=list(temperature=simdata$X.p[index,1:nvisits[k],2]))
        
        # fit model
        out <- PGOcc(occ.formula = ~ forest,
                     det.formula = ~ temperature,
                     data=data.list,
                     n.samples=1500, n.thin=2, n.burn=500, n.chains=3, verbose=FALSE)
        
        # save estimated standard error of occupancy estimate
        results[i,j,k,s] <- sd(out$beta.samples[,1])
      }
    }
  }
}
# saveRDS(results, file='Occupancy_model_survey_design_simulation.rds')

# visualize results
results <- readRDS('Occupancy_model_survey_design_simulation.rds')
dimnames(results)
dimnames(results) <- list(psi, p, nvisits, NULL)
results.df <- array2DF(results) %>% 
  rename_at(1:5, ~ c('psi', 'p', 'K', 'sim', 'SEpsi')) %>% 
  mutate(across(K, as.numeric))
head(results.df)

ggplot(data=results.df, aes(x=K, y=SEpsi, group=p, col=p)) +
  geom_point(size=2) +
  geom_smooth(linewidth=2) +
  facet_grid(~ psi) +
  scale_color_manual(values=c('orange','turquoise3','darkorchid3')) +
  ggtitle(label=expression(paste('True occupancy probability ',Psi))) +
  ylab(label=expression(paste(widehat(SE),' of ',widehat(Psi)))) +
  xlab(label=expression(paste('Number of repeated visits ',italic('K')))) +
  theme_bw() +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14, face='bold'),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14, face='bold'),
        strip.text=element_text(size=12, face='bold'))

# illustrate the possible uncertainty around the occupancy estimate
SEpsi <- c(0.1, 0.35, 0.6)
par(mfrow=c(3,3), mar=c(3,3,4,1))
for(e in 1:length(SEpsi)){
  for(i in 1:length(psi)){
    psi.logit <- rnorm(100, mean=qlogis(psi[i]), sd=SEpsi[e])
    hist(plogis(psi.logit), main=paste('SE =',SEpsi[e]), xlab='', ylab='', las=1,
         xlim=c(0,1), breaks=seq(0,1, by=0.02), col=c('chartreuse4','cyan4','saddlebrown')[e])
  }
}
