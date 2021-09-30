
# runs one model 
run.model <- function(parameters.list, jags.data, model.name, MCMC.params){
  
  if (!file.exists(paste0("RData/jags_out_", model.name, ".rds"))){
    jm <- jags(data = jags.data,
               #inits = inits,
               parameters.to.save= parameters.list,
               model.file = paste0("models/model_DcStranding_", model.name, ".txt"),
               n.chains = MCMC.params$n.chains,
               n.burnin = MCMC.params$n.burnin,
               n.thin = MCMC.params$n.thin,
               n.iter = MCMC.params$n.samples,
               DIC = T, 
               parallel=T)
    
    saveRDS(jm, file = paste0("RData/jags_out_", model.name, ".rds"))
    
  } else {
    
    jm <- readRDS(file = paste0("RData/jags_out_", model.name, ".rds"))
  }
  
  
  jm.MCMC <- MCMC.diag(jm = jm, MCMC.params = MCMC.params)
  
  return(list(jm.out = jm,
              jm.MCMC = jm.MCMC))  
}


# extract MCMC diagnostic statistics, including Rhat, loglikelihood, DIC, and LOOIC
MCMC.diag <- function(jm, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  Rhat <- unlist(lapply(jm$Rhat, FUN = max))
  
  loglik <- cbind(jm$sims.list$loglik[,,1], 
                  jm$sims.list$loglik[,,2])
  
  Reff <- relative_eff(exp(loglik), 
                       chain_id = rep(1:MCMC.params$n.chains, 
                                      each = n.per.chain),
                       cores = 1)
  loo.out <- loo(loglik, r_eff = Reff, cores = 1)
  
  return(list(DIC = jm$DIC,
              loglik.obs = loglik,
              Reff = Reff,
              Rhat = Rhat,
              loo.out = loo.out))
  
}

