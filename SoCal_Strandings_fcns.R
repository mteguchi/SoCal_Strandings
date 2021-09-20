

# extracts loo statistics, including looic and Pareto-k statistics
pareto.k.diag <- function(jm, MCMC.params, jags.data){
  
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.obs <- jm$sims.list$loglik[, !is.na(jags.data$y)]
  # get rid of NA columns - even if data existed (for example the first value) - no likelihood
  # for the first data point
  loglik.obs <- loglik.obs[, colSums(is.na(loglik.obs)) == 0]
  
  #loglik.obs <- jm$sims.list$loglik[, 2:jags.data$T]
  # cores = 1 is needed in the relative_eff function if the number of cores was set to more than
  # 1 with options(mc.cores = parallel::detectCores()) or something similear. See also here:
  # https://discourse.mc-stan.org/t/error-in-loo-relative-eff-related-to-options-mc-cores/5610/2
  
  Reff <- relative_eff(exp(loglik.obs), 
                       chain_id = rep(1:MCMC.params$n.chains, 
                                      each = n.per.chain),
                       cores = 1)
  
  loo.out <- loo(loglik.obs, r_eff = Reff, cores = 1)
  return(list(loglik.obs = loglik.obs,
              Reff = Reff,
              loo.out = loo.out))
}


# computes LOOIC
compute.LOOIC <- function(loglik, data.vector, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.vec <- as.vector(loglik)
  
  # each column corresponds to a data point and rows are MCMC samples
  loglik.mat <- matrix(loglik.vec, nrow = n.per.chain * MCMC.params$n.chains)
  
  # take out the columns that correspond to missing data points
  loglik.mat <- loglik.mat[, !is.na(data.vector)]
  # loglik.mat <- matrix(loglik.vec[!is.na(data.vector)], 
  #                      nrow = MCMC.params$n.chains * n.per.chain)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  #
  loo.out <- loo(loglik.mat, 
                 r_eff = Reff, 
                 cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}


# Extracting posterior samples from jags output:
extract.samples <- function(varname, zm){
  dev <- do.call(rbind, zm[, varname])
  return(dev)
}


# default inputs are for the SCB.
getCoastLine <- function(filename = "data/coast/coast_Epac.txt",
                         lon.limits = c(-121, -116.5),
                         lat.limits = c(32, 35.2)){
  data1 <- read.table(filename, header = F, row.names = NULL, sep = ',')
  # fix longitudes to -180 to 180
  lon.limits[lon.limits > 180] <- lon.limits[lon.limits > 180] - 360
  out.df <- data1[data1[,1] >= min(lon.limits) &
                    data1[,1] <= max(lon.limits) &
                    data1[,2] >= min(lat.limits) &
                    data1[,2] <= max(lat.limits),]
  colnames(out.df) <- c('Longitude', 'Latitude')
  out.df$idx <- NA
  idx.rows <- 1:nrow(out.df)
  na.idx <- is.na(out.df$Longitude)  # T/F for NA
  dif.na.idx <- na.idx[2:length(na.idx)] - na.idx[1:(length(na.idx)-1)]
  idx.neg1 <- idx.rows[dif.na.idx == -1]  # beginning of segments
  idx.pos1 <- idx.rows[dif.na.idx == 1]   # end of segments
  
  for (k in 1:length(idx.neg1)) {
    out.df[idx.neg1[k]:idx.pos1[k], 'idx'] <- k
  }
  
  # change index to a factor variable
  out.df$idx <- as.factor(out.df$idx)
  out.df <- na.omit(out.df)
  # this splits all segments into separate dataframes
  out.list <- split(out.df, out.df$idx)
  
  for (k in 1:length(out.list)){
    line.1 <- out.list[[k]][1,]
    line.end <- out.list[[k]][nrow(out.list[[k]]),]
    if (line.1$Longitude != line.end$Longitude){
      tmp <- out.list[[k]]
      tmp <- rbind(c(max(tmp$Longitude), max(tmp$Latitude), k),
                   tmp, c(max(tmp$Longitude), max(tmp$Latitude), k))
      out.list[[k]] <- tmp
    }
  }
  return(out.list)
  
}