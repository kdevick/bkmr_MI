#################################################
##  source code used to create BKMR plot       ##
##  dataframes when data are multiply imputed  ##
##  on the covariates and/or the exposures     ##
#################################################



###############################################################################
## Function to compile the exposures from all MI datasets                    ##
##               														                                 ##
## this step is needed to so that the contrast used when estimating the      ##
## effects across the MI BKMR fits is consistent. (i.e. the same 25th, 50th, ## 
## and 75th percentile of the metals is used when estimating the effect of a ##
## change in all metals at their 25th to all their 75th on the outcome for   ##
## each MI BKMR fit).                                                        ##
## If the data are not imputed in the exposures, the original Z matrix can   ##                                      
## used as Z.complete.MI                                                     ##
###############################################################################

Z.complete.MI <- function(BKMRfits){
	n <- nrow(BKMRfits[[1]]$Z)
	l <- ncol(BKMRfits[[1]]$Z)
	K <- length(BKMRfits)
	
  	Z.full <- matrix(NA, nrow=n*K, ncol=l)
  	ifelse(is.null(colnames(BKMRfits[[1]]$Z)), colnames(Z.full) <- paste0("z", 1:l), colnames(Z.full) <- colnames(BKMRfits[[1]]$Z))
  	
	for(k in 1:K){
  		fit <- BKMRfits[[k]]
		Z.full[(n*(k-1)+1):(n*k),] <- fit$Z
	}
	Z.full
}


################################################################################
## Function to implment Rubins 1987 method to estimate point estimates and    ##
## std errors when combining information across MI fits (for approx method)   ##
################################################################################


Rubin.MI <- function(mean.vec, variance.vec){
  K    <- length(mean.vec)
  qbar <- mean(mean.vec)
  wbar <- mean(variance.vec)
  B    <- var(mean.vec)
  var  <- wbar+(1+1/K)*B
  c(est=qbar, sd=sqrt(var))
}



################################################################################
## Edited ComputePostmeanHnew.exact function 								                  ##
## 																			                                      ##
## This function now returns the entire mean matrix and variance array        ##
## needed for each of the future functions to obtain an unbiased estimate of  ## 
## the se used to create CI in plots                                          ##
################################################################################


#' ****** EDITED TO OBTAIN THE ENTIRE MEAN MATRIX AND VARIANCE ARRAY ******
#' Compute the posterior mean and variance of \code{h} at a new predictor values
#' Function to estimate the posterior mean and variance by obtaining the posterior mean and variance at particular iterations and then using the iterated mean and variance formulas
#' **** for MI BMKR fits ****

ComputePostmeanHnew.exact.MI <- function(fit, y = NULL, Z = NULL, X = NULL, Znew = NULL, sel = NULL) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) Znew <- matrix(Znew, nrow = 1)
    if (class(Znew) == "data.frame") Znew <- data.matrix(Znew)
    if (ncol(Z) != ncol(Znew)) {
      stop("Znew must have the same number of columns as Z")
    }
  }
  
  if (is.null(dim(X))) X <- matrix(X, ncol=1)
  
  # if (!is.null(fit$Vinv)) {
  #   sel <- attr(fit$Vinv, "sel")
  # }
  
  if (is.null(sel)) {
    sel <- with(fit, seq(floor(iter/2) + 1, iter, 10))
    if (length(sel) < 100) {
      sel <- with(fit, seq(floor(iter/2) + 1, iter, length.out = 100))
    }
    sel <- unique(floor(sel))
  }
  
  family <- fit$family
  data.comps <- fit$data.comps
  post.comps.store <- list(postmean = vector("list", length(sel)),
                           postvar = vector("list", length(sel))
  )
  
  for (i in seq_along(sel)) {
    s <- sel[i]
    beta <- fit$beta[s, ]
    lambda <- fit$lambda[s, ]
    sigsq.eps <- fit$sigsq.eps[s]
    r <- fit$r[s, ]
    
    if (family == "gaussian") {
      ycont <- y
    } else if (family == "binomial") {
      ycont <- fit$ystar[s, ]
    }
    
    Kpart <- makeKpart(r, Z)
    K <- exp(-Kpart)
    Vcomps <- makeVcomps(r = r, lambda = lambda, Z = Z, data.comps = data.comps)
    Vinv <- Vcomps$Vinv
    # if (is.null(fit$Vinv)) {
    # V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*K
    # cholV <- chol(V)
    # Vinv <- chol2inv(cholV)
    # } else {
    #   Vinv <- fit$Vinv[[i]]
    # }
    
    if (!is.null(Znew)) {
      # if(is.null(data.comps$knots)) {
      n0 <- nrow(Z)
      n1 <- nrow(Znew)
      nall <- n0 + n1
      Kpartall <- makeKpart(r, rbind(Z, Znew))
      Kmat <- exp(-Kpartall)
      Kmat0 <- Kmat[1:n0,1:n0 ,drop=FALSE]
      Kmat1 <- Kmat[(n0+1):nall,(n0+1):nall ,drop=FALSE]
      Kmat10 <- Kmat[(n0+1):nall,1:n0 ,drop=FALSE]
      
      lamK10Vinv <- lambda[1]*Kmat10 %*% Vinv
      postvar <- lambda[1]*sigsq.eps*(Kmat1 - lamK10Vinv %*% t(Kmat10))
      postmean <- lamK10Vinv %*% (ycont - X%*%beta)
      # } else {
      # stop("GPP not yet implemented")
      # }
    } else {
      lamKVinv <- lambda[1]*K%*%Vinv
      postvar <- lambda[1]*sigsq.eps*(K - lamKVinv%*%K)
      postmean <- lamKVinv %*% (ycont - X%*%beta)
    }
    
    post.comps.store$postmean[[i]] <- postmean
    post.comps.store$postvar[[i]] <- postvar
    
  }
  
  postmean_mat <- t(do.call("cbind", post.comps.store$postmean))
  m <- colMeans(postmean_mat)
  postvar_arr <- with(post.comps.store, 
                      array(unlist(postvar), 
                            dim = c(nrow(postvar[[1]]), ncol(postvar[[1]]), length(postvar)))
  )
  ve <- var(postmean_mat)
  ev <- apply(postvar_arr, c(1, 2), mean)
  v <- ve + ev
  ret <- list(postmean = m, postvar = v, postmean_mat = postmean_mat, postvar_arr = postvar_arr)
  
  ret
}




#######################################################################
## Edited OverallRiskSummaries function 				        	           ##
#######################################################################


OverallRiskSummaries.MI <- function(BKMRfits, qs = seq(0.25, 0.75, by = 0.05), q.fixed = 0.5, q.alwaysfixed = NULL, index.alwaysfixed = NULL, sel = NULL, method = "approx") {
  
  start.time <- proc.time()["elapsed"]
  cc <- c(-1, 1)
  K <- length(BKMRfits)
  
  Z.MI <- Z.complete.MI(BKMRfits)
  
  toreturn <- data.frame(quantile=qs, est=rep(NA,times=length(qs)), sd=rep(NA,times=length(qs)))
  
  if(method=="exact") {
    print("exact method")
    preds.fun <- function(znew) ComputePostmeanHnew.exact.MI(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    for(i in 1:length(qs)){
      quant <- qs[i]
      
      ## 2 = nrow(newz)  	
      postmean.temp <- matrix(NA, nrow=length(sel)*K, ncol=2)
      postvar.temp  <- array(NA, dim = c(2,2,length(sel)*K)) 
      for(k in 1:K){
        fit <- BKMRfits[[k]]
        y <- fit$y
        Z <- fit$Z
        X <- fit$X
        
        point1 <- apply(Z.MI, 2, quantile, q.fixed)
        point2 <- apply(Z.MI, 2, quantile, quant)
        
        ## if both q.alwaysfixed and index.alwaysfixed are specified,
        ## change the values in the point which we want to keep fixed for all comparisons
        if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
          point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(Z.MI[,index.alwaysfixed, drop=FALSE],2,quantile, q.alwaysfixed)
        }
        
        newz <- rbind(point1, point2)
        
        preds <- preds.fun(newz)
        postmean.temp[((k-1)*length(sel)+1):(length(sel)*k),] <- preds$postmean_mat
        postvar.temp[,,((k-1)*length(sel)+1):(length(sel)*k)] <- preds$postvar_arr
      }
      m  <- colMeans(postmean.temp)
      ve <- var(postmean.temp)
      ev <- apply(postvar.temp, c(1, 2), mean)
      v  <- ve + ev
      
      toreturn[i,"est"] <- drop(cc %*% m)
      toreturn[i,"sd"]  <- drop(sqrt(cc %*% v %*% cc))
      
      end.time <- proc.time()["elapsed"]
      print(paste(i,"out of", length(qs), "complete: ", round((end.time - start.time)/60, digit=2), "min run time" ))				
    }
  } else if(method=="approx") {
    print("approx method")
    preds.fun <- function(znew) ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    
    for(i in 1:length(qs)){
      quant <- qs[i]
      ## 2 = nrow(newz)  	
      
      est.vec <- rep(NA, times=K)
      var.vec <- rep(NA, times=K)
      for(k in 1:K){
        fit <- BKMRfits[[k]]
        y <- fit$y
        Z <- fit$Z
        X <- fit$X
        
        point1 <- apply(Z.MI, 2, quantile, q.fixed)
        point2 <- apply(Z.MI, 2, quantile, quant)
        
        ## if both q.alwaysfixed and index.alwaysfixed are specified,
        ## change the values in the point which we want to keep fixed for all comparisons
        if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
          point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(as.matrix(Z.MI[,index.alwaysfixed]),2,quantile, q.alwaysfixed)
        }
        
        newz <- rbind(point1, point2)
        
        preds <- preds.fun(newz)
        est.vec[k] <- drop(cc %*% preds$postmean)
        var.vec[k] <- drop(cc %*% preds$postvar %*% cc)
      }
      
      if(K==1){MIest <- c(est=est.vec, sd=sqrt(var.vec))}else{MIest <- Rubin.MI(mean.vec = est.vec, variance.vec = var.vec)}
      toreturn[i,"est"] <- MIest["est"]
      toreturn[i,"sd"]  <- MIest["sd"]
      
      end.time <- proc.time()["elapsed"]
      print(paste(i,"out of", length(qs), "complete: ", round((end.time - start.time)/60, digit=2), "min run time" ))				
    }
  } else stop("method must be one of c('approx', 'exact')") 
  toreturn  
}




#######################################################################
## Edited VarRiskSummary and SingVarRiskSummaries functions 	  	   ##
#######################################################################


### ***** a combination function of VarRiskSummary and riskSummary.approx for MI BKMR fits ******
#Compare estimated \code{h} function when a single variable (or a set of variables) is at the 75th versus 25th percentile, when all of the other variables are fixed at a particular percentile

VarRiskSummary.MI <- function(whichz = 1, BKMRfits, Z.MI, qs.diff = c(0.25, 0.75), q.fixed = 0.5, q.alwaysfixed = NULL, index.alwaysfixed = NULL, sel = NULL, method = "approx") {
  
  cc <- c(-1, 1)
  K <- length(BKMRfits)
  
  if(method=="exact") {
    preds.fun <- function(znew) ComputePostmeanHnew.exact.MI(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    
    postmean.temp <- matrix(NA, nrow=length(sel)*K, ncol=2)
    postvar.temp  <- array(NA, dim = c(2,2,length(sel)*K))## 2 = nrow(newz)
    
    for(k in 1:K){
      fit <- BKMRfits[[k]]
      y <- fit$y
      Z <- fit$Z
      X <- fit$X
      
      point1 <- point2 <- apply(Z.MI, 2, quantile, q.fixed)
      point2[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
      point1[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[1])
      
      ## if both q.alwaysfixed and index.alwaysfixed are specified,
      ## change the values in the point which we want to keep fixed for all comparisons
      if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
        point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(Z.MI[,index.alwaysfixed, drop=FALSE],2,quantile, q.alwaysfixed)
      }
      
      newz <- rbind(point1, point2)
      
      preds <- preds.fun(newz)
      
      postmean.temp[((k-1)*length(sel)+1):(length(sel)*k),] <- preds$postmean_mat
      postvar.temp[,,((k-1)*length(sel)+1):(length(sel)*k)] <- preds$postvar_arr
      
    }
    m  <- colMeans(postmean.temp)
    ve <- var(postmean.temp)
    ev <- apply(postvar.temp, c(1, 2), mean)
    v  <- ve + ev
    
    diff     <- drop(cc %*% m)
    diff.sd  <- drop(sqrt(cc %*% v %*% cc))	
    toreturn <- c(est = diff, sd = diff.sd)
    
  } else if(method=="approx") {
    preds.fun <- function(znew) ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, Znew = znew, sel = sel)
    
    est.vec <- rep(NA, times=K)
    var.vec <- rep(NA, times=K)
    
    for(k in 1:K){
      fit <- BKMRfits[[k]]
      y <- fit$y
      Z <- fit$Z
      X <- fit$X
      
      point1 <- point2 <- apply(Z.MI, 2, quantile, q.fixed)
      point2[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[2])
      point1[whichz] <- apply(Z.MI[, whichz, drop = FALSE], 2, quantile, qs.diff[1])
      
      ## if both q.alwaysfixed and index.alwaysfixed are specified,
      ## change the values in the point which we want to keep fixed for all comparisons
      if(!is.null(q.alwaysfixed) & !is.null(index.alwaysfixed)){
        point1[index.alwaysfixed] <- point2[index.alwaysfixed] <- apply(Z.MI[,index.alwaysfixed, drop=FALSE],2,quantile, q.alwaysfixed)
      }
      
      newz <- rbind(point1, point2)
      
      preds <- preds.fun(newz)
      
      est.vec[k] <- drop(cc %*% preds$postmean)
      var.vec[k] <- drop(cc %*% preds$postvar %*% cc)
    }
    if(K==1){toreturn <- c(est=est.vec, sd=sqrt(var.vec))}else{toreturn <- Rubin.MI(mean.vec = est.vec, variance.vec = var.vec)}
  } else stop("method must be one of c('approx', 'exact')")
  toreturn
}



#' Single Variable Risk Summaries
#' 
#' Compute summaries of the risks associated with a change in a single variable in \code{Z} from a single level (quantile) to a second level (quantile), for the other variables in \code{Z} fixed to a specific level (quantile)
#' **** for MI BKMR fits ****

SingVarRiskSummaries.MI <- function(BKMRfits, which.z = 1:ncol(BKMRfits[[1]]$Z), qs.diff = c(0.25, 0.75), q.fixed = c(0.25, 0.50, 0.75), q.alwaysfixed = NULL, index.alwaysfixed = NULL, sel = NULL, z.names = colnames(BKMRfits[[1]]$Z), method="approx",...) {
  
  start.time <- proc.time()["elapsed"]
  if(is.null(z.names)) z.names <- paste0("z", 1:ncol(BKMRfits[[1]]$Z))
  Z.MI <- Z.complete.MI(BKMRfits)
  
  df <- dplyr::data_frame()
  for(i in seq_along(q.fixed)) {
    for(j in seq_along(which.z)) {
      risk <- VarRiskSummary.MI(whichz = which.z[j], BKMRfits=BKMRfits, Z.MI=Z.MI, qs.diff = qs.diff, q.fixed = q.fixed[i], q.alwaysfixed = q.alwaysfixed, index.alwaysfixed = index.alwaysfixed, sel = sel, method = method, ...)
      df0 <- dplyr::data_frame(q.fixed = q.fixed[i], variable = z.names[which.z[j]], est = risk["est"], sd = risk["sd"])
      df <- dplyr::bind_rows(df, df0)
    }
    
    end.time <- proc.time()["elapsed"]
    print(paste(i,"out of", length(q.fixed), "complete: ", round((end.time - start.time)/60, digit=2), "min run time" ))
  }
  df <- dplyr::mutate_(df, variable = ~factor(variable, levels = z.names[which.z]), q.fixed = ~as.factor(q.fixed))
  attr(df, "qs.diff") <- qs.diff
  df
}




#############################################################################
## Edited PredictorResponseUnivarVar and PredictorResponseUnivar functions ##
#############################################################################

#### ****** right now the min.plot.dist is using the Z matrix comprised of data from 
#### all MI Z datasets. If there is truly a min distance you want to observe between
#### observed and predicted, then this will need to be moved inside the K loop so that 
#### it assess this for each dataset individiaully instead of collectively. For datasets 
#### with minimal missing values in Z, this shouldn't make a difference ******

#### **** if the default value of "Inf" is used, then this does not matter ****


PredictorResponseUnivarVar.MI <- function(whichz = 1, BKMRfits, Z.MI, ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE, method = "approx",...) {

    K <- length(BKMRfits)   
    ord <- c(whichz, setdiff(1:ncol(Z.MI), whichz))
    z1 <- seq(min(Z.MI[,ord[1]]), max(Z.MI[,ord[1]]), length = ngrid)
    z.others <- lapply(2:ncol(Z.MI), function(x) quantile(Z.MI[,ord[x]], q.fixed))
    z.all <- c(list(z1), z.others)
    newz.grid <- expand.grid(z.all)
    colnames(newz.grid) <- colnames(Z.MI)[ord]
    newz.grid <- newz.grid[,colnames(Z.MI)]

    if (!is.null(min.plot.dist)) {
        mindists <- rep(NA,nrow(newz.grid))
        for (i in seq_along(mindists)) {
            pt <- as.numeric(newz.grid[i, colnames(Z.MI)[ord[1]]])
            dists <- fields::rdist(matrix(pt, nrow = 1), Z.MI[, colnames(Z.MI)[ord[1]]])
            mindists[i] <- min(dists)
        }
    }
    
    
    
    if(method=="exact") {
      postmean.temp <- matrix(NA, nrow=length(sel)*K, ncol=ngrid)
      postvar.temp  <- array(NA, dim = c(ngrid,ngrid,length(sel)*K))
      
      for(k in 1:K){
      	fit <- BKMRfits[[k]]
    		y <- fit$y
    		Z <- fit$Z
    		X <- fit$X
    
    		preds <- ComputePostmeanHnew.exact.MI(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel)
      	
      	postmean.temp[((k-1)*length(sel)+1):(length(sel)*k),] <- preds$postmean_mat
      	postvar.temp[,,((k-1)*length(sel)+1):(length(sel)*k)] <- preds$postvar_arr
    		
      }
      
      ve <- var(postmean.temp)
      ev <- apply(postvar.temp, c(1, 2), mean)
      v  <- ve + ev
   
      preds.plot <- colMeans(postmean.temp)
      se.plot <- sqrt(diag(v))
  
    } else if(method=="approx") {

      postmean.temp <- matrix(NA, nrow=K, ncol=ngrid)
      postvar.temp  <- matrix(NA, nrow=K, ncol=ngrid)
      
      for(k in 1:K){
        fit <- BKMRfits[[k]]
        y <- fit$y
        Z <- fit$Z
        X <- fit$X
        
        preds <- ComputePostmeanHnew.approx(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel)
        
        postmean.temp[k,] <- preds$postmean
        postvar.temp[k,]  <- diag(preds$postvar)
        
      }
      
      temp <- sapply(1:ngrid, function(x){Rubin.MI(mean.vec = postmean.temp[,x], variance.vec = postvar.temp[,x])})
      
      preds.plot <- temp["est",]
      se.plot    <- temp["sd",]
      
    } else stop("method must be one of c('approx', 'exact')")  
      
    
    if(center) preds.plot <- preds.plot - mean(preds.plot)
    if(!is.null(min.plot.dist)) {
        preds.plot[mindists > min.plot.dist] <- NA
        se.plot[mindists > min.plot.dist] <- NA
    }

    res <- dplyr::data_frame(z = z1, est = preds.plot, se = se.plot)
}




PredictorResponseUnivar.MI <- function(BKMRfits, which.z = 1:ncol(BKMRfits[[1]]$Z), ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = Inf, center = TRUE, method = "approx", ...) {
  
  start.time <- proc.time()["elapsed"]
  Z.MI <- Z.complete.MI(BKMRfits)
  z.names <- colnames(Z.MI)

  df <- dplyr::data_frame()
  for(i in which.z) {
    res <- PredictorResponseUnivarVar.MI(whichz = i, BKMRfits = BKMRfits, Z.MI = Z.MI, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, method = method, ...)
    df0 <- dplyr::mutate(res, variable = z.names[i]) %>%
      dplyr::select_(~variable, ~z, ~est, ~se)
    df <- dplyr::bind_rows(df, df0)
    
    end.time <- proc.time()["elapsed"]
    print(paste(i,"out of", length(which.z), "complete: ", round((end.time - start.time)/60, digit=2), "min run time" ))
  }
  df$variable <- factor(df$variable, levels = z.names[which.z])
  df
}






#############################################################################
## Edited PredictorResponseBivarPair, PredictorResponseBivar, and          ##
## PredictorResponseBivarLevels functions 								                 ##
#############################################################################


#### ****** right now the min.plot.dist is using the Z matrix for each MI BKMR fit individually. ******
#### It will return a NA value if the observed Z matrix is less than 0.5 from any point to be predicted 
#### (or whatever min.plot.dist is set to). The function PredictorResponseBivar.MI used create the dataframe 
#### for plotting the bivariate response curves averages over these values with na.rm=FALSE. Thus, unless the
#### min.plot.dist holds for ALL MI datasets, the point will not be considered (it returns NA).

####  ****** standard errors (SE) are not corrected for MI since the plots do not show these, thus the 
#### SE returned by the function PredictorResponseBivarLevels are biased and not to be trusted ******

PredictorResponseBivarPair.MI <- function(fit, y, Z, X, whichz1 = 1, whichz2 = 2, whichz3 = NULL, method = "approx", prob = 0.5, q.fixed = 0.5, sel = NULL, ngrid = 50, min.plot.dist = 0.5, center = TRUE, Z.MI, ...) {
    if(ncol(Z) < 3) stop("requires there to be at least 3 Z variables")

    if(is.null(colnames(Z))) colnames(Z) <- paste0("z", 1:ncol(Z))

    if(is.null(whichz3)) {
        ord <- c(whichz1, whichz2, setdiff(1:ncol(Z), c(whichz1, whichz2)))
    } else {
        ord <- c(whichz1, whichz2, whichz3, setdiff(1:ncol(Z), c(whichz1, whichz2, whichz3)))
    }
    z1 <- seq(min(Z.MI[,ord[1]]), max(Z.MI[,ord[1]]), length=ngrid)
    z2 <- seq(min(Z.MI[,ord[2]]), max(Z.MI[,ord[2]]), length=ngrid)
    z3 <- quantile(Z.MI[, ord[3]], probs = prob)
    z.all <- c(list(z1), list(z2), list(z3))
    if(ncol(Z) > 3) {
        z.others <- lapply(4:ncol(Z), function(x) quantile(Z.MI[,ord[x]], q.fixed))
        z.all <- c(z.all, z.others)
    }
    newz.grid <- expand.grid(z.all)
    z1save <- newz.grid[, 1]
    z2save <- newz.grid[, 2]
    colnames(newz.grid) <- colnames(Z)[ord]
    newz.grid <- newz.grid[,colnames(Z)]

    if(!is.null(min.plot.dist)) {
        mindists <- rep(NA, nrow(newz.grid))
        for(k in seq_along(mindists)) {
            pt <- as.numeric(newz.grid[k,c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
            dists <- fields::rdist(matrix(pt, nrow = 1), Z[, c(colnames(Z)[ord[1]],colnames(Z)[ord[2]])])
            mindists[k] <- min(dists)
        }
    }

    if (method %in% c("approx", "exact")) {
      preds <- ComputePostmeanHnew(fit = fit, y = y, Z = Z, X = X, Znew = newz.grid, sel = sel, method = method)
      preds.plot <- preds$postmean
      se.plot <- sqrt(diag(preds$postvar))
    } else {
      stop("method must be one of c('approx', 'exact')")
    }
    if(center) preds.plot <- preds.plot - mean(preds.plot)
    if(!is.null(min.plot.dist)) {
        preds.plot[mindists > min.plot.dist] <- NA
        se.plot[mindists > min.plot.dist] <- NA
    }
#     hgrid <- matrix(preds.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))
#     se.grid <- matrix(se.plot, ngrid, ngrid, dimnames=list(z1=round(z1,2), z2=round(z2,2)))

    res <- dplyr::data_frame(z1 = z1save, z2 = z2save, est = preds.plot, se = se.plot)
}

#' Predict the exposure-response function at a new grid of points
#'
#' Predict the exposure-response function at a new grid of points
#'

PredictorResponseBivar.singfit.MI <- function(fit, y = NULL, Z = NULL, X = NULL, z.pairs = NULL, method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(Z), verbose = TRUE, Z.MI,k=1,K=1, ...) {
  
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) y <- fit$y
    if (is.null(Z)) Z <- fit$Z
    if (is.null(X)) X <- fit$X
  }
  
  if (is.null(z.names)) {
    z.names <- colnames(Z.MI)
    if (is.null(z.names)) {
      z.names <- paste0("z", 1:ncol(Z))
    }
  }
  
  if (is.null(z.pairs)) {
    z.pairs <- expand.grid(z1 = 1:ncol(Z), z2 = 1:ncol(Z))
    z.pairs <- z.pairs[z.pairs$z1 < z.pairs$z2, ]
  }
  
  df <- dplyr::data_frame()
  for(i in 1:nrow(z.pairs)) {
    compute <- TRUE
    whichz1 <- z.pairs[i, 1] %>% unlist %>% unname
    whichz2 <- z.pairs[i, 2] %>% unlist %>% unname
    if(whichz1 == whichz2) compute <- FALSE
    z.name1 <- z.names[whichz1]
    z.name2 <- z.names[whichz2]
    names.pair <- c(z.name1, z.name2)
    if(nrow(df) > 0) { ## determine whether the current pair of variables has already been done
      completed.pairs <- df %>%
        dplyr::select_('variable1', 'variable2') %>%
        dplyr::distinct() %>%
        dplyr::transmute(z.pair = paste('variable1', 'variable2', sep = ":")) %>%
        unlist %>% unname
      if(paste(names.pair, collapse = ":") %in% completed.pairs | paste(rev(names.pair), collapse = ":") %in% completed.pairs) compute <- FALSE
    }
    if(compute) {
      if(verbose) message("MI fit ", k," out of ",K, ":  Pair ", i, " out of ", nrow(z.pairs))
      res <- PredictorResponseBivarPair.MI(fit = fit, y = y, Z = Z, X = X, whichz1 = whichz1, whichz2 = whichz2, method = method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, Z.MI=Z.MI, ...)
      df0 <- res
      df0$variable1 <- z.name1
      df0$variable2 <- z.name2
      df0 %<>%
        dplyr::select_(~variable1, ~variable2, ~z1, ~z2, ~est, ~se)
      df <- dplyr::bind_rows(df, df0)
    }
  }
  df$variable1 <- factor(df$variable1, levels = z.names)
  df$variable2 <- factor(df$variable2, levels = z.names)
  df
}


#' Plot cross-sections of the bivariate predictor-response function
#' 
#' Function to plot the \code{h} function of a particular variable at different levels (quantiles) of a second variable
#' 

 
PredictorResponseBivar.MI <- function(BKMRfits, z.pairs = NULL, method = "approx", ngrid = 50, q.fixed = 0.5, sel = NULL, min.plot.dist = 0.5, center = TRUE, z.names = colnames(BKMRfits[[1]]$Z), verbose = TRUE, ...) {
	
  start.time <- proc.time()["elapsed"]
  Z.MI <- Z.complete.MI(BKMRfits)
  l <- ncol(Z.MI)
  z.names <- colnames(Z.MI)

  K <- length(BKMRfits)  
  
  npairs <- length(z.pairs)
  if(is.null(z.pairs)) npairs <- factorial(l)/factorial(l-2)/factorial(2)
  est.matrix <- matrix(NA,nrow=ngrid*ngrid*npairs, ncol=K)
  for(k in 1:K){
    fit <- BKMRfits[[k]]
	  y <- fit$y
  	Z <- fit$Z
  	X <- fit$X	
  	
  	temp <- PredictorResponseBivar.singfit.MI(fit=fit, z.pairs = z.pairs, method = method, ngrid = ngrid, q.fixed = q.fixed, sel = sel, min.plot.dist = min.plot.dist, center = center, z.names = z.names, verbose = TRUE, Z.MI=Z.MI,k=k,K=K, ...)
  	
  	est.matrix[,k] <- temp %>% select(est) %>% unlist(use.names=FALSE)
  	
  	end.time <- proc.time()["elapsed"]
    message(paste("MI fit", k, "out of", K, "complete: ", round((end.time - start.time)/60, digit=2), "min run time" ))	
  }

  est.MI <- apply(est.matrix,1,mean)
  
  
  data.toreturn <- temp ### okay since the grid numbers and variable order are the same
  data.toreturn[,"est"] <- est.MI 

  data.toreturn
}








#################################################
######## other BKMR functions required ##########
#################################################

## these are directly from the BKMR source code
## (no edits made)


# makeKpart <- function(r, Z) {
# Kpart <- as.matrix(dist(sqrt(matrix(r, byrow=TRUE, nrow(Z), ncol(Z)))*Z))^2
# Kpart
# }
makeKpart <- function(r, Z1, Z2 = NULL) {
  Z1r <- sweep(Z1, 2, sqrt(r), "*")
  if (is.null(Z2)) {
    Z2r <- Z1r
  } else {
    Z2r <- sweep(Z2, 2, sqrt(r), "*")
  }
  Kpart <- fields::rdist(Z1r, Z2r)^2
  Kpart
}
makeVcomps <- function(r, lambda, Z, data.comps) {
  if (is.null(data.comps$knots)) {
    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    if (data.comps$nlambda == 2) {
      V <- V + lambda[2]*data.comps$crossTT
    }
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    logdetVinv <- -2*sum(log(diag(cholV)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1
    # Kpartall <- makeKpart(r, rbind(Z, data.comps$knots))
    # Kall <- exp(-Kpartall)
    # K0 <- Kall[1:n0, 1:n0 ,drop=FALSE]
    # K1 <- Kall[(n0+1):nall, (n0+1):nall ,drop=FALSE]
    # K10 <- Kall[(n0+1):nall, 1:n0 ,drop=FALSE]
    K1 <- exp(-makeKpart(r, data.comps$knots))
    K10 <- exp(-makeKpart(r, data.comps$knots, Z))
    Q <- K1 + diag(nugget, n1, n1)
    R <- Q + lambda[1]*tcrossprod(K10)
    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)
    Vinv <- diag(1, n0, n0) - lambda[1]*t(K10) %*% Rinv %*% K10
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}


ComputePostmeanHnew.approx <- function (fit, y = NULL, Z = NULL, X = NULL, Znew = NULL, sel = NULL) {
  if (inherits(fit, "bkmrfit")) {
    if (is.null(y)) 
      y <- fit$y
    if (is.null(Z)) 
      Z <- fit$Z
    if (is.null(X)) 
      X <- fit$X
  }
  if (!is.null(Znew)) {
    if (is.null(dim(Znew))) 
      Znew <- matrix(Znew, nrow = 1)
    if (class(Znew) == "data.frame") 
      Znew <- data.matrix(Znew)
  }
  if (is.null(dim(X))) 
    X <- matrix(X, ncol = 1)
  ests <- ExtractEsts(fit, sel = sel)
  sigsq.eps <- ests$sigsq.eps[, "mean"]
  r <- ests$r[, "mean"]
  beta <- ests$beta[, "mean"]
  lambda <- ests$lambda[, "mean"]
  if (fit$family == "gaussian") {
    ycont <- y
  }
  else if (fit$family == "binomial") {
    ycont <- ests$ystar[, "mean"]
  }
  Kpart <- makeKpart(r, Z)
  K <- exp(-Kpart)
  V <- diag(1, nrow(Z), nrow(Z)) + lambda[1] * K
  cholV <- chol(V)
  Vinv <- chol2inv(cholV)
  if (!is.null(Znew)) {
    n0 <- nrow(Z)
    n1 <- nrow(Znew)
    nall <- n0 + n1
    Kpartall <- makeKpart(r, rbind(Z, Znew))
    Kmat <- exp(-Kpartall)
    Kmat0 <- Kmat[1:n0, 1:n0, drop = FALSE]
    Kmat1 <- Kmat[(n0 + 1):nall, (n0 + 1):nall, drop = FALSE]
    Kmat10 <- Kmat[(n0 + 1):nall, 1:n0, drop = FALSE]
    lamK10Vinv <- lambda[1] * Kmat10 %*% Vinv
    postvar <- lambda[1] * sigsq.eps * (Kmat1 - lamK10Vinv %*% 
                                          t(Kmat10))
    postmean <- lamK10Vinv %*% (ycont - X %*% beta)
  }
  else {
    lamKVinv <- lambda[1] * K %*% Vinv
    postvar <- lambda[1] * sigsq.eps * (K - lamKVinv %*% 
                                          K)
    postmean <- lamKVinv %*% (ycont - X %*% beta)
  }
  ret <- list(postmean = drop(postmean), postvar = postvar)
  ret
}

