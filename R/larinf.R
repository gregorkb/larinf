#' Least angle regression
#'
#' Runs the LAR algorithm.
#'
#' @param X a design matrix.
#' @param y a response vector.
#' @param rescale_y a logical. If `TRUE` then `y` is scaled by one over root n.
#' @returns The `lar()` function returns an object of class \code{lar}. An object of class \code{lar} is a list containing:
#'
#' \describe{
#' \item{`a`}{A matrix containing the sequence of equiangular vectors.}
#' \item{`mu`}{The value of the terminal predictor.}
#' \item{`A`}{A vector containing the sequence of angles.}
#' \item{`C`}{A vector containing the step correlations.}
#' \item{`g`}{A vector containing the sequence of gamma values.}
#' \item{`b`}{A matrix containing the sequence of step coefficient vectors.}
#' \item{`cc`}{A matrix containing the correlations of each column of `X` with the residual across all LAR steps.}
#' \item{`s`}{A vector giving the sign which which the entering variable on each step was signed.}
#' \item{`ord`}{A vector giving the order in which the columns of `X` entered the active set.}
#' \item{`nnew`}{A vector giving the number of the columns of `X` entering on each step.}
#' \item{`varnames`}{The names of the columns of `X` if non-null.}
#' \item{`ehat`}{The residuals.}
#' \item{`rescale_y`}{The argument given for `rescale_y`.}
#' \item{`delta`}{The largest delta for which the separation conditions in Theorem 4.2 of Gregory and Nordman (2025+) are satisfied for this LAR path.}
#' }
#'
#' @details
#' The `lar()` function centers the response vector `y` and centers the columns of `X` and scales them to have
#' unit norm before executing the least angle regression algorithm.  The function can accommodate tied entrances.
#'
#' @references
#' Efron B., Hastie, T., Johnstone, I., and Tibshirani, R. (2004) Least angle regression.
#' *Annals of Statistics*, **32(2)**: 407-499. [doi:10.1214/009053604000000067](https://doi.org/10.1214/009053604000000067)
#'
#' Gregory, K. and Nordman, D. (2025+) Least angle regression inference.
#'
#' @author Karl Gregory
#'
#' @seealso [plot.lar()] for plotting the least angle regression path.
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' b <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate X and y
#' X <- matrix(rnorm(n*p),n,p) %*% chol_Sigma
#' mu <- drop(X %*% b)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e
#'
#' # compute population path and sample path
#' lar_mu <- lar(X,mu)
#' lar_y <- lar(X,y)
#'
#' # plot population path and sample path
#' plot(lar_mu)
#' plot(lar_y)
#' @export
lar <- function(X,y,rescale_y = T){

  tol <- 1e-12

  dm <- dim(X)
  n <- dm[1]
  p <- dm[2]
  K <- min(p,n-1)

  # center and normalize columns of X, center y
  X <- scale(X)/sqrt(n-1)
  y <- y - mean(y)
  if(rescale_y == T) y <- y / sqrt(n)

  # create some matrices to store output
  A <- numeric(K)
  a <- matrix(0,n,K)
  b <- matrix(0,p,K)
  cc <- matrix(0,p,K)
  C <- numeric(K)
  g <- numeric(K)
  s <- integer()
  ord <- integer()
  nnew <- integer(K)
  delta <- Inf # this is going to be my delta from Theorem 1.

  # initialize
  gram <- t(X) %*% X
  c1 <- drop(t(X) %*% y)
  mu <- numeric(n)
  bk <- numeric(p)
  active0 <- integer()
  Lk <- NULL
  for(k in 1:(K + 1)){

    ck <- c1 - drop(t(X) %*% mu)
    Ck <- max(abs(ck))

    # these lines seem really sloppy
    active <- which(abs(ck) > Ck - tol)
    new <- active[which(active %in% active0 == FALSE)] # maybe this is silly
    active <- c(active0,new) # in proper order

    if(Ck < tol * 100 | k == (K+1)) break
    # can sometimes have issues here; Ck is really zero but is not close enough to zero to be under the specified tolerance
    # for this reason I added the condition k == (K + 1)

    delta <- min(delta,Ck - abs(ck[-active])) # update delta

    for(i in new){

      Lk <- updateLk(Lk,Akk = gram[i,i],ak = gram[active0,i])
      active0 <- c(active0,i)

    }

    Gk <- chol2inv(t(Lk))
    s_active <- sign(ck[active])
    Ak <- 1/(sqrt( sum(drop(Gk %*% s_active) * s_active)))
    ak <- Ak * drop(X[,active,drop = FALSE] %*% drop(Gk %*% s_active))
    wk <- drop(t(X) %*% ak)

    # new way
    if(length(active) < p){

      rkj <- sign(ck - Ck/Ak * wk)
      gkj <- (Ck - ck*rkj) / (Ak - wk*rkj)

      sort_gk <- sort(gkj[-active])
      gk <- sort_gk[1]

      if((gk < (Ck / Ak - tol)) & (length(sort_gk) > 2)){
        # if LAR is going to terminate after this update (that is if this is step m), then
        # all the gkj[-active] will be equal to Ck/Ak, so the delta condition does not
        # apply in this step. Here we basically check if this is step m in the prototypical path.
        # Also, we check if there is only one remaining variable, in which case that
        # variable MUST be chosen in the last step, so the delta condition does not apply.

        delta <- min(delta,(sort_gk[2] - gk)*Ak)

      }

    } else if(length(active) == p){

      gk <- Ck / Ak

    }

    # update predictor
    mu <- mu + gk * ak
    bk[active] <- bk[active] + gk * Gk %*% wk[active]

    # store some output
    A[k] <- Ak
    a[,k] <- ak
    C[k] <- Ck
    g[k] <- gk
    b[,k] <- bk
    cc[,k] <- ck

    # keep track of which variables entered on each step
    ord <- c(ord,new)
    nnew[k] <- length(new)

    # store signs of entrance correlations
    s <- c(s,sign(ck[new]))

  }

  m <- k-1

  A <- A[1:m]
  a <- a[,1:m]
  C <- C[1:m]
  g <- g[1:m]
  b <- b[,1:m,drop = FALSE]
  nnew <- nnew[1:m]

  output <- list(a = a,
                 mu = mu,
                 A = A,
                 C = C,
                 g = g,
                 b = b,
                 cc = cc,
                 s = s,
                 ord = ord,
                 nnew = nnew,
                 varnames = colnames(X),
                 ehat = y - mu,
                 rescale_y = rescale_y,
                 delta = delta)

  class(output) <- "lar"

  return(output)

}


#' Compute the lar path based on the gram matrix gram = t(X) %*%X and the vector of marginal correlations cvec = t(X) %*% y. Faster than lar() if gram already computed.
#' @noRd
lar_gram <- function(gram,cvec){

  tol <- 1e-10
  p <- length(cvec)

  # create some matrices to store output
  A <- numeric(p)
  C <- numeric(p)
  g <- numeric(p)
  b <- matrix(0,p,p)
  s <- integer()
  ord <- integer()
  nnew <- integer(p)

  # initialize
  active0 <- integer()
  ord <- integer()
  bk <- numeric(p)
  c1 <- cvec
  mk <- numeric(p)
  Lk <- NULL
  for(k in 1:(p + 1)){

    ck <- c1 - mk
    Ck <- max(abs(ck))
    active <- which(Ck - abs(ck) < tol)
    new <- active[which(active %in% active0 == FALSE)]
    active <- c(active0,new) # in proper order

    if(Ck < tol) break

    for(i in new){

      Lk <- updateLk(Lk,Akk = gram[i,i],ak = gram[active0,i])
      active0 <- c(active0,i)

    }

    s_active <- sign(ck[active])
    fsLks <- forwardsolve(Lk,s_active)

    Ak <- 1/sqrt(sum(fsLks^2))
    wk <- Ak * drop(gram[,active,drop = FALSE] %*% forwardsolve(Lk,fsLks,transpose = TRUE))

    if(length(active) < p){

      rkj <- sign(ck - Ck/Ak * wk)
      gkj <- (Ck - ck*rkj) / (Ak - wk*rkj)

      sort_gk <- sort(gkj[-active])
      gk <- sort_gk[1]

    } else if(length(active) == p){

      gk <- Ck / Ak

    }

    mk <- mk + gk * wk
    bk[active] <- bk[active] + gk * chol2inv(t(Lk)) %*% wk[active]

    # store some output
    A[k] <- Ak
    C[k] <- Ck
    g[k] <- gk
    b[,k] <- bk

    # keep track of which variables entered on each step
    ord <- c(ord,new)
    nnew[k] <- length(new)

    # store signs of entrance correlations
    s <- c(s,sign(ck[new]))

  }

  m <- k - 1

  A <- A[1:m]
  C <- C[1:m]
  g <- g[1:m]
  b <- b[,1:m]
  nnew <- nnew[1:m]

  output <- list(A = A,
                 C = C,
                 g = g,
                 b = b,
                 s = s,
                 ord = ord,
                 nnew = nnew)

  class(output) <- "lar"

  return(output)

}
#' Least angle regression inference
#'
#' Computes bootstrap confidence intervals for least angle regression entrance correlations.
#'
#' @param X a design matrix.
#' @param y a response vector.
#' @param rescale_y a logical. If `TRUE` then `y` is scaled by one over root n.
#' @param B the number of Monte Carlo draws for approximating the bootstrap distribution.
#' @param alpha the significance level.
#' @param m_bar optional argument fixing the number of nonzero entrance correlations at which to threshold. By default the function estimates this with `mest()`.
#' @return The `larinf()` function returns an object of class `larinf`. An object of class `larinf` is a list containing:
#'
#' \describe{
#' \item{`lob`}{a matrix with each column giving the lower bound of the confidence interval for the step coefficients at the corresponding step.}
#' \item{`upb`}{a matrix with each column giving the upper bound of the confidence interval for the step coefficients at the corresponding step.}
#' \item{`b_bar`}{a matrix with each column giving the estimated step coefficients at a step. The number of columns is equal to the estimated number of steps in the population path, and the last column contains the coefficients corresponding to the orthogonal projection of the response on the span of the columns in the active set.}
#' \item{`loC`}{a vector containing the lower limits of the confidence intervals for the step correlations.}
#' \item{`upC`}{a vector containing the upper limits of the confidence intervals for the step correlations.}
#' \item{`sigma_hat`}{the estimate of the error standard deviation.}
#' \item{`lar_out`}{an object of class `lar` returned by a call to the `lar()` function.}
#' \item{`mest_out`}{a list containing information relevant to the estimation of the number of nonzero step correlations.}
#' \item{`actprob`}{a matrix giving the proportion of bootstrap samples on which each variable was active at each step.}
#' \item{`B`}{the number of bootstrap samples drawn.}
#' \item{`alpha`}{the significance level.}
#' }
#'
#' @references
#' Gregory, K. and Nordman, D. (2025+) Least angle regression inference.
#'
#' @details
#' The `larinf()` function centers the response vector `y` and centers and normalizes the columns of `X`.
#'
#' @author Karl Gregory
#'
#' @seealso [plot.larinf()] for plotting and [print.larinf()] for printing the results of least angle regression inference.
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' b <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate X and y
#' X <- matrix(rnorm(n*p),n,p) %*% chol_Sigma
#' mu <- drop(X %*% b)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e
#'
#' # perform least angle regression inference
#' larinf_out<- larinf(X,y)
#' plot(larinf_out)
#' larinf_out
larinf <- function(X,y,rescale_y = T,B = 500, alpha = 0.05, m_bar = NULL){

  dm <- dim(X)
  n <- dm[1]
  p <- dm[2]

  # center and normalize columns of X, center y
  X <- scale(X) / sqrt(n-1)
  y <- y - mean(y)
  if(rescale_y == T) y <- y / sqrt(n)

  lar_out <- lar(X,y,rescale_y = F)
  b_hat <- lar_out$b
  e_hat <- lar_out$ehat
  s_hat <- lar_out$s
  dA_hat <- c(1,sqrt(diff(1/lar_out$A^2)))
  C_hat <- lar_out$C

  sigma_hat <- sqrt(sum(e_hat^2)/(n - p)) * ifelse(rescale_y,sqrt(n),1)

  # estimate the number of nonzero entrance correlations
  mest_out <- mest(lar_out, m_bar = m_bar)
  m_bar <- mest_out$m_bar

  if(m_bar == 0){

    print("Number of non-zero step correlations estimated to be zero.")
    return(NULL)

  }

  # compute the gram matrix
  gram <- t(X) %*% X

  # define b_bar, mu_bar, and C_bar
  b_bar <- b_hat[,1:m_bar,drop = FALSE]
  Am_bar <- lar_out$ord[1:m_bar]
  b_bar[Am_bar,m_bar] <- solve(gram[Am_bar,Am_bar]) %*% t(X[,Am_bar]) %*% y
  mu_bar <- drop(X %*% b_bar[,m_bar])
  C_bar <- c(C_hat[1:m_bar],rep(0,p-m_bar))

  B_boot <- array(0,dim=c(p,m_bar,B))
  T_boot <- matrix(0,B,p)

  ord_boot <- matrix(0,B,p)
  resid_scl <- sqrt(n/(n-p))

  b <- 1
  while(b < B){

    # draw bootstrap residuals
    e_boot <- sample(e_hat,n,replace = TRUE) * resid_scl

    # generate bootstrap data with mu_bar
    y_boot <- mu_bar + e_boot - mean(e_boot)
    c_boot <- t(X) %*% y_boot

    # run lar on the bootstrap data
    lar_boot_out <- lar_gram(gram,c_boot)
    b_hat_boot <- lar_boot_out$b
    s_hat_boot <- lar_boot_out$s
    C_hat_boot <- lar_boot_out$C
    dA_hat_boot <- c(1,sqrt(diff(1/lar_boot_out$A^2)))

    if(length(lar_boot_out$A) < p){

      print("bootstrap lar path had fewer than p steps - discarding")
      next;

    }

    ord_boot[b,] <- lar_boot_out$ord

    # obtain bootstrap version of sigma_hat
    mu_hat_boot <- drop(X %*% lar_boot_out$b[,p])
    e_hat_boot <- y_boot - mu_hat_boot
    sigma_hat_boot <- sqrt(sum(e_hat_boot^2) / (n - p)) * ifelse(rescale_y,sqrt(n),1)

    # get b_bar_boot
    b_bar_boot <- b_hat_boot[,1:m_bar,drop = FALSE]
    Am_bar_boot <- lar_boot_out$ord[1:m_bar]
    b_bar_boot[Am_bar_boot,m_bar] <- solve(gram[Am_bar_boot,Am_bar_boot]) %*% t(X[,Am_bar_boot]) %*% y_boot

    B_boot[,,b] <- (b_bar_boot - b_bar) / sigma_hat_boot
    T_boot[b,] <- s_hat_boot * dA_hat_boot * (C_hat_boot - C_bar) / sigma_hat_boot

    b <- b + 1

  }

  # construct bootstrap confidence intervals for the bk
  Bloalpha2 <- apply(B_boot,c(1,2),function(x) quantile(x,probs=alpha/2))
  Bupalpha2 <- apply(B_boot,c(1,2),function(x) quantile(x,probs=1-alpha/2))

  lob <- b_bar - Bupalpha2 * sigma_hat
  upb <- b_bar - Bloalpha2 * sigma_hat

  # construct bootstrap CIs for the step correlations Ck
  li <- ceiling((alpha/2)*B)
  ui <- ceiling((1 - alpha/2)*B)
  Talpha2 <- apply(T_boot,2,sort)[c(li,ui),,drop = FALSE]

  loC <- numeric(p)
  upC <- numeric(p)
  for(k in 1:p){

    if(s_hat[k] == 1){

      loC[k] <- max(0,C_hat[k] - Talpha2[2,k] * sigma_hat / dA_hat[k] * s_hat[k])
      upC[k] <- C_hat[k] - Talpha2[1,k] * sigma_hat / dA_hat[k] * s_hat[k]

    } else if(s_hat[k] == -1){

      loC[k] <- max(0,C_hat[k] - Talpha2[1,k] * sigma_hat / dA_hat[k] * s_hat[k])
      upC[k] <- C_hat[k] - Talpha2[2,k] * sigma_hat / dA_hat[k] * s_hat[k]

    }

  }

  # keep track of the bootstrap probability of active set membership for each variable
  actprob <- matrix(0,p,p)
  for(j in 1:p)
    for(k in 1:p){

      actprob[j,k] <- mean(apply(ord_boot[,1:k,drop= FALSE] == j,1,any))

    }

  rownames(actprob) <- lar_out$varnames

  output <- list(lob = lob,
                 upb = upb,
                 b_bar = b_bar,
                 loC = loC,
                 upC = upC,
                 sigma_hat = sigma_hat,
                 lar_out = lar_out,
                 mest_out = mest_out,
                 actprob = actprob,
                 B = B,
                 alpha = alpha)

  class(output) <- "larinf"

  return(output)

}


#' Estimate the number of steps in the prototypical LAR path
#' @noRd
mest <- function(lar_out,m_bar = NULL){

  n <- length(lar_out$mu)
  p <- nrow(lar_out$b)
  rescale_y <- lar_out$rescale_y

  sigma_hat <- sqrt(sum(lar_out$ehat^2)/(n - p)) * ifelse(rescale_y,sqrt(n),1)

  A <- c(Inf,lar_out$A)
  W <- diff(1/A^2) * lar_out$C^2 / sigma_hat^2 * ifelse(rescale_y,n,1)

  SW <- sum(W) - c(0,cumsum(W[1:(p-1)]))

  if(is.null(m_bar)){

    if(SW[1] < qchisq(1-1/n,p)){

      m_bar <- 0

    } else {

      cond <- numeric(p)
      for(k in 1:p){

        cond[k] <- all(SW[1:k] > qchisq(1 - 1/n, p - c(1:k) + 1))

      }

      m_bar <- max(which(cond == T))

    }

  }

  if(m_bar == 0){

      mu_bar <- numeric(n)

    } else if((m_bar > 0) & (m_bar < p)){

      # obtain mu_bar
      C <- c(0,lar_out$C)
      a <- cbind(0,lar_out$a)
      mu_bar <- numeric(n)
      for(k in 1:m_bar){

        dk <- a[,k+1]/A[k+1] - a[,k]/A[k]
        mu_bar <- mu_bar + dk*C[k+1]

      }

    } else if(m_bar == p){

      mu_bar <- lar_out$mu

    }

  output <- list(m_bar = m_bar,
                 mu_bar = mu_bar,
                 SW = SW)

  class(output) <- "mest"

  return(output)

}


#' Plot method for class `lar`
#'
#' Plot the steps of the LAR path based on the output of the `lar` function.
#'
#' @param x an object of class `lar`, usually the result of a call to `lar`. The LAR path may not have tied entrances.
#' @param omaadd a vector of length 4 giving amounts to add to each of the four outer margins (which may be necessary to include longer variable names).
#' @param madd a vector of length 4 giving amounts to add to each of the four inner margins (which may be necessary to include longer variable names).
#' @param m the number of steps out to which to plot the LAR path
#' @param mlabel the number of steps out to which to print labels for the variables
#' @param ... additional arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' b <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate X and y
#' X <- matrix(rnorm(n*p),n,p) %*% chol_Sigma
#' mu <- drop(X %*% b)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e
#'
#' # compute population path and sample path
#' lar_mu <- lar(X,mu)
#' lar_y <- lar(X,y)
#'
#' # plot population path and sample path
#' plot(lar_mu)
#' plot(lar_y)
#'
#' @export
plot.lar <- function(x, omaadd = c(0,0,0,0), madd = c(0,0,0,0), m = NULL, mlabel = NULL,...){


  if(any(x$nnew) > 1) stop("LAR path has tied entrances: This plot function is not yet set up to depict paths with tied entrances.")

  C <- x$C
  b <- x$b
  ord <- x$ord
  p <- nrow(b)
  cc <- x$cc

  op <- par(mfrow=c(1,2),
            mar= c(4.1,4.1,1.1,2.1) + madd,
            oma = c(0,0,0,0) + omaadd)

  lab <- x$varnames
  if(is.null(lab)) lab <- paste(1:p)
  if(is.null(m)) m <- length(ord)
  if(is.null(mlabel)) mlabel <- m

  plot(NA,
       xlim = c(1,m),
       ylim = range(C,0),
       xlab = "Step",
       ylab = "Step correlations",
       bty = "l",
       xaxt = "n")

  axis(1,at = 1:m,labels = paste(1:m))

  for(j in ord){

    kj <- which(ord == j)

    lines(abs(cc[j,1:m]),col = ifelse(kj <= m,j,"lightgray"),lty = j)

  }

  for(k in 1:mlabel){

    jk <- ord[k]

    text(x = k,
         y = abs(cc[jk,k]),
         pos = 4,
         labels = lab[jk],
         xpd = NA,
         cex = 0.8)

  }

  plot(NA,
       xlim = c(0,m),
       ylim = range(b),
       xlab = "Step",
       ylab = "Step coefficients",
       bty = "l",
       xaxt = "n")

  axis(1,at = 0:m,labels = paste(0:m))

  for(j in ord){

    kj <- which(ord == j)

    if(kj <= m){

      lines(c(0,b[j,1:m])~c(0:m), lty = j, col = j)

      if(kj <= mlabel){

        text(x = m,
             y = b[j,m],
             pos = 4,
             labels = lab[j],
             xpd = NA,
             cex = 0.8)

      }

    }

  }

  par(op)
  layout(1)

}

#' Plot method for class `larinf`
#'
#' Plot of the inferred LAR path out to the estimated number of steps.
#'
#' @param x an object of class `larinf`, usually the result of a call to `larinf`.
#' @param omaadd a vector of length 4 giving amounts to add to each of the four outer margins (which may be necessary to include longer variable names).
#' @param madd a vector of length 4 giving amounts to add to each of the four inner margins (which may be necessary to include longer variable names).
#' @param text an optional argument giving text to be printed as a title on the plot
#' @param m an optional argument giving the number of steps out to which to plot the results of LAR inference. The default is the estimated number of steps from the `mest()` function.
#' @param which if equal to 1, the inferred LAR path is plotted; if equal to 2 diagnostic plots are plotted, i.e. the tail sum statistics and the bootstrap probability of active set membership at each step for each variable. Equal to 1 by default.
#' @param ... additional arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' b <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate X and y
#' X <- matrix(rnorm(n*p),n,p) %*% chol_Sigma
#' mu <- drop(X %*% b)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e
#'
#' # perform least angle regression inference
#' larinf_out<- larinf(X,y)
#' plot(larinf_out)
#'
#' @export
plot.larinf <- function(x,omaadd=c(0,0,0,0),madd=c(0,0,0,0),text = NULL,m=NULL,which = 1,...){

  if(which == 1){

    # collect values from x
    lab <- x$lar_out$varnames
    th <- x$lar_out$th
    C <- x$lar_out$C
    ord <- x$lar_out$ord
    b_bar <- as.matrix(x$b_bar)
    b <- x$lar_out$b
    p <- nrow(b)
    lob <- cbind(0,x$lob)
    upb <- cbind(0,x$upb)
    m_bar <- x$mest_out$m_bar
    C_bar <- c(C[1:m_bar],0)
    cc <- x$lar_out$cc
    loC <- x$loC
    upC <- x$upC

    if(is.null(lab)) lab <- paste(1:p)
    if(is.null(m)) m <- m_bar

    # set up graphical parameters
    op <- par(mfrow=c(1,2),
              mar=c(4.1,4.1,1.1,2.1) + madd,
              oma = c(0,0,0,0) + omaadd)

    # plot Cks
    plot(NA,
         xlim = c(1,m),
         ylim = range(c(0,loC,upC)),
         xlab = "Step",
         ylab = "Step correlations",
         bty ="l",
         xaxt = "n")

    axis(1,at = 1:m,labels = paste(1:m))

    for(k in 1:m){

      lines(y = c(loC[k],upC[k]),
            x = rep(k,2),
            col = "darkgray")

    }

    xpoly <- c(1:m,m:1)
    ypoly <- c(loC[1:m],upC[m:1])
    polygon(x = xpoly,
            y = ypoly,
            col = rgb(0,0,0,.1),
            border = NA)

    for(k in 1:m){

      jk <- ord[k]

      lines(abs(cc[jk,1:m]),col = jk,lty = jk)

      text(x = k,
           y = abs(cc[jk,k]),
           pos = 4,
           labels = lab[jk],
           xpd = NA,
           cex = 0.8)

    }

    # plot bks
    plot(NA,
         xlim = c(0,m),
         ylim = range(lob[,1:(m+1)],upb[,1:(m+1)]),
         xlab = "Step",
         ylab = "Step coefficients",
         bty = "l",
         xaxt = "n")

    axis(1,at = 0:m,labels = paste(0:m))

    k <- 1
    for(j in ord[1:m]){

      xpoly <- c((k-1):m,m:(k-1))
      ypoly <- c(0,lob[j,(k+1):(m+1)],upb[j,(m+1):(k+1)],0)
      polygon(x = xpoly,
              y = ypoly,
              col = rgb(0,0,0,.1),
              border = NA)

      for(l in k:m){

        lines(x = rep(l-1,2),
              y = c(lob[j,l],upb[j,l]),
              col = "darkgray")

      }

      lines(c(0,b_bar[j,1:m])~c(0:m), lty = j, col = j)

      text(x = m,
           y = b_bar[j,m],
           pos = 4,
           labels = lab[j],
           xpd = NA,
           cex = 0.8)

      k <- k + 1

    }

    } else if(which == 2){

    mest_out <- x$mest_out

    SW <- mest_out$SW
    m_bar <- mest_out$m_bar
    n <- length(mest_out$mu_bar)
    p <- length(SW)

    op <- par(mar = c(4.1,4.1,1.1,1.1), mfrow = c(1,2))

    plot(y = SW[p:1],
         x = p:1,
         xlim = c(1,p),
         ylab = "Tail sum statistics",
         xlab = "Step",
         bty = "l",
         pch = c(rep(1,p - m_bar),rep(19,m_bar)))

    lines(y=qchisq(1-1/n,df=1:p),x=p:1, lty = 3)


    actprob <- x$actprob
    plot(NA,
         xlim = c(0,p),
         ylim = c(0,1),
         xlab = "Step",
         ylab = "Bootstrap active set membership",
         bty = "l")

    for(j in 1:p){

      lines(c(0,actprob[j,])~c(0:p), col = j, lty = j)

    }

  }

  # add optional title text
  if(!is.null(text)){

    mtext(side = 3, line = -1, outer = TRUE,text = text,cex = .8)

  }

  par(op)
  layout(1)



}


#' Update the Cholesky decomposition
#' @noRd
updateLk <- function(Lk,Akk,ak,eps = .Machine$double.eps){

  if(is.null(Lk)){

    Lk <- sqrt(Akk)

  } else {

    lk <- forwardsolve(Lk,ak)
    rad <- Akk - sum(lk^2)

    if(rad <= eps){

      Lkk <- eps

    } else {

      Lkk <- sqrt(rad)

    }

    Lk <- rbind(cbind(Lk,0),cbind(t(lk),Lkk))

  }

  return(Lk)

}


#' Print method for class `larinf`
#'
#' Print tables of output for least angle regression inference.
#'
#' @param x an object of class `larinf`, usually the result of a call to `larinf`.
#' @param ... further arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' b <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate X and y
#' X <- matrix(rnorm(n*p),n,p) %*% chol_Sigma
#' mu <- drop(X %*% b)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e
#'
#' # perform least angle regression inference
#' larinf_out <- larinf(X,y)
#' print(larinf_out)
#'
#' @export
print.larinf <- function(x,...){

  ord <- x$lar_out$ord
  m_bar <- x$mest_out$m_bar
  entered <- ord[1:m_bar]
  C <- x$lar_out$C
  b_bar <- x$b_bar[,m_bar]
  SW <- x$mest_out$SW
  p <- length(SW)
  n <- length(x$lar_out$mu)
  ch <- qchisq(1 - 1/n, p - c(1:p) + 1)

  cat(paste("Estimated number of steps in population path: ",m_bar,"\n\n",sep=""))

  cat("Step correlations in order of entrance:\n\n")

  tabC <- round(cbind(SW,ch,C,x$loC,x$upC),3)
  colnames(tabC) <- c("Sk","Chi^2","Ck",paste(x$alpha/2*100,"%",sep=""),paste((1-x$alpha/2)*100,"%",sep=""))

  tabb <- round(cbind(b_bar[entered],x$lob[entered,m_bar],x$upb[entered,m_bar]),3)
  colnames(tabb) <- c("b",paste(x$alpha/2*100,"%",sep=""),paste((1-x$alpha/2)*100,"%",sep=""))

  if(!is.null(x$lar_out$varnames)){

    rownames(tabC) <- x$lar_out$varnames[ord]
    rownames(tabb) <- x$lar_out$varnames[entered]

  } else {

    rownames(tabC) <- paste("X",ord,sep= "")
    rownames(tabb) <- paste("X",entered,sep= "")

  }

  print(tabC)

  cat("\n\n")

  cat(paste("Step coefficients at step ",m_bar,":\n\n",sep=""))

  print(tabb)

}



