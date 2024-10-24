#' Least angle regression
#'
#' Runs the least angle regression algorithm.
#'
#' @param X a design matrix (typically with columns centered to have mean zero and scaled to have root-n norm).
#' @param y a response vector (typically centered to have mean zero).
#' @returns `lar` returns an object of class \code{lar}. An object of class \code{lar} is a list containing:
#'
#' \describe{
#' \item{`a`}{A matrix containing the sequence of equiangular vectors.}
#' \item{`mu`}{The value of the terminal predictor.}
#' \item{`A`}{A vector containing the sequence of angles.}
#' \item{`C`}{A vector containing the absolute entrance correlations.}
#' \item{`g`}{A vector containing the sequence of gamma values.}
#' \item{`b`}{A matrix containing the sequence of regression coefficient vectors.}
#' \item{`th`}{A vector containing the entrance correlations for the columns of `X`.}
#' \item{`ord`}{A vector giving the order in which the columns of `X` entered the active set.}
#' \item{`nnew`}{A vector giving the number of the columns of `X` entering on each step.}
#' \item{`varnames`}{The names of the columns of `X` if non-null.}
#' }
#'
#' @details
#' The function can accommodate tied entrances.
#'
#' @references
#' Efron B., Hastie, T., Johnstone, I., and Tibshirani, R. (2004) Least angle regression.
#' *Annals of Statistics*, **32(2)**: 407-499. [doi:10.1214/009053604000000067](https://doi.org/10.1214/009053604000000067)
#'
#' Gregory, K. and Nordman, D. (2025+) Least angle regression inference. *Work in progress*
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
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # compute noiseless path and sample path
#' lar_mu <- lar(X,mu)
#' lar_y <- lar(X,y)
#'
#' # plot noiseless path and overlay sample path
#' plot(lar_mu)
#' plot(lar_y,add = TRUE, lty = 3)
#'
#' @export
lar <- function(X,y){

  tol <- 1e-12

  dm <- dim(X)
  n <- dm[1]
  p <- dm[2]
  K <- min(p,n-1)

  # create some matrices to store output
  A <- numeric(K)
  a <- matrix(0,n,K)
  b <- matrix(0,p,K)
  C <- numeric(K)
  g <- numeric(K)
  s <- integer()
  ord <- integer()
  nnew <- integer(K)

  # initialize
  gram <- t(X) %*% X
  c1 <- drop( (1/n) * t(X) %*% y )
  mu <- numeric(n)
  bk <- numeric(p)
  active0 <- integer()
  Lk <- NULL
  for(k in 1:(K + 1)){

    ck <- c1 - drop((1/n) * t(X) %*% mu)
    Ck <- max(abs(ck))

    # these lines seem really sloppy
    active <- which(abs(ck) > Ck - tol)
    new <- active[which(active %in% active0 == FALSE)] # maybe this is silly

    active <- c(active0,new) # in proper order

    if(Ck < tol * 100 | k == (K+1)) break
    # can sometimes have issues here; Ck is really zero but is not close enough to zero to be under the specified tolerance
    # for this reason I added the condition k == (K + 1)

    for(i in new){

      Lk <- updateLk(Lk,Akk = gram[i,i],ak = gram[active0,i])
      active0 <- c(active0,i)

    }

    Gk <- chol2inv(t(Lk))
    s_active <- sign(ck[active])
    Ak <- 1/(sqrt( n * sum(drop(Gk %*% s_active) * s_active)))
    ak <- n * Ak * drop(X[,active,drop = FALSE] %*% drop(Gk %*% s_active))
    wk <- drop(t(X) %*% ak) / n

    # update gamma
    r <- ck[-active]
    d <- wk[-active]

    pl <- r / Ck < d / Ak + tol
    mi <- r / Ck > d / Ak - tol

    if(!any(pl,mi)){

      gk <- Ck / Ak

    } else {

      g1 <- (Ck + r) / (Ak + d)
      g2 <- (Ck - r) / (Ak - d)

      gk <- min(g1[pl],g2[mi])

    }

    mu <- mu + gk * ak
    bk[active] <- bk[active] + gk * Gk %*% wk[active] * n

    # store some output
    A[k] <- Ak
    a[,k] <- ak
    C[k] <- Ck
    g[k] <- gk
    b[,k] <- bk

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
  nnew <- nnew[1:m]

  th <- numeric(p)
  if(m > 0){

    j <- 1
    for(k in 1:m){

      ind <- j:(j + nnew[k] - 1)
      th[ord[ind]] <- s[ind] * C[k]
      j <- max(ind) + 1

    }

  }

  output <- list(a = a,
                 mu = mu,
                 A = A,
                 C = C,
                 g = g,
                 b = b,
                 th = th,
                 ord = ord,
                 nnew = nnew,
                 varnames = colnames(X))

  class(output) <- "lar"

  return(output)

}


# Compute the lar path based on the covariance matrix
# Sigma and the vector of marginal correlations cvec.
# This is a little faster than the lar() function, so it is
# used in the bootstrap loop of larinf().
# It is not exported to the user.
lar_asymp <- function(Sigma,cvec){

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

      Lk <- updateLk(Lk,Akk = Sigma[i,i],ak = Sigma[active0,i])
      active0 <- c(active0,i)

    }

    s_active <- sign(ck[active])
    fsLks <- forwardsolve(Lk,s_active)

    Ak <- 1/sqrt(sum(fsLks^2))
    wk <- Ak * drop(Sigma[,active,drop = FALSE] %*% forwardsolve(Lk,fsLks,transpose = TRUE))

    # update gamma
    r <- ck[-active]
    d <- wk[-active]

    pl <- r / Ck < d / Ak - tol
    mi <- r / Ck > d / Ak + tol

    if(!any(pl,mi)){

      gk <- Ck / Ak

    } else {

      g1 <- (Ck + r) / (Ak + d)
      g2 <- (Ck - r) / (Ak - d)

      gk <- min(g1[pl],g2[mi])

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

  th <- numeric(p)
  j <- 1
  for(k in 1:m){

    ind <- j:(j + nnew[k] - 1)
    th[ord[ind]] <- s[ind] * C[k]
    j <- max(ind) + 1

  }

  output <- list(A = A,
                 C = C,
                 g = g,
                 b = b,
                 th = th,
                 ord = ord,
                 nnew = nnew)

  class(output) <- "lar"

  return(output)

}

#' Least angle regression inference
#'
#' @param X a design matrix with columns centered to have mean zero.
#' @param y a response vector centered to have mean zero.
#' @param B the number of Monte Carlo draws for approximating the bootstrap distribution.
#' @param alpha the significance level.
#' @param thresh an argument passed to the `mest` function.
#' @param ncp an argument passed to the `mest` function.
#' @param m_bar optional argument fixing the number of nonzero entrance correlations at which to threshold.
#' @return `larinf` returns an object of class `larinf`. An object of class `larinf` is a list containing:
#'
#' \describe{
#'  \item{`lo`}{The lower bounds of the entrance correlation confidence intervals.}
#'  \item{`up`}{The upper bounds of the entrance correlation confidence intervals.}
#'  \item{`rej`}{An indicator of whether the entrance correlation confidence intervals exclude zero.}
#'  \item{`b_lo`}{The lower bounds of the regression coefficient confidence intervals.}
#'  \item{`b_up`}{The upper bounds of the regression coefficient confidence intervals.}
#'  \item{`b_rej`}{An indicator of whether the regression coefficient confidence intervals exclude zero.}
#'  \item{`sigma_hat`}{The estimator of the error term standard deviation.}
#'  \item{`lar_out`}{The object returned by `lar(X,y)`.}
#'  \item{`mest_out`}{The object returned by `mest(X,y,thresh,ncp,m_bar)`.}
#'  \item{`B`}{The number of Monte Carlo draws used to approximate the bootstrap distribution.}
#'  \item{`alpha`}{The significance level.}
#' }
#'
#' @references
#' Gregory, K. and Nordman, D. (2025+) Least angle regression inference. *In progress*
#'
#' @author Karl Gregory
#'
#' @seealso [plot.larinf()] for plotting and [print.larinf()] for printing the results of least angle regression inference and [mest()] for selecting the number of nonzero entrance correlations at which to threshold for the bootstrap.
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate some data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # perform least angle regression inference
#' larinf_out <- larinf(X,y)
#' plot(larinf_out)
#' print(larinf_out)
larinf <- function(X,y,B = 500, alpha = 0.05, thresh = 0.80, ncp = 0.1, m_bar = NULL){

  # check dimensions of X
  n <- nrow(X)
  p <- ncol(X)

  if(p >= (n-1)) stop("Design matrix X has p >= n-1")

  # compute the lar path and store the full output
  lar_out <- lar(X,y)

  # compute sample lar path
  th_hat <- lar_out$th

  if(any(th_hat == 0)) stop("Sample path has entrance correlations equal to zero")

  mu_hat <- lar_out$mu
  e_hat <- y - mu_hat
  sigma_hat <- sqrt(sum(e_hat^2)/(n - p))

  mest_out <- mest(X, y, thresh = thresh, ncp = ncp, m_bar = m_bar)

  m_bar <- mest_out$m_bar
  th_bar <- mest_out$th_bar
  mu_bar <- mest_out$mu_bar

  gram <- t(X) %*% X / n

  R_boot <- matrix(0,B,p)
  resid_scl <- sqrt(n/(n-p))

  b <- 1
  while(b < B){

    # draw bootstrap residuals
    e_boot <- sample(e_hat,n,replace = TRUE) * resid_scl

    # generate bootstrap data with mu_bar
    y_boot <- mu_bar + e_boot - mean(e_boot)
    c_boot <- t(X) %*% y_boot / n

    # run lar on the bootstrap data
    lar_boot_out <- lar_asymp(gram,c_boot)
    th_hat_boot <- lar_boot_out$th

    if(length(lar_boot_out$A) < p){

      print("bootstrap lar path had fewer than p steps - discarding")
      next;

    }

    # obtain bootstrap version of sigma_hat
    mu_hat_boot <- drop(X %*% lar_boot_out$b[,p])
    e_hat_boot <- y_boot - mu_hat_boot
    sigma_hat_boot <- sqrt(sum(e_hat_boot^2) / (n - p))

    R_boot[b,] <- (th_hat_boot - th_bar) / sigma_hat_boot

    b <- b + 1

  }

  # construct bootstrap confidence intervals
  li <- ceiling((alpha/2)*B)
  ui <- ceiling((1 - alpha/2)*B)
  Ralpha2 <- apply(R_boot,2,sort)[c(li,ui),]

  lo <- th_hat - Ralpha2[2,] * sigma_hat
  up <- th_hat - Ralpha2[1,] * sigma_hat
  rej <- sign(up) == sign(lo)

  # classical inference on slope coefficients
  Om <- solve(gram)
  b_hat <- as.numeric(Om %*% t(X) %*% y / n)
  se_b_hat <- sigma_hat * sqrt(diag(Om)) / sqrt(n)

  t_val <- qt(1-alpha/2,df = n-p)
  b_lo <- b_hat - t_val * se_b_hat
  b_up <- b_hat + t_val * se_b_hat
  b_rej <- sign(b_up) == sign(b_lo)

  output <- list(lo = lo,
                 up = up,
                 rej = rej,
                 b_lo = b_lo,
                 b_up = b_up,
                 b_rej = b_rej,
                 sigma_hat = sigma_hat,
                 lar_out = lar_out,
                 mest_out = mest_out,
                 B = B,
                 alpha = alpha)

  class(output) <- "larinf"

  return(output)

}


#' Estimate the number of nonzero entrance correlations
#'
#' Estimate the number of nonzero entrance correlations in the least angle regression path.
#'
#' @param X a design matrix with columns centered to have mean zero
#' @param y a response vector centered to have mean zero
#' @param thresh rejection threshold
#' @param ncp explain later
#' @param m_bar a user-specified number of nonzero entrance correlations.
#'
#' @return `mest` returns an object of class `mest`. An object of class `mest` is a list containing:
#'
#' \describe{
#'  \item{`m_bar`}{The selected number of nonzero entrance correlations.}
#'  \item{`th_bar`}{A vector containing the thresholded entrance correlations.}
#'  \item{`Wk`}{A vector containing the Wk values.}
#'  \item{`thresh`}{The value of the `thresh` argument.}
#'  \item{`ncp`}{The value of the `ncp` argument.}
#'  \item{`Fval`}{The quantile of the F distribution to which the Wk are compared.}
#' }
#'
#' @references
#' Gregory, K. and Nordman, D. (2025+) Least angle regression inference. *In progress*
#'
#' @author Karl Gregory
#'
#' @seealso [plot.mest()] for producing plots related to the selection of the number of nonzero entrance correlations.
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # select the number of nonzero entrance correlations at which to threshold
#' mest_out <- mest(X,y)
#' plot(mest_out)
#' @export
mest <- function(X,y,thresh = 0.80,ncp = 0.1,m_bar = NULL){

  lar_out <- lar(X,y)
  n <- nrow(X)
  p <- ncol(X)

  sigma_hat <- sqrt(sum((y - lar_out$mu)^2)/(n-p))
  A <- c(Inf,lar_out$A)
  Wk <- n * diff(1/A^2) * lar_out$C^2 / sigma_hat^2
  Fval <- qf(thresh,df1=1,df2=n-p,ncp = ncp*log(n))

  # plot(Wk)
  if(thresh == 0){

    mu_bar <- lar_out$mu
    th_bar <- lar_out$th
    m_bar <- length(lar_out$ord)

  } else {

    if(is.null(m_bar)){

      # m_bar <- min(which((Wk > Fval) == FALSE)) - 1
      m_bar <- sum(Wk > Fval)

    }

  if(m_bar == 0){

      mu_bar <- numeric(n)
      th_bar <- numeric(p)

    } else {

      # obtain mu_bar
      C <- c(0,lar_out$C)
      a <- cbind(0,lar_out$a)
      mu_bar <- numeric(n)
      for(k in 1:m_bar){

        dk <- a[,k+1]/A[k+1] - a[,k]/A[k]
        mu_bar <- mu_bar + dk*C[k+1]

      }

      # obtain th_bar
      th_bar <- numeric(p)
      first_m_bar <- lar_out$ord[1:m_bar]
      th_bar[first_m_bar] <- lar_out$th[first_m_bar]

    }

  }

  output <- list(m_bar = m_bar,
                 mu_bar = mu_bar,
                 th_bar = th_bar,
                 Wk = Wk,
                 thresh = thresh,
                 ncp = ncp,
                 Fval = Fval)

  class(output) <- "mest"

  return(output)

}



#' Plot method for class `lar`
#'
#' Make a cobweb plot of the least angle regression path based on the output of the `lar` function.
#'
#' @param x an object of class `lar`, usually the result of a call to `lar`.
#' @param main the main title of the plot.
#' @param add whether the cobweb plot should be added to the existing plot.
#' @param lty the line type.
#' @param xlim the limits of the horizontal axis if `add == FALSE`.
#' @param ylim the limits of the vertical axis if `add == FALSE`.
#' @param which_lab optional vector giving the indices of variables to be labeled.
#' @param ... additional arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # compute noiseless path and sample path
#' lar_mu <- lar(X,mu)
#' lar_y <- lar(X,y)
#'
#' # plot noiseless path and overlay sample path
#' plot(lar_mu)
#' plot(lar_y,add = TRUE, lty = 3)
#'
#' @export
plot.lar <- function(x,main = "", add = FALSE, lty = 1, xlim = NULL, ylim = NULL, which_lab = NULL,...){

  op <- par(mar=c(4.1,4.1,1.1,1.1))

  th <- x$th
  C <- x$C
  ord <- x$ord
  b <- x$b

  K <- length(C)
  p <- length(th)
  s <- sign(th)

  if(is.null(xlim)) xlim <- range(th)
  if(is.null(ylim)) ylim <- range(b)

  if(add == FALSE){

    plot(NA,
         xlim = xlim,
         ylim = ylim,
         xlab = "Entrance correlations",
         ylab = "Regression coefficients",
         main = main)

    abline(v = 0)
    abline(h = 0)

  }

  hundy <- (grconvertY(1,from="nfc", to ="user") - grconvertY(0,from="nfc", to ="user"))/100

  if(!is.null(which_lab)){

    for(j in which_lab){

      text(labels = paste(x$varnames[j]),
           x = th[j],
           y = ifelse(s[j] > 0,-hundy,hundy),
           pos = ifelse(s[j] > 0,2,4),
           offset = 0,
           cex = .8,
           xpd = NA,
           srt = 90)

    }

    for(j in 1:p){


      lines(x = c(C[1:K]*s[j],0), # need to add zero
          y = c(0,b[j,1:K]),
          col = ifelse(j %in% which_lab,j,"gray"),
          lty = lty)

    }

  } else {

    for(j in 1:p){

      lines(x = c(C[1:K]*s[j],0), # need to add zero
            y = c(0,b[j,1:K]),
            col = j,
            lty = lty)

    }
  }

  abline(h = 0)
  abline(v = 0)

  par(op)

}

#' Plot method for class `mest`
#'
#' Make plots related to the selection of the number of nonzero entrance correlations at which to threshold.
#'
#' @param x an object of class `mest`, usually the result of a call to `mest`.
#' @param ... further arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # select the number of nonzero entrance correlations at which to threshold
#' mest_out <- mest(X,y)
#' plot(mest_out)
#'
#' @export
plot.mest <- function(x,...){

  op <- par(mar = c(4.1,4.1,1.1,1.1))
  layout( matrix(c(1,2),nrow = 1, byrow = TRUE))

    plot(x$Wk,
         ylab = "Wk",
         xlab = "",
         pch = c(rep(19,x$m_bar),rep(1,length(x$th_bar) - x$m_bar)))
    abline(h = x$Fval,lty = 3)

    plot(y = cumsum(x$Wk[length(x$Wk):1]),
         x = length(x$Wk):1,
         xlim = c(length(x$Wk),1),
         ylab = "Cumulative sum of Wk",
         xlab = "",
         pch = c(rep(1,length(x$th_bar) - x$m_bar),rep(19,x$m_bar)))

    mtext(outer = TRUE, side = 1, line = -1.5 , text = "Least angle regression step")

  par(op)
  layout(1)

}


#' Plot method for class `larinf`
#'
#' Make a cobweb plot of the least angle regression path and display confidence intervals for the entrance correlations as well as for the regression coefficients.
#'
#' @param x an object of class `larinf`, usuall the result of a call to `larinf`.
#' @param omaadd a vector of length 4 giving amounts to add to each of the four outer margins (which may be necessary to include longer variable names).
#' @param ... additional arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate some data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # perform least angle regression inference
#' larinf_out <- larinf(X,y)
#' plot(larinf_out)
#'
#' @export
plot.larinf <- function(x,omaadd=c(0,0,0,0),...){

  p <- nrow(x$lar_out$b)

  lab <- x$lar_out$varnames
  if(is.null(lab)) lab <- paste(1:p)

  op <- par(mar=c(0,0,0,0), oma = c(4.1,5.1,1.1,2.1) + omaadd)
  layout( matrix(c(3,4,1,2),nrow = 2, byrow = FALSE))

  plot(NA,
       xlim = range(x$lo,x$up),
       ylim = range(x$b_lo,x$b_up),
       xlab = "",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       xpd = NA,
       bty = "n")

  for(j in 1:p){

    lines(x = c(x$lar_out$C * sign(x$lar_out$th[j]),0),
          y = c(0,x$lar_out$b[j,]),
          col = ifelse(x$rej[j],j,"gray"))

  }

  abline(v = 0, lty = 1)
  abline(h = 0, lty = 1)

  # plot CIs for the thetas
  plot(NA,
       xlim = range(x$lo,x$up),
       ylim = c(0,1),
       xlab = "Entrance correlations",
       ylab = "",
       yaxt = "n",
       xpd = NA,
       bty = "n")

  hundy <- (grconvertY(1,from="ndc", to ="user") - grconvertY(0,from="ndc", to ="user"))/300

  abline(v = 0, lty = 1)


  k <- 1
  for(j in order(x$lar_out$th)){

    curvy(x = c(x$lo[j],x$up[j]), y = k/p, col = ifelse(x$rej[j],j,"darkgray"))

    if(j %in% which(x$rej )){

      text(labels = lab[j],
           x = x$up[j] + hundy,
           y = rep(k/p,2),
           adj = 0,
           cex = .8,
           xpd = NA)

    }

    k <- k + 1

  }

  # plot CIs for the betas
  plot(NA,
       xlim = c(0,1),
       ylim = range(x$b_lo,x$b_up),
       xlab = "",
       ylab = "Regression coefficients",
       xaxt = "n",
       xpd = NA,
       bty = "n")

  abline(h = 0, lty = 1)

  k <- 1
  for(j in order(x$lar_out$th)){

    curvy(x = k/p, y = c(x$b_lo[j],x$b_up[j]),col = ifelse(x$b_rej[j],j,"darkgray"))

    if(j %in% which(x$b_rej)){

      text(labels = lab[j],
           x = rep(k/p,2),
           y = x$b_up[j] + hundy,
           cex = .8,
           adj = 0,
           xpd = NA,
           srt = 90)

    }

    k <- k + 1

  }

  par(op)
  layout(1)

}


# Draw curvy horizontal or vertical lines
curvy <- function(x,y,amp = 1,frq = 100, col = col){

  if( (length(x) == 2) & (length(y) == 1)){

    x.ndc <- grconvertX(c(x[1],x[2]),from="user",to="ndc")
    y.ndc <- grconvertY(y,from="user",to="ndc")

    x.seq.ndc <- seq(x.ndc[1], x.ndc[2], length = 1000 )
    y.seq.ndc <- y.ndc + .25*amp / 100 * sin(  x.seq.ndc*frq*2*pi )

  } else if( (length(x) == 1) & (length(y) == 2)){

    x.ndc <- grconvertX(x,from="user",to="ndc")
    y.ndc <- grconvertY(c(y[1],y[2]),from="user",to="ndc")

    y.seq.ndc <- seq(y.ndc[1], y.ndc[2], length = 1000 )
    x.seq.ndc <- x.ndc + .25*amp / 100 * sin(  y.seq.ndc*frq*2*pi )

  } else{

    print("must input a horizontal or vertical line")

  }

  x.seq <- grconvertX(x.seq.ndc,from="ndc",to="user")
  y.seq <- grconvertY(y.seq.ndc,from="ndc",to="user")

  lines(x.seq,y.seq, col = col)

}


# Update the Cholesky decomposition
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
#' Print a table of estimated entrance correlations with bootstrap confidence intervals.
#'
#' @param x an object of class `larinf`, usuall the result of a call to `larinf`.
#' @param ... further arguments.
#'
#' @author Karl Gregory
#'
#' @examples
#' # set parameters for generating data
#' n <- 100
#' p <- 10
#' sigma <- 1
#' beta <- c(1,-1,2,-1.5,rep(0,p-4))
#' Sigma <- (1/2)^abs(outer(1:p,1:p,FUN = "-"))
#' chol_Sigma <- chol(Sigma)
#'
#' # generate some data
#' X <- scale(matrix(rnorm(n*p),n,p) %*% chol_Sigma)*sqrt(n/(n-1))
#' mu <- drop(X %*% beta)
#' e <- rnorm(n,0,sigma)
#' y <- mu + e - mean(e)
#'
#' # perform least angle regression inference
#' larinf_out <- larinf(X,y)
#' print(larinf_out)
#'
#' @export
print.larinf <- function(x,...){

  cat("In order of entrance:\n\n")

  ord <- x$lar_out$ord
  k <- length(ord)
  tab <- round(cbind(x$lar_out$th[ord],x$lo[ord],x$up[ord]),3)

  colnames(tab) <- c("Entrance correlation",paste(x$alpha/2*100,"%",sep=""),paste((1-x$alpha/2)*100,"%",sep=""))

  if(!is.null(x$lar_out$varnames)){

    rownames(tab) <- x$lar_out$varnames[ord]

  } else {

    rownames(tab) <- paste("X",ord,sep= "")
  }

  print(tab)

}



