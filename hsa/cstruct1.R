# Publications that use results obtained from this software please include a citation of the paper:
# Zou, C., Zhang, Y., Ouyang, Z. (2016) HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure. Submitted.

library(MASS)
colSds <- matrixStats::colSds

## sumsq(x) is a faster version of sum(x^2)
sumsq <- inline::cfunction(sig = methods::signature(x = "numeric"), language = "C", body = '
  SEXP res;
  R_xlen_t n = xlength(x);
  double *xx = REAL(x);
  double s = 0;
  for (R_xlen_t ii = 0; ii < n; ++ii) s += xx[ii] * xx[ii];
  return ScalarReal(s);
')

## sumprod(x, y) is a faster version of sum(x * y)
sumprod <- inline::cfunction(sig = methods::signature(x = "numeric", y = "numeric"), language = "C", body = '
  SEXP res;
  R_xlen_t n = xlength(x);
  double *xx = REAL(x);
  double *yy = REAL(y);
  double s = 0;
  if (xlength(y) != n) error("Argument \'x\' and \'y\' are of different lengths");
  for (R_xlen_t ii = 0; ii < n; ++ii) s += xx[ii] * yy[ii];
  return ScalarReal(s);
')

## PERFORMANCE: dist_matrix(x) is a faster version of as.matrix(dist(x)),
## because it avoids the overhead from S3 method dispatching, handling of
## non-needed attributes etc.  Moreover, dist_matrix(x, square = TRUE) is
## avoids internal duplication of the distance matrix.
dist_matrix <- function(x, square = FALSE) {
  x <- dist(x)
  size <- attr(x, "Size")
  df <- matrix(0, nrow = size, ncol = size)
  if (square) x <- x ^ 2
  df[row(df) > col(df)] <- x
  df + t(df)
}

## PERFORMANCE: log_det(x) is a faster version of log(det(x)), because
## it avoids overhead from S3 method dispatching and skips an internal
## log(exp(t)) step.
log_det <- function(x) {
  z <- determinant.matrix(x, logarithm = TRUE)
#  d <- c(z$sign * exp(z$modulus))
#  log(d)
  if (z$sign < 0) stop("Log-determinant: NaN")
  c(z$modulus)
}

## PERFORMANCE: Avoid overhead from S3 dispatch on solve()
solve <- local({
  bs <- list()
  function(a, b = NULL) {
    if (is.null(b)) {
      n <- nrow(a)
      ## Memoization of 'b'
      if (n <= length(bs)) b <- bs[[n]]
      if (is.null(b)) {
        b <- diag(1, nrow = n)
        bs[[n]] <<- b
      }
    }
    solve.default(a, b)
  }
})

## PERFORMANCE: Avoid overhead from setting dimnames in cbind() and rbind()
cbind <- local({
  base_cbind <- base::cbind
  function(...) base_cbind(..., deparse.level = 0L)
})

rbind <- local({
  base_rbind <- base::rbind
  function(...) base_rbind(..., deparse.level = 0L)
})

normP <- function(P) {
  P1 <- t(t(P) - colMeans(P))
  5 / sqrt(max(rowSums(P1^2))) * P1
}

kinetic0 <- function(p0) {
  sumsq(p0) / dim(p0)[1]
}

kinetic0_1 <- function(p0) {
  sumsq(p0)
}

kinetic <- function(p0, N, rho) {
  (sumsq(p0) - sumprod(p0[-1, ], p0[-N, ]) * rho * 2) / N
}

kinetic_1 <- function(p0, N, rho) {
  sumsq(p0) - sumprod(p0[-1, ], p0[-N, ]) * rho * 2
}

momentum <- local({
  zeros_3 <- rep(0, times = 3)
  function(p0, N, rho) {
    p0 - rbind(p0[-1, ] * rho, zeros_3) - rbind(zeros_3, p0[-N, ] * rho)
  }
})

momentum0 <- function(p0) p0

finistructure <- function(S0, bin) {
  n <- dim(S0)[1]
  N <- nrow(bin)
  if (n == N) {
    if (dim(S0)[2] == 3) {
      S <- S0
    } else {
      S <- S0[, 3:5]
    }
  } else {
    pts <- c(S0[, 1], S0[n, 2])
    S0_c35 <- S0[, 3:5]
    Y <- as_rbind(S0_c35, rnorm(3, mean = S0_c35[n, ], sd = colSds(S0_c35[-1, ] - S0_c35[-n, ])))
    bin_c1 <- bin[, 1]
    S <- normP(sapply(1:3, FUN = function(x) splinefun(pts, Y[, x])(bin_c1)))
    S <- S + matrix(rnorm(3 * N, mean = 0, sd = sqrt(5 / N)), nrow = N, ncol = 3L)
  }
  S
}

## PERFORMANCE: fmkorder_temp(), fmkorder(), and fmkorder2() now returns a list of
## parameter estimates that can be used "as is" in the likelihood calculations.
## Previously these functions would stack these parameters up as columns in a matrix
## which then had to be extracted again.  Returning the parameters as is also makes
## it much clear what is calculated.
fmkorder_temp <- function(m, A, b, sigma, S) {
  if (m < 2) {
    mu <- A %*% S + b
    Sigma <- sigma
  } else {
    tmp <- fmkorder_temp(m - 1, A, b, sigma, S)
    mu <- A %*% tmp$mu + b
    temp <- tmp$A %*% A
    Sigma <- tmp$Sigma + t(temp) %*% sigma %*% temp
    A <- temp
  }
  list(mu = mu, Sigma = Sigma, A = A)
}

fmkorder <- function(m, A, b, sigma, S) {
  temp <- fmkorder_temp(m, A, b, sigma, S)
  Sigma <- temp$Sigma
  Sigma <- solve(Sigma)
  temp$mu <- c(temp$mu)
  temp$Sigma <- Sigma
  temp
}

fmkorder2 <- local({
  zeros_3 <- rep(0, times = 3)
  function(m, A, b, sigma) fmkorder(m, A, b, sigma, S = zeros_3)
})

fnormvec <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
}

fangle <- function(a, b) {
  acos(sumprod(a, b) / sqrt(sumsq(a) * sumsq(b)))
}

frotanyvec <- function(x, v, theta) {
  n <- v / sqrt(sumsq(v))
  t <- sumprod(x, n) * n
  cos(theta) * (x - t) + sin(theta) * fnormvec(n, x) + t
}

fbead <- function(S1, S2) {
  m <- dim(S1)[1]
  # S=S1*sqrt(sum((S2[1,]-S2[2,])^2)/sum((S1[1,]-S1[m,])^2))
  S <- S1
  S_r1 <- S[1, ]
  S_rm <- S[m, ]
  S2_r1 <- S2[1, ]
  S2_r2 <- S2[2, ]
  n <- fnormvec(S_rm - S_r1, S2_r2 - S2_r1)
  theta <- fangle(S_rm - S_r1, S2_r2 - S2_r1)
  S_t <- t(S) - S_r1
  S_t <- apply(S_t, MARGIN = 2L, FUN = frotanyvec, v = n, theta = theta)
  S <- t(S_t - S_t[, 1] + S1[1, ])
  S
}

fmirror <- local({
  diag_3 <- diag(3)
  function(v) {
    if (any(v)) {
      diag_3 - v %*% t(v) / sumsq(v)
    } else {
      diag_3
    }
  }
})

tranS <- local({
  zeros_9 <- rep(0, times = 9)
  function(S1, S2, I_scale = TRUE) {
    tmp <- cbind(1, S1)
    tmp_t <- t(tmp)
    beta <- ginv(tmp_t %*% tmp) %*% (tmp_t %*% S2)
    # beta=solve(t(tmp)%*%tmp,t(tmp)%*%S2)
    s <- svd(beta[-1, ])
    tmp2 <- s$u %*% (diag(sign(s$d))) %*% t(s$v)
    beta[-1, ] <- tmp2
    S <- tmp %*% beta
    if (I_scale) {
      beta[-1, ] <- mean(abs(s$d)) * tmp2
      tmp <- optim(c(1, 0, 0, 0, 0, 0, 0), fn = function(x) {
        sum(sqrt(rowSums(
          (t(t(x[1] * S %*% angle2mat(x[2:4])) + x[5:7]) - S2)^2
        )))
      })
      tmp <- tmp$par
      S <- t(t(tmp[1] * S %*% angle2mat(tmp[2:4])) + tmp[5:7])
    } else {
      tmp <- optim(zeros_9, fn = function(x) {
        sum(sqrt(rowSums(
          (t(t(S %*% angle2mat(x[4:6]) %*% fmirror(x[1:3])) + x[7:9]) - S2)^2
        )))
      })
      tmp <- tmp$par
      S <- t(t(S %*% angle2mat(tmp[4:6]) %*% fmirror(tmp[1:3])) + tmp[7:9])
    }
    S
  }
})

rmol <- function(loci, P) {
  n <- dim(P)[1]
  m <- dim(P)[2]
  P1 <- P
  v <- !is.finite(P1)
  if (any(v)) {
    outlier <- v[, 1] | v[, 2] | v[, 3]
    spf <- splinefun
    P2 <- P[!outlier, ]
    loci_outlier <- loci[outlier]
    loci_nonoutlier <- loci[!outlier]
    tmp <- sapply(1:m, FUN = function(x) spf(loci_nonoutlier, P2[, x])(loci_outlier))
    # print(tmp)
    P1[outlier, ] <- tmp
  }
  
  d1 <- sqrt(rowSums((P1[-1, ] - P1[-n, ])^2))
  cutoff <- 10 * median(d1)
  outlier <- (d1 >= cutoff)
  if (any(outliers)) {
    v <- which(outlier)
    ## HB: The below does lots of t(x) over and over.
    ## Not optimized because not covered by the test.
    for (i in 1:length(v)) {
      v_i <- v[i]
      idxs <- 1:v[i]
      P1[idxs, ] <- t(t(P1[idxs, ]) + (P[v_i + 1L, ] - P[v_i, ]) * 0.8)
    }
    P1 <- tranS(P1, S2 = P)
  }
  # plot3d(P1[,1],P1[,2],P1[,3],type="l")
  # points3d(P1[!outlier,1],P1[!outlier,2],P1[!outlier,3],col='green')
  # points3d(P[outlier,1],P[outlier,2],P[outlier,3],col='blue')
  # points3d(P1[outlier,1],P1[outlier,2],P1[outlier,3],col='red')
  P1
}

avsmth <- local({
  thirds_3 <- rep(1/3, times = 3L)
  function(bin, P) {
    N <- length(bin)
    # dP=sqrt(rowSums((P[-1,]-P[-N,])^2))
    P0 <- filter(P, thirds_3, sides = 2L)
    P0[1, ] <- (P[1, ] + P[2, ]) / 2
    P0[N, ] <- (P[N, ] + P[N - 1L, ]) / 2
    matrix(P0, nrow = N, ncol = 3L)
  }
})

loglikelihood0 <- function(P0, A, b, invSigma, beta, cx, mat, pos = NULL, v = NULL, mak = NULL) {
  L <- 0
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.na(C)) {
    C <- 1L
    cx <- array(cx, dim = c(dim(cx), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  if (is.null(pos)) {
    pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
    if (!is.list(pos)) {
      pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
    }
  }
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    cx_i <- cx[pos_i, pos_i, i]
    ## HB: The following calculation is one of the most expensive ones
    ##     in the whole program. /HB 2018-04-02
    temp <- -cx_i * distmat_i^beta_i
    temp[v_i] <- temp[v_i] + mat[pos_i, pos_i + 1L, i][v_i] * (beta_i * log(distmat_i[v_i]) + log(cx_i[v_i]))
    ## sum2(..., idxs)?
    L <- L + sum(temp[upper.tri(temp)]) / (3 * N) # +sum(temp[lower.tri(temp))
  }
  L
}

dloglikelihood0 <- function(P0, A, b, invSigma, beta, cx, mat, pos, v = NULL, mak = NULL) {
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  # pos=apply(mat,3,function(x) which(!is.na(x[,1])))
  # if(!is.list(pos)){pos=lapply(1:nrow(t(pos)),function(i) t(pos)[i,])}
  dL <- matrix(0, nrow = N, ncol = 3L)
  distmat <- dist_matrix(P, square = TRUE) # apply(P*P,1,sum)%*%t(rep(1,N))+rep(1,N)%*%t(apply(P*P,1,sum))-2*P%*%t(P)
  temp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] - beta_i * cx[pos_i, pos_i, i] * (distmat_i^(beta_i / 2 - 1))
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dL[, 1] <- colSums(t(tmp) - tmp) / (3 * N)
  tmp <- temp * P[, 2]
  dL[, 2] <- colSums(t(tmp) - tmp) / (3 * N)
  tmp <- temp * P[, 3]
  dL[, 3] <- colSums(t(tmp) - tmp) / (3 * N)
  dL
}

loglikelihood <- function(P0, A, b, invSigma, beta, cx, mat, pos = NULL, v = NULL, mak = NULL) {
  L <- 0
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.na(C)) {
    C <- 1L
    cx <- array(cx, dim = c(dim(cx), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  if (is.null(pos)) {
    pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
    if (!is.list(pos)) {
      pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
    }
  }
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    cx_i <- cx[pos_i, pos_i, i]
    ## HB: The following calculation is one of the most expensive ones
    ##     in the whole program. /HB 2018-04-02
    temp <- -cx_i * distmat_i^beta_i
    temp[v_i] <- temp[v_i] + mat[pos_i, pos_i + 1L, i][v_i] * (beta_i * log(distmat_i[v_i]) + log(cx_i[v_i]))
    ## sum2(..., idxs)?
    L <- L + sum(temp[upper.tri(temp)]) / (3 * N)
    # temp=-cx[pos[[i]],pos[[i]],i]*distmat[pos[[i]],pos[[i]]]^beta[i]/N/3+mat[pos[[i]],pos[[i]]+1,i]*(beta[i]*log(distmat[pos[[i]],pos[[i]]])+log(cx[pos[[i]],pos[[i]],i]))/N/3
    # L=L+sum(temp[upper.tri(temp)])# +sum(temp[lower.tri(temp))
  }
  if (is.null(mak)) {
    L <- L + sum(sapply(2:N, FUN = function(ii) {
      nmp <- fmkorder(P0[ii, 1] - P0[ii - 1L, 1], A = A, b = b, sigma = sigma, S = P[ii - 1L, ])
      R_ii <- P[ii, ] - nmp$mu
      Sigma <- nmp$Sigma
      -R_ii %*% Sigma %*% R_ii / 2 + log_det(Sigma) / 2
    })) / (3 * N)
  } else {
    L <- L + sum(sapply(2:N, FUN = function(ii) {
      mak_ii1 <- mak[[ii - 1L]]
      mu <- mak_ii1$mu
      Sigma <- mak_ii1$Sigma
      A <- mak_ii1$A
      mu <- mu + A %*% P[ii - 1L, ]
      R_ii <- P[ii, ] - mu
      -t(R_ii) %*% Sigma %*% R_ii / 2 + log_det(Sigma) / 2
    })) / (3 * N)
  }
  L
}

dloglikelihood <- function(P0, A, b, invSigma, beta, cx, mat, pos, v = NULL, mak = NULL) {
  P <- P0[, -1]
  sigma <- solve(invSigma)
  # pos=apply(mat,3,function(x) which(!is.na(x[,1])))
  # if(!is.list(pos)){pos=lapply(1:nrow(t(pos)),function(i) t(pos)[i,])}
  N <- dim(P)[1]
  C <- dim(cx)[3]
  dL <- matrix(0, nrow = N, ncol = 3L)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P, square = TRUE)
  # distmat=apply(P*P,1,sum)%*%t(rep(1,N))+rep(1,N)%*%t(apply(P*P,1,sum))-2*P%*%t(P)
  temp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] - beta_i * cx[pos_i, pos_i, i] * (distmat_i^(beta_i / 2 - 1))
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dL[, 1] <- colSums(t(tmp) - tmp) / (3 * N)
  tmp <- temp * P[, 2]
  dL[, 2] <- colSums(t(tmp) - tmp) / (3 * N)
  tmp <- temp * P[, 3]
  dL[, 3] <- colSums(t(tmp) - tmp) / (3 * N)
  if (is.null(mak)) {
    nmp <- lapply(2:N, FUN = function(ii) {
      fmkorder(P0[ii, 1] - P0[ii - 1L, 1], A = A, b = b, sigma = sigma, S = P[ii - 1L, ])
    })
    
    temp1 <- sapply(2:N, FUN = function(ii) {
      nmp_iim1 <- nmp[[ii - 1L]]
      mu <- nmp_iim1$mu
      Sigma <- nmp_iim1$Sigma
      -t(P[ii, ] - mu) %*% Sigma
    }) / (3 * N)
    dL[2:N, ] <- dL[2:N, ] + t(temp1)
    
    temp2 <- sapply(1:(N - 1L), FUN = function(ii) {
      nmp_ii <- nmp[[ii]]
      mu <- nmp_ii$mu
      Sigma <- nmp_ii$Sigma
      A <- nmp_ii$A
      A_t <- t(A)
      -A_t %*% Sigma %*% A %*% P[ii, ] - A_t %*% Sigma %*% (mu - A %*% P[ii, ] - P[ii + 1L, ])
    }) / (3 * N)
    dL[1:(N - 1L), ] <- dL[1:(N - 1L), ] + t(temp2)
  } else {
    temp1 <- sapply(2:N, FUN = function(ii) {
      mak_ii1 <- mak[[ii - 1L]]
      mu <- mak_ii1$mu
      Sigma <- mak_ii1$Sigma
      A <- mak_ii1$A
      mu <- mu + A %*% P[ii - 1L, ]
      -t(P[ii, ] - mu) %*% Sigma
    }) / (3 * N)
    dL[2:N, ] <- dL[2:N, ] + t(temp1)
    temp2 <- sapply(1:(N - 1L), FUN = function(ii) {
      mak_ii <- mak[[ii]]
      mu <- mak_ii$mu
      Sigma <- mak_ii$Sigma
      A <- mak_ii$A
      A_t <- t(A)
      -A_t %*% Sigma %*% A %*% P[ii, ] - A_t %*% Sigma %*% (mu - P[ii + 1L, ])
    }) / (3 * N)
    dL[1:(N - 1L), ] <- dL[1:(N - 1L), ] + t(temp2)
  }
  dL
}

mkcloglikelihood <- function(theta, P0) {
  P <- P0[, -1]
  N <- dim(P)[1]
  A <- matrix(theta[1:9], nrow = 3L, ncol = 3L)
  b <- theta[10:12]
  invSigma <- matrix(0, nrow = 3L, ncol = 3L)
  invSigma[upper.tri(invSigma, diag = TRUE)] <- theta[-(1:12)]
  invSigma <- invSigma + t(invSigma) - diag(diag(invSigma))
  sigma <- solve(invSigma)
  temp <- sapply(2:N, FUN = function(ii) {
    nmp <- fmkorder(P0[ii, 1] - P0[ii - 1L, 1], A = A, b = b, sigma = sigma, S = P[ii - 1L, ])
    R_ii <- P[ii, ] - nmp$mu
    Sigma <- nmp$Sigma
    -R_ii %*% Sigma %*% R_ii / 2 + log_det(Sigma) / 2
  })
  mean(temp) / 3
}

dhllk <- function(index, theta, P0, A, b, invSigma, beta, cx, mat, pos, v = NULL) {
  P <- rbind(P0[index[1, 1]:index[1, 2], -1], do.call(rbind, args = sapply(2:dim(index)[1], FUN = function(x) {
    t(t(P0[index[x, 1]:index[x, 2], -1] %*% matrix(theta[x - 1L, 1:9], nrow = 3L, ncol = 3L)) + theta[x - 1L, 10:12])
  })))
##HB:  P <- as_matrix(P)
  # sigma=solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  # cat(dim(theta),";")
  dL <- matrix(0, nrow = dim(theta)[1], ncol = 12L)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P, square = TRUE)
  temp <- matrix(0, nrow = N, ncol = N)
  dD <- array(0, dim = c(N, N, 3L))
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] - beta_i * cx[pos_i, pos_i, i] * (distmat_i^(beta_i / 2 - 1))
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }

  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dD[, , 1] <- t(tmp) - tmp
  tmp <- temp * P[, 2]
  dD[, , 2] <- t(tmp) - tmp
  tmp <- temp * P[, 3]
  dD[, , 3] <- t(tmp) - tmp
  index_rn1 <- index[-1, ]
  tmp <- t(apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]    
    P0t <- P0[idxs, -1]
    c(colSums(dD[-idxs, idxs, 1] %*% P0t),
      colSums(dD[-idxs, idxs, 2] %*% P0t),
      colSums(dD[-idxs, idxs, 3] %*% P0t))
  })) / (3 * N)
  dL[, 1:9] <- tmp # t(apply(index[-1,],1,function(x) c(colSums(dD[-(x[1]:x[2]),x[1]:x[2],1]%*%P0[x[1]:x[2],-1]),colSums(dD[-(x[1]:x[2]),x[1]:x[2],2]%*%P0[x[1]:x[2],-1]),colSums(dD[-(x[1]:x[2]),x[1]:x[2],3]%*%P0[x[1]:x[2],-1]))))/N/3

  dL[, 1:9] <- dL[, 1:9] + t(apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sapply(2:4, FUN = function(k) {
      P0t <- P0[idxs, k]
      T_1 <- dD[idxs, idxs, 1] * P0t
      T_2 <- dD[idxs, idxs, 2] * P0t
      T_3 <- dD[idxs, idxs, 3] * P0t
      c(sum(upper.tri(t(T_1)) - upper.tri(T_1)),
        sum(upper.tri(t(T_2)) - upper.tri(T_2)),
        sum(upper.tri(t(T_3)) - upper.tri(T_3)))
    })
  })) / (3 * N)
  
  dL[, 10] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 1])
  }) / (3 * N)
  dL[, 11] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 2])
  }) / (3 * N)
  dL[, 12] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 3])
  }) / (3 * N)
  dL
}

rotamat <- function(theta, full = TRUE) {
  sin_theta <- sin(theta)
  cos_theta <- cos(theta)
  R <- matrix(c(
    cos_theta[1] * cos_theta[3] - cos_theta[2] * sin_theta[1] * sin_theta[3],
   -cos_theta[2] * cos_theta[3] * sin_theta[1] - cos_theta[1] * sin_theta[3],
    sin_theta[1] * sin_theta[2],
    cos_theta[3] * sin_theta[1] + cos_theta[1] * cos_theta[2] * sin_theta[3],
    cos_theta[1] * cos_theta[2] * cos_theta[3] - sin_theta[1] * sin_theta[3],
   -cos_theta[1] * sin_theta[2],
    sin_theta[2] * sin_theta[3],
    cos_theta[3] * sin_theta[2],
    cos_theta[2]
  ), nrow = 3L, ncol = 3L)
  
  if (full) R <- cbind(R, sin_theta, cos_theta)

  R
}

dDtotheta <- function(p, b, temp) {
  l <- matrix(0, nrow = dim(b)[1], ncol = 3L)

  temp_11 <- temp[1, 1]
  temp_12 <- temp[1, 2]
  temp_13 <- temp[1, 3]
  temp_14 <- temp[1, 4]
  temp_15 <- temp[1, 5]
  temp_21 <- temp[2, 1]
  temp_22 <- temp[2, 2]
  temp_23 <- temp[2, 3]
  temp_31 <- temp[3, 1]
  temp_32 <- temp[3, 2]
  temp_33 <- temp[3, 3]
  temp_35 <- temp[3, 5]

  b_c1 <- b[, 1]
  b_c2 <- b[, 2]
  b_c3 <- b[, 3]
  
  l[, 1] <- -2 * (p[2] * -temp_22 - p[1] * temp_12 - temp_32 * p[3]) * (b_c1 - p[1] * temp_11 + p[2] * -temp_21 - temp_31 * p[3]) - 2 * (p[1] * temp_11 - p[2] * -temp_21 + temp_31 * p[3]) * (b_c2 - p[1] * temp_12 + p[2] * -temp_22 - temp_32 * p[3])

  l[, 2] <- 2 * (temp_15 * temp_33 * p[3] + temp_15 * temp_23 * p[2] + temp_15 * temp_13 * p[1]) * (b_c2 - p[1] * temp_12 + p[2] * -temp_22 - temp_32 * p[3]) - 2 * (temp_33 * temp_14 * p[3] + temp_35 * temp_31 * p[2] + temp_14 * temp_13 * p[1]) * (b_c1 - p[1] * temp_11 + p[2] * -temp_21 - temp_31 * p[3]) + 2 * (temp_33 * temp_35 * p[2] - temp[2, 4] * p[3] + temp_33 * temp[3, 4] * p[1]) * (temp_33 * p[3] - b_c3 + temp_23 * p[2] + temp_13 * p[1])

  l[, 3] <- 2 * (temp_23 * p[1] - temp_13 * p[2]) * (temp_33 * p[3] - b_c3 + temp_23 * p[2] + temp_13 * p[1]) + 2 * (p[1] * -temp_22 + p[2] * temp_12) * (b_c2 - p[1] * temp_12 + p[2] * -temp_22 - temp_32 * p[3]) + 2 * (p[1] * -temp_21 + p[2] * temp_11) * (b_c1 - p[1] * temp_11 + p[2] * -temp_21 - temp_31 * p[3])

  # l[1]=- 2*(p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + cos(theta[1])*sin(theta[2])*p[3])*(b[1] - p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) + p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) - sin(theta[1])*sin(theta[2])*p[3]) - 2*(p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) - p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) + sin(theta[1])*sin(theta[2])*p[3])*(b[2] - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + cos(theta[1])*sin(theta[2])*p[3])

  # l[2]=2*(cos(theta[1])*cos(theta[2])*p[3] + cos(theta[1])*cos(theta[3])*sin(theta[2])*p[2] + cos(theta[1])*sin(theta[2])*sin(theta[3])*p[1])*(b[2] - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + cos(theta[1])*sin(theta[2])*p[3]) - 2*(cos(theta[2])*sin(theta[1])*p[3] + cos(theta[3])*sin(theta[1])*sin(theta[2])*p[2] + sin(theta[1])*sin(theta[2])*sin(theta[3])*p[1])*(b[1] - p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) + p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) - sin(theta[1])*sin(theta[2])*p[3]) + 2*(cos(theta[2])*cos(theta[3])*p[2] - sin(theta[2])*p[3] + cos(theta[2])*sin(theta[3])*p[1])*(cos(theta[2])*p[3] - b[3] + cos(theta[3])*sin(theta[2])*p[2] + sin(theta[2])*sin(theta[3])*p[1])

  # l[3]=2*(cos(theta[3])*sin(theta[2])*p[1] - sin(theta[2])*sin(theta[3])*p[2])*(cos(theta[2])*p[3] - b[3] + cos(theta[3])*sin(theta[2])*p[2] + sin(theta[2])*sin(theta[3])*p[1]) + 2*(p[1]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + p[2]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])))*(b[2] - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + cos(theta[1])*sin(theta[2])*p[3]) + 2*(p[1]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) + p[2]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])))*(b[1] - p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) + p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) - sin(theta[1])*sin(theta[2])*p[3])
  l
}

dhllk1 <- function(index, theta, P0, A, b, invSigma, beta, cx, mat, pos, v = NULL) {
  matheta <- lapply(2:dim(index)[1], FUN = function(x) rotamat(theta[x - 1L, ]))
  P <- rbind(P0[index[1, 1]:index[1, 2], -1], do.call(rbind, args = sapply(2:dim(index)[1], FUN = function(x) {
    t(t(P0[index[x, 1]:index[x, 2], -1] %*% matheta[[x - 1L]][, 1:3]) + theta[x - 1L, 4:6])
  })))
## HB:  P <- as_matrix(P)
  # sigma=solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  # cat(dim(theta),";")
  dL <- matrix(0, nrow = dim(theta)[1], ncol = 6L)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  }
  distmat <- dist_matrix(P, square = TRUE)
  temp <- matrix(0, nrow = N, ncol = N)
  dD <- array(0, dim = c(N, N, 3L))
  for (i in 1:C) {
    pos_i <- pos[[i]]
    v_i <- v[[i]]
    beta_i <- beta[i]
    distmat_i <- distmat[pos_i, pos_i]
    temp[pos_i, pos_i] <- temp[pos_i, pos_i] - beta_i * cx[pos_i, pos_i, i] * (distmat_i^(beta_i / 2 - 1))
    temp[pos_i, pos_i][v_i] <- temp[pos_i, pos_i][v_i] + beta_i * mat[pos_i, pos_i + 1L, i][v_i] / distmat_i[v_i]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dD[, , 1] <- t(tmp) - tmp
  tmp <- temp * P[, 2]
  dD[, , 2] <- t(tmp) - tmp
  tmp <- temp * P[, 3]
  dD[, , 3] <- t(tmp) - tmp
  dL[, 1:3] <- t(sapply(2:dim(index)[1], FUN = function(x) {
    idxs <- index[x, 1]:index[x, 2]
    rowSums(sapply(idxs, FUN = function(i) colSums(dDtotheta(P0[i, -1], b = t(t(P[-idxs, ]) - theta[x - 1, 4:6]), temp = matheta[[x - 1]]) * temp[-idxs, i])))
  }))
  index_rn1 <- index[-1, ]
  dL[, 4] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 1])
  })
  dL[, 5] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 2])
  })
  dL[, 6] <- apply(index_rn1, MARGIN = 1L, FUN = function(x) {
    idxs <- x[1]:x[2]
    sum(dD[-idxs, idxs, 3])
  })
  dL / (3 * N)
}

angle2mat <- function(theta) {
  rotamat(theta, full = FALSE)
}

v2mat <- function(theta) {
  matrix(theta[1:9], nrow = 3L, ncol = 3L)
}

Leapfrog <- function(grad_U, L, epsilon, p0, q0, fM) {
  N <- dim(p0)[1]
  if (L < 2) {
    # q = (q0 + epsilon * p0)
    q <- q0 + epsilon * fM(p0)
    p <- p0 - epsilon * grad_U(q) / 2
  } else {
    temp <- Leapfrog(grad_U, L - 1, epsilon, p0, q0, fM)
    # q=temp[[2]]+epsilon*(temp[[1]])
    q <- temp[[2]] + epsilon * fM(temp[[1]])
    p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
  }
  list(p, q)
}

HMC <- function(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE) {
  N <- dim(current_q0)[1]
  m <- dim(current_q0)[2]
  current_q <- current_q0
  q <- current_q
  if (I_trans) {
    q <- rmol(1:N, P = current_q)
    # q=tranS(q,current_q)
  }
  # p = array(0,dim(q))
  p <- array(rnorm(N * m, mean = 0, sd = 1 / sqrt(N)), dim = dim(q)) # independent standard normal variates
  current_p <- p
  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  temp <- Leapfrog(grad_U, L, epsilon, p, q, fM)
  # Make a half step for momentum at the end.
  q <- temp[[2]]
  p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # p[1:2,]=0
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- fK(current_p) / 2
  proposed_U <- U(q)
  proposed_K <- fK(p) / 2
  const <- (current_U - proposed_U + current_K - proposed_K)
  # if(N>1000){
  # cat(current_U,current_K,proposed_U,proposed_K,const,";\n")
  # plot3d(q[,1],q[,2],q[,3],type="l")
  # points3d(q[,1],q[,2],q[,3],col='red')
  # }

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (any(is.infinite(const) + is.na(const))) {
    cat("NA or inf produced!!", "\n")
    # return (matrix(rnorm(N*m,rmol(1:N,current_q),0.005),N,m))
    return(rmol(1:N, P = current_q))
  } else {
    if (runif(1) < exp((const) / T0)) {
      if (I_trans) {
        q <- rmol(1:N, P = q)
        q <- tranS(q, S2 = current_q)
      }
      return(q) # accept
    } else {
      return(current_q) # reject
    }
  }
}

HMC1 <- function(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans = FALSE) {
  N <- dim(current_q0)[1]
  m <- dim(current_q0)[2]
  current_q <- current_q0
  q <- current_q
  # p = array(0,dim(q))
  p <- array(rnorm(N * m, mean = 0, sd = 1 / sqrt(N)), dim = dim(q)) # independent standard normal variates
  current_p <- p
  # Make a half step for momentum at the beginning
  p <- p - epsilon * grad_U(q) / 2
  # Alternate full steps for position and momentum
  temp <- Leapfrog(grad_U, L, epsilon, p, q, fM)
  # Make a half step for momentum at the end.
  q <- temp[[2]]
  p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
  # Negate momentum at end of trajectory to make the proposal symmetric
  p <- -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U <- U(current_q)
  current_K <- fK(current_p) / 2
  proposed_U <- U(q)
  proposed_K <- fK(p) / 2
  const <- (current_U - proposed_U) / log(runif(1) * 0.5)
  const <- ifelse(const > 0, const, 1)

  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  if (any(c(is.infinite(const), is.na(const)))) {
    cat("NA or inf produced!!", "\n")
    return(rmol(1:N, P = current_q))
    # return (matrix(rnorm(N*m,current_q,abs(current_q)/15),N,m))
  } else {
    if (runif(1) < exp((const) / T0 / const)) {
      if (I_trans) {
        q <- rmol(1:N, P = q)
        q <- tranS(q, S2 = current_q)
      }
      return(q) # accept
    } else {
      return(current_q) # reject
    }
  }
}


Qsis <- function(N, pbin, A, b, invSigma, beta, cx, mat, q0, fL) {
  # cat(c(length(pbin[1:N]),dim(mat[1:N,1:(N+1),]),dim(cx[1:N,1:N,])),"~")
  if (N == 2) {
    pbin_t <- pbin[1:2]
    cx_t <- cx[1:2, 1:2, ]
    mat_t <- mat[1:2, 1:3, ]
    Padd <- optim(b + rnorm(3, sd = 1 / 5) + q0, fn = function(q) {
      fL(cbind(pbin_t, rbind(q0, q)), A, b, invSigma, beta, cx_t, mat_t)
    })
    # cat(Padd$par,"\n")
    return(rbind(q0, t(Padd$par)))
  } else {
    # cat(pbin,"\n")
    temp <- Qsis(N - 1, pbin, A, b, invSigma, beta, cx, mat, q0, fL)
    pbin_t <- pbin[1:N]
    cx_t <- cx[1:N, 1:N, ]
    mat_t <- mat[1:N, 1:(N + 1L), ]
    Padd <- optim(rnorm(3, sd = 1 / 5) + A %*% temp[dim(temp)[1], ] + b, fn = function(q) {
      fL(cbind(pbin_t, rbind(temp, q)), A, b, invSigma, beta, cx_t, mat_t)
    })
    # cat(Padd$par,"\n")
    return(rbind(temp, t(Padd$par)))
  }
}

Sis <- function(d, pbin, A, b, invSigma, beta, cx0, mat0, q0, fL) {
  N <- dim(mat0)[1]
  # cat(length(pbin),",")
  if (is.na(dim(cx0)[3])) {
    cx <- array(cx0, dim = c(dim(cx0), 1L))
    mat <- array(mat0, dim = c(dim(mat0), 1L))
  } else {
    cx <- cx0
    mat <- mat0
  }
  if (N == d) {
    # cat(c(length(pbin),dim(mat),dim(cx)),"~")
    return(Qsis(d, pbin, A, b, invSigma, beta, cx, mat, q0, fL))
  } else {
    # cat(c(length(pbin[-N]),dim(mat[-N,-(N+1),])),"\n")
    temp <- Sis(d, pbin[-N], A, b, invSigma, beta, cx[-N, -N, ], mat[-N, -(N + 1L), ], q0, fL)
    # cat(dim(temp),"\n")
    idxs <- (N - d):N
    pbin_t <- pbin[idxs]
    nrow <- nrow(temp)
    temp_t <- temp[(nrow - d + 1L):nrow, ]
    cx_t <- cx[idxs, idxs, ]
    mat_t <- mat[idxs, c(1, (N - d + 1L):(N + 1L)), ]
    addq <- optim(rnorm(3, sd = 1 / 5) + A %*% temp[nrow, ] + b, fn = function(x) {
      fL(cbind(pbin_t, rbind(temp_t, x)), A, b, invSigma, beta, cx_t, mat_t)
    })
    addq <- addq$par
    return(rbind(temp, t(addq)))
  }
}

subinitial <- function(pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike) {
  N <- length(pbin)
  if (is.na(dim(covmat0)[3])) {
    covmat0 <- array(covmat0, dim = c(dim(covmat0), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
  if (!is.list(pos)) {
    pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
  }
  P0 <- Sis(4, pbin, A0, b0, invSigma0, beta1, covmat0, mat, c(0, 0, 0), function(x, ...) -loglikelihood(x, ...))
  u <- 0
  while (u < 100) {
    P <- HMC(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos), 0.002, 20, P0, 10 * exp(-u / 20), function(p) kinetic(p, N = N, rho = 0.1), function(p) momentum(p, N = N, rho = 0.1), 1)
    P0 <- P
    u <- u + 1
  }
  P
}

suboptimz <- function(pbin, P0, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike, I_mle = TRUE) {
  epslon <- 0.0005
  stp <- 35
  M <- 100
  if (is.na(dim(covmat0)[3])) {
    covmat0 <- array(covmat0, dim = c(dim(covmat0), 1L))
    mat <- array(mat, dim = c(dim(mat), 1L))
  }
  pos <- apply(mat, MARGIN = 3L, FUN = function(x) which(!is.na(x[, 1])))
  if (!is.list(pos)) {
    pos <- lapply(1:ncol(pos), FUN = function(i) pos[, i])
  }
  u <- 0
  Po <- P0
  N <- dim(P0)[1]
  Lo <- floglike(cbind(pbin, P0), A0, b0, invSigma0, beta1, covmat0, mat)
  P <- HMC1(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos), epslon, stp, P0, exp(-u * 3.5 / M), kinetic0_1, momentum0)
  P0 <- P
  if (I_mle) {
    while (u < M) {
      P <- HMC1(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos), epslon, stp, P0, exp(-u * 3.5 / M), kinetic0_1, momentum0)
      P0 <- P
      L <- floglike(cbind(pbin, P), A0, b0, invSigma0, beta1, covmat0, mat)
      if (L >= Lo) {
        Po <- P
        Lo <- L
      }
      u <- u + 1
    }
  } else {
    while (u < M) {
      P <- HMC1(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos), epslon, stp, P0, exp(-u * 3.5 / M), kinetic0_1, momentum0)
      P0 <- P
      u <- u + 1
    }
    Po <- P0
  }
  Po
}

piece <- function(theta0, P1, P2, pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike) {
  theta <- matrix(theta0, nrow = 4L, ncol = 3L)
  -floglike(cbind(pbin, rbind(P1, t(t(P2 %*% theta[-4, ]) + theta[4, ]))), A0, b0, invSigma0, beta1, covmat0, mat)
}

finital <- function(pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike) {
  N <- length(pbin)
  if (N <= 500) {
    P <- Sis(5, pbin, A0, b0, invSigma0, beta1, covmat0, mat, c(0, 0, 0), function(x, ...) -loglikelihood(x, ...))
  } else {
    if (N > 500 && N <= 2000) {
      m <- 100
    } else {
      m <- 200
    }
    index <- cbind(c(1, seq(from = m, to = N - m, by = m)), c(seq(from = m, to = N - m, by = m), N))
    index[-1, 1] <- index[-1, 1] + 1L
    lP <- lapply(1:dim(index)[1], FUN = function(x) {
      idxs <- index[x, 1]:index[x, 2]
      subinitial(pbin[idxs], A0, b0, invSigma0, beta1, covmat0[idxs, idxs, ], mat[idxs, c(1, idxs + 1L), ], floglike, fdloglike)
    })
    P <- matrix(0, nrow = N, ncol = 3L)
    P[index[1, 1]:index[1, 2], ] <- lP[[1]]
    for (i in 2:dim(index)[1]) {
      index_ri_c1 <- index[i, 1]
      index_ri_c2 <- index[i, 2]
      idxs <- 1:index_ri_c2
      P_i <- P[1:(index[i - 1, 2]), ]
      pbin_i <- pbin[idxs]
      lP_i <- lP[[i]]
      covmat0_i <- covmat0[idxs, idxs, ]
      mat_i <- mat[idxs, 1:(index_ri_c2 + 1L), ]
      theta <- optim(as.vector(rbind(diag(3), P[index[i - 1, 2], ] - P[index_ri_c1, ] + rnorm(3, sd = 1 / 100))), fn = function(x) {
        piece(x, P_i, lP_i, pbin_i, A0, b0, invSigma0, beta1, covmat0_i, mat_i, floglike)
      })
      theta <- matrix(theta$par, nrow = 4L, ncol = 3L)
      # cat(dim(theta),"\t")
      P[index_ri_c1:index_ri_c2, ] <- lP_i %*% theta[-4, ] + rep(1, times = index_ri_c2 - index_ri_c1 + 1L) %*% t(theta[4, ])
    }
  }
  avsmth(pbin, P = P)
  # return(P)
}

fmain <- function(lsmap0, lscov0, outfile, Maxiter, submaxiter, lambda, Leapfrog, epslon, mkfix = 0, rho = 0, mk, initialS = NULL, coarsefit = TRUE, rmoutlier = FALSE, fitmode = 0) {
  floglike <- loglikelihood0
  fdloglike <- dloglikelihood0
  fcorrect <- rmol
  recodllk <- rep(0, times = 7)
  if (mk) {
    floglike <- loglikelihood
    fdloglike <- dloglikelihood
  }
  C <- length(lsmap0)
  bin <- NULL
  for (i in 1:C) {
    bin <- rbind(bin, lsmap0[[i]][, 1:2])
  }
  bin <- unique(bin, MARGIN = 1L)
  bin <- bin[order(bin[, 1], decreasing = FALSE), ]
  N <- dim(bin)[1]
  bin_c1 <- bin[, 1]
  mbin <- mean(bin[, 2] - bin_c1)
  tmp <- c(1, pmax(floor((bin_c1[-1] - bin_c1[-N]) / mbin), 1))
  neigdc <-  max(tmp)
  pbin <- cumsum(tmp)
  gldata <- vector("list", length = C)
  A0 <- diag(3) #+0*rbind(c(0,-runif(1)/5,0),c(runif(1)/5,0,0),rep(0,3))	
  b0 <- c(0, 0, 0)
  t <- 0:(N - 1)
  invSigma0 <- diag(3)
  invS <- diag(3) * sqrt(N)
  mat <- array(NA_real_, dim = c(N, N + 1L, C))
  mat0 <- mat
  beta1 <- rep(-1.3, times = C)
  alpha <- seq(from = 0.5, to = 1.5, by = 0.5)
  covmat0 <- array(1, dim = c(N, N, C))
  Beta <- vector("list", length = C)
  impute <- 0
  pos <- vector("list", length = C)
  if (is.numeric(lscov0)) {
    for (c in 1:C) {
      pos_c <- which(bin_c1 %in% lsmap0[[c]][, 1])
      mat[pos_c, 1, c] <- pos_c
      pos[[c]] <- pos_c
      # cat(temp,"\n")
      temp <- as.matrix(lsmap0[[c]][, -(1:2)])
      if (isSymmetric(temp)) {
        mat[pos_c, pos_c + 1L, c] <- temp
        gldata[[c]] <- temp[upper.tri(temp)]
      } else {
        temp_t <- t(temp)
        mat[pos_c, pos_c + 1L, c] <- temp + temp_t
        temp <- temp + temp_t
        gldata[[c]] <- temp[upper.tri(temp)]
      }
      mat0[, , c] <- mat[, , c]
      gldata[[c]] <- data.frame(cbind(temp[upper.tri(temp)], rep(1, times = length(gldata[[c]]))))
    }

    # P10=finital(pbin,A0,b0,invSigma0,beta1,covmat0,mat,floglike,fdloglike)
    if (is.null(initialS)) {
      P10 <- finital(pbin, A0, b0, invS, beta1, covmat0, mat, floglike, fdloglike)
    } else {
      P10 <- finistructure(initialS, bin = bin)
    }
    P10 <- normP(P10)
    P01 <- P10
    dmat <- dist_matrix(P10)
    for (c in 1:C) {
      pos_c <- pos[[c]]
      dmat1 <- dmat[pos_c, pos_c]
      dmat1 <- dmat1[upper.tri(dmat1)]
      gldata[[c]][, dim(gldata[[c]])[2]] <- log(dmat1)
      colnames(gldata[[c]]) <- paste0("V", 1:dim(gldata[[c]])[2])
      betaest <- glm(V1 ~ ., family = poisson(), data = gldata[[c]])
      Beta[[c]] <- as.vector(betaest$coefficients)
      cat(Beta[[c]], "\n")
      beta1[c] <- min(-abs(Beta[[c]][length(Beta[[c]])]), beta1[c])
      covmat0[, , c] <- covmat0[, , c] * exp(Beta[[c]][1])
    }
  } else {
    beta11 <- beta1
    # covmat00=array(1,c(N,N,C))
    lscov <- lscov0
    for (c in 1:C) {
      pos_c <- which(bin_c1 %in% lsmap0[[c]][, 1])
      # cat(temp)
      mat[pos_c, 1, c] <- pos_c
      pos[[c]] <- pos_c
      temp <- as.matrix(lsmap0[[c]][, -(1:2)])
      if (isSymmetric(temp)) {
        mat[pos_c, pos_c + 1L, c] <- temp
      } else {
        temp_t <- t(temp)
        mat[pos_c, pos_c + 1L, c] <- temp + temp_t
        temp <- temp + temp_t
      }
      mat0[, , c] <- mat[, , c]
      if (is.numeric(lscov[[c]])) {
        gldata[[c]] <- temp[upper.tri(temp)]
      } else {
        gldata[[c]] <- cbind(temp[upper.tri(temp)], sapply(lscov[[c]], FUN = function(x) log(x[upper.tri(x)])))
      }
      gldata[[c]] <- data.frame(gldata[[c]])
      colnames(gldata[[c]]) <- paste0("V", 1:dim(gldata[[c]])[2])
      betaest <- glm(V1 ~ ., family = poisson(), data = gldata[[c]])
      gldata[[c]] <- cbind(gldata[[c]], rep(1, times = dim(gldata[[c]])[1]))
      Beta[[c]] <- c(as.vector(betaest$coefficients), beta1[c])
      cat(Beta[[c]], "\n")
      covmat0[, , c] <- covmat0[, , c] * exp(Beta[[c]][1])
      # covmat00[,,c]=covmat00[,,c]*exp(Beta[[c]][1])
      if (!is.numeric(lscov[[c]])) {
        pos_c <- pos[[c]]
        for (k in 2:length(Beta[[c]][-1])) {
          covmat0[pos_c, pos_c, c] <- covmat0[pos_c, pos_c, c] * lscov[[c]][[k - 1]]^Beta[[c]][k]
        }
      }
    }

    if (is.null(initialS)) {
      P10 <- finital(pbin, A0, b0, invS, beta1, covmat0, mat, floglike, fdloglike)
    } else {
      P10 <- finistructure(initialS, bin = bin)
    }
    P01 <- P10
    dmat <- dist_matrix(P10)
    covmat0 <- array(1, dim = c(N, N, C))
    for (c in 1:C) {
      pos_c <- pos[[c]]
      dmat1 <- dmat[pos_c, pos_c]
      dmat1 <- dmat1[upper.tri(dmat1)]
      gldata[[c]][, dim(gldata[[c]])[2]] <- log(dmat1)
      colnames(gldata[[c]]) <- paste0("V", 1:dim(gldata[[c]])[2])
      betaest <- glm(V1 ~ ., family = poisson(), data = gldata[[c]])
      Beta[[c]] <- as.vector(betaest$coefficients)
      cat(Beta[[c]], "\n")
      beta1[c] <- min(-abs(Beta[[c]][length(Beta[[c]])]), beta1[c])
      covmat0[, , c] <- covmat0[, , c] * exp(Beta[[c]][1])
      if (!is.numeric(lscov)) {
        if (!is.numeric(lscov[[c]])) {
          for (k in 2:length(Beta[[c]][-1])) {
            covmat0[pos_c, pos_c, c] <- covmat0[pos_c, pos_c, c] * lscov[[c]][[k - 1]]^Beta[[c]][k]
          }
        }
      }
    }
  }
  Loglike0 <- floglike(cbind(pbin, P10), A0, b0, invSigma0, beta1, covmat0, mat)
  cat("LLK0: ", Loglike0, "\n")
  cat("number of nodes:", N, "\n")
  v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1L, i] > 0))
  if (mk) {
    sigma <- solve(invSigma0)
    mak <- lapply(2:N, FUN = function(ii) fmkorder2(pbin[ii] - pbin[ii - 1], A = A0, b = b0, sigma = sigma))
    # lapply(mak,function(x) cat(dim(x),";"))
  } else {
    mak <- NULL
  }
  if (N > 1000) {
    m <- 100
    if (N > 2000) {
      m <- 200
    }
    index <- cbind(c(1, seq(from = m, to = N - m, by = m)), c(seq(from = m, to = N - m, by = m), N))
    index[-1, 1] <- index[-1, 1] + 1L
    m2 <- floor(N / 100)
    index2 <- seq(from = 1, to = N - m2, by = m2)
    index2 <- c(index2, N)
    Psub <- P01
    Psub <- do.call(rbind, args = lapply(1:dim(index)[1], FUN = function(x) {
      idxs <- index[x, 1]:index[x, 2]
      tranS(suboptimz(pbin[idxs], P01[idxs, ], A0, b0, invSigma0, beta1, covmat0[idxs, idxs, ], mat[idxs, c(1L, idxs + 1L), ], floglike, fdloglike, 0), S2 = P01[idxs, ])
    }))
    Logsub <- floglike(cbind(pbin, Psub), A0, b0, invSigma0, beta1, covmat0, mat)
    cat("LLK after suboptimization: ", Logsub, "\n")
    if (Loglike0 <= Logsub) {
      P01 <- Psub
      Loglike0 <- Logsub
    }
    Pframe <- tranS(suboptimz(pbin[index2], P01[index2, ], A0, b0, invSigma0, beta1, covmat0[index2, index2, ], mat[index2, c(1L, index2 + 1L), ], floglike, fdloglike, 0), S2 = P01[index2, ])
    Psub2 <- Psub
    for (i_sub in 1:(length(index2) - 1)) {
      idxs2 <- index2[i_sub]:index2[i_sub + 1L]
      P_tmp <- P01[idxs2, ]
      if (i_sub > 1) {
        P_tmp <- t(t(P_tmp) - P_tmp[1, ] + Psub2[index2[i_sub], ])
      }
      Psub2[idxs2, ] <- fbead(P_tmp, S2 = Pframe[i_sub:(i_sub + 1L), ])
    }
    Psub2 <- tranS(Psub2, S2 = P01)
    # plot3d(Psub2[,1],Psub2[,2],Psub2[,3],type="l")
    # points3d(Psub2[,1],Psub2[,2],Psub2[,3],col='red')
    # lines3d(P01[,1],P01[,2],P01[,3],col='green')
    Logsub <- floglike(cbind(pbin, Psub2), A0, b0, invSigma0, beta1, covmat0, mat)
    cat("LLK after suboptimization2: ", Logsub, "\n")
    if (Loglike0 <= Logsub) {
      P01 <- Psub2
      Loglike0 <- Logsub
    }
  }

  Ao <- A0
  bo <- b0
  invSigmao <- invSigma0
  Po <- P01
  betao <- Beta
  iternum <- 0
  P <- P01
  write.table(cbind(bin, normP(Po)), file = paste0(outfile, ".txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
  write.table(unlist(betao), file = paste0(outfile, "beta.txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
  if (fitmode) {
    fHMC <- HMC
    fknt <- kinetic
    e_c <- 0.02
    eps_coarse <- ifelse(min(sapply(lsmap0, FUN = function(x) sum(x[, -(1:2)] > 0) / dim(x)[1] / dim(x)[1])) < 0.15, 0.005, e_c)
    num_coarse <- min(50, floor(Maxiter / 4))
  } else {
    fHMC <- HMC1
    fknt <- kinetic_1
    e_c <- 0.01
    eps_coarse <- ifelse(min(sapply(lsmap0, FUN = function(x) sum(x[, -(1:2)] > 0) / dim(x)[1] / dim(x)[1])) < 0.15, 0.005, min(e_c, 10 * epslon))
    num_coarse <- min(20, floor(Maxiter / 4))
  }
  while (iternum < Maxiter) {
    u <- 0
    if (iternum < 20) {
      submaxiter1 <- max(submaxiter + 50, 150)
      lambda1 <- max(submaxiter + 50, 50)
    } else {
      if (iternum < 0.8 * Maxiter) {
        submaxiter1 <- submaxiter
        lambda1 <- lambda
      } else {
        submaxiter1 <- min(submaxiter * 1.1, 300)
        lambda1 <- min(lambda * 1.1, 100)
      }
    }

    current_epslon <- ifelse(iternum < num_coarse && coarsefit, eps_coarse, epslon)

    while (u < submaxiter && (N < 1000 || (iternum %% 10 == 1 && iternum > 5) || u < 3)) {
      P <- fHMC(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat, pos, v, mak), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos, v, mak), current_epslon, Leapfrog, P01, exp(-u / lambda), function(p) fknt(p, N, rho), function(p) momentum(p, N = N, rho = rho))
      P01 <- P
      # plot3d(P01[,1],P01[,2],P01[,3],type="l")
      # points3d(P01[,1],P01[,2],P01[,3],col='red')
      u <- u + 1
    }

    if (!mk) {
      P <- normP(P)
    }

    if (mkfix && iternum >= 5) {
      P0 <- cbind(pbin, P)
      thetam <- optim(c(as.vector(A0), b0, as.vector(invSigma0[upper.tri(invSigma0, diag = TRUE)])), fn = function(theta) -mkcloglikelihood(theta, P0 = P0))
      thetam <- thetam$par
      # cat(thetam,"\n")
      A <- matrix(thetam[1:9], nrow = 3L, ncol = 3L)
      b <- thetam[10:12]
      invSigma <- matrix(0, nrow = 3L, ncol = 3L)
      invSigma[upper.tri(invSigma, diag = TRUE)] <- thetam[-(1:12)]
      invSigma <- invSigma + t(invSigma) - diag(diag(invSigma))
      A0 <- A
      b0 <- b
      invSigma0 <- invSigma
    } else {
      A <- A0
      b <- b0
      invSigma <- invSigma0
    }
    # cat(loglikelihood(P,A,b,invSigma,beta1,covmat0,mat),"\n")
    covmat0 <- array(1, dim = c(N, N, C))
    dmat <- dist_matrix(P)
    for (c in 1:C) {
      pos_c <- pos[[c]]
      dmat1 <- dmat[pos_c, pos_c]
      dmat1 <- dmat1[upper.tri(dmat1)]
      gldata[[c]][, dim(gldata[[c]])[2]] <- log(dmat1)
      colnames(gldata[[c]]) <- paste0("V", 1:dim(gldata[[c]])[2])
      betaest <- glm(V1 ~ ., family = poisson(), data = gldata[[c]])
      Beta[[c]] <- as.vector(betaest$coefficients)
      cat(iternum, ": ", Beta[[c]], "\n")
      beta1[c] <- ifelse(iternum <= 5, min(-abs(Beta[[c]][length(Beta[[c]])]), -1.3), min(-abs(Beta[[c]][length(Beta[[c]])]), -0.6))
      covmat0[, , c] <- covmat0[, , c] * exp(Beta[[c]][1])
      if (!is.numeric(lscov0)) {
        if (!is.numeric(lscov0[[c]])) {
          for (k in 2:(length(Beta[[c]]) - 1)) {
            covmat0[pos_c, pos_c, c] <- covmat0[pos_c, pos_c, c] * lscov[[c]][[k - 1]]^Beta[[c]][k]
          }
        }
      }
      # temp=covmat0[,,c]*(dmat^beta1[c])
      # mat[,,c][!mat0[,,c]]=temp[!mat0[,,c]]
    }
    usLoglike <- floglike(cbind(pbin, P), A, b, invSigma, beta1, covmat0, mat)
    cat("LLK after GLM: ", usLoglike, "\n")
    if (N > 1000 && iternum < (Maxiter - 5)) {
      Psub <- do.call(rbind, args = lapply(1:dim(index)[1], FUN = function(x) {
        idxs <- index[x, 1]:index[x, 2]
        tranS(suboptimz(pbin[idxs], P[idxs, ], A0, b0, invSigma0, beta1, covmat0[idxs, idxs, ], mat[idxs, c(1L, idxs + 1L), ], floglike, fdloglike, 0), S2 = P[idxs, ])
      }))
      Logsub <- floglike(cbind(pbin, Psub), A0, b0, invSigma0, beta1, covmat0, mat)
      cat("LLK after suboptimization: ", Logsub, "\n")
      if (Logsub >= usLoglike) {
        P <- Psub
        usLoglike <- Logsub
      }
      Pframe <- tranS(suboptimz(pbin[index2], P[index2, ], A0, b0, invSigma0, beta1, covmat0[index2, index2, ], mat[index2, c(1L, index2 + 1L), ], floglike, fdloglike, 0), S2 = P[index2, ])
      # Psub2=do.call(rbind,lapply(1:(length(index2)-1),function(x)fbead(Psub[index2[x]:index2[x+1],],Pframe[x:(x+1),])[-1,]))
      # Psub2=tranS(rbind(Pframe[1,],Psub2),P)
      Psub2 <- P
      for (i_sub in 1:(length(index2) - 1)) {
        idxs2 <- index2[i_sub]:index2[i_sub + 1L]
        P_tmp <- P[idxs2, ]
        if (i_sub > 1) {
          P_tmp <- t(t(P_tmp) - P_tmp[1, ] + Psub2[index2[i_sub], ])
        }
        Psub2[idxs2, ] <- fbead(P_tmp, S2 = Pframe[i_sub:(i_sub + 1L), ])
      }
      Psub2 <- tranS(Psub2, S2 = P)
      Logsub <- floglike(cbind(pbin, Psub2), A0, b0, invSigma0, beta1, covmat0, mat)
      cat("LLK after suboptimization2: ", Logsub, "\n")
      if (Logsub >= usLoglike) {
        P <- Psub2
        usLoglike <- Logsub
      }
    }

    Pf <- fcorrect(pbin, P)
    Loglike <- usLoglike
    if (any(Pf != P)) {
      tmp_opt <- optim(c(1, beta1), fn = function(x) -floglike(cbind(pbin, x[1] * Pf), A, b, invSigma, x[-1], covmat0, mat, pos, v, mak))
      sLoglike <- -tmp_opt$value
      cat("LLK before and after outlier removal:", c(usLoglike, sLoglike), "(", tmp_opt$par, ")", "\n")
      if (usLoglike <= sLoglike) {
        Loglike <- sLoglike
        beta1 <- tmp_opt$par
        P <- beta1[1] * matrix(Pf, nrow = N, ncol = 3L)
        beta1 <- beta1[-1]
      }
    }
    Pf <- avsmth(pbin, P = P)
    beta2 <- sapply(Beta, FUN = function(x) x[length(x)])
    if (N < 1000) {
      tmp_opt <- optim(c(1, beta1), fn = function(x) -floglike(cbind(pbin, x[1] * Pf), A, b, invSigma, x[-1], covmat0, mat, pos, v, mak))
      sLoglike <- -tmp_opt$value
      cat("LLK before and after smoothing:", c(Loglike, sLoglike), "(", tmp_opt$par, ")", "\n")
      if (!(Loglike > sLoglike || (any(beta2 > -1) && runif(1) < 0.5))) {
        Loglike <- sLoglike
        beta1 <- tmp_opt$par
        P <- beta1[1] * matrix(Pf, nrow = N, ncol = 3L)
        beta1 <- beta1[-1]
      }
    } else {
      sLoglike <- floglike(cbind(pbin, Pf), A, b, invSigma, beta1, covmat0, mat, pos, v, mak)
      cat("LLK before and after smoothing:", c(Loglike, sLoglike), "\n")
      if (!(Loglike > sLoglike || (any(beta2 > -1) && runif(1) < 0.5))) {
        Loglike <- sLoglike
        P <- matrix(Pf, nrow = N, ncol = 3L)
      }
    }

    beta1 <- pmin(beta1, ifelse(iternum > 10, -0.6, -1))
    P01 <- P
    if (iternum == (num_coarse - 1) && coarsefit) {
      P01 <- matrix(Pf, nrow = N, ncol = 3L)
    }
    write.table(cbind(bin, P), file = paste0(outfile, "temp.txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
    recodllk[(iternum %% 7) + 1L] <- Loglike > Loglike0
    if (Loglike > Loglike0) {
      Loglike0 <- Loglike
      Ao <- A
      bo <- b
      invSigmao <- invSigma
      Po <- P
      if (rmoutlier) {
        Po <- fcorrect(pbin, Po)
      }
      betao <- Beta
      cat(bo, Ao, invSigmao, Loglike0, "\n")
      write.table(cbind(bin, Po), file = paste0(outfile, ".txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
      write.table(unlist(betao), file = paste0(outfile, "beta.txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
    }
    # if(iternum>6&& !sum(recodllk)){P01=matrix(Pf,N,3)+matrix(rnorm(3*N,0,sqrt(5/N)),N,3)}

    iternum <- iternum + 1
  }
  list(Ao, bo, invSigmao, Po, betao)
}
