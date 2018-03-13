# Publications that use results obtained from this software please include a citation of the paper:
# Zou, C., Zhang, Y., Ouyang, Z. (2016) HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure. Submitted.

library(stats) # filter(), glm(), median(), rnorm(), splinefun()
library(MASS)  # ginv()

normP <- function(P) {
  P1 <- t(t(P) - colMeans(P))
  P1 <- 5 * P1 / sqrt(max(rowSums(P1^2)))
  P1
}

kinetic0 <- function(p0) {
  sum(p0^2) / dim(p0)[1]
}

kinetic0_1 <- function(p0) {
  sum(p0^2)
}

kinetic <- function(p0, N, rho) {
  (sum(p0^2) - sum(p0[-1, ] * p0[-N, ]) * rho * 2) / N
}

kinetic_1 <- function(p0, N, rho) {
  (sum(p0^2) - sum(p0[-1, ] * p0[-N, ]) * rho * 2)
}

momentum <- function(p0, N, rho) {
  p0 - rbind(p0[-1, ] * rho, rep(0, times = 3)) - rbind(rep(0, times = 3), p0[-N, ] * rho)
}

momentum0 <- function(p0) {
  p0
}

#' @importFrom stats rnorm splinefun
finistructure <- function(S0, bin) {
  n <- dim(S0)[1]
  N <- length(bin[, 1])
  if (n == N) {
    if (dim(S0)[2] == 3) {
      S <- S0
    } else {
      S <- S0[, 3:5]
    }
  } else {
    pts <- c(S0[, 1], S0[n, 2])
    Y <- as.matrix(rbind(S0[, 3:5], rnorm(3, mean = as.numeric(S0[n, 3:5]), sd = as.numeric(apply(S0[-1, 3:5] - S0[-n, 3:5], MARGIN = 2, FUN = sd)))))
    S <- normP(sapply(1:3, FUN = function(x) splinefun(pts, y = Y[, x])(bin[, 1])))
    S <- S + matrix(rnorm(3 * N, mean = 0, sd = sqrt(5 / N)), nrow = N, ncol = 3)
  }
  S
}

fmkorder_temp <- function(m, A, b, sigma, S) {
  if (m < 2) {
    mu <- A %*% S + b
    Sigma <- sigma
    return(as.matrix(cbind(mu, Sigma, A)))
  } else {
    tmp <- Recall(m - 1, A, b, sigma, S)
    mu <- A %*% tmp[, 1] + b
    temp <- tmp[, -c(1:4)] %*% A
    Sigma <- tmp[, 2:4] + (t(temp)) %*% sigma %*% (temp)
    return(as.matrix(cbind(mu, Sigma, temp)))
  }
}

fmkorder <- function(m, A, b, sigma, S) {
  temp <- fmkorder_temp(m, A, b, sigma, S)
  temp[, 2:4] <- solve(temp[, 2:4])
  temp
}

fmkorder2 <- function(m, A, b, sigma) {
  fmkorder(m, A, b, sigma, rep(0, times = 3))
}

fnormvec <- function(a, b) {
  c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
}

fangle <- function(a, b) {
  acos(sum(a * b) / sqrt(sum(a^2)) / sqrt(sum(b^2)))
}

frotanyvec <- function(x, v, theta) {
  n <- v / sqrt(sum(v^2))
  cos(theta) * (x - sum(x * n) * n) + sin(theta) * fnormvec(n, x) + sum(x * n) * n
}

fbead <- function(S1, S2) {
  m <- dim(S1)[1]
  # S=S1*sqrt(sum((S2[1,]-S2[2,])^2)/sum((S1[1,]-S1[m,])^2))
  S <- S1
  n <- fnormvec(S[m, ] - S[1, ], S2[2, ] - S2[1, ])
  theta <- fangle(S[m, ] - S[1, ], S2[2, ] - S2[1, ])
  S <- t(t(S) - S[1, ])
  S <- t(apply(S, MARGIN = 1, FUN = function(x) frotanyvec(x, n, theta)))
  S <- t(t(S) - S[1, ] + S1[1, ])
  S
}

fmirror <- function(v) {
  y <- diag(3)
  if (any(v)) {
    y <- diag(3) - v %*% t(v) / sum(v^2)
  }
  y
}

#' @importFrom MASS ginv
tranS <- function(S1, S2, I_scale = TRUE) {
  tmp <- cbind(1, S1)
  beta <- ginv(t(tmp) %*% tmp) %*% (t(tmp) %*% S2)
  # beta=solve(t(tmp)%*%tmp,t(tmp)%*%S2)
  s <- svd(beta[-1, ])
  beta[-1, ] <- s$u %*% (diag(sign(s$d))) %*% t(s$v)
  S <- tmp %*% beta
  if (I_scale) {
    beta[-1, ] <- mean(abs(s$d)) * s$u %*% (diag(sign(s$d))) %*% t(s$v)
    tmp <- optim(c(1, 0, 0, 0, 0, 0, 0), fn = function(x) sum(sqrt(rowSums((t(t(x[1] * S %*% angle2mat(x[2:4])) + x[5:7]) - S2)^2))))
    tmp <- tmp$par
    S <- t(t(tmp[1] * S %*% angle2mat(tmp[2:4])) + tmp[5:7])
  } else {
    tmp <- optim(rep(0, times = 9), fn = function(x) sum(sqrt(rowSums((t(t(S %*% angle2mat(x[4:6]) %*% fmirror(x[1:3])) + x[7:9]) - S2)^2))))
    tmp <- tmp$par
    S <- t(t(S %*% angle2mat(tmp[4:6]) %*% fmirror(tmp[1:3])) + tmp[7:9])
  }
  S
}

#' @importFrom stats median splinefun
rmol <- function(loci, P) {
  n <- dim(P)[1]
  m <- dim(P)[2]
  P1 <- P
  v <- sapply(1:m, FUN = function(x) is.na(P1[, x]) | is.infinite(P1[, x]))
  if (any(v)) {
    outlier <- v[, 1] | v[, 2] | v[, 3]
    spf <- splinefun
    P2 <- P[!outlier, ]
    tmp <- sapply(1:m, FUN = function(x) spf(loci[!outlier], y = P2[, x])(loci[outlier]))
    # print(tmp)
    P1[outlier, ] <- tmp
  }
  d1 <- sqrt(rowSums((P1[-1, ] - P1[-n, ])^2))
  cutoff <- 10 * median(d1)
  if (any(d1 >= cutoff)) {
    outlier <- d1 >= cutoff
    v <- which(outlier)
    for (i in 1:length(v)) {
      P1[1:v[i], ] <- t(t(P1[1:v[i], ]) + (P[v[i] + 1, ] - P[v[i], ]) * 0.8)
    }
    P1 <- tranS(P1, P)
  }
  # plot3d(P1[,1],P1[,2],P1[,3],type="l")
  # points3d(P1[!outlier,1],P1[!outlier,2],P1[!outlier,3],col='green')
  # points3d(P[outlier,1],P[outlier,2],P[outlier,3],col='blue')
  # points3d(P1[outlier,1],P1[outlier,2],P1[outlier,3],col='red')
  P1
}

#' @importFrom stats filter
avsmth <- function(bin, P) {
  N <- length(bin)
  # dP=sqrt(rowSums((P[-1,]-P[-N,])^2))
  P0 <- filter(P, filter = rep(1, times = 3) / 3, sides = 2)
  P0[1, ] <- (P[1, ] + P[2, ]) / 2
  P0[N, ] <- (P[N, ] + P[N - 1, ]) / 2
  P0 <- matrix(P0, nrow = N, ncol = 3)
}

loglikelihood0 <- function(P0, A, b, invSigma, beta, cx, mat, pos=NULL, v=NULL, mak=NULL) {
  L <- 0
  P <- P0[, -1]
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.na(C)) {
    C <- 1
    cx <- array(cx, dim = c(dim(cx), 1))
    mat <- array(mat, dim = c(dim(mat), 1))
  }
  if (is.null(pos)) {
    pos <- apply(mat, MARGIN = 3, FUN = function(x) which(!is.na(x[, 1])))
    if (!is.list(pos)) {
      pos <- lapply(1:nrow(t(pos)), FUN = function(i) t(pos)[i, ])
    }
  }
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0))
  }
  distmat <- as.matrix(dist(P))
  for (i in 1:C) {
    temp <- -cx[pos[[i]], pos[[i]], i] * distmat[pos[[i]], pos[[i]]]^beta[i]
    temp[v[[i]]] <- temp[v[[i]]] + mat[pos[[i]], pos[[i]] + 1, i][v[[i]]] * (beta[i] * log(distmat[pos[[i]], pos[[i]]][v[[i]]]) + log(cx[pos[[i]], pos[[i]], i][v[[i]]]))
    L <- L + sum(temp[upper.tri(temp)]) / N / 3 # +sum(temp[lower.tri(temp))
  }
  L
}

dloglikelihood0 <- function(P0, A, b, invSigma, beta, cx, mat, pos, v=NULL, mak=NULL) {
  P <- as.matrix(P0[, -1])
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0))
  }
  # pos=apply(mat,3,function(x) which(!is.na(x[,1])))
  # if(!is.list(pos)){pos=lapply(1:nrow(t(pos)),function(i) t(pos)[i,])}
  dL <- matrix(0, nrow = N, ncol = 3)
  distmat <- as.matrix(dist(P))^2 # apply(P*P,1,sum)%*%t(rep(1,N))+rep(1,N)%*%t(apply(P*P,1,sum))-2*P%*%t(P)
  temp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:C) {
    temp[pos[[i]], pos[[i]]] <- temp[pos[[i]], pos[[i]]] - beta[i] * cx[pos[[i]], pos[[i]], i] * (distmat[pos[[i]], pos[[i]]]^(beta[i] / 2 - 1))
    temp[pos[[i]], pos[[i]]][v[[i]]] <- temp[pos[[i]], pos[[i]]][v[[i]]] + beta[i] * mat[pos[[i]], pos[[i]] + 1, i][v[[i]]] / distmat[pos[[i]], pos[[i]]][v[[i]]]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dL[, 1] <- colSums(t(tmp) - tmp) / N / 3
  tmp <- temp * P[, 2]
  dL[, 2] <- colSums(t(tmp) - tmp) / N / 3
  tmp <- temp * P[, 3]
  dL[, 3] <- colSums(t(tmp) - tmp) / N / 3
  dL
}

loglikelihood <- function(P0, A, b, invSigma, beta, cx, mat, pos=NULL, v=NULL, mak=NULL) {
  L <- 0
  P <- as.matrix(P0[, -1])
  sigma <- solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  if (is.na(C)) {
    C <- 1
    cx <- array(cx, dim = c(dim(cx), 1))
    mat <- array(mat, dim = c(dim(mat), 1))
  }
  if (is.null(pos)) {
    pos <- apply(mat, MARGIN = 3, FUN = function(x) which(!is.na(x[, 1])))
    if (!is.list(pos)) {
      pos <- lapply(1:nrow(t(pos)), FUN = function(i) t(pos)[i, ])
    }
  }
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0))
  }
  distmat <- as.matrix(dist(P))
  for (i in 1:C) {
    temp <- -cx[pos[[i]], pos[[i]], i] * distmat[pos[[i]], pos[[i]]]^beta[i]
    temp[v[[i]]] <- temp[v[[i]]] + mat[pos[[i]], pos[[i]] + 1, i][v[[i]]] * (beta[i] * log(distmat[pos[[i]], pos[[i]]][v[[i]]]) + log(cx[pos[[i]], pos[[i]], i][v[[i]]]))
    L <- L + sum(temp[upper.tri(temp)]) / N / 3
    # temp=-cx[pos[[i]],pos[[i]],i]*distmat[pos[[i]],pos[[i]]]^beta[i]/N/3+mat[pos[[i]],pos[[i]]+1,i]*(beta[i]*log(distmat[pos[[i]],pos[[i]]])+log(cx[pos[[i]],pos[[i]],i]))/N/3
    # L=L+sum(temp[upper.tri(temp)])# +sum(temp[lower.tri(temp))
  }
  if (is.null(mak)) {
    L <- L + sum(sapply(2:N, FUN = function(ii) {
      nmp <- fmkorder(P0[ii, 1] - P0[ii - 1, 1], A, b, sigma, P[ii - 1, ])
      -(P[ii, ] - nmp[, 1]) %*% (nmp[, 2:4]) %*% (P[ii, ] - nmp[, 1]) / 2 + log(det(nmp[, 2:4])) / 2
    }) / N / 3)
  } else {
    L <- L + sum(sapply(2:N, FUN = function(ii) {
      mu <- mak[[ii - 1]][, 1] + mak[[ii - 1]][, 5:7] %*% P[ii - 1, ]
      -t(P[ii, ] - mu) %*% mak[[ii - 1]][, 2:4] %*% (P[ii, ] - mu) / 2 + log(det(mak[[ii - 1]][, 2:4])) / 2
    })) / N / 3
  }
  L
}

dloglikelihood <- function(P0, A, b, invSigma, beta, cx, mat, pos, v=NULL, mak=NULL) {
  P <- P0[, -1]
  sigma <- solve(invSigma)
  # pos=apply(mat,3,function(x) which(!is.na(x[,1])))
  # if(!is.list(pos)){pos=lapply(1:nrow(t(pos)),function(i) t(pos)[i,])}
  N <- dim(P)[1]
  C <- dim(cx)[3]
  dL <- matrix(0, nrow = N, ncol = 3)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0))
  }
  distmat <- as.matrix(dist(P))^2
  # distmat=apply(P*P,1,sum)%*%t(rep(1,N))+rep(1,N)%*%t(apply(P*P,1,sum))-2*P%*%t(P)
  temp <- matrix(0, nrow = N, ncol = N)
  for (i in 1:C) {
    temp[pos[[i]], pos[[i]]] <- temp[pos[[i]], pos[[i]]] - beta[i] * cx[pos[[i]], pos[[i]], i] * (distmat[pos[[i]], pos[[i]]]^(beta[i] / 2 - 1))
    temp[pos[[i]], pos[[i]]][v[[i]]] <- temp[pos[[i]], pos[[i]]][v[[i]]] + beta[i] * mat[pos[[i]], pos[[i]] + 1, i][v[[i]]] / distmat[pos[[i]], pos[[i]]][v[[i]]]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dL[, 1] <- colSums(t(tmp) - tmp) / N / 3
  tmp <- temp * P[, 2]
  dL[, 2] <- colSums(t(tmp) - tmp) / N / 3
  tmp <- temp * P[, 3]
  dL[, 3] <- colSums(t(tmp) - tmp) / N / 3
  if (is.null(mak)) {
    nmp <- lapply(2:N, FUN = function(ii) fmkorder(P0[ii, 1] - P0[ii - 1, 1], A, b, sigma, P[ii - 1, ]))
    temp1 <- sapply(2:N, FUN = function(ii) -t(P[ii, ] - nmp[[ii - 1]][, 1]) %*% nmp[[ii - 1]][, 2:4]) / N / 3
    dL[2:N, ] <- dL[2:N, ] + t(temp1)
    temp2 <- sapply(1:(N - 1), FUN = function(ii) -t(nmp[[ii]][, 5:7]) %*% nmp[[ii]][, 2:4] %*% nmp[[ii]][, 5:7] %*% P[ii, ] - t(nmp[[ii]][, 5:7]) %*% nmp[[ii]][, 2:4] %*% (nmp[[ii]][, 1] - nmp[[ii]][, 5:7] %*% P[ii, ] - P[ii + 1, ])) / N / 3
    dL[1:(N - 1), ] <- dL[1:(N - 1), ] + t(temp2)
  } else {
    temp1 <- sapply(2:N, FUN = function(ii) {
      mu <- mak[[ii - 1]][, 1] + mak[[ii - 1]][, 5:7] %*% P[ii - 1, ]
      -t(P[ii, ] - mu) %*% mak[[ii - 1]][, 2:4]
    }) / N / 3
    dL[2:N, ] <- dL[2:N, ] + t(temp1)
    temp2 <- sapply(1:(N - 1), FUN = function(ii) -t(mak[[ii]][, 5:7]) %*% mak[[ii]][, 2:4] %*% mak[[ii]][, 5:7] %*% P[ii, ] - t(mak[[ii]][, 5:7]) %*% mak[[ii]][, 2:4] %*% (mak[[ii]][, 1] - P[ii + 1, ])) / N / 3
    dL[1:(N - 1), ] <- dL[1:(N - 1), ] + t(temp2)
  }
  dL
}

mkcloglikelihood <- function(theta, P0) {
  P <- P0[, -1]
  N <- dim(P)[1]
  A <- matrix(theta[1:9], nrow = 3, ncol = 3)
  b <- theta[10:12]
  invSigma <- matrix(0, nrow = 3, ncol = 3)
  invSigma[upper.tri(invSigma, diag = TRUE)] <- theta[-c(1:12)]
  invSigma <- invSigma + t(invSigma) - diag(diag(invSigma))
  sigma <- solve(invSigma)
  temp <- sapply(2:N, FUN = function(ii) {
    nmp <- fmkorder(P0[ii, 1] - P0[ii - 1, 1], A, b, sigma, P[ii - 1, ])
    -(P[ii, ] - nmp[, 1]) %*% (nmp[, 2:4]) %*% (P[ii, ] - nmp[, 1]) / 2 + log(det(nmp[, 2:4])) / 2
  })
  mean(temp) / 3
}

dhllk <- function(index, theta, P0, A, b, invSigma, beta, cx, mat, pos, v=NULL) {
  P <- rbind(P0[index[1, 1]:index[1, 2], -1], do.call(rbind, args = sapply(2:dim(index)[1], FUN = function(x) t(t(P0[index[x, 1]:index[x, 2], -1] %*% matrix(theta[x - 1, 1:9], nrow = 3, ncol = 3)) + theta[x - 1, 10:12]))))
  P <- as.matrix(P)
  # sigma=solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  # cat(dim(theta),";")
  dL <- matrix(0, nrow = dim(theta)[1], ncol = 12)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0))
  }
  distmat <- as.matrix(dist(P))^2
  temp <- matrix(0, nrow = N, ncol = N)
  dD <- array(0, dim = c(N, N, 3))
  for (i in 1:C) {
    temp[pos[[i]], pos[[i]]] <- temp[pos[[i]], pos[[i]]] - beta[i] * cx[pos[[i]], pos[[i]], i] * (distmat[pos[[i]], pos[[i]]]^(beta[i] / 2 - 1))
    temp[pos[[i]], pos[[i]]][v[[i]]] <- temp[pos[[i]], pos[[i]]][v[[i]]] + beta[i] * mat[pos[[i]], pos[[i]] + 1, i][v[[i]]] / distmat[pos[[i]], pos[[i]]][v[[i]]]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }

  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dD[, , 1] <- t(tmp) - tmp
  tmp <- temp * P[, 2]
  dD[, , 2] <- t(tmp) - tmp
  tmp <- temp * P[, 3]
  dD[, , 3] <- t(tmp) - tmp
  tmp <- t(apply(index[-1, ], MARGIN = 1, FUN = function(x) c(colSums(dD[-(x[1]:x[2]), x[1]:x[2], 1] %*% P0[x[1]:x[2], -1]), colSums(dD[-(x[1]:x[2]), x[1]:x[2], 2] %*% P0[x[1]:x[2], -1]), colSums(dD[-(x[1]:x[2]), x[1]:x[2], 3] %*% P0[x[1]:x[2], -1])))) / N / 3
  dL[, 1:9] <- tmp # t(apply(index[-1,],1,function(x) c(colSums(dD[-(x[1]:x[2]),x[1]:x[2],1]%*%P0[x[1]:x[2],-1]),colSums(dD[-(x[1]:x[2]),x[1]:x[2],2]%*%P0[x[1]:x[2],-1]),colSums(dD[-(x[1]:x[2]),x[1]:x[2],3]%*%P0[x[1]:x[2],-1]))))/N/3
  dL[, 1:9] <- dL[, 1:9] + t(apply(index[-1, ], MARGIN = 1, FUN = function(x) c(sapply(2:4, FUN = function(k) c(sum(upper.tri(t(dD[x[1]:x[2], x[1]:x[2], 1] * P0[x[1]:x[2], k])) - upper.tri(dD[x[1]:x[2], x[1]:x[2], 1] * P0[x[1]:x[2], k])), sum(upper.tri(t(dD[x[1]:x[2], x[1]:x[2], 2] * P0[x[1]:x[2], k])) - upper.tri(dD[x[1]:x[2], x[1]:x[2], 2] * P0[x[1]:x[2], k])), sum(upper.tri(t(dD[x[1]:x[2], x[1]:x[2], 3] * P0[x[1]:x[2], k])) - upper.tri(dD[x[1]:x[2], x[1]:x[2], 3] * P0[x[1]:x[2], k]))))))) / N / 3

  dL[, 10] <- apply(index[-1, ], MARGIN = 1, FUN = function(x) sum(dD[-(x[1]:x[2]), x[1]:x[2], 1])) / N / 3
  dL[, 11] <- apply(index[-1, ], MARGIN = 1, FUN = function(x) sum(dD[-(x[1]:x[2]), x[1]:x[2], 2])) / N / 3
  dL[, 12] <- apply(index[-1, ], MARGIN = 1, FUN = function(x) sum(dD[-(x[1]:x[2]), x[1]:x[2], 3])) / N / 3
  dL
}

rotamat <- function(theta) {
  temp <- matrix(0, nrow = 3, ncol = 5)
  temp[, 4:5] <- matrix(c(sin(theta[1:3]), cos(theta[1:3])), nrow = 3, ncol = 2)
  temp[, 1:3] <- matrix(c(
     temp[1, 5] * temp[3, 5] - temp[2, 5] * temp[1, 4] * temp[3, 4],
    -temp[2, 5] * temp[3, 5] * temp[1, 4] - temp[1, 5] * temp[3, 4],
     temp[1, 4] * temp[2, 4],
     temp[3, 5] * temp[1, 4] + temp[1, 5] * temp[2, 5] * temp[3, 4],
     temp[1, 5] * temp[2, 5] * temp[3, 5] - temp[1, 4] * temp[3, 4],
    -temp[1, 5] * temp[2, 4],
     temp[2, 4] * temp[3, 4],
     temp[3, 5] * temp[2, 4],
     temp[2, 5]
  ), nrow = 3, ncol = 3)
  temp
}

dDtotheta <- function(p, b, temp) {
  l <- matrix(0, nrow = dim(b)[1], ncol = 3)
  l[, 1] <- -2 * (p[2] * (-temp[2, 2]) - p[1] * temp[1, 2] - temp[3, 2] * p[3]) * (b[, 1] - p[1] * temp[1, 1] + p[2] * (-temp[2, 1]) - temp[3, 1] * p[3]) - 2 * (p[1] * temp[1, 1] - p[2] * (-temp[2, 1]) + temp[3, 1] * p[3]) * (b[, 2] - p[1] * temp[1, 2] + p[2] * (-temp[2, 2]) - temp[3, 2] * p[3])

  l[, 2] <- 2 * (temp[1, 5] * temp[3, 3] * p[3] + temp[1, 5] * temp[2, 3] * p[2] + temp[1, 5] * temp[1, 3] * p[1]) * (b[, 2] - p[1] * temp[1, 2] + p[2] * (-temp[2, 2]) - temp[3, 2] * p[3]) - 2 * (temp[3, 3] * temp[1, 4] * p[3] + temp[3, 5] * temp[3, 1] * p[2] + temp[1, 4] * temp[1, 3] * p[1]) * (b[, 1] - p[1] * temp[1, 1] + p[2] * (-temp[2, 1]) - temp[3, 1] * p[3]) + 2 * (temp[3, 3] * temp[3, 5] * p[2] - temp[2, 4] * p[3] + temp[3, 3] * temp[3, 4] * p[1]) * (temp[3, 3] * p[3] - b[, 3] + temp[2, 3] * p[2] + temp[1, 3] * p[1])

  l[, 3] <- 2 * (temp[2, 3] * p[1] - temp[1, 3] * p[2]) * (temp[3, 3] * p[3] - b[, 3] + temp[2, 3] * p[2] + temp[1, 3] * p[1]) + 2 * (p[1] * (-temp[2, 2]) + p[2] * temp[1, 2]) * (b[, 2] - p[1] * temp[1, 2] + p[2] * (-temp[2, 2]) - temp[3, 2] * p[3]) + 2 * (p[1] * (-temp[2, 1]) + p[2] * temp[1, 1]) * (b[, 1] - p[1] * temp[1, 1] + p[2] * (-temp[2, 1]) - temp[3, 1] * p[3])

  # l[1]=- 2*(p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + cos(theta[1])*sin(theta[2])*p[3])*(b[1] - p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) + p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) - sin(theta[1])*sin(theta[2])*p[3]) - 2*(p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) - p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) + sin(theta[1])*sin(theta[2])*p[3])*(b[2] - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + cos(theta[1])*sin(theta[2])*p[3])

  # l[2]=2*(cos(theta[1])*cos(theta[2])*p[3] + cos(theta[1])*cos(theta[3])*sin(theta[2])*p[2] + cos(theta[1])*sin(theta[2])*sin(theta[3])*p[1])*(b[2] - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + cos(theta[1])*sin(theta[2])*p[3]) - 2*(cos(theta[2])*sin(theta[1])*p[3] + cos(theta[3])*sin(theta[1])*sin(theta[2])*p[2] + sin(theta[1])*sin(theta[2])*sin(theta[3])*p[1])*(b[1] - p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) + p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) - sin(theta[1])*sin(theta[2])*p[3]) + 2*(cos(theta[2])*cos(theta[3])*p[2] - sin(theta[2])*p[3] + cos(theta[2])*sin(theta[3])*p[1])*(cos(theta[2])*p[3] - b[3] + cos(theta[3])*sin(theta[2])*p[2] + sin(theta[2])*sin(theta[3])*p[1])

  # l[3]=2*(cos(theta[3])*sin(theta[2])*p[1] - sin(theta[2])*sin(theta[3])*p[2])*(cos(theta[2])*p[3] - b[3] + cos(theta[3])*sin(theta[2])*p[2] + sin(theta[2])*sin(theta[3])*p[1]) + 2*(p[1]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + p[2]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])))*(b[2] - p[1]*(cos(theta[3])*sin(theta[1]) + cos(theta[1])*cos(theta[2])*sin(theta[3])) + p[2]*(sin(theta[1])*sin(theta[3]) - cos(theta[1])*cos(theta[2])*cos(theta[3])) + cos(theta[1])*sin(theta[2])*p[3]) + 2*(p[1]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) + p[2]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])))*(b[1] - p[1]*(cos(theta[1])*cos(theta[3]) - cos(theta[2])*sin(theta[1])*sin(theta[3])) + p[2]*(cos(theta[1])*sin(theta[3]) + cos(theta[2])*cos(theta[3])*sin(theta[1])) - sin(theta[1])*sin(theta[2])*p[3])
  l
}

dhllk1 <- function(index, theta, P0, A, b, invSigma, beta, cx, mat, pos, v=NULL) {
  matheta <- lapply(2:dim(index)[1], FUN = function(x) rotamat(theta[x - 1, ]))
  P <- rbind(P0[index[1, 1]:index[1, 2], -1], do.call(rbind, args = sapply(2:dim(index)[1], FUN = function(x) t(t(P0[index[x, 1]:index[x, 2], -1] %*% matheta[[x - 1]][, 1:3]) + theta[x - 1, 4:6]))))
  P <- as.matrix(P)
  # sigma=solve(invSigma)
  N <- dim(P)[1]
  C <- dim(cx)[3]
  # cat(dim(theta),";")
  dL <- matrix(0, nrow = dim(theta)[1], ncol = 6)
  if (is.null(v)) {
    v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0))
  }
  distmat <- as.matrix(dist(P))^2
  temp <- matrix(0, nrow = N, ncol = N)
  dD <- array(0, dim = c(N, N, 3))
  for (i in 1:C) {
    temp[pos[[i]], pos[[i]]] <- temp[pos[[i]], pos[[i]]] - beta[i] * cx[pos[[i]], pos[[i]], i] * (distmat[pos[[i]], pos[[i]]]^(beta[i] / 2 - 1))
    temp[pos[[i]], pos[[i]]][v[[i]]] <- temp[pos[[i]], pos[[i]]][v[[i]]] + beta[i] * mat[pos[[i]], pos[[i]] + 1, i][v[[i]]] / distmat[pos[[i]], pos[[i]]][v[[i]]]
    # temp[pos[[i]],pos[[i]]]=temp[pos[[i]],pos[[i]]]-beta[i]*cx[pos[[i]],pos[[i]],i]*(distmat[pos[[i]],pos[[i]]]^(beta[i]/2-1))+beta[i]*mat[pos[[i]],pos[[i]]+1,i]/distmat[pos[[i]],pos[[i]]]
  }
  diag(temp) <- 0
  tmp <- temp * P[, 1]
  dD[, , 1] <- t(tmp) - tmp
  tmp <- temp * P[, 2]
  dD[, , 2] <- t(tmp) - tmp
  tmp <- temp * P[, 3]
  dD[, , 3] <- t(tmp) - tmp
  dL[, 1:3] <- t(sapply(2:dim(index)[1], FUN = function(x) rowSums(sapply(index[x, 1]:index[x, 2], FUN = function(i) colSums(dDtotheta(P0[i, -1], t(t(P[-(index[x, 1]:index[x, 2]), ]) - theta[x - 1, 4:6]), matheta[[x - 1]]) * temp[-(index[x, 1]:index[x, 2]), i])))))
  dL[, 4] <- apply(index[-1, ], MARGIN = 1, FUN = function(x) sum(dD[-(x[1]:x[2]), x[1]:x[2], 1]))
  dL[, 5] <- apply(index[-1, ], MARGIN = 1, FUN = function(x) sum(dD[-(x[1]:x[2]), x[1]:x[2], 2]))
  dL[, 6] <- apply(index[-1, ], MARGIN = 1, FUN = function(x) sum(dD[-(x[1]:x[2]), x[1]:x[2], 3]))
  dL / N / 3
}

angle2mat <- function(theta) {
  rotamat(theta)[, 1:3]
}

v2mat <- function(theta) {
  matrix(theta[1:9], nrow = 3, ncol = 3)
}

Leapfrog <- function(grad_U, L, epsilon, p0, q0, fM) {
  N <- dim(p0)[1]
  if (L < 2) {
    # q = (q0 + epsilon * p0)
    q <- (q0 + epsilon * (fM(p0)))
    p <- p0 - epsilon * grad_U(q) / 2

    return(list(p, q))
  } else {
    temp <- Recall(grad_U, L - 1, epsilon, p0, q0, fM)
    # q=temp[[2]]+epsilon*(temp[[1]])
    q <- (temp[[2]] + epsilon * (fM(temp[[1]])))
    p <- temp[[1]] - epsilon * grad_U(temp[[2]]) / 2
    return(list(p, q))
  }
}


#' @importFrom stats rnorm
HMC <- function(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans=FALSE) {
  N <- dim(current_q0)[1]
  m <- dim(current_q0)[2]
  current_q <- current_q0
  q <- current_q
  if (I_trans) {
    q <- rmol(1:N, current_q)
    # q=tranS(q,current_q)
  }
  # p = array(0,dim(q))
  p <- array(rnorm(N * m, mean = 0) / sqrt(N), dim = dim(q)) # independent standard normal variates
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
    return(rmol(1:N, current_q))
  } else {
    if (runif(1) < exp((const) / T0)) {
      if (I_trans) {
        q <- rmol(1:N, q)
        q <- tranS(q, current_q)
      }
      return(q) # accept
    } else {
      return(current_q) # reject
    }
  }
}


#' @importFrom stats rnorm
HMC1 <- function(U, grad_U, epsilon, L, current_q0, T0, fK, fM, I_trans=FALSE) {
  N <- dim(current_q0)[1]
  m <- dim(current_q0)[2]
  current_q <- current_q0
  q <- current_q
  # p = array(0,dim(q))
  p <- array(rnorm(N * m, mean = 0) / sqrt(N), dim = dim(q)) # independent standard normal variates
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
    return(rmol(1:N, current_q))
    # return (matrix(rnorm(N*m,current_q,abs(current_q)/15),N,m))
  } else {
    if (runif(1) < exp((const) / T0 / const)) {
      if (I_trans) {
        q <- rmol(1:N, q)
        q <- tranS(q, current_q)
      }
      return(q) # accept
    } else {
      return(current_q) # reject
    }
  }
}


#' @importFrom stats rnorm
Qsis <- function(N, pbin, A, b, invSigma, beta, cx, mat, q0, fL) {
  # cat(c(length(pbin[1:N]),dim(mat[1:N,1:(N+1),]),dim(cx[1:N,1:N,])),"~")
  if (N == 2) {
    Padd <- optim(b + rnorm(3) / 5 + q0, fn = function(q) fL(cbind(pbin[1:N], rbind(q0, q)), A, b, invSigma, beta, cx[1:N, 1:N, ], mat[1:N, 1:(N + 1), ]))
    # cat(Padd$par,"\n")
    return(rbind(q0, t(Padd$par)))
  } else {
    # cat(pbin,"\n")
    temp <- Recall(N - 1, pbin, A, b, invSigma, beta, cx, mat, q0, fL)
    Padd <- optim(rnorm(3) / 5 + A %*% temp[dim(temp)[1], ] + b, fn = function(q) fL(cbind(pbin[1:N], rbind(temp, q)), A, b, invSigma, beta, cx[1:N, 1:N, ], mat[1:N, 1:(N + 1), ]))
    # cat(Padd$par,"\n")
    return(rbind(temp, t(Padd$par)))
  }
}

Sis <- function(d, pbin, A, b, invSigma, beta, cx0, mat0, q0, fL) {
  N <- dim(mat0)[1]
  # cat(length(pbin),",")
  if (is.na(dim(cx0)[3])) {
    cx <- array(cx0, dim = c(dim(cx0), 1))
    mat <- array(mat0, dim = c(dim(mat0), 1))
  } else {
    cx <- cx0
    mat <- mat0
  }
  if (N == d) {
    # cat(c(length(pbin),dim(mat),dim(cx)),"~")
    return(as.matrix(Qsis(d, pbin, A, b, invSigma, beta, cx, mat, q0, fL)))
  } else {
    # cat(c(length(pbin[-N]),dim(mat[-N,-(N+1),])),"\n")
    temp <- Recall(d, pbin[-N], A, b, invSigma, beta, cx[-N, -N, ], mat[-N, -(N + 1), ], q0, fL)
    # cat(dim(temp),"\n")
    addq <- optim(rnorm(3) / 5 + A %*% temp[dim(temp)[1], ] + b, fn = function(x) fL(cbind(pbin[(N - d):N], rbind(temp[(nrow(temp) - d + 1):nrow(temp), ], x)), A, b, invSigma, beta, cx[(N - d):N, (N - d):N, ], mat[(N - d):N, c(1, (N - d + 1):(N + 1)), ]))
    addq <- addq$par
    return(as.matrix(rbind(temp, t(addq))))
  }
}

subinitial <- function(pbin, A0, b0, invSigma0, beta1, covmat0, mat, floglike, fdloglike) {
  N <- length(pbin)
  if (is.na(dim(covmat0)[3])) {
    covmat0 <- array(covmat0, dim = c(dim(covmat0), 1))
    mat <- array(mat, dim = c(dim(mat), 1))
  }
  pos <- apply(mat, MARGIN = 3, FUN = function(x) which(!is.na(x[, 1])))
  if (!is.list(pos)) {
    pos <- lapply(1:nrow(t(pos)), FUN = function(i) t(pos)[i, ])
  }
  P0 <- Sis(4, pbin, A0, b0, invSigma0, beta1, covmat0, mat, c(0, 0, 0), function(x, ...) -loglikelihood(x, ...))
  u <- 0
  while (u < 100) {
    P <- HMC(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos), 0.002, 20, P0, 10 * exp(-u / 20), function(p) kinetic(p, N, 0.1), function(p) momentum(p, N, 0.1), 1)
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
    covmat0 <- array(covmat0, dim = c(dim(covmat0), 1))
    mat <- array(mat, dim = c(dim(mat), 1))
  }
  pos <- apply(mat, MARGIN = 3, FUN = function(x) which(!is.na(x[, 1])))
  if (!is.list(pos)) {
    pos <- lapply(1:nrow(t(pos)), FUN = function(i) t(pos)[i, ])
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
  theta <- matrix(theta0, nrow = 4, ncol = 3)
  -floglike(cbind(pbin, rbind(P1, t(t(P2 %*% theta[-4, ]) + theta[4, ]))), A0, b0, invSigma0, beta1, covmat0, mat)
}


#' @importFrom stats rnorm
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
    index[-1, 1] <- index[-1, 1] + 1
    lP <- lapply(1:dim(index)[1], FUN = function(x) subinitial(pbin[index[x, 1]:index[x, 2]], A0, b0, invSigma0, beta1, covmat0[index[x, 1]:index[x, 2], index[x, 1]:index[x, 2], ], mat[index[x, 1]:index[x, 2], c(1, index[x, 1]:index[x, 2] + 1), ], floglike, fdloglike))
    P <- matrix(0, nrow = N, ncol = 3)
    P[index[1, 1]:index[1, 2], ] <- lP[[1]]
    for (i in 2:dim(index)[1])
    {
      theta <- optim(as.vector(rbind(diag(3), P[index[i - 1, 2], ] - P[index[i, 1], ] + rnorm(3) / 100)), fn = function(x) piece(x, P[1:(index[i - 1, 2]), ], lP[[i]], pbin[1:index[i, 2]], A0, b0, invSigma0, beta1, covmat0[1:index[i, 2], 1:index[i, 2], ], mat[1:index[i, 2], c(1:(index[i, 2] + 1)), ], floglike))
      theta <- matrix(theta$par, nrow = 4, ncol = 3)
      # cat(dim(theta),"\t")
      P[index[i, 1]:index[i, 2], ] <- lP[[i]] %*% theta[-4, ] + rep(1, times = index[i, 2] - index[i, 1] + 1) %*% t(theta[4, ])
    }
  }
  as.matrix(avsmth(pbin, P))
  # return(P)
}

## lsmap0 contact map
## lscov0 covariance file
## output output name
## Maxiter control parameter
## submaxiter control parameter
## lambda control parameter
## Leapfrom control parameter
## epslon control parameter
## mkfix control parameter defaulted to 0
## rho control parameter defaulted to 0
## mk Markov parameter 0 or 1
## initialS initial S
## coarsefit control parameter defaulted to TRUE
## rmoutfiler remove outlier, defaulted to FALSE
## fitmode control parameter defaulted to 0

#' @importFrom stats glm
fmain <- function(lsmap0, lscov0, output, Maxiter, submaxiter, lamda, Leapfrog, epslon, mkfix=0, rho=0, mk, initialS=NULL, coarsefit=TRUE, rmoutlier=FALSE, fitmode=0) {
  floglike <- loglikelihood0 ## loglikelihood0 and dloglikelihood0 are functions
  fdloglike <- dloglikelihood0
  fcorrect <- rmol ## rmol is a function, appears to be related to outliers
  recodllk <- rep(0, times = 7)
  if (mk) {
    floglike <- loglikelihood ## if Markov then different likehood functions called, minus 0s
    fdloglike <- dloglikelihood
  }
  C <- length(lsmap0) ## how many contact matrices?
  bin <- NULL
  for (i in 1:C) { ## stack the contract map bins on top of each other
    bin <- rbind(bin, lsmap0[[i]][, 1:2])
  }
  bin <- unique(bin, MARGIN = 1) ## finds unique rows
  bin <- bin[order(bin[, 1], decreasing = FALSE), ] ## orders the bins in ascending order
  N <- dim(bin)[1] ## how many rows, that is, bins
  mbin <- mean(bin[, 2] - bin[, 1]) ## average bin size
  neigdc <- max(c(1, pmax(floor((bin[-1, 1] - bin[-N, 1]) / mbin), 1))) ## biggest jump between bins (scaled by average jump)
  pbin <- cumsum(c(1, pmax(floor((bin[-1, 1] - bin[-N, 1]) / mbin), 1))) ## sum of jumps between bins (scaled by average jump)
  gldata <- vector("list", length = C) ## list of size of contact matrices
  A0 <- diag(3) ## 3-dimensional identity matrix
  b0 <- c(0, 0, 0)
  t <- c(0:(N - 1))
  invSigma0 <- diag(3) ## 3-dimensional identity matrix
  invS <- diag(3) * sqrt(N) ## diagonal matrix consisting of the square root of the number of bins
  mat <- array(NA, dim = c(N, N + 1, C)) ## 3-dimensional array of N by N+1 by C all filled with NAs
  mat0 <- mat ## Same NA matrix as mat
  beta1 <- -rep(1.3, times = C) ## the number 1.3 repeated C times
  alpha <- seq(from = 0.5, to = 1.5, by = 0.5) ## the squence 0.5, 1, 1.5
  covmat0 <- array(1, dim = c(N, N, C)) ## a N by N by C matrix of 1s
  Beta <- vector("list", length = C) ## list same size as contact matrix
  impute <- 0
  pos <- vector("list", length = C) ## list same size as contact matrix
  if (is.numeric(lscov0)) {
    for (c in 1:C) { ## iterate over contact matrices
      temp <- which(bin[, 1] %in% lsmap0[[c]][, 1]) ## which bins are in the particular contact matrix
      mat[temp, 1, c] <- temp ## update mat with whether the bin is contained in contact matrix c
      pos[[c]] <- temp ## update pos with whether the bin is contained in contact matrix c
      # cat(temp,"\n")
      temp <- as.matrix(lsmap0[[c]][, -c(1, 2)]) ## temp updated to contain the contact matrices without the bins positions
      if (isSymmetric(temp)) { ## If symmetric contact matrix
        mat[pos[[c]], pos[[c]] + 1, c] <- temp ## update mat
        gldata[[c]] <- temp[upper.tri(temp)] ## gldata is upper triangular version of temp
      } else { ## If not symmetric contact matrix
        mat[pos[[c]], pos[[c]] + 1, c] <- temp + t(temp) ## make mat symmetric if not already
        temp <- temp + t(temp) ## make temp symmetric
        gldata[[c]] <- temp[upper.tri(temp)] ## gldata is upper triangular version of temp, same so should come out of loop
      }
      mat0[, , c] <- mat[, , c] ## third dimension of mat0 is same as mat
      gldata[[c]] <- data.frame(cbind(temp[upper.tri(temp)], rep(1, times = length(gldata[[c]])))) ## overwrite previous version of gldata, so maybe should be updated
    }

    # P10=finital(pbin,A0,b0,invSigma0,beta1,covmat0,mat,floglike,fdloglike)
    if (is.null(initialS)) { ## if not initial S
      P10 <- finital(pbin, A0, b0, invS, beta1, covmat0, mat, floglike, fdloglike) ## P10 is created
    } else {
      P10 <- finistructure(initialS, bin) ## If initial S P10 is created using bin info
    }
    P10 <- normP(P10) ## P10 is normalized
    P01 <- P10 ## P01 is same as P10
    dmat <- as.matrix(dist(P10)) ## distance matrix of P10, which is formed by initial S
    for (c in 1:C) { ## iterate over contact matrices
      dmat1 <- dmat[pos[[c]], pos[[c]] ] ## distance between bins that actually exists within contact matrices
      dmat1 <- dmat1[upper.tri(dmat1)] ## upper triangular version of dmat1
      gldata[[c]][, dim(gldata[[c]])[2]] <- log(dmat1) ## gldata gets log of distance matrix
      colnames(gldata[[c]]) <- paste0("V", 1:dim(gldata[[c]])[2]) ## gldata columns get named
      betaest <- glm(V1 ~ ., family = poisson(), data = gldata[[c]]) ## first column of gldata is regressed on everything else using Poisson regression
      Beta[[c]] <- as.vector(betaest$coefficients) ## Beta contains coefficients
      cat(Beta[[c]], "\n") ## print betas
      beta1[c] <- min(-abs(Beta[[c]][length(Beta[[c]])]), beta1[c]) ## beta1 contains smallest betas
      covmat0[, , c] <- covmat0[, , c] * exp(Beta[[c]][1]) ## covmat0 is adjusted with new Betas
    }
  } else { ## if covariance matrix is not nummeric
    beta11 <- beta1 ## beta11 is set to 1.3, why?
    # covmat00=array(1,c(N,N,C))
    lscov <- lscov0 ## set cov to initial cov
    for (c in 1:C) { ## loop through contact matrix loops similar, perhaps exactly the same
      temp <- which(bin[, 1] %in% lsmap0[[c]][, 1])
      # cat(temp)
      mat[temp, 1, c] <- temp
      pos[[c]] <- temp
      temp <- as.matrix(lsmap0[[c]][, -c(1, 2)])
      if (isSymmetric(temp)) {
        mat[pos[[c]], pos[[c]] + 1, c] <- temp
      } else {
        mat[pos[[c]], pos[[c]] + 1, c] <- temp + t(temp)
        temp <- temp + t(temp)
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
        for (k in 2:length(Beta[[c]][-1])) {
          covmat0[pos[[c]], pos[[c]], c] <- covmat0[pos[[c]], pos[[c]], c] * lscov[[c]][[k - 1]]^Beta[[c]][k]
        }
      }
    }

    if (is.null(initialS)) {
      P10 <- finital(pbin, A0, b0, invS, beta1, covmat0, mat, floglike, fdloglike)
    } else {
      P10 <- finistructure(initialS, bin)
    }
    P01 <- P10
    dmat <- as.matrix(dist(P10))
    covmat0 <- array(1, dim = c(N, N, C)) ## covmat0 re-written as 3-dimensional array of 1s
    for (c in 1:C) { ## why a second round of Poisson regression? This loop same as above
      dmat1 <- dmat[pos[[c]], pos[[c]]] ## distance between bins that actually exists within contact matrices
      dmat1 <- dmat1[upper.tri(dmat1)] ## upper triangular version of dmat1
      gldata[[c]][, dim(gldata[[c]])[2]] <- log(dmat1) ## gldata gets log of distance matrix
      colnames(gldata[[c]]) <- paste0("V", 1:dim(gldata[[c]])[2]) ## gldata columns get named
      betaest <- glm(V1 ~ ., family = poisson(), data = gldata[[c]]) ## first column of gldata is regressed on everything else using Poisson regression
      Beta[[c]] <- as.vector(betaest$coefficients) ## Beta contains coefficients
      cat(Beta[[c]], "\n") ## print betas
      beta1[c] <- min(-abs(Beta[[c]][length(Beta[[c]])]), beta1[c]) ## beta1 contains smallest betas
      covmat0[, , c] <- covmat0[, , c] * exp(Beta[[c]][1]) ## covmat0 is adjusted with new Betas
      if (!is.numeric(lscov)) { ## updating covariance matrix but only is not numeric lscov and lscov[[c]], why?
        if (!is.numeric(lscov[[c]])) {
          for (k in 2:length(Beta[[c]][-1])) {
            covmat0[pos[[c]], pos[[c]], c] <- covmat0[pos[[c]], pos[[c]], c] * lscov[[c]][[k - 1]]^Beta[[c]][k]
          }
        }
      }
    }
  }
  Loglike0 <- floglike(cbind(pbin, P10), A0, b0, invSigma0, beta1, covmat0, mat) ## update likelihood
  cat("LLK0: ", Loglike0, "\n") ## print likelihood
  cat("number of nodes:", N, "\n") ## print number of bins
  v <- lapply(1:C, FUN = function(i) which(mat[pos[[i]], pos[[i]] + 1, i] > 0)) ## list of positions greater than 0 between adjacent bins
  if (mk) { ## if Markov property
    mak <- lapply(2:N, FUN = function(ii) fmkorder2(pbin[ii] - pbin[ii - 1], A0, b0, solve(invSigma0))) ## Update Markov matrix
    # lapply(mak,function(x) cat(dim(x),";"))
  } else {
    mak <- NULL ## Markov matrix is null if not Markov property
  }
  if (N > 1000) { ## hard-coded step if over 1000 bins
    m <- 100
    if (N > 2000) {
      m <- 200
    }
    index <- cbind(c(1, seq(from = m, to = N - m, by = m)), c(seq(from = m, to = N - m, by = m), N))
    index[-1, 1] <- index[-1, 1] + 1
    m2 <- floor(N / 100)
    index2 <- seq(from = 1, to = N - m2, by = m2)
    index2 <- c(index2, N)
    Psub <- P01
    Psub <- do.call(rbind, args = lapply(1:dim(index)[1], FUN = function(x) tranS(suboptimz(pbin[index[x, 1]:index[x, 2]], P01[index[x, 1]:index[x, 2], ], A0, b0, invSigma0, beta1, covmat0[index[x, 1]:index[x, 2], index[x, 1]:index[x, 2], ], mat[index[x, 1]:index[x, 2], c(1, index[x, 1]:index[x, 2] + 1), ], floglike, fdloglike, 0), P01[index[x, 1]:index[x, 2], ]))) ## complex rbind
    Logsub <- floglike(cbind(pbin, Psub), A0, b0, invSigma0, beta1, covmat0, mat) ## updated likelihood
    cat("LLK after suboptimization: ", Logsub, "\n") ## print update likelihood
    if (Loglike0 <= Logsub) { ## If first sub-optimization is an improvement update
      P01 <- Psub
      Loglike0 <- Logsub
    }
    Pframe <- tranS(suboptimz(pbin[index2], P01[index2, ], A0, b0, invSigma0, beta1, covmat0[index2, index2, ], mat[index2, c(1, index2 + 1), ], floglike, fdloglike, 0), P01[index2, ])
    Psub2 <- Psub
    for (i_sub in 1:(length(index2) - 1)) {
      P_tmp <- P01[index2[i_sub]:index2[i_sub + 1], ]
      if (i_sub > 1) {
        P_tmp <- t(t(P_tmp) - P_tmp[1, ] + Psub2[index2[i_sub], ])
      }
      Psub2[index2[i_sub]:index2[i_sub + 1], ] <- fbead(P_tmp, Pframe[i_sub:(i_sub + 1), ])
    }
    Psub2 <- tranS(Psub2, P01)
    # plot3d(Psub2[,1],Psub2[,2],Psub2[,3],type="l")
    # points3d(Psub2[,1],Psub2[,2],Psub2[,3],col='red')
    # lines3d(P01[,1],P01[,2],P01[,3],col='green')
    Logsub <- floglike(cbind(pbin, Psub2), A0, b0, invSigma0, beta1, covmat0, mat)
    cat("LLK after suboptimization2: ", Logsub, "\n")
    if (Loglike0 <= Logsub) { ## If second sub-optimization is an improvement update
      P01 <- Psub2
      Loglike0 <- Logsub
    }
  } ## End of update loop if greater than 1000 bins

  Ao <- A0 ## Update parameters
  bo <- b0
  invSigmao <- invSigma0
  Po <- P01
  betao <- Beta
  iternum <- 0
  P <- P01
  write.table(cbind(bin, normP(Po)), file = paste0(output, ".txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE) ## Write out parameters
  write.table(unlist(betao), file = paste0(output, "beta.txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE) ## Write out betas
  if (fitmode) { ## default to not fit mode
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
  while (iternum < Maxiter) { ## Optimization loop, everything to here is setup
    u <- 0
    if (iternum < 20) { ## first twenty do one thing
      submaxiter1 <- max(submaxiter + 50, 150)
      lamda1 <- max(submaxiter + 50, 50)
    } else {
      if (iternum < 0.8 * Maxiter) { ## from 20 to 80% of Maxiter -1 do another
        submaxiter1 <- submaxiter
        lamda1 <- lamda
      } else { ## from 80% of Maxiter to end do another
        submaxiter1 <- min(submaxiter * 1.1, 300)
        lamda1 <- min(lamda * 1.1, 100)
      }
    }

    current_epslon <- ifelse(iternum < num_coarse && coarsefit, eps_coarse, epslon)

    while (u < submaxiter && (N < 1000 || (iternum %% 10 == 1 && iternum > 5) || u < 3)) { ## likelihood optimization not totally understood
      P <- fHMC(function(x) -floglike(cbind(pbin, x), A0, b0, invSigma0, beta1, covmat0, mat, pos, v, mak), function(y) -fdloglike(cbind(pbin, y), A0, b0, invSigma0, beta1, covmat0, mat, pos, v, mak), current_epslon, Leapfrog, P01, exp(-u / lamda), function(p) fknt(p, N, rho), function(p) momentum(p, N, rho))
      P01 <- P
      # plot3d(P01[,1],P01[,2],P01[,3],type="l")
      # points3d(P01[,1],P01[,2],P01[,3],col='red')
      u <- u + 1
    }

    if (!mk) {
      P <- normP(P)
    } ## If not Markov property normalize P

    if (mkfix & iternum >= 5) { ## Furth optimization not totally understood
      thetam <- optim(c(as.vector(A0), b0, as.vector(invSigma0[upper.tri(invSigma0, diag = TRUE)])), fn = function(theta) -mkcloglikelihood(theta, cbind(pbin, P)))
      thetam <- thetam$par
      # cat(thetam,"\n")
      A <- matrix(thetam[1:9], nrow = 3, ncol = 3)
      b <- thetam[10:12]
      invSigma <- matrix(0, nrow = 3, ncol = 3)
      invSigma[upper.tri(invSigma, diag = TRUE)] <- thetam[-c(1:12)]
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
    dmat <- as.matrix(dist(P))
    for (c in 1:C) { ## yet another Poisson regression fit separately by contact matrix
      dmat1 <- dmat[pos[[c]], pos[[c]]]
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
            covmat0[pos[[c]], pos[[c]], c] <- covmat0[pos[[c]], pos[[c]], c] * lscov[[c]][[k - 1]]^Beta[[c]][k]
          }
        }
      }
      # temp=covmat0[,,c]*(dmat^beta1[c])
      # mat[,,c][!mat0[,,c]]=temp[!mat0[,,c]]
    }
    usLoglike <- floglike(cbind(pbin, P), A, b, invSigma, beta1, covmat0, mat)
    cat("LLK after GLM: ", usLoglike, "\n")
    if (N > 1000 && iternum < (Maxiter - 5)) {
      Psub <- do.call(rbind, args = lapply(1:dim(index)[1], FUN = function(x) tranS(suboptimz(pbin[index[x, 1]:index[x, 2]], P[index[x, 1]:index[x, 2], ], A0, b0, invSigma0, beta1, covmat0[index[x, 1]:index[x, 2], index[x, 1]:index[x, 2], ], mat[index[x, 1]:index[x, 2], c(1, index[x, 1]:index[x, 2] + 1), ], floglike, fdloglike, 0), P[index[x, 1]:index[x, 2], ])))
      Logsub <- floglike(cbind(pbin, Psub), A0, b0, invSigma0, beta1, covmat0, mat)
      cat("LLK after suboptimization: ", Logsub, "\n")
      if (Logsub >= usLoglike) {
        P <- Psub
        usLoglike <- Logsub
      }
      Pframe <- tranS(suboptimz(pbin[index2], P[index2, ], A0, b0, invSigma0, beta1, covmat0[index2, index2, ], mat[index2, c(1, index2 + 1), ], floglike, fdloglike, 0), P[index2, ])
      # Psub2=do.call(rbind,lapply(1:(length(index2)-1),function(x)fbead(Psub[index2[x]:index2[x+1],],Pframe[x:(x+1),])[-1,]))
      # Psub2=tranS(rbind(Pframe[1,],Psub2),P)
      Psub2 <- P
      for (i_sub in 1:(length(index2) - 1)) {
        P_tmp <- P[index2[i_sub]:index2[i_sub + 1], ]
        if (i_sub > 1) {
          P_tmp <- t(t(P_tmp) - P_tmp[1, ] + Psub2[index2[i_sub], ])
        }
        Psub2[index2[i_sub]:index2[i_sub + 1], ] <- fbead(P_tmp, Pframe[i_sub:(i_sub + 1), ])
      }
      Psub2 <- tranS(Psub2, P)
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
        P <- beta1[1] * matrix(Pf, nrow = N, ncol = 3)
        beta1 <- beta1[-1]
      }
    }
    Pf <- avsmth(pbin, P)
    beta2 <- sapply(Beta, FUN = function(x) x[length(x)])
    if (N < 1000) { ## fewer than 1000 bins do one optimization
      tmp_opt <- optim(c(1, beta1), fn = function(x) -floglike(cbind(pbin, x[1] * Pf), A, b, invSigma, x[-1], covmat0, mat, pos, v, mak))
      sLoglike <- -tmp_opt$value
      cat("LLK before and after smoothing:", c(Loglike, sLoglike), "(", tmp_opt$par, ")", "\n")
      if (!(Loglike > sLoglike || (any(beta2 > -1) && runif(1) < 0.5))) {
        Loglike <- sLoglike
        beta1 <- tmp_opt$par
        P <- beta1[1] * matrix(Pf, nrow = N, ncol = 3)
        beta1 <- beta1[-1]
      }
    } else { ## greater than 1000 bins do another optimization
      sLoglike <- floglike(cbind(pbin, Pf), A, b, invSigma, beta1, covmat0, mat, pos, v, mak)
      cat("LLK before and after smoothing:", c(Loglike, sLoglike), "\n")
      if (!(Loglike > sLoglike || (any(beta2 > -1) && runif(1) < 0.5))) {
        Loglike <- sLoglike
        P <- matrix(Pf, nrow = N, ncol = 3)
      }
    }

    beta1 <- pmin(beta1, ifelse(iternum > 10, -0.6, -1))
    P01 <- P
    if (iternum == (num_coarse - 1) && coarsefit) {
      P01 <- matrix(Pf, nrow = N, ncol = 3)
    }
    write.table(cbind(bin, P), file = paste0(output, "temp.txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
    recodllk[(iternum %% 7) + 1] <- Loglike > Loglike0
    if (Loglike > Loglike0) { ## If Loglike is greater than Loglike0, update
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
      write.table(cbind(bin, Po), file = paste0(output, ".txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
      write.table(unlist(betao), file = paste0(output, "beta.txt"), sep = "\t", eol = "\n", row.names = FALSE, col.names = FALSE)
    }
    # if(iternum>6&& !sum(recodllk)){P01=matrix(Pf,N,3)+matrix(rnorm(3*N,0,sqrt(5/N)),N,3)}

    iternum <- iternum + 1
  } ## end of maximization
  list(Ao, bo, invSigmao, Po, betao) ## parameters that are output of function but not written as output
}
