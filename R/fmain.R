# Publications that use results obtained from this software please include a citation of the paper:
# Zou, C., Zhang, Y., Ouyang, Z. (2016) HSA: integrating multi-track Hi-C data for genome-scale reconstruction of 3D chromatin structure. Submitted.

#' @importFrom stats glm optim poisson runif
#' @export
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
  A0 <- diag(3)
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
    lscov <- lscov0
    for (c in 1:C) {
      pos_c <- which(bin_c1 %in% lsmap0[[c]][, 1])
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
  write_tsv(cbind(bin, normP(Po)), file = paste0(outfile, ".txt"))
  write_tsv(unlist(betao), file = paste0(outfile, "beta.txt"))
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
      u <- u + 1
    }

    if (!mk) {
      P <- normP(P)
    }

    if (mkfix && iternum >= 5) {
      P0 <- cbind(pbin, P)
      thetam <- optim(c(as.vector(A0), b0, as.vector(invSigma0[upper.tri(invSigma0, diag = TRUE)])), fn = function(theta) -mkcloglikelihood(theta, P0 = P0))
      thetam <- thetam$par

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
    write_tsv(cbind(bin, P), file = paste0(outfile, "temp.txt"))
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
      write_tsv(cbind(bin, Po), file = paste0(outfile, ".txt"))
      write_tsv(unlist(betao), file = paste0(outfile, "beta.txt"))
    }

    iternum <- iternum + 1
  }
  list(Ao, bo, invSigmao, Po, betao)
}